% Retrieve the missing pixels of an image.
% The goal of this script is to inpaint the missing pixels of an image. It
% does so by constructing a patch graph of the known pixels. According to
% some priority, it then iteratively connect unknown patches to the graph.
% A global optimization is run in the end.

close all; clear; clc;
gsp_start();

imtype = 'lena3'; % Type of line.
imsize = 100;         % Image size.
holesize = 20;

gparam.psize = 5;  % Patch size.
gparam.knn = 10;  % Patch graph minimum number of connections (KNN).
gparam.sigma = 0.1; % Variance of the distance kernel.
gparam.loc = 0.001;  % Importance of local information. (default 0.001, 0.1)
gparam.priority_threshold = 1e-3; % Threshold when creating priority from diffused energy.
gparam.cheb_order = 30; % Order of the Chebyshev approximation (number of hopes).
gparam.heat_scale = 500; % 1000 for lena
gparam.max_unknown_pixels = gparam.psize; % Maximum number of unknown pixels to connect a patch.
gparam.inpainting_retrieve = 'copy'; % Average connected patches or copy the strongest.
gparam.inpainting_compose = 'overwrite'; % Keep known pixels or overwrite everything.
gparam.optim_prior = 'thikonov'; % Global optimization constraint : thikonov or tv.
gparam.optim_maxit = 100; % Maximum number of iterations.
gparam.optim_sigma = 0; % Noise level.

plot = false;
savefig = false;

%% Image

[img, vertices] = giin_image(imtype, imsize, gparam);

% Unknown pixels are negative (known ones are in [0,1]).
% Negative enough such that they won't connect to anything other than other unknown patches.
imsize = length(img);
bordersize = (imsize-holesize)/2;
xyrange = bordersize+1:imsize-bordersize;
obsimg = img;
obsimg(xyrange,xyrange) = -1e3;
clear imtype holesize bordersize xyrange

%% Patch graph
% Unknown patches are disconnected.

tic;

% We want some local information for the spatially close patches to be
% connected in case of large uniform surfaces.
param.rho = gparam.loc;

param.patch_size = gparam.psize;
param.nnparan.center = 0;
param.nnparam.resize = 0;
param.nnparam.k = gparam.knn;
param.nnparam.sigma = gparam.sigma;
[G, pixels, patches] = gsp_patch_graph(obsimg, param);

% Execution time.
fprintf('Time to create graph : %f seconds\n', toc);

clear param

% Visualize a patch : reshape(patches(4380,1:end-2),psize,psize)

%% Visualize graph with image signal

if plot
    % Signal graph.
    fig = figure();

    % Any negative value is a missing pixel --> red.
    cmap = [1,0,0;gray];
    colormap(fig, cmap);
    param.climits = [-1/(length(cmap)-1),1];
    
    param.colorbar = 0;
%     param.vertex_highlight = connected; % draw by hand in different colors instead
%     param.show_edges = true; % very slow

    gsp_plot_signal(G, pixels, param);
    if savefig, saveas(gcf,[imtype,'_patch_graph.png']); end

    clear param fig cmap
end

%% Iterative inpainting

tic;

% Each unknown pixel has a value of -1e3. A patch with 4 unknown pixels
% will end up with a value of -4. The minimum is -psize^2.
unknowns = (patches<0) .* patches;
unknowns = sum(unknowns,2) / 1000;

% List of new vertices considered for inpainting.
news = find(unknowns<0).';
currents = [];
inpainted = [];

% Patches which contain no other information than their position cannot be
% connected in the non-local graph.

fprintf('There is %d incomplete patches :\n', sum(unknowns<0));
fprintf('  %d without any information\n', sum(unknowns==-gparam.psize^2));
% fprintf('  %d considered for inpainting\n', length(news));

% List of fully known patches to which we can connect.
knowns = find(unknowns==0);
if sum(unknowns<0)+length(knowns) ~= G.N
    error('Missing vertices !');
end

priorities = nan(size(pixels));

% currents  : vertices to be inpainted
% news      : vertices to be connected, i.e. newly considered pixels
% inpainted : vertices already inpainted, to avoid an infinite loop

% Until no more vertices to inpaint.
first = true;
while ~isempty(currents) || first
    first = false;
    
    % Each unknown pixel has a value of -1e3. A patch with 4 unknown pixels
    % will end up with a value of -4. The minimum is -psize^2.
    unknowns = (patches<0) .* patches;
    unknowns = sum(unknowns,2) / 1000;
    
    news = find(unknowns<0).';

    % We only consider patches with less than some number of missing pixels.
    news = news(unknowns(news)>=-gparam.max_unknown_pixels);
    % Which are not already connected.
    news = news(~ismember(news, currents));
    % Neither already visited (to prevent infinite loop and reconnections).
    news = news(~ismember(news, inpainted));
    currents = [currents, news];

    if any(ismember(currents, inpainted))
        error('A vertex could be visited again !');
    end
    
    % What if some patch has no more unknown pixels ?
    % Do we always inpaint over ?

    % Connect the newly reachable vertices.
    G = giin_connect(G, news, knowns, patches, gparam);

    % Compute their priorities, i.e. update the priority signal.
    priorities = giin_priorities(news, priorities, G, gparam);

    % TODO: we also need to take into account the data priority, and normalize
    % the two to give them the same weight.

    % Highest priority patch.
    [~,vertex] = max(priorities);
    priorities(vertex) = -1-priorities(vertex);

    % Update pixels and patches.
    [pixels, patches, news] = giin_inpaint(vertex, G, pixels, patches, gparam);

    % Remove the currently impainted vertex from the lists.
%     news = news(news~=vertex);
    currents = currents(currents~=vertex);
    inpainted = [inpainted, vertex];
    
    % Live plot.
    if plot
        figure(10);
        width = max(G.coords(:,1));
        height = max(G.coords(:,2));
        imshow(reshape(pixels,height,width), 'InitialMagnification',600);
        drawnow;
    end

    fprintf('Inpainted vertices : %d (%d waiting)\n', length(inpainted), length(currents));
end

% Execution time
fprintf('Iterative inpainting : %f\n', toc);

clear unknowns news currents vertex first knowns

%% Visualize priorities

if plot
    vertices = [2238,4370,3493,3589,4380,5703]; % 2238,4370,3600
    giin_plot_priorities(vertices, priorities, G, gparam, savefig);
end

% Restore priorities.
priorities = -1-priorities;
% figure();
% imshow(imadjust(reshape(priorities,imsize,imsize)));
% colormap(hot);

%% Global stage by convex optimization
% Now we inpaint again the image using the created non-local graph.

init_unlocbox();
verbose = 1;

% Observed signal (image).
M = reshape(obsimg>=0, [], 1);
y = M .* reshape(img, [], 1);

% Data term.
% fdata.grad = @(x) 2*M.*(M.*x-y);
% fdata.eval = @(x) norm(M.*x-y)^2;
param_b2.verbose = verbose -1;
param_b2.y = y;
param_b2.A = @(x) M.*x;
param_b2.At = @(x) M.*x;
param_b2.tight = 1;
param_b2.epsilon = gparam.optim_sigma*sqrt(sum(M(:)));
fdata.prox = @(x,T) proj_b2(x,T,param_b2);
fdata.eval = @(x) eps;

% Prior.
param_prior.verbose = verbose-1;
switch(gparam.optim_prior)
    
    % Thikonov prior.
    case 'thikonov'
        fprior.prox = @(x,T) gsp_prox_tik(x,T,G,param_prior);
        fprior.eval = @(x) gsp_norm_tik(G,x);
    
    % TV prior.
    case 'tv'
        G = gsp_adj2vec(G);
        G = gsp_estimate_lmax(G);
        fprior.prox = @(x,T) gsp_prox_tv(x,T,G,param_prior);
        fprior.eval = @(x) gsp_norm_tv(G,x);
end

% Solve the convex optimization problem.
param_solver.verbose = verbose;
param_solver.tol = 1e-7;
param_solver.maxit = gparam.optim_maxit;
tic;
[sol, info] = douglas_rachford(y,fprior,fdata,param_solver);

% Execution time.
fprintf('Global optimization : %f (%d iterations)\n', toc, info.iter);

clear verbose M fdata fprior param_b2 param_prior param_solver info

%% Results

% Images.
subplot(2,2,1);
imshow(img);
title('Original');
subplot(2,2,2);
imshow(reshape(y,imsize,imsize));
title('Masked');
subplot(2,2,3);
imshow(reshape(pixels,imsize,imsize));
title('Inpainted');
subplot(2,2,4);
imshow(reshape(sol,imsize,imsize));
title(['Globally optimized (',gparam.optim_prior,')']);
saveas(gcf,'inpainting_last.png');

% Reconstruction errors.
fprintf('Observed image error (L2-norm) : %f\n', norm(reshape(img,[],1) - y));
fprintf('Inpainting reconstruction error : %f\n', norm(reshape(img,[],1) - pixels));
fprintf('Globally optimized reconstruction error : %f\n', norm(reshape(img,[],1) - sol));