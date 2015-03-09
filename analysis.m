% Regularity analysis of the perfect graph, the reconstructed graph and a
% dumb graph.

close all; clear; clc;
gsp_start();

% Experiment parameters.
imtype = 'lena3';
plot = false;
savefig = false;

addpath('./lib');
addpath('./data');

gparam = giin_default_parameters();
[img, obsimg, imsize] = giin_image(imtype, true);

%% Perfect graph

[Gperfect, pixels, patches]  = giin_patch_graph(img, gparam, false);

%% Disconnected graph

[Gdisc, obspixels, obspatches] = giin_patch_graph(obsimg, gparam, false);

%% Reconstructed graph

Grec = giin_inpaint(Gdisc, obspixels, obspatches, gparam, plot);

%% Dumb graph

W = Gdisc.W;

% Connect the unconnected vertices in a grid.
weight = 1;  % Has very limited impact.
for k = 1:size(obspatches,1)
    if any(obspatches(k,:) < 0)
        W(k,k+1) = weight;      % bottom
        W(k,k-1) = weight;      % top
        W(k,k+imsize) = weight; % right
        W(k,k-imsize) = weight; % left
    end
end

% Symmetrize the weights.
W = (W + W.') / 2;

% Construct the graph object.
Gdumb = gsp_graph(W, Gdisc.coords, Gdisc.plotting.limits);

% Sanity check : all(all(Gdisc.W==Gdumb.W))

%% Analysis

% Nathanaël example :
% N = 500;
% X = rand(N,25); % 25 = 5*5 pixels
% G = gsp_sensor(N);
% G = gsp_create_laplacian( G,'normalized' );
% sum(gsp_norm_tik(G,X))

% Mutli-dimensional signal : sum of the individual norms.
% sum(gsp_norm_tik(G,X))
% sum(sum(X .* (G.L* X) ))
% reg=0; for k=1:size(X,2), reg=reg + X(:,k).'*G.L*X(:,k); end; reg

% Normalize Laplacians.
graphs = {Gperfect, Gdisc, Grec, Gdumb};
graphs = gsp_create_laplacian( graphs,'normalized' );

% Regularity measure : x L x^T
regularity = zeros(length(graphs),2);
for k = 1:length(graphs)
    regularity(k,1) = sum(gsp_norm_tik(graphs{k}, pixels));
    regularity(k,2) = sum(gsp_norm_tik(graphs{k}, patches));
end
% Similar : regularity(:,1) ./ regularity(:,2)

titles = {sprintf('Perfect graph : %.2f (%.2d)',regularity(1,1),regularity(1,2)), ...
    sprintf('Disconnected graph : %.2f (%.2d)',regularity(2,1),regularity(2,2)), ...
    sprintf('Reconstructed graph : %.2f (%.2d)',regularity(3,1),regularity(3,2)), ...
    sprintf('Dumb graph : %.2f (%.2d)',regularity(4,1),regularity(4,2))};

for str = titles
    disp(str{1})
end

%% Graphs visualization
% Show the connections around the center of the graph.

% Vertex numbers.
vertices = 1:imsize^2;
vertices = reshape(vertices,imsize,imsize);

% Vertex numbers of interest.
holesize = round(imsize / 4);
vizsize = holesize + 10;
xyrange = floor((imsize-vizsize)/2)+1 : imsize-ceil((imsize-vizsize)/2);
vertices = vertices(xyrange,xyrange);
vertices = vertices(:);

% New coordinates.
% coords = Gdumb.coords(vertices,:);
coords = zeros(vizsize^2, 2);  % gsp_plot_signal takes [x,y]
coords(:,2) = repmat((vizsize:-1:1).', vizsize, 1);
coords(:,1) = reshape(repmat((1:+1:vizsize), vizsize, 1), [], 1);

% Reduced signal.
sig = obsimg(:);
sig = sig(vertices);

fig = figure();
for k = 1:length(graphs)
    G = graphs{k};
    
    % Reduce the graph.
    W = G.W(vertices,vertices);
    Gviz = gsp_graph(W, coords, [1,vizsize,1,vizsize]);

    % Plot the connections.
    subplot(2,2,k);
    giin_plot_signal(Gviz, sig, true);
    title(titles{k});
end
% saveas(fig,'results/graphs.png');