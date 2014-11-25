function [ G, pixels, patches ] = giin_patch_graph( img, gparam, plot )
%GIIN_PATCH_GRAPH Create a patch graph from an image.
%   Unknown patches are disconnected.
%
%   Examples :
%       [G, pixels, patches] = giin_patch_graph(obsimg, gparam, plot);
%       reshape(patches(4380,1:end-2), gparam.graph.psize, gparam.graph.psize)

tstart = tic;

% We want some local information for the spatially close patches to be
% connected in case of large uniform surfaces.
param.rho = gparam.graph.loc;

param.patch_size = gparam.graph.psize;
param.nnparan.center = 0;
param.nnparam.resize = 0;
param.nnparam.k = gparam.graph.knn;
param.nnparam.sigma = gparam.graph.sigma;
[G, pixels, patches] = gsp_patch_graph(img, param);

% Estimate the maximum eigenvalue of the Lapalacian (for filters).
% Done only once as it does not change much when adding new vertices.
G = gsp_estimate_lmax(G);

% Execution time.
fprintf('Time to create graph : %f seconds\n', toc(tstart));

if plot
    figure();
    hist(G.W(:), -.05:.1:1);
    xlim([eps,1]);
    title(['Graph weights distribution, \sigma=',num2str(gparam.graph.sigma)]);
end

end