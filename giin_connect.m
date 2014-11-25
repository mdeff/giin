function [ Gnew ] = giin_connect( Gold, vertices, knowns, patches, gparam )
%GIIN_CONNECT Connect patches in a graph
%   Create a new graph given a list of patches to insert in the old one.

tic;

knn = gparam.graph.knn;
spi = zeros(knn*length(vertices),1);
spj = zeros(knn*length(vertices),1);
spv = zeros(knn*length(vertices),1);

% Connect the considered patches.
for ii = 1:length(vertices)
    vertex = vertices(ii);
    
    % Custom-made comparison mask.
    M = patches(vertex,:)>=0;
    
    % Find the K closest patches. One at a time because we need a custom
    % comparison function or a custom dictionary.
    dict = repmat(M,size(knowns)) .* patches(knowns,:);
    [idx,dist] = knnsearch(dict, M.*patches(vertex,:), 'K',knn, 'Distance','euclidean');
    neighbors = knowns(idx);
    
    % Find the K closest patches.
%     e = nan(Gold.N, 1);
%     for known = knowns.'
%         e(known) = norm(M.*patches(vertex,:) - M.*patches(known,:), 2);
%     end
%     [~,idx] = sort(e);
%     neighbors = idx(1:knn);
%     dist = e(neighbors);
    
    % Fill the 3-col values with [i, j, exp(-d(i,j)^2 / sigma)]
    start = (ii-1)*knn+1;
    spi(start:start+knn-1) = vertex;
    spj(start:start+knn-1) = neighbors;
    spv(start:start+knn-1) = exp(-dist.^2 / gparam.graph.sigma);
end

% The new connections.
W = sparse(spi, spj, spv, size(Gold.W,1), size(Gold.W,2));

% Do not symmetrize the graph because we don't have tested the known
% patches agains the considered : would be unfair.
W = (W + W.') / 2;

% G.W = G.W + W;
% G = gsp_graph_default_parameters(G);
% G = gsp_create_laplacian(G); % Laplacian was not updated
Gnew = gsp_graph(Gold.W+W, Gold.coords, Gold.plotting.limits);

if Gold.Ne + nnz(W) ~= Gnew.Ne
    error('Some of the new connections were already there !');
end

% Execution time.
% fprintf('giin_connect : %f seconds\n', toc);

end