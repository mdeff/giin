% Regularity analysis of the perfect graph, the reconstructed graph and a
% dumb graph.

close all; clear; clc;
gsp_start();

% Experiment parameters.
imtype = 'lena3';
plot = false;
savefig = false;

gparam = giin_default_parameters();
[img, obsimg, imsize] = giin_image(imtype, true);

%% Perfect graph

Gperfect = giin_patch_graph(img, gparam, false);

%% Disconnected graph

[Gdisc, pixels, patches] = giin_patch_graph(obsimg, gparam, false);

%% Reconstructed graph

Grec = giin_inpaint(Gdisc, pixels, patches, gparam, plot);

%% Dumb graph

W = Gdisc.W;

% Connect the unconnected vertices in a grid.
weight = 0.5;
for k = 1:size(patches,1)
    if any(patches(k,:) < 0)
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

%% Analysis

% Regularity measure : x L x^T
regPerfect = img(:).' * Gperfect.L * img(:);
regDisc = img(:).' * Gdisc.L * img(:);
regRec = img(:).' * Grec.L * img(:);
regDumb = img(:).' * Gdumb.L * img(:);

titles = {sprintf('Perfect graph : %f',regPerfect), ...
    sprintf('Disconnected graph : %f',regDisc), ...
    sprintf('Reconstructed graph : %f',regRec), ...
    sprintf('Dumb graph : %f',regDumb)};

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

graphs = {Gperfect, Gdisc, Grec, Gdumb};

fig = figure();
for k = 1:4
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