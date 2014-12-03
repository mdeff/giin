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