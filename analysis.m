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

%% Reconstructed graph

[Grec, pixels, patches] = giin_patch_graph(obsimg, gparam, false);
Grec = giin_inpaint(Grec, pixels, patches, gparam, plot);

%% Dumb graph

[Gdumb, ~, patches] = giin_patch_graph(obsimg, gparam, false);

% Connect the unconnected vertices in a grid.
for k = 1:size(patches,1)
    if any(patches(k,:) < 0)
        Gdumb.W(k,k+1) = 1;      % bottom
        Gdumb.W(k,k-1) = 1;      % top
        Gdumb.W(k,k+imsize) = 1; % right
        Gdumb.W(k,k-imsize) = 1; % left
    end
end

% Symmetrize the weights.
Gdumb.W = (Gdumb.W + Gdumb.W.') / 2;

%% Analysis

% Regularity measure : x L x^T
regPerfect = img(:).' * Gperfect.L * img(:);
regRec = img(:).' * Grec.L * img(:);
regDumb = img(:).' * Gdumb.L * img(:);

fprintf('Perfect graph : %f\n', regPerfect);
fprintf('Reconstructed graph : %f\n', regRec);
fprintf('Dumb graph : %f\n', regDumb);

%% Graphs visualization

giin_plot_signal(Gperfect, obsimg(:), true)
