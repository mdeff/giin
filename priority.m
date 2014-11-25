% Assign structure priority to every pixel of an image.

close all; clear; clc;
gsp_start();

% Experiment parameters.
imtype = 'lenafull'; % Type of line.
imsize = 100; % Image size.
plot = false;
savefig = false;

gparam = giin_default_parameters();

%% Image and graph

[img, ~, imsize, vertices] = giin_image(imtype, imsize);
G = giin_patch_graph(img, gparam, plot);

%% Priorities

tstart = tic;
Pstructure = nan(G.N, 1);
N = G.N / 1000; % Chunks of 1000 to save runtime memory.
for n = 0:N-1;
    Pstructure = giin_priorities((1:1000)+n*1000, Pstructure, G, gparam);
end
fprintf('Priorities : %f seconds\n', toc(tstart));

if any(isnan(Pstructure))
    error('Some pixels have no assigned priority !');
end

%% Visualize priorities

% Show some vertices of interest.
giin_plot_priorities(vertices, G, gparam, savefig);

% Plot the priority.
figure();
subplot(1,2,1);
imshow(img);
title('Original image');
subplot(1,2,2);
imshow(reshape(Pstructure,imsize,imsize) / max(Pstructure));
title('Structure priority');
% colormap(hot);
saveas(gcf,'results/priority.png');