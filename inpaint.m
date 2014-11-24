% Retrieve the missing pixels of an image.
% The goal of this script is to inpaint the missing pixels of an image. It
% does so by constructing a patch graph of the known pixels. According to
% some priority, it then iteratively connect unknown patches to the graph.
% A global optimization is run in the end.

close all; clear; clc;
gsp_start();

% Experiment parameters.
imtype = 'lena3'; % Type of line.
imsize = 100; % Image size.
holesize = 20; % Hole size.
plot = true;
savefig = false;

%% Inpainting algorithm

gparam = giin_default_parameters();
[img, obsimg, vertices] = giin_image(imtype, imsize, holesize);
[G, pixels, patches] = giin_patch_graph(obsimg, gparam, plot);
[G, pixels, Pstructure, Pinformation] = giin_inpaint(G, pixels, patches, gparam, plot);
sol = giin_global(G, obsimg, gparam);

%% Visualize graph with image signal

if plot
    % Signal graph.
    fig = figure(); %#ok<UNRCH>

    % Any negative value is a missing pixel --> red.
    cmap = [1,0,0;gray];
    colormap(fig, cmap);
    param.climits = [-1/(length(cmap)-1),1];
    
    param.colorbar = 0;
%     param.vertex_highlight = connected; % draw by hand in different colors instead
%     param.show_edges = true; % very slow

    gsp_plot_signal(G, pixels, param);
    if savefig, saveas(gcf,['results/',imtype,'_patch_graph.png']); end

    clear param fig cmap
end

%% Visualize priorities

if plot
    % Show some vertices of interest.
    giin_plot_priorities(vertices, G, gparam, savefig);
    
    % Plot the various priorities.
    figure();
    subplot(2,2,1);
%     imshow(imadjust(reshape(Pstructure,imsize,imsize)));
    imshow(reshape(Pstructure,imsize,imsize) / max(Pstructure));
    title('Structure priority');
    subplot(2,2,2);
    imshow(reshape(Pinformation(:,1), imsize, imsize));
    title('Pixel infomation priority');
    subplot(2,2,4);
    imshow(reshape(Pinformation(:,2), imsize, imsize));
    title('Patch infomation priority');
    subplot(2,2,3);
    imshow(reshape(Pstructure .* Pinformation(:,2), imsize, imsize) / max(Pstructure .* Pinformation(:,2)));
    title('Global priority');
    colormap(hot);
end

%% Results

% Images.
figure();
subplot(2,2,1);
imshow(img);
title('Original');
subplot(2,2,2);
imshow(reshape(obsimg,imsize,imsize));
title('Masked');
subplot(2,2,3);
imshow(reshape(pixels,imsize,imsize));
title('Inpainted');
subplot(2,2,4);
imshow(reshape(sol,imsize,imsize));
title(['Globally optimized (',gparam.optim.prior,')']);
saveas(gcf,'results/inpainting_last.png');

% Reconstruction errors.
% fprintf('Observed image error (L2-norm) : %f\n', norm(reshape(img,[],1) - y));
fprintf('Inpainting reconstruction error : %f\n', norm(reshape(img,[],1) - pixels));
fprintf('Globally optimized reconstruction error : %f\n', norm(reshape(img,[],1) - sol));