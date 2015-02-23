function [ ] = visualize( imname )
%VISUALIZE Visualize the results of an inpaining.
%   Show a comparison of the original image with the inpainted one. Show
%   the priorities.

addpath('./lib');
addpath('./data');
filename = ['results/',imname];
load([filename,'.mat']);

%% Visualize graph with image signal

% for ii = 1:size(pixels,2)
%     figure;
%     giin_plot_signal(G, pixels(:,ii), false);
%     if savefig, saveas(gcf,['results/',imtype,...
%             '_',num2str(ii),'_patch_graph.png']); end
% end

%% Visualize priorities

% Show some vertices of interest.
giin_plot_priorities(vertices, G, gparam, filename);

% Plot the various priorities.
figure();
subplot(2,2,1);
imshow(imadjust(reshape(Pstructure,imsize,imsize)));
% imshow(reshape(Pstructure,imsize,imsize) / max(Pstructure));
title('Structure priority');
subplot(2,2,2);
imshow(reshape(Pinformation(:,1), imsize, imsize));
title('Pixel infomation priority');
subplot(2,2,4);
imshow(reshape(Pinformation(:,2), imsize, imsize));
title('Patch infomation priority');
subplot(2,2,3);
imshow(imadjust(reshape(Pstructure .* Pinformation(:,2),imsize,imsize)));
% imshow(reshape(Pstructure .* Pinformation(:,2), imsize, imsize) / max(Pstructure .* Pinformation(:,2)));
title('Global priority');
colormap(hot);

saveas(gcf,[filename,'_priorities.png']);

%% Results

% Images.
figure();
subplot(2,2,1);
imshow(img);
title('Original');
subplot(2,2,2);
imshow(reshape(obsimg,imsize,imsize,Nc));
title('Masked');
subplot(2,2,3);
imshow(reshape(pixels,imsize,imsize,Nc));
title('Inpainted');
subplot(2,2,4);
imshow(reshape(sol,imsize,imsize,Nc));
title(['Globally optimized (',gparam.optim.prior,')']);
saveas(gcf,[filename,'_inpainting.png']);

% Reconstruction errors.
% fprintf('Observed image error (L2-norm) : %f\n', norm(reshape(img,[],1) - y));
fprintf('Inpainting reconstruction error : %f\n', norm(reshape(img,[],Nc) - pixels));
fprintf('Globally optimized reconstruction error : %f\n', norm(reshape(img,[],Nc) - sol));

end