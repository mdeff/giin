function [ sol ] = inpaint( imname )
%INPAINT Retrieve the missing pixels of an image.
%   Usage :
%       vertices = giin_image('vertical'); 
%       inpaint('vertical');
%
%   Input parameters :
%       imname : name of the image file
%
%   Output parameters :
%       sol    : the inpainted image
%
% The goal of this script is to inpaint the missing pixels of an image. It
% does so by constructing a patch graph of the known pixels. According to
% some priority, it then iteratively connect unknown patches to the graph.
% A global optimization is run in the end.

% Author: Michael Defferrard
% Date: February 2015

addpath('./lib');
addpath('./data');
gsp_start();

img = double(imread([imname,'.png'])) / 255;
Nx = size(img,1);
Ny = size(img,2);
Nc = size(img,3);

% Extract the mask.
if Nc == 1
    mask = img == 1;
else
    mask = img(:,:,1)==0 & img(:,:,2)==1 & img(:,:,3)==0;
    mask = repmat(mask, [1,1,3]);
end

% Unknown pixels are negative (known ones are in [0,1]). Negative enough
% such that they don't connect to anything else than other unknown patches.
unknown = nan;
img(mask) = unknown;

%% Inpainting algorithm

gparam = giin_default_parameters();
[G, pixels, patches] = giin_patch_graph(img, gparam, false);
[G, pixels, Pstructure, Pinformation] = giin_inpaint(G, pixels, patches, gparam, false);

% Global optimization.
sol = zeros(size(pixels));
G = gsp_estimate_lmax(G);
for ii = 1:Nc
    sol(:,ii) = giin_global(G, img(:,:,ii), reshape(pixels(:,ii),Nx,Ny), gparam);
end

%% Results saving

filename = ['results/',imname];
save([filename,'.mat']);

imwrite(reshape(pixels,size(img)), ...
    [filename,'_inpainted.png'], 'png');
imwrite(reshape(sol,size(img)), ...
    [filename,'_inpainted_global.png'], 'png');
imwrite(imadjust(reshape(Pstructure,Nx,Ny))*255, ...
    hot(64), [filename,'_priority_structure.png'], 'png');
imwrite(reshape(Pinformation(:,1),Nx,Ny)*255, ...
    hot(64), [filename,'_priority_information_pixel.png'], 'png');
imwrite(reshape(Pinformation(:,2),Nx,Ny)*255, ...
    hot(64), [filename,'_priority_information_patch.png'], 'png');
imwrite(imadjust(reshape(Pstructure .* Pinformation(:,2),Nx,Ny))*255, ...
    hot(64), [filename,'_priority_global.png'], 'png');

%% Visualization

% load([filename,'.mat']);
% giin_plot_signal(G, pixels(:,1), false);
% giin_plot_priorities(vertices, G, gparam, filename);