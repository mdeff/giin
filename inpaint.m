function [ sol ] = inpaint( imname, fillcolor )
%INPAINT Retrieve the missing pixels of an image.
%   Usage :
%       inpaint('bungee', [0,255,0]); 
%       [img, vertices] = giin_image('horizontal', 5, 50); 
%
%   Input parameters :
%       imname    : name of the image file
%       fillcolor : color of the fill region
%
%   Output parameters :
%       sol       : the inpainted image
%
% The goal of this script is to inpaint the missing pixels of an image. It
% does so by constructing a patch graph of the known pixels. According to
% some priority, it then iteratively connect unknown patches to the graph.
% A global optimization is run in the end.

addpath('./lib');
addpath('./data');
gsp_start();

%% Inpainting algorithm

gparam = giin_default_parameters();
[img, obsimg, imsize, vertices] = giin_image(imname, true);
[G, pixels, patches] = giin_patch_graph(obsimg, gparam, false);


%%
Nc = size(pixels,2);

[G, pixels, Pstructure, Pinformation] = giin_inpaint(G, pixels, patches, gparam, false);
%%
sol = zeros(size(pixels));
G = gsp_estimate_lmax(G);

for ii = 1:Nc
    sol(:,ii) = giin_global(G, obsimg(:,:,ii),reshape(pixels(:,ii),size(img,1),size(img,2)), gparam);
end

%% Save the results

filename = ['results/',imname];
save([filename,'.mat']);
imwrite(reshape(pixels,size(img)), [filename,'_inpainted.png'], 'png');
imwrite(reshape(sol,size(img)), [filename,'_inpainted_global.png'], 'png');