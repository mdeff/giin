function [ Pstructure ] = priority( imname )
%INPAINT Assign structure priority to every pixel of an image.
%   Usage :
%       giin_image('vertical'); 
%       priority('vertical_original');
%
%   Input parameters :
%       imname     : name of the image file
%
%   Output parameters :
%       Pstructure : structure priority
%

% Author: Michael Defferrard
% Date: February 2015

addpath('./lib');
addpath('./data');
gsp_start();

%% Image loading

img = double(imread([imname,'.png'])) / 255;
Nx = size(img,1);
Ny = size(img,2);

%% Patch graph

gparam = giin_default_parameters();
G = giin_patch_graph(img, gparam, false);

%% Priorities

tstart = tic;
Pstructure = nan(G.N, 1);
N = ceil(G.N / 1000); % Chunks of 1000 to save runtime memory.
for n = 0:N-1;
    range = n*1000+1 : min((n+1)*1000, G.N);
    Pstructure = giin_priorities(range, Pstructure, G, gparam);
end
fprintf('Priorities : %f seconds\n', toc(tstart));

if any(isnan(Pstructure))
    error('Some pixels have no assigned priority !');
end

%% Results saving

filename = ['results/',imname,'_priority'];
save([filename,'.mat']);
imwrite(imadjust(reshape(Pstructure,Nx,Ny)), [filename,'.png']);

%% Visualization

% load([filename,'.mat']);
% imshow(imadjust(reshape(Pstructure,Nx,Ny)));
% giin_plot_priorities(vertices, G, gparam, filename);