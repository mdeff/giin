function [ img, obsimg, vertices ] = giin_image( imname, fillcolor )
%GIIN_IMAGE Load or generate an image.
%   Usage :
%       img = giin_image('horizontal'); 
%       [img, obsimg, vertices] = giin_image('bungee', [0,255,0]); 
%
%   Input parameters :
%       imname    : name of the image file
%       fillcolor : color of the fill region
%
%   Output parameters :
%       img      : matrix representing the original image
%       obsimg   : matrix representing the masked image
%       vertices : some vertices of interest for this image
%

% Author: Michael Defferrard
% Date: February 2015

if ~exist('fillcolor', 'var')
    fillcolor = [0,255,0];
end

% Default parameters.
imsize = 30;
img = zeros(imsize);
contrast = 0.8;
vertices = [];

% Unknown pixels are negative (known ones are in [0,1]). Negative enough
% such that they don't connect to anything else than other unknown patches.
unknown = -1e3;

switch(imname)
    
    % Synthetic images.
    case 'horizontal'
        img(imsize/2+1:end,:) = contrast;
        vertices = [];
    case 'vertical'
        img(:,imsize/2+1:end) = contrast;
        vertices = [-6,-3 ; -2,-6 ; 1,2 ; 6,-5];
    case 'diagonal'
        img = repmat(1:imsize,imsize,1);
        img = img + repmat((1:imsize).',1,imsize);
        img = (img > imsize) * contrast;
        vertices = [-8,-8 ; -2,-6 ; 1,2 ; 4,-5];
    case 'cross'
        img(imsize/2+1:end,1:imsize/2) = contrast;
        img(1:imsize/2,imsize/2+1:end) = contrast;
        vertices = [-9,-9 ; -2,-2 ; 1,2 ; 6,-3];
        
    % Variations on Lena.
    case 'lena1'
        imsize = 30;
        img = rgb2gray(imread('lena.png'));
        img = imcrop(img, [120,100,imsize-1,imsize-1]);
    case 'lena2'
        imsize = 100;
        img = rgb2gray(imread('lena.png'));
        img = imcrop(img, [100,100,imsize-1,imsize-1]);
        vertices = [-40,40 ; -26,-20 ; -6,-3 ; 10,1 ; 16,10 ; 20,23];
    case 'lena3'
        imsize = 100;
        img = rgb2gray(imread('lena.png'));
        img = imcrop(img, [70,100,imsize-1,imsize-1]);
        vertices = [-40,40 ; 5,-20 ; -6,-3 ; 10,1 ; 16,7];
    case 'lena3_color'
        imsize = 100;
        img = imread('lena.png');
        img = imcrop(img, [70,100,imsize-1,imsize-1]);
        vertices = [-40,40 ; 5,-20 ; -6,-3 ; 10,1 ; 16,7];
    case 'lena4'
        imsize = 200;
        img = rgb2gray(imread('lena.png'));
        img = imcrop(img, [200,1,imsize-1,imsize-1]);
	case 'lena4_color'
        imsize = 200;
        img = imread('lena.png');
        img = imcrop(img, [200,1,imsize-1,imsize-1]);
    case 'lenafull'
        img = imread('lena.png');
        
    % General case: load an image along with a mask.
    otherwise
        img = imread([imname,'.png']);
        obsimg = double(imread([imname,'_masked.png']));
        mask = obsimg(:,:,1)==fillcolor(1) & ...
            obsimg(:,:,2)==fillcolor(2) & obsimg(:,:,3)==fillcolor(3);
        mask = repmat(mask, [1,1,3]);
        obsimg(mask) = unknown;
end

% The priority computation of the vertices in this list will be displayed.
if numel(vertices) > 0
    vertices = vertices + floor(imsize/2);
    vertices = (vertices(:,1)-1)*imsize + vertices(:,2);
    vertices = vertices.';
    if any(vertices<0)
        error('Image size too small to show the vertices of interest.');
    end
end

% Create a square hole if we didn't load a mask.
if ~exist('obsimg', 'var')
    holesize = round(imsize / 4);
    xyrange = floor((imsize-holesize)/2)+1 : imsize-ceil((imsize-holesize)/2);
    obsimg = double(img);
    obsimg(xyrange,xyrange,:) = unknown;
end

% Normalize the image in [0,1].
if any(img(:) > 1)
	img = double(img) / 255;
end
if any(obsimg(:) > 1)
	obsimg = double(obsimg) / 255;
    obsimg(obsimg<0) = unknown;
end

end