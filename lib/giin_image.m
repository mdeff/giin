function [ vertices ] = giin_image( imname )
%GIIN_IMAGE Generate an image. Result saved in the data sub-folder.
%   Usage :
%       vertices = giin_image('vertical'); 
%       inpaint('vertical'); 
%
%   Input parameters :
%       imname   : name of the image file
%
%   Output parameters :
%       vertices : some vertices of interest for this image
%

% Author: Michael Defferrard
% Date: February 2015

% Default parameters.
imsize = 50;
img = zeros(imsize);
contrast = 0.8;
vertices = [];

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
    case 'lena3c'
        imsize = 100;
        img = imread('lena.png');
        img = imcrop(img, [70,100,imsize-1,imsize-1]);
        vertices = [-40,40 ; 5,-20 ; -6,-3 ; 10,1 ; 16,7];
    case 'lena4'
        imsize = 200;
        img = rgb2gray(imread('lena.png'));
        img = imcrop(img, [200,1,imsize-1,imsize-1]);
	case 'lena4c'
        imsize = 200;
        img = imread('lena.png');
        img = imcrop(img, [200,1,imsize-1,imsize-1]);
    case 'lenafull'
        img = imread('lena.png');
        
    otherwise
        error('Unknown image.');
end

% The priority computation of the vertices in this list will be displayed.
if numel(vertices) > 0
    vertices = vertices + floor(imsize/2);
    vertices = (vertices(:,1)-1)*imsize + vertices(:,2);
    vertices = vertices.';
    if any(vertices < 0)
        error('Image size too small to show the vertices of interest.');
    end
end

% Normalize the image in [0,1].
if any(img(:) > 1)
	img = double(img) / 255;
end

% Create a square hole.
holesize = round(imsize / 4);
xyrange = floor((imsize-holesize)/2)+1 : imsize-ceil((imsize-holesize)/2);
obsimg = double(img);
if size(obsimg,3) == 1
    obsimg(xyrange,xyrange) = 1;
else
    obsimg(xyrange,xyrange,1) = 0;
    obsimg(xyrange,xyrange,2) = 1;
    obsimg(xyrange,xyrange,3) = 0;
end

imwrite(img, ['data/',imname,'_original.png']);
imwrite(obsimg, ['data/',imname,'.png']);

end