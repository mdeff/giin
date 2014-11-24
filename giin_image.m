function [ img, vertices ] = giin_image( imtype, imsize, gparam )
%GIIN_IMAGE Images to inpaint.
%   Usage :
%       img = giin_image('horizontal', 5); 
%       [img, vertices] = giin_image('horizontal', 5, 50); 
%
%   Input parameters :
%       imtype : image to generate
%       psize  : size of a patch (psize x psize)
%       imsize : size of the image (imsize x imsize)
%
%   Output parameters :
%       img      : matrix representing the image
%       vertices : some vertices of interest for this image
%

% Author: Michaël Defferrard
% Date: November 2014

if nargin < 3
    imsize = 30;
end

img = zeros(imsize);
contrast = 0.8;

if strcmp(imtype, 'horizontal')
    img(imsize/2+1:end,:) = contrast;
    vertices = [0,0];
elseif strcmp(imtype, 'vertical')
    img(:,imsize/2+1:end) = contrast;
    vertices = [-6,-3 ; -2,-6 ; 1,2 ; 6,-5];
elseif strcmp(imtype, 'diagonal')
    img = repmat(1:imsize,imsize,1);
    img = img + repmat((1:imsize).',1,imsize);
    img = (img > imsize) * contrast;
    vertices = [-8,-8 ; -2,-6 ; 1,2 ; 4,-5];
elseif strcmp(imtype, 'cross')
    img(imsize/2+1:end,1:imsize/2) = contrast;
    img(1:imsize/2,imsize/2+1:end) = contrast;
    vertices = [-9,-9 ; -2,-2 ; 1,2 ; 6,-3];
elseif strcmp(imtype, 'lena1')
    imsize = 30;
    img = imread('lena.png');
    img = imcrop(img, [120,100,imsize-1,imsize-1]);
    img = double(img) / 255;
    vertices = [0,0];
elseif strcmp(imtype, 'lena2')
    imsize = 100;
    img = imread('../lena.png');
    img = imcrop(img, [100,100,imsize-1,imsize-1]);
    img = double(img) / 255;
    vertices = [-40,40 ; -26,-20 ; -6,-3 ; 10,1 ; 16,7];
elseif strcmp(imtype, 'lena3')
    imsize = 100;
    img = imread('../lena.png');
    img = imcrop(img, [70,100,imsize-1,imsize-1]);
    img = double(img) / 255;
    vertices = [0,0];
elseif strcmp(imtype, 'lenafull')
    img = imread('../lena.png');
    img = double(img) / 255;
    vertices = [0,0];
else
    error('Unknown image type.');
end

height = imsize - 2*floor(gparam.psize/2);
vertices = vertices + floor(imsize/2);
vertices = (vertices(:,1)-1)*height + vertices(:,2);
vertices = vertices.';

if any(vertices(:)<0)
    error('Image size too small to show the vertices of interest.');
end

% Corner markers. Can break KNN 
% if strcmp(imtype, 'horizontal') || strcmp(imtype, 'vertical')
%     margin = floor(psize / 2);
%     img(1+margin,1+margin) = 0.2;
%     img(imsize-margin,imsize-margin) = 1.0;
% end

end