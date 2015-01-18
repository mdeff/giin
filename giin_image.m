function [ img, obsimg, imsize, vertices ] = giin_image( imtype, hole )
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

% Author: Michael Defferrard
% Date: November 2014

if nargin < 2
    hole = false;
end

imsize = 30;
img = zeros(imsize);
contrast = 0.8;

switch(imtype)
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
    case 'lena1'
        imsize = 30;
        img = lena*256;
        img = imcrop(img, [120,100,imsize-1,imsize-1]);
        img = double(img) / 255;
        vertices = [];
    case 'lena2'
        imsize = 100;
        img = lena*256;
        img = imcrop(img, [100,100,imsize-1,imsize-1]);
        img = double(img) / 255;
        vertices = [-40,40 ; -26,-20 ; -6,-3 ; 10,1 ; 16,10 ; 20,23];
    case 'lena3'
        imsize = 100;
        img = lena*256;
        img = imcrop(img, [70,100,imsize-1,imsize-1]);
        img = double(img) / 255;
        vertices = [-40,40 ; 5,-20 ; -6,-3 ; 10,1 ; 16,7];
    case 'lena3_color'
        imsize = 100;
        img = lena(1);
        img = imcrop(img, [70,100,imsize-1,imsize-1]);
        vertices = [-40,40 ; 5,-20 ; -6,-3 ; 10,1 ; 16,7];
    case 'lena4'
        imsize = 200;
        img = lena*256;
        img = imcrop(img, [200,1,imsize-1,imsize-1]);
        img = double(img) / 255;
        vertices = [];
    case 'lenafull'
        img = lena;
        vertices = [];
    case 'bungee'
        [img, mbungee] = extract_bungee();
        % fix pixels in the border
        mbungee(:,65:76)=0;
        %img(end-4:end,65:76) = mean(img(end,[64,77]));
        vertices = [];
    case 'bungee_color'
        [img, mbungee] = extract_bungee(1);
        % fix pixels in the border
        mbungee(:,65:76)=0;
        %img(end-4:end,65:76,:) = repmat(mean(img(end,[64,77],:),2),5,76-65+1,1);
        vertices = [];
	otherwise
        error('Unknown image type.');
end

if numel(vertices) > 0
    vertices = vertices + floor(imsize/2);
    vertices = (vertices(:,1)-1)*imsize + vertices(:,2);
    vertices = vertices.';
end

if any(vertices<0)
    error('Image size too small to show the vertices of interest.');
end

% Corner markers (could break KNN, it was a bug).
if strcmp(imtype, 'horizontal') || strcmp(imtype, 'vertical')
    margin = floor(psize / 2);
    img(1+margin,1+margin) = 0.2;
    img(imsize-margin,imsize-margin) = 1.0;
end

% Unknown pixels are negative (known ones are in [0,1]). Negative enough
% such that they don't connect to anything else than other unknown patches.
if hole
    holesize = round(imsize / 4);
else
    holesize = 0;
end

if strcmp(imtype,'bungee')
    obsimg = img;
    obsimg(mbungee) = -1e3;
elseif strcmp(imtype,'bungee_color')
    obsimg = img;
    obsimg(mbungee) = -1e3;
else   
    xyrange = floor((imsize-holesize)/2)+1 : imsize-ceil((imsize-holesize)/2);
    obsimg = img;
    obsimg(xyrange,xyrange,:) = -1e3;
end

end