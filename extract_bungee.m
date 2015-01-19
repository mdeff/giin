function [img,mask] =  extract_bungee(color)
%% Bungee signal

if nargin<1
    color = 0;
end
    

bungee = imread('bungee02.png');

if color
    bungee = double(bungee)/256;
else
    bungee = double(rgb2gray(bungee))/256;
end
% figure;
% imagesc(bungee);
% if ~color
%     colormap(gray);
% end

mask = imread('bungee12.png');

mask = mask(:,:,1)==0 & mask(:,:,2)==255 & mask(:,:,3)==0;
% figure;
% imagesc(mask)

if color
    mask = cat(3, mask,mask,mask);
end

img = bungee;

end