function [img,mask] =  extract_bungee()
%% Bungee signal

bungee = imread('bungee0.png');

bungee = double(rgb2gray(bungee))/256;
figure;
imagesc(bungee);
colormap(gray);

mask = imread('bungee1.png');

mask = mask(:,:,1)==0 & mask(:,:,2)==255 & mask(:,:,3)==0;
% figure;
% imagesc(mask)

img = bungee;

end