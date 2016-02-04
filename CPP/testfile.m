% 
% nPx = 50000;
% 
% A = speye(nPx);
% b = rand(nPx,1);
% tic
% flexboxCPP(A,b);
% toc;
% %clear all;
% %%

clear all;close all;clc;
%% denoising
clc
addpath(genpath('../'));

image = im2double((imread('test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end

imageNoisy = image + randn(size(image)) * 0.05;


%figure(1);imagesc(image);axis image;colormap(gray)
%figure(2);imagesc(imageNoisy);axis image;colormap(gray)

%%L2L2 denoising
tic;
resultL2L2 = denoising(imageNoisy,0.1,'L2','IsoTV');
toc

figure(5);imagesc(resultL2L2);axis image;colormap(gray)