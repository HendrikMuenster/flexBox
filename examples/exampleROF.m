%%clear all;close all;clc;

%% read data
addpath(genpath('..'));

image = im2double((imread('data/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end

imageNoisy = image + randn(size(image)) * 0.05;

figure(1);imagesc(image);axis image;colormap(gray);title('Input Image')
figure(2);imagesc(imageNoisy);axis image;colormap(gray);title('Noisy Image')
%% ROF denoising
main = flexBox;

%add primal var u
numberU = main.addPrimalVar(size(image));

%add data-fidelity: 1/2\|u-f\|_2^2
main.addTerm(L2dataTerm(1,imageNoisy),numberU);

%add regularizer: 0.08*\|\nabla u\|_1
main.addTerm(L1gradientIso(0.08,size(image)),numberU);

main.runAlgorithm;

%get result
result = main.getPrimal(numberU);

figure(3);imagesc(result);axis image;colormap(gray);;title('Denoised Image');