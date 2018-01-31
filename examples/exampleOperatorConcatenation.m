clear all;close all;clc;

%% read data
addpath(genpath('..'));

image = im2double((imread('data/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end

imageNoisy = image + randn(size(image)) * 0.05;

figure(1);imagesc(image);axis image;colormap(gray);title('Input Image')

%% create blurry + noisy image

kernelSize = 11;
sigma = 2;

%create convolution operator with gaussian kernel and identity operator
blurOperator = blurOperator(size(image),kernelSize,sigma);
identityOp = identityOperator(numel(image));

%concat both operators using the concatOperator class
concatOp = concatOperator(blurOperator,identityOp,'composition');


imageBlurred = reshape(concatOp * imageNoisy(:),size(image));

figure(2);imagesc(real(imageBlurred));axis image;colormap(gray);title('Blurred and Noisy Image')
%% ROF deblurring
main = flexBox;

%add primal var u
numberU = main.addPrimalVar(size(image));

%add data-fidelity: 1/2\|u-f\|_2^2
main.addTerm(L2dataTermOperator(1,concatOp,imageBlurred),numberU);

%add regularizer: 0.1*\|\nabla u\|_1
main.addTerm(L1gradientIso(0.1,size(image)),numberU);


main.runAlgorithm();

%get result
result = main.getPrimal(numberU);

figure(3);imagesc(real(result));axis image;colormap(gray);title('Deblurred and Denoised Image');