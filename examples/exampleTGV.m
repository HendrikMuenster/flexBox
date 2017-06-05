clear all;close all;clc;

%% read data
addpath(genpath('..'));

image = im2double((imread('data/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end

imageNoisy = image + randn(size(image)) * 0.05;

figure(1);imagesc(image);axis image;colormap(gray);title('Input Image')
figure(2);imagesc(imageNoisy);axis image;colormap(gray);title('Noisy Image')
%% denoising with tgv
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled
main.params.tol = 1e-5;

%add primal vars u,w
numberU = main.addPrimalVar(size(image));
numberW1 = main.addPrimalVar(size(image));
numberW2 = main.addPrimalVar(size(image));

%add data fit 1/2\|u-f\|_2^2
main.addTerm(L2dataTerm(1,imageNoisy),numberU);

%add regularizer 0.05*\|\nabla (u)-w\|_1
main.addTerm(L1secondOrderGradientIso(0.05,size(image)),[numberU,numberW1,numberW2]);

%create symmetrized gradients
gradientBackX = generateBackwardGradND( size(image),[1,1], 1);
gradientBackY = generateBackwardGradND( size(image),[1,1], 2);

gradientSymmetrized = 0.5 * (gradientBackX + gradientBackY);

main.addTerm(frobeniusOperator(10,2,{gradientBackX,gradientSymmetrized,gradientSymmetrized,gradientBackY}),[numberW1,numberW2]);

tic;main.runAlgorithm;toc;
%% get result
result = main.getPrimal(numberU);

figure(3);imagesc(result);axis image;colormap(gray);title('Denoised Image');

figure(4);plot(1:size(image),result(150,:));title('1D Result - Linear behavior');
