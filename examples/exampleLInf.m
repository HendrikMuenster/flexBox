clear all;close all;clc;

%% read data
addpath(genpath('..'));

image = im2double((imread('data/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end
%image = imresize(image,3);
%imageNoisy = image + randn(size(image)) * 0.05;
imageNoisy = image;

%show clean input and noisy image
figure(1);imagesc(image);axis image;colormap(gray);title('Input Image')
figure(2);imagesc(imageNoisy);axis image;colormap(gray);title('Noisy Image')
%%
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal var u
numberU = main.addPrimalVar(size(image));

%main.addTerm(L1InfdataTerm(1,imageNoisy),numberU);
main.addTerm(L2dataTerm(1,imageNoisy),numberU);

%main.addTerm(L1gradientIso(0.1,size(image)),numberU);
main.addTerm(L1Infgradient(0.01,size(image)),numberU);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
result = main.getPrimal(numberU);

%% show result
figure(3);imagesc(result);axis image;colormap(gray);title('Result');
