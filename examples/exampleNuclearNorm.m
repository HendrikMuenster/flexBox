clear all;close all;clc;

%% read data
addpath(genpath('..'));

image = im2double((imread('data/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end
imageNoisy = imnoise(image, 'salt & pepper', 0.01);

%show clean input and noisy image
figure(1);imagesc(image);axis image;colormap(gray);title('Input Image')
figure(2);imagesc(imageNoisy);axis image;colormap(gray);title('Noisy Image')
%% FlexBox
main = flexBox;
main.params.tryCPP = 0;

numberL = main.addPrimalVar(size(image));
numberS = main.addPrimalVar(size(image));
weightL2 = 3;
weightNuclear = 0.3;
weightL1 = 0.01;

nPx = numel(image);

%data term: 1/2 * || id*(L+S) - f||
main.addTerm(L2dataTermOperator(weightL2,{identityOperator(nPx), identityOperator(nPx)},imageNoisy),[numberL, numberS]);

%nuclear norm regularizer
main.addTerm(nuclearIdentity(weightNuclear,size(image)), numberL);

%1-norm regularizer (doesnt matter if L1
main.addTerm(L1identity(weightL1,size(image)), numberS);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
minLFlex = main.getPrimal(numberL);
minSFlex = main.getPrimal(numberS);


%% print results
figure(3);imagesc(minLFlex);axis image;colormap(gray);title('L');
figure(4);imagesc(minSFlex);axis image;colormap(gray);title('S');

%% print diffs
%figure(5);imagesc(minLFlex + minSFlex);axis image;colormap(gray);title('L+S');
%figure(6);imagesc(minLFlex + minSFlex - imageNoisy);axis image;colormap(gray);title('L+S - data');
%figure(7);imagesc(minLFlex - image);axis image;colormap(gray);title('L - original');

