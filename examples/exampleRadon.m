clear all;close all;clc;

%% read data
addpath(genpath('..'));

sizeImage = 128;
angles = 1:3:180;

%generate Shepp-Logan phantom
image = phantom('Modified Shepp-Logan', sizeImage);

%read MATLAB MRI example
%image = im2double(imread('mri.tif'));
%%
%show clean input image
figure(1);imagesc(image);axis image;colormap(gray);title('Input image');

%generate matrix representative of radon transform using FlexBox
[radonMatrix,sizeSinogram] = generateRadonMatrix ( size(image),angles );

imageRadonTransform = reshape(radonMatrix*image(:),sizeSinogram);

figure(2);imagesc(imageRadonTransform);axis image;colormap(gray);title('Radon Transform of Input Image');

imageRadonTransformNoise = imageRadonTransform + randn(size(imageRadonTransform))*0.05;

figure(3);imagesc(imageRadonTransformNoise);axis image;colormap(gray);title('Radon Transform of Input Image with Noise');

filteredBackprojection = iradon(imresize(imageRadonTransformNoise,sizeSinogram),angles);
figure(4);imagesc(filteredBackprojection);axis image;colormap(gray);title('Iradon Reconstruction');drawnow

%% Reconstruction and denoising using ROF model
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal var u
numberU = main.addPrimalVar(size(image));

%add data-fidelity: 1/2\|Au-f\|_2^2
main.addTerm(L2dataTermOperator(1,radonMatrix,imageRadonTransformNoise),numberU);

%add regularizer
%main.addTerm(L2identity(1,size(image)),numberU); %Tikhonov-regularization
%main.addTerm(L2gradient(2,size(image)),numberU); %Smoothness
main.addTerm(L1gradientIso(0.08,size(image)),numberU); %TV-regularization


% %BEGIN TGV REGULARIZATION
% numberW1 = main.addPrimalVar(size(image));
% numberW2 = main.addPrimalVar(size(image));
% main.addTerm(L1secondOrderGradientIso(0.05,size(image)),[numberU,numberW1,numberW2]);
% 
% %create symmetrized gradients
% gradientBackX = generateBackwardGradND( size(image),[1,1], 1);
% gradientBackY = generateBackwardGradND( size(image),[1,1], 2);
% gradientSymmetrized = 0.5 * (gradientBackX + gradientBackY);
% 
% main.addTerm(frobeniusOperator(10,2,{gradientBackX,gradientSymmetrized,gradientSymmetrized,gradientBackY}),[numberW1,numberW2]);
% %END TGV REGULARIZATION

%box constraint ensures the result to stay in [0,1]
main.addTerm(boxConstraint(0,1,size(image)),numberU);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%get result
result = main.getPrimal(numberU);

%show result
figure(5);imagesc(result);axis image;colormap(gray);title('Reconstructed Image');
