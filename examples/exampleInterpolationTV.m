clear all;close all;clc;

%% read data
addpath(genpath('..'));
f1 = imread('data/hendrik.jpg');

if (size(f1,3) > 1)
    f1 = rgb2gray(f1);
end

f1 = im2double(f1);

f1 = imresize(f1,0.1); %shrink hendrik

figure(1);imagesc(f1);axis image;colormap(gray);title('Image 1')
%%

%parameters
upsamplingFactor = 4;

main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal var u
numberU = main.addPrimalVar(size(f1)*upsamplingFactor);

%init downsampling operator K
downsamplingOp = superpixelOperator(size(f1),upsamplingFactor);

%add data term Ku-f where K downsamples u
main.addTerm(L1dataTermOperator(1,downsamplingOp,f1),numberU);

%add regularizers - one for each component
main.addTerm(L1gradientIso(0.01,upsamplingFactor*size(f1)),numberU);

%run minimization algorithm
tic;main.runAlgorithm;toc;

% get result
result = main.getPrimal(numberU);

figure(3);clf;imagesc(result);colormap(gray);axis image;title('Result')