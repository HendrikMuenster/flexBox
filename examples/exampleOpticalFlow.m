clear all;close all;clc;

%% read data
addpath(genpath('..'));
f1 = imread('data/a.png');
f2 = imread('data/b.png');

if (size(f1,3) > 1)
    f1 = rgb2gray(f1);
end
if (size(f2,3) > 1)
    f2 = rgb2gray(f2);
end

f1 = im2double(f1);
f2 = im2double(f2);

%f1 = imresize(f1,2);
%f2 = imresize(f2,2);

figure(1);imagesc(f1);axis image;colormap(gray);title('Image 1')
figure(2);imagesc(f2);axis image;colormap(gray);title('Image 2')
%%
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled
main.params.verbose = 1; %change, if C++ module is compiled

%add primal vars v_1,v_2
numberV1 = main.addPrimalVar(size(f1));
numberV2 = main.addPrimalVar(size(f1));

%add optical flow data term
main.addTerm(L1opticalFlowTerm(1,f1,f2),[numberV1,numberV2]);

%uncomment for gradient constancy
%main.addTerm(L1opticalFlowTerm(1,f1,f2,'termType','gradientConstancy','constancyDimension',1),[numberV1,numberV2]);
%main.addTerm(L1opticalFlowTerm(1,f1,f2,'termType','gradientConstancy','constancyDimension',2),[numberV1,numberV2]);

%add regularizers - one for each component
main.addTerm(huberGradient(0.1,size(f1),0.01),numberV1);
main.addTerm(huberGradient(0.1,size(f1),0.01),numberV2);

%run minimization algorithm
tic;main.runAlgorithm;toc;

% get result
% switch components because MATLAB switches axes
resultV2 = main.getPrimal(numberV1);
resultV1 = main.getPrimal(numberV2);

figure(3);clf;imagesc(flowToColorV2(cat(3,resultV1,resultV2),5));title('Color-coded flow field')