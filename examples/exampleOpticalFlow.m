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

figure(1);imagesc(f1);axis image;colormap(gray);title('Image 1')
figure(2);imagesc(f2);axis image;colormap(gray);title('Image 2')
%%
main = flexBox;

%add primal vars v_1,v_2
numberV1 = main.addPrimalVar(size(f1));
numberV2 = main.addPrimalVar(size(f2));

%add optical flow data term
main.addTerm(L1opticalFlowTerm(1,f1,f2),[numberV1,numberV2]);

%add regularizers - one for each component
main.addTerm(L1gradientIso(0.05,size(f1)),numberV1);
main.addTerm(L1gradientIso(0.05,size(f1)),numberV2);

main.runAlgorithm;

% get result
resultV1 = main.getPrimal(numberV1);
resultV2 = main.getPrimal(numberV2);

figure(3);clf;imagesc(flowToColorV2(cat(3,resultV1,resultV2),5));title('Color-coded flow field')