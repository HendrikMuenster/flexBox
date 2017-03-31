clear all;close all;clc;

%% read data
addpath(genpath('..'));
image = imread('data/hendrik.jpg');
if (size(image,3) > 1)
    image = rgb2gray(image);
end

image = im2double(image);
image = imresize(image,0.25);

figure(1);imagesc(image);axis image;colormap(gray);colorbar;title('Image to cluster')
%%
numberOfLabels = 5;
dims = size(image);
labels = [0.1,0.3,0.5,0.7,0.9];

main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

for i=1:numberOfLabels
    %add one primal var for each label
    main.addPrimalVar(size(image));
end

%init data term
main.addTerm(labelingTerm(1,image,labels),1:numberOfLabels);

for i=1:numberOfLabels
    %add regularizer for each label
    main.addTerm(L1gradientIso(0.05,size(image)),i);
end

%run minimization algorithm
tic;main.runAlgorithm;toc;

for i=1:numberOfLabels
    labelMatrix(:,:,i) = main.getPrimal(i);
end

%%
result = 0;
for i=1:numberOfLabels
    result = result + main.getPrimal(i) * labels(i);
end

figure(2);imagesc(result);axis image;colormap(gray);colorbar;title('Result')