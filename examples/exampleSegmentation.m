clear all;close all;clc;

%% read data
image = imread('data/a.png');
if (size(image,3) > 1)
    image = rgb2gray(image);
end

image = im2double(image);

figure(1);imagesc(image);axis image;colormap(gray);colorbar;title('Image to cluster')
%%

numberOfLabels = 3;
dims = size(image);
labels = [0.1,0.5,0.75];


main = flexBox;

for i=1:numberOfLabels
    main.addPrimalVar(size(image));
end

%init data term
main.addTerm(labelingTerm(1,image,labels),1:numberOfLabels);

for i=1:numberOfLabels
    main.addTerm(L1gradientIso(0.01,size(image)),i);
end

main.runAlgorithm;

for i=1:numberOfLabels
    labelMatrix(:,:,i) = main.getPrimal(i);
end

%
for i=1:numberOfLabels
    figure(1+i);imagesc(main.getPrimal(i));axis image;colormap(gray);colorbar;title(['Region ',num2str(i)])
    completeLabels(:,:,i) = main.getPrimal(i);
end

[~,indexMap] = max(completeLabels,[],3);

figure(2+numberOfLabels);imagesc(indexMap);axis image;colormap(gray);colorbar;title('Segmentation')