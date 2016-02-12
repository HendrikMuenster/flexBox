clear all;close all;clc;
%% read data
addpath(genpath(cd));

image = im2double((imread('examples/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end

imageNoisy = image + randn(size(image)) * 0.05;

figure(1);imagesc(image);axis image;colormap(gray)
figure(2);imagesc(imageNoisy);axis image;colormap(gray)
%% denoising
main = flexbox;

%add primal var u
numberU = main.addPrimalVar(size(image));

%add data-fidelity: 1/2\|u-f\|_2^2
main.addTerm(L2dataTerm(1,imageNoisy),numberU);

%add regularizer: 0.08*\|\nabla u\|_1
main.addTerm(L1gradientIso(0.08,size(image)),numberU);

main.params.tol = 1e-5;
main.params.tryCPP = 1;
tic;
main.runAlgorithm;
toc;



%get result
result = main.getPrimal(numberU);

figure(3);imagesc(result);axis image;colormap(gray);
%% denoising with tgv
main = flexbox;

%add primal vars u,w
numberU = main.addPrimalVar(size(image));
numberW1 = main.addPrimalVar(size(image));
numberW2 = main.addPrimalVar(size(image));

%add data fit 1/2\|u-f\|_2^2
main.addTerm(L2dataTerm(1,imageNoisy),numberU);

%add regularizer 0.05*\|\nabla (u)-w\|_1
main.addTerm(L1secondOrderGradientIso(0.05,size(image)),[numberU,numberW1,numberW2]);

%add regularizer 0.05*\|\nabla (w)\|_1
main.addTerm(L1gradientIso(20,size(image),'discretization','backward'),numberW1);
main.addTerm(L1gradientIso(20,size(image),'discretization','backward'),numberW2);

main.params.showPrimals = 500;
main.params.tryCPP = 0;
tic
main.runAlgorithm;
toc
%% get result
result = main.getPrimal(numberU);
result2 = main.getPrimal(numberW1);

figure(3);imagesc(result);axis image;colormap(gray)
figure(4);imagesc(result2);axis image;colormap(gray)

%% optical flow
f1 = imread('examples/a.png');
f2 = imread('examples/b.png');

if (size(f1,3) > 1)
    f1 = rgb2gray(f1);
end
if (size(f2,3) > 1)
    f2 = rgb2gray(f2);
end

f1 = im2double(f1);
f2 = im2double(f2);

main = flexbox;

%add primal vars v_1,v_2
numberV1 = main.addPrimalVar(size(f1));
numberV2 = main.addPrimalVar(size(f2));

%add optical flow data term
main.addTerm(L1opticalFlowTerm(1,f1,f2),[numberV1,numberV2]);

%add regularizers - one for each component
main.addTerm(L1gradientIso(0.05,size(f1)),numberV1);
main.addTerm(L1gradientIso(0.05,size(f1)),numberV2);

main.runAlgorithm;

resultV1 = main.getPrimal(numberV1);
resultV2 = main.getPrimal(numberV2);

figure(1);imagesc(f1);axis image;colormap(gray)
figure(2);imagesc(f2);axis image;colormap(gray)

figure(3);clf;imagesc(flowToColorV2(cat(3,resultV1,resultV2),5));drawnow
%% segmentation
image = imread('examples/a.png');
if (size(image,3) > 1)
    image = rgb2gray(image);
end

image = im2double(image);

numberOfLabels = 3;
dims = size(image);
labels = rand(numberOfLabels,1);


main = flexbox;

for i=1:numberOfLabels
    main.addPrimalVar(size(image));
end

%init data term
main.addTerm(labelingTerm(1,image,labels),1:numberOfLabels);

for i=1:numberOfLabels
    main.addTerm(L1gradientIso(0.5,size(image)),i);
end

main.runAlgorithm;

for i=1:numberOfLabels
    labelMatrix(:,:,i) = main.getPrimal(i);
end

%%

figure(1);imagesc(image);axis image;colormap(gray);colorbar;
for i=1:numberOfLabels
    figure(1+i);imagesc(main.getPrimal(i));axis image;colormap(gray);colorbar;
    completeLabels(:,:,i) = main.getPrimal(i);
end

[~,indexMap] = max(completeLabels,[],3);



figure(2+numberOfLabels);imagesc(indexMap);axis image;colormap(gray);colorbar;