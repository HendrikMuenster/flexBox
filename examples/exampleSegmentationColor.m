clear all;close all;clc;

%% read data
addpath(genpath('..'));
image = imread('data/hendrik.jpg');


image = im2double(image);
image = imresize(image,0.25);

figure(1);imagesc(image);axis image;title('Image to cluster')
%%
numberOfLabels = 10;
dims = size(image);

labels = rand(numberOfLabels,3); %each label is a 1x3 vector for RGB images


for k = 1:10
    main = flexBox;
    main.params.tryCPP = 0; %change, if C++ module is compiled

    for i=1:numberOfLabels
        %add one primal var for each label
        main.addPrimalVar(dims(1:2));
    end

    %init data term
    main.addTerm(labelingTerm(1,image,labels,'lastDim','feature'),1:numberOfLabels); %tell flexBox that the last dimension is a feature

    for i=1:numberOfLabels
        %add regularizer for each label
        main.addTerm(L1gradientIso(0.05,dims(1:2)),i);
    end

    %run minimization algorithm
    tic;main.runAlgorithm;toc;

    %% extract result and update labels
    result = zeros(size(image));
    for i=1:numberOfLabels
        for j=1:3
            result(:,:,j) = result(:,:,j) + main.getPrimal(i) * labels(i,j);

            imageTmp = image(:,:,j);

            set = imageTmp(main.getPrimal(i) > 0.5);
            if (numel(set)>0)
                labels(i,j) = mean(set(:));
            end
        end
    end

    figure(20);imagesc(result);axis image;title(['Colored segmentation after ',num2str(k),'th iteration.'])
end
