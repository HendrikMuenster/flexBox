clear all;close all;clc;
%% read data
addpath(genpath(cd));

image = im2double((imread('examples/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end

imageNoisy = image + randn(size(image)) * 0.05;

figure(1);imagesc(image,[0,1]);axis image;colormap(gray)
figure(2);imagesc(imageNoisy,[0,1]);axis image;colormap(gray)


%% denoising
main = flexBox;

%add primal var u
numberU = main.addPrimalVar(size(image));


%main.addTerm(L1dataTerm(1,imageNoisy),numberU);
main.addTerm(L2dataTerm(1,imageNoisy),numberU);
%main.addTerm(KLdataTermOperator(1,speye(numel(image)),imageNoisy),numberU);
%main.addTerm(L1dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);
%main.addTerm(L2dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);

%grad
main.addTerm(L1gradientAniso(1,size(image)),numberU);
%main.addTerm(L1gradientIso(0.1,size(image)),numberU);
%main.addTerm(L2gradient(0.1,size(image)),numberU);
%main.addTerm(huberGradient(0.1,size(image),0.001),numberU);
%main.addTerm(frobeniusGradient(10,size(image)),numberU);

%generalOperator
%main.addTerm(L1operatorAniso(0.5,1,speye(numel(image))),numberU);
%main.addTerm(L1operatorIso(0.5,1,speye(numel(image))),numberU);
%main.addTerm(L2operator(0.5,1,speye(numel(image))),numberU);
%main.addTerm(frobeniusOperator(10,1,speye(numel(image))),numberU);


%vector field
%main.addTerm(L1divergence(0.5,size(image),'usedims',[0,1]),numberU);
%main.addTerm(L2divergence(0.5,size(image),'usedims',[0,1]),numberU);
%main.addTerm(L1curl(0.5,size(image)),[numberU,numberU]);
%main.addTerm(L2curl(0.5,size(image)),[numberU,numberU]);

%other
%main.addTerm(L1identity(0.5,size(image)),numberU);
%main.addTerm(L2identity(0.5,size(image)),numberU);

main.params.tryCPP = 1;
tic;
main.runAlgorithm;
toc;

%%


%get result
result = main.getPrimal(numberU);

figure(3);imagesc(result,[0,1]);axis image;colormap(gray);