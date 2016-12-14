clear all;close all;clc;
%% read data
clc;
addpath(genpath(cd));

image = im2double((imread('examples/data/test.png')));
if (size(image,3) > 1)
    image = rgb2gray(image);
end


image = imresize(image,1);

imageNoisy = image + randn(size(image)) * 0.05;

%imageNoisy = image;
%figure(1);imagesc(image,[0,1]);axis image;colormap(gray)
%figure(2);imagesc(imageNoisy,[0,1]);axis image;colormap(gray)


% denoising
main = flexBox;

%add primal var u
numberU = main.addPrimalVar(size(image));
numberU2 = main.addPrimalVar(size(image));

oneM = rand(numel(image),1);
oneM1 = rand(numel(image),1);
oneM2 = rand(numel(image),1);
opEye = spdiags([oneM,-oneM1,oneM2],[0,1,2],numel(image),numel(image));
opEye = spdiags([oneM],[1],numel(image),numel(image));

%opEye = speye(numel(image));
%opEye(end-30000:end,:) = [];


main.addTerm(L1dataTermOperator(1,diagonalOperator(ones(numel(image),1)),imageNoisy),numberU);
main.addTerm(L1dataTermOperator(1,opEye,imageNoisy),numberU);
main.addTerm(L1dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);
main.addTerm(L2dataTerm(1,imageNoisy),numberU);
main.addTerm(L2dataTermOperator(1,opEye,imageNoisy),numberU);
main.addTerm(L1gradientAniso(1,size(image)),numberU);
main.addTerm(L1operatorAniso(1,1,opEye),numberU);

main.addTerm(L2dataTerm(1,imageNoisy),numberU2);
main.addTerm(L1gradientAniso(1,size(image)),numberU2);

main.addTerm(L1dataTerm(1,imageNoisy),numberU);
main.addTerm(KLdataTermOperator(1,15*speye(numel(image)),imageNoisy),numberU);
main.addTerm(KLdataTermOperator(1,identityOperator(numel(image)),imageNoisy),numberU);
main.addTerm(L1dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);
main.addTerm(L2dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);

%grad
main.addTerm(L1gradientIso(0.1,size(image)),numberU);
main.addTerm(L2gradient(10,size(image)),numberU);
main.addTerm(huberGradient(0.8,size(image),0.01),numberU);
main.addTerm(frobeniusGradient(10,size(image)),numberU);

%generalOperator
main.addTerm(L1operatorAniso(0.5,1,speye(numel(image))),numberU);
main.addTerm(L1operatorIso(0.5,1,speye(numel(image))),numberU);
main.addTerm(L2operator(0.5,1,speye(numel(image))),numberU);
main.addTerm(frobeniusOperator(10,1,speye(numel(image))),numberU);


%vector field
main.addTerm(L1divergence(0.5,size(image),'usedims',[0,1]),numberU);
main.addTerm(L2divergence(0.5,size(image),'usedims',[0,1]),numberU);
main.addTerm(L1curl(0.5,size(image)),[numberU,numberU]);
main.addTerm(L2curl(0.5,size(image)),[numberU,numberU]);

%other
main.addTerm(L1identity(0.5,size(image)),numberU);
main.addTerm(L2identity(0.5,size(image)),numberU);

%main.params.tryCPP = 0;
%tic;main.runAlgorithm;toc;

%main.params.tryCPP = 0;
%tic;main.runAlgorithm;toc;
%%
%main.params.tryCPP = 0;
%main.params.tol = 1e-7;
main.params.verbose = 2;
%for i=1:100
tic;main.runAlgorithm;toc;
%end
figure(3);imagesc(main.getPrimal(1));colormap(gray);colorbar
%%

main2 = flexBox;

numberU = main2.addPrimalVar(size(image));

main2.addTerm(L2dataTerm(1,imageNoisy),numberU);
%main2.addTerm(L1operatorAniso(1,1,opEye),numberU);
main2.addTerm(L1gradientAniso(1,size(image)),numberU);

%add primal var u

%main2.addTerm(L1dataTerm(1,imageNoisy),numberU);
%main2.addTerm(L2dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);
%main2.addTerm(KLdataTermOperator(1,15*speye(numel(image)),imageNoisy),numberU);
%main2.addTerm(L2dataTermOperator(1,speye(numel(image)),imageNoisy),numberU);
%main2.addTerm(L1gradientAniso(1,size(image)),numberU);
%main2.addTerm(L2gradient(10,size(image)),numberU);
%main2.addTerm(frobeniusGradient(10,size(image)),numberU);
%main2.params.maxIt = 1000;
tic;main2.runAlgorithm;toc;

%%


%get result
result = main.getPrimal(numberU);
result2 = main2.getPrimal(numberU);

figure(4);imagesc(abs(result-result2),[0,1])

figure(6);imagesc(result,[0,1]);axis image;colormap(gray);
figure(7);imagesc(result2,[0,1]);axis image;colormap(gray);
