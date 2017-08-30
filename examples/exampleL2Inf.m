clear all;close all;clc;

%% read data
addpath(genpath('..'));

dimension = 500000;
weightL2 = 3;
weightL2Inf = 0.1;
dataPart = [5*rand(1,dimension);5*rand(1,dimension)];
%dataPart = [1:10; 10:-1:1];
%% FlexBox
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

numberT = main.addPrimalVar(size(dataPart));


Dx = gradientOperator(size(dataPart), 1);
Dy = gradientOperator(size(dataPart), 2);

func = @(x) weightL2/2 * norm(x - dataPart, 2)^2 + weightL2Inf * sqrt(max(Dx*x(:).^2 + Dy*x(:).^2));

main.addTerm(L2dataTerm(weightL2,dataPart),numberT);
main.addTerm(L2Infgradient(weightL2Inf,size(dataPart)),numberT);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
minXFlex = main.getPrimal(numberT);

%% print result
%disp('FlexBox min: ')
%disp(num2str(minXFlex));
disp(['with value of: ', num2str(func(minXFlex))]);
%tic;resMatlab = fminsearch(func, dataPart);toc;
%disp(['rel erro: ', num2str(norm(resMatlab - minXFlex))]);
