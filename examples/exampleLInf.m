clear all;close all;clc;

%% read data
addpath(genpath('..'));

dimension = 10;
weightL2 = 2;
weightLInf = 0.5;
dataPart = 5*rand(1,dimension);
%func receives n dimensional x
func = @(x) weightL2/2 * norm(x - dataPart, 2)^2 + weightLInf * norm(x, Inf);
%%
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

numberT = main.addPrimalVar(size(dataPart));
main.addTerm(L2dataTerm(weightL2,dataPart),numberT);
main.addTerm(LInfidentity(weightLInf,size(dataPart)),numberT);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
minXFlex = main.getPrimal(numberT);

%% print result
disp(['FlexBox min: ', num2str(minXFlex)]);
disp(['with value of: ', num2str(func(minXFlex))]);
