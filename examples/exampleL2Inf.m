clear all;close all;clc;

%% read data
addpath(genpath('..'));

dimension = 4;
weightL2 = 0.1;
weightL2Inf = 0.00001;
dataPart = [5*rand(1,dimension);5*rand(1,dimension)]
%func receives n dimensional x
func = @(x) weightL2/2 * norm(x - dataPart, 2)^2 + weightL2Inf * max(sum(x.^2));
minMatlab = fminunc(func, dataPart);
disp('Matlab min: ')
disp(num2str(minMatlab));
disp(['with value of: ', num2str(func(minMatlab))]);
%%
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

numberT = main.addPrimalVar(size(dataPart));
main.addTerm(L2dataTerm(weightL2,dataPart),numberT);
main.addTerm(L2Infidentity(weightL2Inf,size(dataPart)),numberT);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
minXFlex = main.getPrimal(numberT);

%% print result
disp('FlexBox min: ')
disp(num2str(minXFlex));
disp(['with value of: ', num2str(func(minXFlex))]);
