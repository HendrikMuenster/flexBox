clear all;close all;clc;

%% read data
addpath(genpath('..'));

dimension = 10;
weightL2 = 2;
weightLInf = 0.2;
dataPart = linspace(0,10,dimension);
%func receives n dimensional x
func = @(x) weightL2/2 * norm(x - dataPart, 2)^2 + weightLInf * norm(x, Inf);
minXMatlab = fminunc(func, dataPart);
disp(["Matlab min: ", minXMatlab]);
disp(["with value of: ", func(minXMatlab)]);
%%
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal var u
numberT = main.addPrimalVar(size(dataPart));

main.addTerm(L2dataTerm(weightL2,dataPart),numberT);

main.addTerm(LInfidentity(weightLInf,size(dataPart)),numberT);
%main.addTerm(LInfgradient(weightLInf,size(dataPart)),numberT);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
minXFlex = main.getPrimal(numberT);
disp(["FlexBox min: ", minXFlex]);
disp(["with value of: ", func(minXFlex)]);
%% show result
%figure(2);
%plot(result);
