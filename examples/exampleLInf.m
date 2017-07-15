clear all;close all;clc;

%% read data
addpath(genpath('..'));

weightL2 = 0.1;
weightLInf = 2.5;
dataPart = [10, 12];
[X,Y] = meshgrid(-20:0.5:20,-20:0.5:20);
fx = @(X,Y) weightL2/2*((X - dataPart(1)).^2 + (Y - dataPart(2)).^2) + weightLInf*(max(abs(X),abs(Y)));
Z = fx(X,Y);
surf(X,Y,Z);rotate3d on;

[minX,fval] = fminsearch(@(x) fx(x(1),x(2)), [1,1]);
disp('Calculating minimal function value using MATLAB')
disp(['argmin: ', num2str(minX)]);
disp(['min: ', num2str(fval)]);

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
result = main.getPrimal(numberT);

disp('Calculating minimal function value using FlexBox')
disp(['argmin: ', num2str([result(1),result(2)])]);
disp(['min: ', num2str(fx(result(1),result(2)))]);
