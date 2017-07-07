clear all;close all;clc;

%% read data
addpath(genpath('..'));

weightL2 = 3;
weightLInf = 25.5;
dataPart = [10, 12];
[X,Y] = meshgrid(-20:0.5:20,-20:0.5:20);
%fx = @(X,Y) weightL2/2*((X - dataPart(1)).^2 + (Y - dataPart(2)).^2) + weightLInf*(max(abs(X),abs(Y)));
fx = @(X,Y) weightLInf*max(abs(X - dataPart(1)),abs(Y - dataPart(2))) + weightL2/2*(X.^2 + Y.^2);
Z = fx(X,Y);
surf(X,Y,Z);
minX = fminunc(@(x) fx(x(1),x(2)), [0,0]);
disp(["MATLAB: ", minX]);
%%
main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal var u
numberT = main.addPrimalVar(size(dataPart));

main.addTerm(LInfdataTerm(weightLInf,dataPart),numberT);

main.addTerm(L2identity(weightL2,size(dataPart)),numberT);

%run minimization algorithm
tic;main.runAlgorithm;toc;

%% get result
result = main.getPrimal(numberT);
disp(result);
fx(minX(1),minX(2))
fx(result(1),result(2))
%% show result
%figure(2);
%plot(result);