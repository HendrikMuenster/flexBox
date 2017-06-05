clear all;close all;clc;

%% read data
addpath(genpath('..'));
f1 = imread('data/a.png');

if (size(f1,3) > 1)
    f1 = rgb2gray(f1);
end

f1 = im2double(f1);

figure(1);imagesc(f1,[0,1]);axis image;colormap(gray);title('Image');colorbar;
%% definition of required parts for function handle operator
%define a function handle, here the -indentity
myHandle = @(x) -x;
%transposed is also identity
myHandleT = @(x) -x;
%size of the input argument is size of the image
argumentSize = size(f1);
%operator norm is 1
operatorNorm = 1;
%create operator
myOperator = functionHandleOperator(myHandle,myHandleT,argumentSize,operatorNorm);

%%

% let us minimize the functional
% min_{u} |u-f1|_2^2 + |myHandle(u)|_2^2

main = flexBox;
%add primal vars v_1,v_2

numberU = main.addPrimalVar(size(f1));

%add data-fidelities: 
main.addTerm(L2dataTerm(1,f1),numberU);
main.addTerm(L2operator(1,1,myOperator),numberU);

main.runAlgorithm;

%get result
result = main.getPrimal(numberU);

figure(3);imagesc(result,[0,1]);axis image;colormap(gray);colorbar;title('Result Image');