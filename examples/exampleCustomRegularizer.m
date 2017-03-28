%clear all;close all;clc;

%% read data
addpath(genpath('..'));
f1 = imread('data/a.png');
f2 = imread('data/b.png');

if (size(f1,3) > 1)
    f1 = rgb2gray(f1);
end
if (size(f2,3) > 1)
    f2 = rgb2gray(f2);
end

f1 = im2double(f1);
f2 = im2double(f2);

figure(1);imagesc(f1);axis image;colormap(gray);title('Image 1')
figure(2);imagesc(f2);axis image;colormap(gray);title('Image 2')
%%
% let us minimize the functional
% min_{u_1,u_2} |u_1-f1|_2^2 + |u_2-f2|_2^2 + 0.5*| sqrt(  (u_1)^2 + (u_1 -
% u_2)^2  ) |_1

main = flexBox;
main.params.tryCPP = 0; %change, if C++ module is compiled

%add primal vars v_1,v_2
numberU1 = main.addPrimalVar(size(f1));
numberU2 = main.addPrimalVar(size(f2));

%add data-fidelities:
main.addTerm(L2dataTerm(1,f1),numberU1);
main.addTerm(L2dataTerm(1,f2),numberU2);

% the regularizer can be written as |K(u_1,u_2)^T|_{1,2} with K = [I & 0\\
% I & -I]. We use the L1operatorIso class which has the input (alpha,numVars,operator)
% where the operator is given row-wise as cell-array of blocks.
% Here, the blocks could be specified as regular sparse matrices, but flexBox
% contains the more efficient equivalent classes identityOperator and zeroOperator
nPx = numel(f1);

main.addTerm(L1operatorIso(1,2,{identityOperator(nPx),zeroOperator(nPx),identityOperator(nPx),-identityOperator(nPx)}),[numberU1,numberU2]);

tic;main.runAlgorithm;toc;

%get result
result1 = main.getPrimal(numberU1);
result2 = main.getPrimal(numberU2);

figure(3);imagesc(result1);axis image;colormap(gray);title('Image 1');
figure(4);imagesc(result2);axis image;colormap(gray);title('Image 2');
