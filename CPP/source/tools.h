/*
% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.1
% Date: 2015-06-23
%
% History:
% 1.1: Added Upwind and dt3

% All Rights Reserved
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for any purpose other than its incorporation into a
% commercial product is hereby granted without fee, provided that the
% above copyright notice appear in all copies and that both that
% copyright notice and this permission notice appear in supporting
% documentation, and that the name of the author and University of Muenster not be used in
% advertising or publicity pertaining to distribution of the software
% without specific, written prior permission.
*/
#include "mex.h"

float dxp(const float *data, const mwSize *sizeMat, int i, int j);
float dyp(const float *data, const mwSize *sizeMat, int i, int j);
float dxm(const float *data, const mwSize *sizeMat, int i, int j);
float dym(const float *data, const mwSize *sizeMat, int i, int j);
int index2DtoLinear(const mwSize *sizeMat, int i, int j);
int index3DtoLinear(const mwSize *sizeMat, int i, int j, int k);
float myAbs(float x);
float myMin(float a, float b);
float myMax(float a, float b);
float myPow2(float x);
float dxc(const float *data,const float *u, const mwSize *sizeMat, int i, int j);
float dyc(const float *data,const float *u, const mwSize *sizeMat, int i, int j);
float dxcT(const float *data,const float *u, const mwSize *sizeMat, int i, int j);
float dycT(const float *data,const float *u, const mwSize *sizeMat, int i, int j);
float dxp3(const float *data, const mwSize *sizeMat, int i, int j, int k);
float dyp3(const float *data, const mwSize *sizeMat, int i, int j, int k);
float dtp3(const float *data, const mwSize *sizeMat, int i, int j, int k);
float dxm3(const float *data, const mwSize *sizeMat, int i, int j, int k);
float dym3(const float *data, const mwSize *sizeMat, int i, int j, int k);
float dtm3(const float *data, const mwSize *sizeMat, int i, int j, int k);

float dxc3(const float *data, const float *u, const mwSize *sizeMat, int i, int j, int k);
float dyc3(const float *data, const float *u, const mwSize *sizeMat, int i, int j, int k);
float dxcT3(const float *data, const float *u, const mwSize *sizeMat, int i, int j, int k);
float dycT3(const float *data, const float *u, const mwSize *sizeMat, int i, int j, int k);

float dxUpwind(const float *data,const float *v, const mwSize *sizeMat, int i, int j, int k);
float dyUpwind(const float *data,const float *v, const mwSize *sizeMat, int i, int j, int k);
float dxUpwindT(const float *data,const float *v, const mwSize *sizeMat, int i, int j, int k);
float dyUpwindT(const float *data,const float *v, const mwSize *sizeMat, int i, int j, int k);

int linearTo3Di(const mwSize *sizeMat, int index);
int linearTo3Dj(const mwSize *sizeMat, int index);
int linearTo3Dk(const mwSize *sizeMat, int index);