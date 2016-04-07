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

float myAbs(float x)
{
	return x > 0 ? x : -x;
}

float myMin(float a, float b)
{
	return a > b ? b : a;
}

float myMax(float a, float b)
{
	return a < b ? b : a;
}



float myPow2(float x)
{
	return x * x;
}