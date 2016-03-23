/*
% Author Information: 
% Hendrik Dirks
% Institute for Computational and Applied Mathematics
% University of Muenster, Germany
%
% Contact: hendrik.dirks@wwu.de
%
%
% Version 1.0
% Date: 2015-06-17

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
#include "math.h"
#include <omp.h>
#include "tools.h"
#include <stdio.h>
#include <sys/types.h>



#include <cstddef>
#include <ctime>

#include "flexMatrix.h"
#include "flexBox.h"

//primal
#include "flexL1dataTerm.h"
#include "flexL2dataTerm.h"
#include "flexEmptyDataTerm.h"

//dual
#include "flexL1dualizedOperatorAniso.h"
#include "flexL1dualizedOperatorIso.h"
#include "flexL2dualizedOperator.h"
#include "flexFrobeniusDualizedOperator.h"


#include "flexL1DualizedDataTerm.h"
#include "flexL2dualizedDataTerm.h"
#include "flexKLdualizedDataTerm.h"

using namespace std;

void copyMatlabToFlexmatrix(const mxArray *input, flexMatrix<float> &output);
void copyMatlabToFlexmatrixNew(const mxArray *input, flexMatrix<float> &output);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L1TVOpticalFlow(image1,image2,tol,lambda,maxIterations,norm,inputV,inputY,stepsize,discretization,numberOfWarps)
{
	char *func_type, *class_name;

	//List of primal terms
	//flexVector<flexTermPrimal<float>* > termsPrimal;
	//List of dual terms
	//flexVector<flexTermDual<float>* > termsDual;

	// Initialize main flexBox object
	flexBox<float> mainObject;
	mainObject.isMATLAB = true;

	int entry = 0;
	while (entry < nrhs-2)
	{
		int oldEntry = entry;
		//printf("%d\n", entry);
		//check type and name
		func_type = mxArrayToString(prhs[entry]);
		
		//printf("%s\n", func_type);
		if (strcmp(func_type, "primalVar") == 0)
		{
			//extract vector dims from second argument
			flexVector<int> _dims;
			double *input_dims = mxGetPr(prhs[entry + 1]);
			for (int i = 0; i < mxGetN(prhs[entry + 1]) * mxGetM(prhs[entry + 1]); ++i)
			{
				_dims.push_back((int)input_dims[i]);
			}
			//printf("%d\n", _dims[0]);
			//printf("%d\n", _dims[1]);

			//add primal variable
			mainObject.addPrimalVar(_dims);

			//jump two entries
			entry += 2;
		}
		else if (strcmp(func_type, "primal") == 0)
		{
			class_name = mxArrayToString(prhs[entry + 1]);
			//printf("%s\n", class_name);
			
			//all primal term have alpha as first arguments

			//first argument is factor:
			float alpha = (float)mxGetScalar(prhs[entry + 2]);

			//second argument is list of corresponding primals:
			flexVector<int> _correspondingPrimals;
			double *input_correspondingPrimals = mxGetPr(prhs[entry + 3]);
			for (int i = 0; i < mxGetN(prhs[entry + 3]) * mxGetM(prhs[entry + 3]); ++i)
			{
				//decrease number by 1 because C++ internal counter starts at 0
				_correspondingPrimals.push_back((int)input_correspondingPrimals[i] - 1);
			}

			if (strcmp(class_name, "L2dataTerm") == 0)
			{
				//third is data array
				double *input_f = mxGetPr(prhs[entry + 4]);
				int size_f = mxGetN(prhs[entry + 4]) * mxGetM(prhs[entry + 4]);
				flexVector<float> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = (float)input_f[i];
				}

				//append primal term to list of primal terms
				mainObject.addPrimal(new flexL2DataTerm<float>(alpha, f), _correspondingPrimals);

				//increase entry counter by number of arguments
				entry += 5;
			}
			else if (strcmp(class_name, "L1dataTerm") == 0)
			{
				//third is data array
				double *input_f = mxGetPr(prhs[entry + 4]);
				int size_f = mxGetN(prhs[entry + 4]) * mxGetM(prhs[entry + 4]);
				flexVector<float> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = (float)input_f[i];
				}

				//append primal term to list of primal terms
				mainObject.addPrimal(new flexL1DataTerm<float>(alpha, f),_correspondingPrimals);

				//increase entry counter by number of arguments
				entry += 5;
			}
			else if (strcmp(class_name, "emptyDataTerm") == 0)
			{
				//append primal term to list of primal terms
				mainObject.addPrimal(new flexEmptyDataTerm<float>(), _correspondingPrimals);

				//increase entry counter by number of arguments
				entry += 4;
			}
		}
		else if (strcmp(func_type, "dual") == 0)
		{
			class_name = mxArrayToString(prhs[entry + 1]);
			//printf("Dual: %s\n", class_name);

			//first argument is factor:
			float alpha = (float)mxGetScalar(prhs[entry + 2]);

			//printf("Factor is %f\n", alpha);

			int nfields = mxGetNumberOfElements(prhs[entry + 3]);
			//printf("Number of Operators is %d\n", nfields);

			flexVector<flexMatrix<float> > operatorList;
			for (int i = 0; i < nfields; ++i)
			{
				const mxArray *pointerA = mxGetCell(prhs[entry + 3], i);

				int m = mxGetM(pointerA);
				int n = mxGetN(pointerA);

				flexMatrix<float> A(mxGetM(pointerA), mxGetN(pointerA));

				copyMatlabToFlexmatrixNew(pointerA, A);

				operatorList.push_back(A);
			}

			double *input_corrPrimal = mxGetPr(prhs[entry + 4]);
			int size_corrPrimal = mxGetN(prhs[entry + 4]) * mxGetM(prhs[entry + 4]);
			flexVector<int> corrPrimal(size_corrPrimal, 0);
			for (int i = 0; i < size_corrPrimal; ++i)
			{
				corrPrimal[i] = (int)input_corrPrimal[i] - 1;
			}

			if (strcmp(class_name, "L1dualizedOperatorIso") == 0)
			{
				mainObject.addDual(new flexL1dualizedOperatorIso<float>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "L1dualizedOperatorAniso") == 0)
			{
				mainObject.addDual(new flexL1dualizedOperatorAniso<float>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "L2dualizedOperator") == 0)
			{
				mainObject.addDual(new flexL2dualizedOperator<float>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "FrobeniusDualizedOperator") == 0)
			{
				mainObject.addDual(new flexFrobeniusDualizedOperator<float>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}


			
			else if (strcmp(class_name, "L1dualizedDataTerm") == 0 || strcmp(class_name, "L2dualizedDataTerm") == 0 || strcmp(class_name, "KLdualizedDataTerm") == 0)
			{
				//6th is data array
				double *input_f = mxGetPr(prhs[entry + 5]);
				int size_f = mxGetN(prhs[entry + 5]) * mxGetM(prhs[entry + 5]);
				flexVector<float> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = (float)input_f[i];
				}

				if (strcmp(class_name, "L1dualizedDataTerm") == 0)
				{
					mainObject.addDual(new flexL1dualizedDataTerm<float>(alpha, operatorList, f), corrPrimal);
				}
				else if (strcmp(class_name, "L2dualizedDataTerm") == 0)
				{
					mainObject.addDual(new flexL2dualizedDataTerm<float>(alpha, operatorList, f), corrPrimal);
				}
				else if (strcmp(class_name, "KLdualizedDataTerm") == 0)
				{
					mainObject.addDual(new flexKLdualizedDataTerm<float>(alpha, operatorList, f), corrPrimal);
				}
				entry += 6;
			}

			//increase entry counter by number of arguments
			
		}

		if (oldEntry == entry)
		{
			mexErrMsgTxt("Problem finding C++ class!\n");
		}
	}

	mainObject.runAlgorithm();

	//retrieve results from FlexBox
	int numPrimalVars = mainObject.getNumPrimalVars();

	printf("Number of primal vars is %d\n", numPrimalVars);

	for (int i = 0; i < numPrimalVars; ++i)
	{
		int numberDims = mainObject.getDims(i).size();

		int numberElements = 1;
		size_t *resultSize = new size_t[numberDims];
		for (int j = 0; j < numberDims; ++j)
		{
			resultSize[j] = mainObject.getDims(i)[j];
			numberElements *= resultSize[j];
		}

		plhs[i] = mxCreateNumericArray(numberDims, resultSize, mxDOUBLE_CLASS, mxREAL);
		double *ptrResult = mxGetPr(plhs[i]);

		flexVector<float> flexResult = mainObject.getPrimal(i);

		for (int j = 0; j < numberElements; ++j)
		{
			ptrResult[j] = flexResult[j];
		}
	}

	//input: Operator A , data b, tau, sigma
	//

	/*int rows = mxGetN(prhs[0]); 
	int cols = mxGetM(prhs[0]);

	flexMatrix<float> A(rows, cols);

	copyMatlabToFlexmatrix(prhs[0], A);

	double *input_b = mxGetPr(prhs[1]);
	flexVector<float> b(cols,0.0f);
	for (int i = 0; i < cols; ++i)
	{
		b[i] = input_b[i];
	}

	flexBox<float> mainObject;

	flexVector<int> _dims; _dims.push_back(cols); _dims.push_back(1);

	flexEmptyDataterm<float> dataTerm(_dims);

	flexL1DualizedDataterm<float> reg1(1.0f, A, b, 4.0f, 4.0f);

	//b.print();
	//A.printMatrix();
	
	mainObject.addPrimal(&dataTerm);

	flexVector<int> corrDual; corrDual.push_back(0);
	

	mainObject.addDual(&reg1, corrDual);

	mainObject.runAlgorithm();

	//write output
	flexVector<float> flexResult = mainObject.getPrimal(0);

	const mwSize resultSize[2] = { cols, 1 };

	//flexResult.print();

	// Output v1
	plhs[0] = mxCreateNumericArray(2, resultSize, mxDOUBLE_CLASS, mxREAL);
	double *result = mxGetPr(plhs[0]);
	for (int i = 0; i < rows; ++i)
	{
		result[i] = flexResult[i];
	}

	//int m = 
	//int n = 

	/*
	const char* infile = "a.png";

	double *u = mxGetPr(prhs[0]);

	const mwSize *sizeImage = mxGetDimensions(prhs[0]);

	int nPx = (int)(sizeImage[0] * sizeImage[1]);

	flexBox<float> mainObject;

	flexVector<int> _dims;_dims.push_back(sizeImage[0]);_dims.push_back(sizeImage[1]);

	float _alpha = 1;

	flexVector<float> _f(nPx);

	//define gradient matrix
	flexMatrix<float> Dx(nPx, nPx), DxT(nPx, nPx), Dy(nPx, nPx), DyT(nPx, nPx);
	for (int j = 0; j < sizeImage[1] - 1; ++j)
	{
		for (int i = 0; i < sizeImage[0] - 1; ++i)
		{
			_f[index2DtoLinear(sizeImage, i, j)] = (float)u[index2DtoLinear(sizeImage, i, j)];


			Dx.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i + 1, j), 1.0);
			Dx.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j), -1.0);

			Dy.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j + 1), 1.0);
			Dy.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j), -1.0);
		}
	}


	flexL2DataTerm<float> dataTerm(_dims, _alpha, _f);

	mainObject.addPrimal(dataTerm);

	flexVector<flexMatrix<float> > _operatorList(2);
	_operatorList[0] = Dx;
	_operatorList[1] = Dy;

	flexVector<float> sigmaTau;
	sigmaTau.push_back(4);
	sigmaTau.push_back(4);

	flexL2Gradient<float> regularizer(0.2, _operatorList, sigmaTau, sigmaTau);
	mainObject.addDual(regularizer,0);


	
	double *u1 = mxGetPr(prhs[0]);
	double *u2 = mxGetPr(prhs[1]);

	float tol = (float)mxGetScalar(prhs[2]);
	float lambda = (float)mxGetScalar(prhs[3]);

	int maxIterations = (int)mxGetScalar(prhs[4]);

	const mwSize *sizeImage = mxGetDimensions(prhs[0]);

	int nPx = (int)(sizeImage[0] * sizeImage[1]);

	const mwSize sizeY[2] = { 4 * nPx,1 };

	// Output v1
	plhs[0] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
	double *Outv1 = mxGetPr(plhs[0]);

	// Output v2
	plhs[1] = mxCreateNumericArray(2, sizeImage, mxDOUBLE_CLASS, mxREAL);
	double *Outv2 = mxGetPr(plhs[1]);

	flexVector<float> ux(nPx), uy(nPx), ut(nPx), uxut(nPx), uyut(nPx);

	//float* c1 = new float[nPx];
	//float* c2 = new float[nPx];
	//float* c3 = new float[nPx];
	//float* teiler = new float[nPx];

	float tau = 0.25;
	float sigma = 0.5;

	//residuals
	float p = 0;
	float d = 0;
	float err = 1.0;

	float ssl = 1.0f - sigma / (sigma + lambda);

	flexVector<float> c1(nPx), c2(nPx), c3(nPx), teiler(nPx), b1(nPx), b2(nPx);

	flexVector<float> v1(nPx), v2(nPx), v1Old(nPx), v2Old(nPx);

	flexVector<float> y11(nPx), y12(nPx), y21(nPx), y22(nPx), y11Old(nPx), y12Old(nPx), y21Old(nPx), y22Old(nPx);
	flexVector<float> Kx11(nPx), Kx12(nPx), Kx21(nPx), Kx22(nPx), Kx11Old(nPx), Kx12Old(nPx), Kx21Old(nPx), Kx22Old(nPx);
	flexVector<float> Kty1(nPx), Kty2(nPx), Kty1Old(nPx), Kty2Old(nPx);

	flexVector<float> p1(nPx), p2(nPx);
	flexVector<float> d1(nPx), d2(nPx), d3(nPx), d4(nPx);

	
	//define gradient matrix
	fleMatrix<float> Dx(nPx, nPx), DxT(nPx, nPx), Dy(nPx, nPx), DyT(nPx, nPx);
	for (int j = 0; j < sizeImage[1]-1; ++j)
	{
		for (int i = 0; i < sizeImage[0]-1; ++i)
		{
			Dx.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i + 1, j), 1.0);
			Dx.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i , j), -1.0);

			Dy.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j + 1), 1.0);
			Dy.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j), -1.0);
		}
	}


	
	DxT = Dx; DxT.transpose();
	DyT = Dy; DyT.transpose();

	//Dx.printMatrix();
	//DxT.printMatrix();
	
	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			//Index for gradients
			ut[tmpIndex] = (float)(u2[tmpIndex] - u1[tmpIndex]);
			
			if (i>0 && i < sizeImage[0] - 1)
			{
				uy[tmpIndex] = (float)(0.5f * (u1[index2DtoLinear(sizeImage, i + 1, j)] - u1[index2DtoLinear(sizeImage, i - 1, j)]));
			}
			else
			{
				uy[tmpIndex] = 0.0f;
			}

			if (j>0 && j < sizeImage[1] - 1)
			{
				ux[tmpIndex] = (float)(0.5f * (u1[index2DtoLinear(sizeImage, i, j + 1)] - u1[index2DtoLinear(sizeImage, i, j - 1)]));
			}
			else
			{
				ux[tmpIndex] = 0.0f;
			}

			uxut[tmpIndex] = ux[tmpIndex] * ut[tmpIndex];
			uyut[tmpIndex] = uy[tmpIndex] * ut[tmpIndex];

			c1[tmpIndex] = 1.0f + tau * ux[tmpIndex] * ux[tmpIndex];
			c2[tmpIndex] = tau * ux[tmpIndex] * uy[tmpIndex];
			c3[tmpIndex] = 1.0f + tau * uy[tmpIndex] * uy[tmpIndex];

			teiler[tmpIndex] = 1.0f / (c1[tmpIndex] * c3[tmpIndex] - c2[tmpIndex] * c2[tmpIndex]);
		}
	}
	
	int iterations = 0;

	double elapsed_secs = 0.0;
	clock_t begin, end;

	while (err > tol && iterations <= maxIterations)
	{
		++iterations;

		if (iterations % 50 == 0)
		{
			p = 0.0f;
			d = 0.0f;
		}
			Kty1Old = Kty1;
			Kty2Old = Kty2;

			v1Old = v1;
			v2Old = v2;

			y11Old = y11;
			y12Old = y12;
			y21Old = y21;
			y22Old = y22;


		//begin = clock();
		DxT.times(y11, Kty1);
		DyT.timesPlus(y12, Kty1);
		DxT.times(y21, Kty2);
		DyT.timesPlus(y22, Kty2);
		//end = clock();
		//elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;

		b1 = v1 - (Kty1 + uxut) * tau;
		b2 = v2 - (Kty2 + uyut) * tau;

		v1 = (b1 * c3 - c2 * b2) * teiler;
		v2 = (b2 * c1 - c2 * b1) * teiler;

		//primal residual
		if (iterations % 50 == 0)
		{
			p1 = (v1Old - v1) / tau - Kty1Old + Kty1; p1.abs();
			p2 = (v2Old - v2) / tau - Kty2Old + Kty2; p2.abs();

			p = p1.sum() + p2.sum();
		}

		Kx11Old = Kx11;
		Kx12Old = Kx12;
		Kx21Old = Kx21;
		Kx22Old = Kx22;

		//begin = clock();
		Dx.times(v1, Kx11);
		Dy.times(v1, Kx12);
		Dx.times(v2, Kx21);
		Dy.times(v2, Kx22);
		//end = clock();
		//elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
//
		y11 = (y11 + (Kx11 + Kx11 - Kx11Old) * sigma) * ssl;
		y12 = (y12 + (Kx12 + Kx12 - Kx12Old) * sigma) * ssl;
		y21 = (y21 + (Kx21 + Kx21 - Kx21Old) * sigma) * ssl;
		y22 = (y22 + (Kx22 + Kx22 - Kx22Old) * sigma) * ssl;

		if (iterations % 50 == 0)
		{
			d1 = (y11Old - y11) / sigma - Kx11Old + Kx11; d1.abs();
			d2 = (y12Old - y12) / sigma - Kx12Old + Kx12; d2.abs();
			d3 = (y21Old - y21) / sigma - Kx21Old + Kx21; d3.abs();
			d4 = (y22Old - y22) / sigma - Kx22Old + Kx22; d4.abs();

			d = d1.sum() + d2.sum() + d3.sum() + d4.sum();

			err = (d*d + p*p) / nPx;
		}

		if (iterations % 1000 == 0)
		{
			mexPrintf("Iteration %d,Residual %e\n", iterations, err);
			mexEvalString("drawnow;");
		}
	}

	//write output
	#pragma omp parallel for
	for (int j = 0; j < sizeImage[1]; ++j)
	{
		for (int i = 0; i < sizeImage[0]; ++i)
		{
			int tmpIndex = index2DtoLinear(sizeImage, i, j);

			Outv1[tmpIndex] = (double)v1[tmpIndex];
			Outv2[tmpIndex] = (double)v2[tmpIndex];
		}
	}

	mexPrintf("Multi Time %f\n", elapsed_secs);*/
}

void copyMatlabToFlexmatrixNew(const mxArray *input, flexMatrix<float> &output)
{
	double  *pr, *pi;
	mwIndex  *ir, *jc;
	mwSize      col, total = 0;
	mwIndex   starting_row_index, stopping_row_index, current_row_index;
	mwSize      n;

	flexVector<int> indexI(0, 0.0f);
	flexVector<int> indexJ(0, 0.0f);
	flexVector<float> indexVal(0, 0.0f);

	/* Get the starting positions of all four data arrays. */
	pr = mxGetPr(input);
	pi = mxGetPi(input);
	ir = mxGetIr(input);
	jc = mxGetJc(input);

	/* Display the nonzero elements of the sparse array. */
	n = mxGetN(input);
	for (col = 0; col<n; col++)
	{
		starting_row_index = jc[col];
		stopping_row_index = jc[col + 1];
		if (starting_row_index == stopping_row_index)
			continue;
		else
		{
			for (current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++)
			{
				//mexPrintf("\t(%"FMT_SIZE_T"u,%"FMT_SIZE_T"u) = %g\n", ir[current_row_index] + 1, col + 1, pr[total]);
				// matalab uses +1 notation for col and row
				//int indexI = ir[current_row_index];
				//int indexJ = col;
				//float value = pr[total];

				indexI.push_back(ir[current_row_index]);
				indexJ.push_back(col);
				indexVal.push_back(pr[total]);

				//output.insertElement(indexI, indexJ, value);

				total++;
			}
		}
	}

	output.blockInsert(indexI, indexJ, indexVal);
}

void copyMatlabToFlexmatrix(const mxArray *input, flexMatrix<float> &output)
{
	double  *pr, *pi;
	mwIndex  *ir, *jc;
	mwSize      col, total = 0;
	mwIndex   starting_row_index, stopping_row_index, current_row_index;
	mwSize      n;

	/* Get the starting positions of all four data arrays. */
	pr = mxGetPr(input);
	pi = mxGetPi(input);
	ir = mxGetIr(input);
	jc = mxGetJc(input);

	/* Display the nonzero elements of the sparse array. */
	n = mxGetN(input);
	for (col = 0; col<n; col++)
	{
		starting_row_index = jc[col];
		stopping_row_index = jc[col + 1];
		if (starting_row_index == stopping_row_index)
			continue;
		else
		{
			for (current_row_index = starting_row_index; current_row_index < stopping_row_index; current_row_index++)
			{
				//mexPrintf("\t(%"FMT_SIZE_T"u,%"FMT_SIZE_T"u) = %g\n", ir[current_row_index] + 1, col + 1, pr[total]);
				// matalab uses +1 notation for col and row
				int indexI = ir[current_row_index];
				int indexJ = col;
				float value = pr[total];
				output.insertElement(indexI, indexJ, value);

				total++;
			}
		}
	}
}