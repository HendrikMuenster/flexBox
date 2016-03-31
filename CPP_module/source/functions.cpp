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

#include <stdio.h>
#include <sys/types.h>



#include <cstddef>
#include <ctime>
#include <cmath>

#include "tools.h"

#include "flexLinearOperator.h"
#include "flexIdentityOperator.h"
#include "flexZeroOperator.h"
#include "flexDiagonalOperator.h"
#include "flexMatrix.h"
#include "flexGradientOperator.h"


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

typedef float floatingType;

bool checkClassType(mxArray *object, char *className);
bool checkSparse(mxArray *object);
void copyMatlabToFlexmatrix(const mxArray *input, flexMatrix<floatingType> *output);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
//L1TVOpticalFlow(image1,image2,tol,lambda,maxIterations,norm,inputV,inputY,stepsize,discretization,numberOfWarps)
{
	char *func_type, *class_name;
	int verbose = 0;

	// Initialize main flexBox object
	flexBox<floatingType> mainObject;
	mainObject.isMATLAB = true;

	int entry = 0;
	while (entry < nrhs-2)
	{
		int oldEntry = entry;

		//check type and name
		func_type = mxArrayToString(prhs[entry]);
		
		if (strcmp(func_type, "primalVar") == 0)
		{
			//extract vector dims from second argument
			flexVector<int> _dims;
			double *input_dims = mxGetPr(prhs[entry + 1]);
			for (int i = 0; i < mxGetN(prhs[entry + 1]) * mxGetM(prhs[entry + 1]); ++i)
			{
				_dims.push_back((int)input_dims[i]);
			}

			//add primal variable
			mainObject.addPrimalVar(_dims);

			//jump two entries
			entry += 2;
		}
		else if (strcmp(func_type, "primal") == 0)
		{
			class_name = mxArrayToString(prhs[entry + 1]);

			//all primal term have alpha as first argument
			//first argument is factor:
			floatingType alpha = (floatingType)mxGetScalar(prhs[entry + 2]);

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
				int size_f = static_cast<int>(mxGetN(prhs[entry + 4]) * mxGetM(prhs[entry + 4]));
				flexVector<floatingType> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = static_cast<floatingType>(input_f[i]);
				}

				//append primal term to list of primal terms
				mainObject.addPrimal(new flexL2DataTerm<floatingType>(alpha, f), _correspondingPrimals);

				//increase entry counter by number of arguments
				entry += 5;
			}
			else if (strcmp(class_name, "L1dataTerm") == 0)
			{
				//third is data array
				double *input_f = mxGetPr(prhs[entry + 4]);
				int size_f = static_cast<int>(mxGetN(prhs[entry + 4]) * mxGetM(prhs[entry + 4]));
				flexVector<floatingType> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = (floatingType)input_f[i];
				}

				//append primal term to list of primal terms
				mainObject.addPrimal(new flexL1DataTerm<floatingType>(alpha, f),_correspondingPrimals);

				//increase entry counter by number of arguments
				entry += 5;
			}
			else if (strcmp(class_name, "emptyDataTerm") == 0)
			{
				//append primal term to list of primal terms
				mainObject.addPrimal(new flexEmptyDataTerm<floatingType>(), _correspondingPrimals);

				//increase entry counter by number of arguments
				entry += 4;
			}
		}
		else if (strcmp(func_type, "dual") == 0)
		{
			//type of term comes first:
			class_name = mxArrayToString(prhs[entry + 1]);

			//second argument is weight:
			floatingType alpha = (floatingType)mxGetScalar(prhs[entry + 2]);

			//argument 3+4: string describing operator and operator itself
			int nfields = static_cast<int>(mxGetNumberOfElements(prhs[entry + 3]));
			//printf("Number of Operators is %d\n", nfields);

			//extract the corresponding primals first, because they're needed for operators
			double *input_corrPrimal = mxGetPr(prhs[entry + 4]);
			int size_corrPrimal = static_cast<int>(mxGetN(prhs[entry + 4]) * mxGetM(prhs[entry + 4]));
			flexVector<int> corrPrimal(size_corrPrimal, 0);
			for (int i = 0; i < size_corrPrimal; ++i)
			{
				corrPrimal[i] = (int)input_corrPrimal[i] - 1;
			}

			//create list of operators
			if (verbose > 1)
			{
				printf("Transfering operators:\n");
			}

			flexVector<flexLinearOperator<floatingType>*> operatorList;
			for (int i = 0; i < nfields; ++i)
			{
				int correspondingNumberPrimalVar = i%corrPrimal.size();

				int numElementsPrimalVar = mainObject.getDims(corrPrimal[correspondingNumberPrimalVar]).product();

				mxArray *pointerA = mxGetCell(prhs[entry + 3], i);

				if (checkClassType(pointerA, "gradientOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <gradientOperator>\n",i);
					}

					char *gradientType = mxArrayToString(mxGetProperty(pointerA, 0, "type"));
					int gradientDirection = static_cast<int>(mxGetScalar(mxGetProperty(pointerA, 0, "gradDirection"))) - 1; //substract one!

					int gradientTypeInt = 0;
					if (strcmp(gradientType, "backward") == 0)
					{
						gradientTypeInt = 1;
					}

					flexGradientOperator<floatingType>*A = new flexGradientOperator<floatingType>(mainObject.getDims(corrPrimal[correspondingNumberPrimalVar]), gradientDirection, gradientTypeInt);
					operatorList.push_back(A);
				}
				else if (checkClassType(pointerA, "identityOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <identityOperator>\n", i);
					}

					bool isMinus = mxGetScalar(mxGetProperty(pointerA, 0, "minus")) > 0;

					flexIdentityOperator<floatingType>*A = new flexIdentityOperator<floatingType>(numElementsPrimalVar, numElementsPrimalVar, isMinus);
					operatorList.push_back(A);
				}
				else if (checkClassType(pointerA, "zeroOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <zeroOperator>\n", i);
					}

					flexZeroOperator<floatingType>*A = new flexZeroOperator<floatingType>(numElementsPrimalVar, numElementsPrimalVar);
					operatorList.push_back(A);
				}
				else if (checkClassType(pointerA, "diagonalOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <diagonalOperator>\n", i);
					}

					//copy diagonal vector
					double *tmpDiagonalVector = mxGetPr(mxGetProperty(pointerA, 0, "diagonalElements"));
					flexVector<floatingType> tmpDiagonalFlexVector(numElementsPrimalVar, static_cast<floatingType>(0));

					for (int i = 0; i < numElementsPrimalVar; ++i)
					{
						tmpDiagonalFlexVector[i] = static_cast<floatingType>(tmpDiagonalVector[i]);
					}
					
					flexDiagonalOperator<floatingType>*A = new flexDiagonalOperator<floatingType>(tmpDiagonalFlexVector);
					operatorList.push_back(A);
				}
				else if (checkSparse(pointerA))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <matrix>\n", i);
					}

					flexMatrix<floatingType>*A = new flexMatrix<floatingType>(static_cast<int>(mxGetM(pointerA)), static_cast<int>(mxGetN(pointerA)));
					copyMatlabToFlexmatrix(pointerA, A);
					operatorList.push_back(A);
				}
				else
				{
					mexErrMsgTxt("Operator type not supported!\n");
				}
			}

			if (verbose > 1)
			{
				printf("Creating dual term\n");
			}

			if (strcmp(class_name, "L1dualizedOperatorIso") == 0)
			{
				mainObject.addDual(new flexL1dualizedOperatorIso<floatingType>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "L1dualizedOperatorAniso") == 0)
			{
				mainObject.addDual(new flexL1dualizedOperatorAniso<floatingType>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "L2dualizedOperator") == 0)
			{
				mainObject.addDual(new flexL2dualizedOperator<floatingType>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "FrobeniusDualizedOperator") == 0)
			{
				mainObject.addDual(new flexFrobeniusDualizedOperator<floatingType>(alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "L1dualizedDataTerm") == 0 || strcmp(class_name, "L2dualizedDataTerm") == 0 || strcmp(class_name, "KLdualizedDataTerm") == 0)
			{
				//6th is data array
				double *input_f = mxGetPr(prhs[entry + 5]);
				int size_f = static_cast<int>(mxGetN(prhs[entry + 5]) * mxGetM(prhs[entry + 5]));
				flexVector<floatingType> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = (floatingType)input_f[i];
				}

				if (strcmp(class_name, "L1dualizedDataTerm") == 0)
				{
					mainObject.addDual(new flexL1dualizedDataTerm<floatingType>(alpha, operatorList, f), corrPrimal);
				}
				else if (strcmp(class_name, "L2dualizedDataTerm") == 0)
				{
					mainObject.addDual(new flexL2dualizedDataTerm<floatingType>(alpha, operatorList, f), corrPrimal);
				}
				else if (strcmp(class_name, "KLdualizedDataTerm") == 0)
				{
					mainObject.addDual(new flexKLdualizedDataTerm<floatingType>(alpha, operatorList, f), corrPrimal);
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
			numberElements *= static_cast<int>(resultSize[j]);
		}

		plhs[i] = mxCreateNumericArray(numberDims, resultSize, mxDOUBLE_CLASS, mxREAL);
		double *ptrResult = mxGetPr(plhs[i]);

		flexVector<floatingType> flexResult = mainObject.getPrimal(i);

		for (int j = 0; j < numberElements; ++j)
		{
			ptrResult[j] = flexResult[j];
		}
	}
}

bool checkClassType(mxArray *object,char *className)
{

	mxArray *output[1], *input[2];

	input[0] = object;
	input[1] = mxCreateString(className);
	
	mexCallMATLAB(1, output, 2, input, "isa");

	if (mxGetScalar(output[0]) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
	
}

bool checkSparse(mxArray *object)
{

	mxArray *output[1], *input[1];

	input[0] = object;

	mexCallMATLAB(1, output, 1, input, "issparse");

	if (mxGetScalar(output[0]) > 0)
	{
		return true;
	}
	else
	{
		return false;
	}

}

void copyMatlabToFlexmatrix(const mxArray *input, flexMatrix<floatingType> *output)
{
	double  *pr, *pi;
	mwIndex  *ir, *jc;
	mwSize      col, total = 0;
	mwIndex   starting_row_index, stopping_row_index, current_row_index;
	mwSize      n;

	flexVector<int> indexI(0, 0);
	flexVector<int> indexJ(0, 0);
	flexVector<floatingType> indexVal(0, 0.0f);

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
				indexI.push_back(static_cast<int>(ir[current_row_index]));
				indexJ.push_back(static_cast<int>(col));
				indexVal.push_back(static_cast<floatingType>(pr[total]));

				total++;
			}
		}
	}

	output->blockInsert(indexI, indexJ, indexVal);
}