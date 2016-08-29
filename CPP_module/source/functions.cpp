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

//uncomment later and compiler directive
//#define __CUDACC__ 1
//#define DO_CUDA_CHECK 1

#include "mex.h"
#include "math.h"
#include <omp.h>
#include <iostream>

#include <stdio.h>
#include <sys/types.h>
#include <string.h>
#include <cstddef>
#include <ctime>

#include "tools.h"

#include "flexLinearOperator.h"
#include "flexIdentityOperator.h"
#include "flexZeroOperator.h"
#include "flexDiagonalOperator.h"
#include "flexMatrix.h"
#include "flexGradientOperator.h"


#include "flexBox.h"

//primal
#include "flexTermPrimal.h"

//dual
#include "flexDualizedOperator.h"
#include "flexDualizedDataTerm.h"

using namespace std;

typedef float floatingType;

#if __CUDACC__
	#include "flexMatrixGPU.h"

	typedef thrust::device_vector<floatingType> vectorData;
#else
	typedef std::vector<floatingType> vectorData;
#endif



void copyToVector(std::vector<floatingType> &vector, const double *input, int numElements);
bool checkClassType(mxArray *object, char *className);
bool checkSparse(mxArray *object);

void copyMatlabToFlexmatrix(const mxArray *input, flexMatrix<floatingType,vectorData> *output);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	#if __CUDACC__
		printf("Using CUDA\n");
	#endif

	char *func_type, *class_name;
	int verbose = 0;

	// Initialize main flexBox object
	flexBox<floatingType,vectorData> mainObject;
	mainObject.isMATLAB = true;

	int entry = 0;
	while (entry < nrhs-2)
	{
		int oldEntry = entry;

		//check type and name
		func_type = mxArrayToString(prhs[entry]);
		
		if (strcmp(func_type, "parameter") == 0)
		{
			char *parameterName = mxArrayToString(prhs[entry + 1]);

			if (strcmp(parameterName, "maxIt") == 0)
			{
				mainObject.maxIterations = (int)mxGetScalar(prhs[entry + 2]);
			}

			//jump three entries
			entry += 3;
		}
		else if (strcmp(func_type, "primalVar") == 0)
		{
			//extract vector dims from second argument
			std::vector<int> _dims;
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
			std::vector<int> _correspondingPrimals;
			double *input_correspondingPrimals = mxGetPr(prhs[entry + 3]);
			for (int i = 0; i < mxGetN(prhs[entry + 3]) * mxGetM(prhs[entry + 3]); ++i)
			{
				//decrease number by 1 because C++ internal counter starts at 0
				_correspondingPrimals.push_back((int)input_correspondingPrimals[i] - 1);
			}

			if (strcmp(class_name, "emptyDataTerm") == 0)
			{
				//append primal term to list of primal terms
				mainObject.addPrimal(new flexTermPrimal<floatingType, vectorData>(1, alpha, primalEmptyProx), _correspondingPrimals);

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
			std::vector<int> corrPrimal(size_corrPrimal, 0);
			for (int i = 0; i < size_corrPrimal; ++i)
			{
				corrPrimal[i] = (int)input_corrPrimal[i] - 1;
			}

			//create list of operators
			if (verbose > 1)
			{
				printf("Transfering operators:\n");
			}

			std::vector<flexLinearOperator<floatingType, vectorData>*> operatorList;
			for (int i = 0; i < nfields; ++i)
			{
				int correspondingNumberPrimalVar = i%corrPrimal.size();

				int numElementsPrimalVar = vectorProduct(mainObject.getDims(corrPrimal[correspondingNumberPrimalVar]));
				
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
						//gradientTypeInt = 1;
					}

					operatorList.push_back(new flexGradientOperator<floatingType, vectorData>(mainObject.getDims(corrPrimal[correspondingNumberPrimalVar]), gradientDirection, gradientTypeInt));
				}
				else if (checkClassType(pointerA, "identityOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <identityOperator>\n", i);
					}

					bool isMinus = mxGetScalar(mxGetProperty(pointerA, 0, "minus")) > 0;

					operatorList.push_back(new flexIdentityOperator<floatingType, vectorData>(numElementsPrimalVar, numElementsPrimalVar, isMinus));
				}
				else if (checkClassType(pointerA, "zeroOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <zeroOperator>\n", i);
					}

					operatorList.push_back(new flexZeroOperator<floatingType, vectorData>(numElementsPrimalVar, numElementsPrimalVar));
				}
				else if (checkClassType(pointerA, "diagonalOperator"))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <diagonalOperator>\n", i);
					}

					//copy diagonal vector
					double *tmpDiagonalVector = mxGetPr(mxGetProperty(pointerA, 0, "diagonalElements"));
					std::vector<floatingType> tmpDiagonal(numElementsPrimalVar, static_cast<floatingType>(0));

					for (int i = 0; i < numElementsPrimalVar; ++i)
					{
						tmpDiagonal[i] = static_cast<floatingType>(tmpDiagonalVector[i]);
					}
					
					operatorList.push_back(new flexDiagonalOperator<floatingType, vectorData>(tmpDiagonal));
				}
				else if (checkSparse(pointerA))
				{
					if (verbose > 1)
					{
						printf("Operator %d is type <matrix>\n", i);
					}

					#if __CUDACC__
						mwIndex  *ir, *jc;

						jc = mxGetJc(pointerA);
						ir = mxGetIr(pointerA);
						double * pr = mxGetPr(pointerA);

						//matlab stores in compressed column format
						int numCols = mxGetN(pointerA);
						int* colList = new int[numCols + 1];
						for (int i = 0; i <= numCols; ++i)
						{
							colList[i] = jc[i];
						}

						int nnz = colList[numCols];

						int* rowList = new int[nnz];
						float* valList = new float[nnz];
						for (int i = 0; i < nnz; ++i)
						{
							rowList[i] = ir[i];
							valList[i] = pr[i];
						}

						//printf("%d,%d,%d\n%d,%d,%d,%d\n%f,%f,%f", colList[0], colList[1], colList[2], rowList[0], rowList[1], rowList[2], rowList[3], valList[0], valList[1], valList[2]);

						operatorList.push_back(new flexMatrixGPU<floatingType, vectorData>((int)mxGetM(pointerA), (int)mxGetN(pointerA), rowList, colList, valList,false));
					#else
						flexMatrix<floatingType, vectorData>*A = new flexMatrix<floatingType, vectorData>(static_cast<int>(mxGetM(pointerA)), static_cast<int>(mxGetN(pointerA)));
						copyMatlabToFlexmatrix(pointerA, A);
						operatorList.push_back(A);
					#endif

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

			if (strcmp(class_name, "L1dualizedOperatorIso") == 0 || strcmp(class_name, "L1dualizedOperatorAniso") == 0 || strcmp(class_name, "L2dualizedOperator") == 0 || strcmp(class_name, "FrobeniusDualizedOperator") == 0 || strcmp(class_name, "HuberDualizedOperator") == 0  )
			{
				prox proxName;

				if (strcmp(class_name, "L1dualizedOperatorIso") == 0)
				{
					proxName = dualL1IsoProx;
				}
				else if (strcmp(class_name, "L1dualizedOperatorAniso") == 0)
				{
					proxName = dualL1AnisoProx;
				}
				else if (strcmp(class_name, "L2dualizedOperator") == 0)
				{
					proxName = dualL2Prox;
				}
				else if (strcmp(class_name, "HuberDualizedOperator") == 0)
				{
					proxName = dualHuberProx;
				}
				else if (strcmp(class_name, "FrobeniusDualizedOperator") == 0)
				{
					proxName = dualFrobeniusProx;
				}
				else if (strcmp(class_name, "FrobeniusDualizedOperator") == 0)
				{
					proxName = dualFrobeniusProx;
				}

				mainObject.addDual(new flexDualizedOperator<floatingType, vectorData>(proxName, alpha, corrPrimal.size(), operatorList), corrPrimal);
				entry += 5;
			}
			else if (strcmp(class_name, "L1dualizedDataTerm") == 0 || strcmp(class_name, "L2dualizedDataTerm") == 0 || strcmp(class_name, "KLdualizedDataTerm") == 0)
			{
				//6th is data array
				double *input_f = mxGetPr(prhs[entry + 5]);
				int size_f = static_cast<int>(mxGetN(prhs[entry + 5]) * mxGetM(prhs[entry + 5]));
				std::vector<floatingType> f(size_f, 0.0f);
				for (int i = 0; i < size_f; ++i)
				{
					f[i] = (floatingType)input_f[i];
				}

				prox proxName;

				if (strcmp(class_name, "L1dualizedDataTerm") == 0)
				{
					proxName = dualL1DataProx;
				}
				else if (strcmp(class_name, "L2dualizedDataTerm") == 0)
				{
					proxName = dualL2DataProx;
				}
				else if (strcmp(class_name, "KLdualizedDataTerm") == 0)
				{
					proxName = dualKLDataProx;
				}

				mainObject.addDual(new flexDualizedDataTerm<floatingType, vectorData>(proxName, alpha, operatorList, f), corrPrimal);

				entry += 6;
			}

			//increase entry counter by number of arguments
			
		}
		else if (strcmp(func_type, "x") == 0)
		{
			int numberPrimalVars = mainObject.getNumPrimalVars();
			
			for (int i = 0; i < numberPrimalVars; ++i)
			{
				mxArray *xInput = mxGetCell(prhs[entry + 1], i);
				double *xPointerInput = mxGetPr(xInput);

				int size_input = static_cast<int>(mxGetN(xInput) * mxGetM(xInput));

				std::vector<floatingType> tmpVector(size_input, 0.0);

				copyToVector(tmpVector, xPointerInput, size_input);
				
				mainObject.setPrimal(i, tmpVector);
			}

			entry += 2;
		}
		else if (strcmp(func_type, "y") == 0)
		{
			int numberDualVars = mainObject.getNumDualVars();

			for (int i = 0; i < numberDualVars; ++i)
			{
				mxArray *yInput = mxGetCell(prhs[entry + 1], i);
				double *yPointerInput = mxGetPr(yInput);

				int size_input = static_cast<int>(mxGetN(yInput) * mxGetM(yInput));

				std::vector<floatingType> tmpVector(size_input, 0.0);

				copyToVector(tmpVector, yPointerInput, size_input);

				mainObject.setDual(i, tmpVector);
			}

			entry += 2;
		}

		if (oldEntry == entry)
		{
			mexErrMsgTxt("Problem finding C++ class!\n");
		}
	}

	mainObject.runAlgorithm();

	//send content of primal vars

	//retrieve results from FlexBox
	int numPrimalVars = mainObject.getNumPrimalVars();
	for (int i = 0; i < numPrimalVars; ++i)
	{
		std::vector<floatingType> flexResult = mainObject.getPrimal(i); 

		size_t *resultSize = new size_t[2];
		resultSize[0] = flexResult.size();
		resultSize[1] = 1;

		plhs[i] = mxCreateNumericArray(2, resultSize, mxDOUBLE_CLASS, mxREAL);
		double *ptrResult = mxGetPr(plhs[i]);
		for (int j = 0; j < resultSize[0]; ++j)
		{
			ptrResult[j] = flexResult[j];
		}
	}

	//send content of dual vars
	//retrieve results from FlexBox
	int numDualVars = mainObject.getNumDualVars();
	for (int i = 0; i < numDualVars; ++i)
	{
		std::vector<floatingType> flexResult = mainObject.getDual(i); 
		
		size_t *resultSize = new size_t[2];
		resultSize[0] = flexResult.size();
		resultSize[1] = 1;

		plhs[numPrimalVars + i] = mxCreateNumericArray(2, resultSize, mxDOUBLE_CLASS, mxREAL);
		double *ptrResult = mxGetPr(plhs[numPrimalVars+i]);

		for (int j = 0; j < resultSize[0]; ++j)
		{
			ptrResult[j] = flexResult[j];
		}
	}
}

void copyToVector(std::vector<floatingType> &vector,const double *input, int numElements)
{
	for (int j = 0; j < numElements; ++j)
	{
		vector[j] = (floatingType)input[j];
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


void copyMatlabToFlexmatrix(const mxArray *input, flexMatrix<floatingType, vectorData> *output)
{
	double  *pr;
	mwIndex  *ir, *jc;
	mwSize      col, total = 0;
	mwIndex   starting_row_index, stopping_row_index, current_row_index;
	mwSize      n;

	std::vector<int> indexI(0, 0);
	std::vector<int> indexJ(0, 0);
	std::vector<floatingType> indexVal(0, 0.0f);

	pr = mxGetPr(input);
	ir = mxGetIr(input);
	jc = mxGetJc(input);


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
