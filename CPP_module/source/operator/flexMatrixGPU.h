#ifndef flexMatrixGPU_H
#define flexMatrixGPU_H



#include <cuda_runtime.h>
#include <cusparse_v2.h>


#include "flexLinearOperator.h"

#include <vector>

template < typename T, typename Tvector >
class flexMatrixGPU : public flexLinearOperator<T, Tvector>
{
private:
	cusparseHandle_t handle;    
	cusparseMatDescr_t descrA;

	int* listRowEntries;
	int* listColIndices;
	T* listValues;

	int nnz;

	bool transposed;

public:
	flexMatrixGPU(int  _numRows, int  _numCols, int* rowList, int *colList, T* indexVal, bool formatCRS) : flexLinearOperator<T, Tvector>(_numRows, _numCols, matrixGPUOp)
	{
		transposed = false;

		//create sparse matrix
		cusparseCreate(&this->handle);
		cusparseCreateMatDescr(&this->descrA);

		cusparseSetMatType(this->descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(this->descrA, CUSPARSE_INDEX_BASE_ZERO);

		//if formatCRS is true then the input data is already in compressed row storage format, otherwise we have to convert it
		if (formatCRS)
		{
			this->nnz = rowList[_numRows]; //access last entry
		}
		else
		{
			this->nnz = colList[_numCols]; //access last entry
		}

        cudaMalloc(&this->listValues, this->nnz * sizeof(T)); CUDA_CHECK;
		cudaMalloc(&this->listColIndices, this->nnz * sizeof(int)); CUDA_CHECK;
		cudaMalloc(&this->listRowEntries, (_numRows + 1) * sizeof(int)); CUDA_CHECK;

		if (formatCRS == false)
		{
			//copy input to device memory
			T* listValuesTmp;
			int* listColIndicesTmp;
			int* listRowEntriesTmp;

			cudaMalloc(&listValuesTmp, this->nnz * sizeof(T));	CUDA_CHECK;
			cudaMemcpy(listValuesTmp, indexVal, this->nnz * sizeof(T), cudaMemcpyHostToDevice); CUDA_CHECK;
			cudaMalloc(&listColIndicesTmp, (_numCols + 1) * sizeof(int));	CUDA_CHECK;
			cudaMemcpy(listColIndicesTmp, colList, (_numCols + 1) * sizeof(int), cudaMemcpyHostToDevice); CUDA_CHECK;
			cudaMalloc(&listRowEntriesTmp, this->nnz  * sizeof(int));	CUDA_CHECK;
			cudaMemcpy(listRowEntriesTmp, rowList, this->nnz * sizeof(int), cudaMemcpyHostToDevice); CUDA_CHECK;


			cudaDeviceSynchronize();
			cusparseStatus_t status = cusparseScsr2csc(this->handle, _numCols, _numRows, this->nnz, listValuesTmp, listColIndicesTmp, listRowEntriesTmp, this->listValues, this->listColIndices, this->listRowEntries, CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO); CUDA_CHECK;
			cudaDeviceSynchronize();

			cudaFree(listValuesTmp); CUDA_CHECK;
			cudaFree(listColIndicesTmp); CUDA_CHECK;
			cudaFree(listRowEntriesTmp); CUDA_CHECK;

			if (VERBOSE > 0)
			{
				switch (status)
				{
				case CUSPARSE_STATUS_SUCCESS:
				{
					printf("Copy was successfull\n");
					break;
				}
				case CUSPARSE_STATUS_NOT_INITIALIZED:
				{
					printf("the library was not initialized\n");
					break;
				}
				case CUSPARSE_STATUS_ALLOC_FAILED:
				{
					printf("the resources could not be allocated\n");
					break;
				}
				case CUSPARSE_STATUS_INVALID_VALUE:
				{
					printf("invalid parameters were passed(m, n, nnz<0)\n");
					break;
				}
				case CUSPARSE_STATUS_ARCH_MISMATCH:
				{
					printf("the device does not support double precision\n");
					break;
				}
				case CUSPARSE_STATUS_EXECUTION_FAILED:
				{
					printf("the function failed to launch on the GPU\n");
					break;
				}
				case CUSPARSE_STATUS_INTERNAL_ERROR:
				{
					printf("the function failed to launch on the GPU\n");
					break;
				}
				default:
				{
					printf("Error Copy!");
					break;
				}
				}
			}
		}
		else
		{
			cudaMemcpy(this->listValues, indexVal, this->nnz * sizeof(T), cudaMemcpyHostToDevice);
			cudaMemcpy(this->listColIndices, colList, this->nnz * sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(this->listRowEntries, rowList, (_numRows + 1) * sizeof(int), cudaMemcpyHostToDevice);
		}
	};

	~flexMatrixGPU()
	{
		if (VERBOSE > 0) printf("MatrixGPU destructor!");
		//free cuda memory
		cudaFree(this->listValues); CUDA_CHECK;
		cudaFree(this->listColIndices); CUDA_CHECK;
		cudaFree(this->listRowEntries); CUDA_CHECK;
	}

	flexMatrixGPU<T, Tvector>* copy()
	{
		//copy matrix data to host

		//allocate memory
		T *hostValues = (T *)malloc(this->nnz * sizeof(T));
		int *hostRowIndices = (int *)malloc((this->getNumRows() + 1) * sizeof(int));
		int *hostColIndices = (int *)malloc(this->nnz * sizeof(int));

		cudaMemcpy(hostValues, this->listValues, this->nnz * sizeof(T), cudaMemcpyDeviceToHost);
		cudaMemcpy(hostRowIndices, this->listRowEntries, (this->getNumRows() + 1) * sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(hostColIndices, this->listColIndices, this->nnz * sizeof(int), cudaMemcpyDeviceToHost);

		flexMatrixGPU<T, Tvector>* A = new flexMatrixGPU<T, Tvector>(this->getNumRows(), this->getNumCols(), hostRowIndices, hostColIndices, hostValues,true);

		free(hostValues);
		free(hostRowIndices);
		free(hostColIndices);

		return A;
	}

	void times(const Tvector &input, Tvector &output)
	{
		const T alpha = (T)1;
		const T beta = (T)0;

		T* ptrOutput = thrust::raw_pointer_cast(output.data());
		const T* ptrInput = thrust::raw_pointer_cast(input.data());

		if (transposed == false)
		{
			cusparseScsrmv(this->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, this->getNumRows(), this->getNumCols(), nnz, &alpha, this->descrA, this->listValues, this->listRowEntries, this->listColIndices, ptrInput, &beta, ptrOutput);
		}
		else
		{
			cusparseScsrmv(this->handle, CUSPARSE_OPERATION_TRANSPOSE, this->getNumRows(), this->getNumCols(), nnz, &alpha, this->descrA, this->listValues, this->listRowEntries, this->listColIndices, ptrInput, &beta, ptrOutput);
		}
	}

	void timesPlus(const Tvector &input, Tvector &output)
	{
		const T alpha = (T)1;
		const T beta = (T)1;

		T* ptrOutput = thrust::raw_pointer_cast(output.data());
		const T* ptrInput = thrust::raw_pointer_cast(input.data());

		if (transposed == false)
		{
			cusparseScsrmv(this->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, this->getNumRows(), this->getNumCols(), nnz, &alpha, this->descrA, this->listValues, this->listRowEntries, this->listColIndices, ptrInput, &beta, ptrOutput);
		}
		else
		{
			cusparseScsrmv(this->handle, CUSPARSE_OPERATION_TRANSPOSE, this->getNumRows(), this->getNumCols(), nnz, &alpha, this->descrA, this->listValues, this->listRowEntries, this->listColIndices, ptrInput, &beta, ptrOutput);
		}
	}

	void timesMinus(const Tvector &input, Tvector &output)
	{
		const T alpha = -(T)1;
		const T beta = (T)1;

		T* ptrOutput = thrust::raw_pointer_cast(output.data());
		const T* ptrInput = thrust::raw_pointer_cast(input.data());

		if (transposed == false)
		{
			cusparseScsrmv(this->handle, CUSPARSE_OPERATION_NON_TRANSPOSE, this->getNumRows(), this->getNumCols(), nnz, &alpha, this->descrA, this->listValues, this->listRowEntries, this->listColIndices, ptrInput, &beta, ptrOutput);
		}
		else
		{
			cusparseScsrmv(this->handle, CUSPARSE_OPERATION_TRANSPOSE, this->getNumRows(), this->getNumCols(), nnz, &alpha, this->descrA, this->listValues, this->listRowEntries, this->listColIndices, ptrInput, &beta, ptrOutput);
		}
	}

	T getMaxRowSumAbs()
	{
		//todo

		return 1;
	}

	std::vector<T> getAbsRowSum()
	{
		std::vector<T> result(this->getNumRows(),1);

		return result;
	}

	void printRow(int i)
	{

	}
	void printMatrix()
	{
		T *hostValues = (T *)malloc(this->nnz * sizeof(T));
		int *hostRowIndices = (int *)malloc((this->getNumRows() + 1) * sizeof(int));
		int *hostColIndices = (int *)malloc(this->nnz * sizeof(int));

		cudaMemcpy(hostValues, this->listValues, this->nnz * sizeof(T), cudaMemcpyDeviceToHost);
		cudaMemcpy(hostRowIndices, this->listRowEntries, (this->getNumRows() + 1) * sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(hostColIndices, this->listColIndices, this->nnz * sizeof(int), cudaMemcpyDeviceToHost);

		cudaDeviceSynchronize();

		free(hostValues);
		free(hostRowIndices);
		free(hostColIndices);
	}

	//transpose current matrix
	void transpose()
	{
		if (transposed == false)
		{
            transposed = true;
		}
		else
		{
			transposed = false;
		}
	}

	T timesElement(int index, const T* input)
	{
		return (T)0;
	}

#if __CUDACC__
	__device__ T timesElementCUDA(int index, const T* input)
	{
		T rowsum = (T)0;
		// initialize result
		int indexNext = this->listRowEntries[index+1];
		for (int elementIndex = this->listRowEntries[index]; elementIndex < indexNext; ++elementIndex)
		{
			rowsum += input[this->listColIndices[elementIndex]] * this->listValues[elementIndex];
		}
		
		return rowsum;
	}

	__device__ T getRowsumElementCUDA(int index)
	{
		T rowsum = (T)0;
		// initialize result
		int indexNext = listRowEntries[index + 1];
		for (int elementIndex = listRowEntries[index]; elementIndex < indexNext; ++elementIndex)
		{
			rowsum += fabs(listValues[elementIndex]);
		}

		return rowsum;
	}
	
	thrust::device_vector<T> getAbsRowSumCUDA()
	{
		Tvector result(this->getNumRows(),(T)1);

		return result;
	}
#endif
};

#endif