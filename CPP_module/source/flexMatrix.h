#ifndef flexMatrix_H
#define flexMatrix_H

#include "flexLinearOperator.h"

#include <vector>

template < typename T, typename Tvector >
class flexMatrix : public flexLinearOperator<T, Tvector>
{
private:
	std::vector<int> rowToIndexList;
	std::vector<int> indexList;
	Tvector valueList;

public:
	flexMatrix(void) : indexList(), valueList(), rowToIndexList(), flexLinearOperator<T, Tvector>(0, 0, matrixGPUOp){};

	flexMatrix(int  _numRows, int  _numCols) : rowToIndexList(_numRows + 1, static_cast<int>(0)), indexList(0, 0), valueList(0, 0), flexLinearOperator<T, Tvector>(_numRows, _numCols, matrixGPUOp){};

	flexMatrix<T, Tvector>* copy()
	{
		flexMatrix<T, Tvector>* A = new flexMatrix<T, Tvector>(this->getNumRows(), this->getNumCols());
		
		A->rowToIndexList = rowToIndexList;
		A->indexList = indexList;
		A->valueList = valueList;

		return A;
	}


	void times(const Tvector &input, Tvector &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < this->getNumRows(); ++i)
		{
			T rowsum = (T)0;
			// initialize result
			int indexNext = rowToIndexList[i + 1];
			for (int index = rowToIndexList[i]; index < indexNext; ++index)
			{
				rowsum += input[indexList[index]] * valueList[index];
			}
			output[i] = rowsum;
		}
	}

	void timesPlus(const Tvector &input, Tvector &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < this->getNumRows(); ++i)
		{
			T rowsum = (T)0;
			// initialize result
			int indexNext = rowToIndexList[i + 1];
			for (int index = rowToIndexList[i]; index < indexNext; ++index)
			{
				rowsum += input[indexList[index]] * valueList[index];
			}
			output[i] += rowsum;
		}
	}

	void timesMinus(const Tvector &input, Tvector &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < this->getNumRows(); i++)
		{
			T rowsum = (T)0;
			// initialize result
			int indexNext = rowToIndexList[i + 1];
			for (int index = rowToIndexList[i]; index < indexNext; ++index)
			{
				rowsum += input[indexList[index]] * valueList[index];
			}
			output[i] -= rowsum;
		}
	}

	T timesElement(int index, const T* input)
	{
		T rowsum = (T)0;
		// initialize result
		int indexNext = this->rowToIndexList[index + 1];
		for (int elementIndex = this->rowToIndexList[index]; elementIndex < indexNext; ++elementIndex)
		{
			rowsum += input[this->indexList[elementIndex]] * this->valueList[elementIndex];
		}

		return rowsum;
	}

	//this is the fast way to fill flexMatrix
	void blockInsert(std::vector<int> &indexI,const  std::vector<int> &indexJ,const Tvector &indexVal)
	{
		//clear matrix
		//clear();

		int numberListElements = (int)indexI.size();

		//initialize vecvector
		std::vector<int> emptyBucket(0, 0);
		std::vector < std::vector<int> > buckets(this->getNumRows(), emptyBucket);

		//add elements to buckets
		for (int indexInput = 0; indexInput < numberListElements; indexInput++)
		{
			int bucketIndex = indexI[indexInput];
			buckets[bucketIndex].push_back(indexInput);
		}

		//go trough all rows:
		for (int indexRow = 0; indexRow < this->getNumRows(); indexRow++)
		{
			int numElements = 0;

			//go through bucket
			for (int indexBucket = 0; indexBucket < (int)buckets[indexRow].size(); indexBucket++)
			{
				int tmpIndex = buckets[indexRow][indexBucket];

				indexList.push_back(indexJ[tmpIndex]);
				valueList.push_back(indexVal[tmpIndex]);
				++numElements;
			}

			//update rowToIndexList
			rowToIndexList[indexRow + 1] = rowToIndexList[indexRow] + numElements;
		}

	}

	/*
	//inserts new matrix element val at position [i][j]. This is SLOW!
	void insertElement(int i, int j, T val)
	{
		//get start position of next row
		int startIndexNextRow = rowToIndexList[i + 1];

		int numElt = indexList.size();

		//increment size of index and value list by 1
		indexList.push_back(0);
		valueList.push_back(static_cast<T>(0));
		//indexList.resize(indexList.size() + 1,static_cast<T>(0));
		//valueList.resize(valueList.size() + 1,static_cast<T>(0));

		//shift all elements starting with startIndexNextRow to next position
		for (int index = indexList.size()-1; index > startIndexNextRow; index--)
		{
			indexList[index] = indexList[index - 1];
			valueList[index] = valueList[index - 1];
		}

		//update indexList and valueList at current position
		indexList[startIndexNextRow] = j;
		valueList[startIndexNextRow] = val;

		//increase all elemets above i in rowToIndexList
		for (int index = i + 1; index < numRows+1; index++)
		{
			++rowToIndexList[index];
		}
	}*/

	//clear current matrix
	void clear()
	{
		indexList.clear();
		valueList.clear();

		//set all rowToIndexList to 0
		for (int i = 0; i < this->getNumRows() + 1; i++)
		{
			rowToIndexList[i] = 0;
		}
	}

	T getMaxRowSumAbs()
	{
		T maxSum = static_cast<T>(0);

		for (int i = 0; i < this->getNumRows(); ++i)
		{
			T tmpSum = static_cast<T>(0);
			for (int index = rowToIndexList[i]; index < rowToIndexList[i + 1]; ++index)
			{
				tmpSum += abs(valueList[index]);
			}

			if (tmpSum > maxSum)
			{
				maxSum = tmpSum;
			}
		}

		return maxSum;
	}


	std::vector<T> getAbsRowSum()
	{
		std::vector<T> result(this->getNumRows());
		#pragma omp parallel for
		for (int k = 0; k < this->getNumRows(); ++k)
		{
			T tmpSum = static_cast<T>(0);
			for (int index = rowToIndexList[k]; index < rowToIndexList[k + 1]; ++index)
			{
				tmpSum += std::abs(valueList[index]);
			}


			result[k] = tmpSum;
		}

		return result;
	}


	void printRow(int i)
	{
		for (int index = rowToIndexList[i]; index < rowToIndexList[i+1]; ++index)
		{
			printf("(%d,%d,%f)|", i, indexList[index], valueList[index]);
		}

		printf("\n");

	}
	void printMatrix()
	{
		for (int i = 0; i < this->getNumRows(); i++)
		{
			printRow(i);
		}
	}

	//transpose current matrix
	void transpose()
	{
		std::vector<int> tmpindexI(0, 0);
		std::vector<int> tmpindexJ(0, 0);
		Tvector tmpindexVal(0,(T)0);

		for (int i = 0; i < this->getNumRows(); i++)
		{
			for (int index = rowToIndexList[i]; index < rowToIndexList[i + 1]; ++index)
			{
				tmpindexI.push_back(indexList[index]);
				tmpindexJ.push_back(i);
				tmpindexVal.push_back(valueList[index]);
			}
		}

		this->clear();

		int numRowsTmp = this->getNumRows();
		this->setNumRows(this->getNumCols());
		this->setNumCols(numRowsTmp);
		this->rowToIndexList.resize(this->getNumRows() + 1);

		this->blockInsert(tmpindexI, tmpindexJ, tmpindexVal);
	}
};

#endif