#ifndef flexMatrix_H
#define flexMatrix_H

#include "flexVector.h"

template < typename T >
class flexMatrix
{
private:
	flexVector<int> rowToIndexList;
	flexVector<int> indexList;
	flexVector<T> valueList;

	int numCols;
	int numRows;

public:
	flexMatrix(void) : indexList(), valueList(), rowToIndexList()
	{
		numCols = 0;
		numRows = 0;
	}

	flexMatrix(const int  _numRows, const int  _numCols) : rowToIndexList(_numRows + 1, static_cast<int>(0)), indexList(0, 0), valueList(0, 0)
	{
		numCols = _numCols;
		numRows = _numRows;

		//rowToIndexList.print();
	};

	int getNumCols()
	{
		return numCols;
	}

	int getNumRows()
	{
		return numRows;
	}

	void times(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			output[i] = rowMulti(i,input);
		}
	}

	/*void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			output[i] += rowMulti(i, input);
		}
	}*/

	void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			T rowsum = static_cast<T>(0);
			// initialize result
			int indexNext = rowToIndexList[i + 1];
			for (int index = rowToIndexList[i]; index < indexNext; ++index)
			{
				rowsum += input[indexList[index]] * valueList[index];
			}
			output[i] += rowsum;
		}
	}

	void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			T rowsum = static_cast<T>(0);
			// initialize result
			int indexNext = rowToIndexList[i + 1];
			for (int index = rowToIndexList[i]; index < indexNext; ++index)
			{
				rowsum += input[indexList[index]] * valueList[index];
			}
			output[i] -= rowsum;
		}
	}



	/*void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			output[i] -= rowMulti(i, input);
		}
	}*/

	T rowMulti(const int i, const flexVector<T> &input)
	{
		// initialize result
		T result = static_cast<T>(0);

		int indexNext = rowToIndexList[i + 1];
		for (int index = rowToIndexList[i]; index < indexNext; ++index)
		{
			result += input[indexList[index]] * valueList[index];
		}

		return result;
	}

	/*T rowMulti(const flexVector<T> &input, const flexVector<int> &inputI, const flexVector<T> &inputV)
	{
		// initialize result
		T result = static_cast<T>(0);

		//int elements = inputI._size;

		if (inputI._size > 0)
		{
			for (int j = 0; j < inputI._size; j++)
			{
				result += input[inputI[j]] * inputV[j];
			}
		}
		return result;
	}*/

	void blockInsert(flexVector<int> indexI, flexVector<int> indexJ, flexVector<T> indexVal)
	{
		//clear matrix
		clear();

		int numberListElements = indexI.size();

		//initialize vecvector
		flexVector<int> emptyBucket(0, 0);
		flexVector < flexVector<int> > buckets(numRows, emptyBucket);

		//add elements to buckets
		for (int indexInput = 0; indexInput < numberListElements; indexInput++)
		{
			int bucketIndex = indexI[indexInput];
			buckets[bucketIndex].push_back(indexInput);
		}

		//go trough all rows:
		for (int indexRow = 0; indexRow < numRows; indexRow++)
		{
			int numElements = 0;

			//go through bucket
			for (int indexBucket = 0; indexBucket < buckets[indexRow].size(); indexBucket++)
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

	//inserts new matrix element val at position [i][j] this is SLOW!
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

		//rowToIndexList.print();
		//indexList.print();
	}

	//clear current matrix
	void clear()
	{
		indexList.clear();
		valueList.clear();

		//set all rowToIndexList to 0
		for (int i = 0; i < numRows+1; i++)
		{
			rowToIndexList[i] = 0;
		}
	}

	T getMaxRowSumAbs()
	{
		T maxSum = static_cast<T>(0);

		for (int i = 0; i < numRows; ++i)
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

	void printRow(int i)
	{
		//printf("|");
		for (int index = rowToIndexList[i]; index < rowToIndexList[i+1]; ++index)
		{
			printf("(%d,%d,%f)|", i, indexList[index], valueList[index]);
		}

		printf("\n");

	}
	void printMatrix()
	{
		for (int i = 0; i < numRows; i++)
		{
			//printf("Row: %d|", i);
			printRow(i);
			//printf("\n");
		}
	}

	//transpose current matrix
	void transpose()
	{
		flexVector<int> tmpindexI(0, 0);
		flexVector<int> tmpindexJ(0, 0);
		flexVector<T> tmpindexVal(0, static_cast<T>(0));

		for (int i = 0; i < numRows; i++)
		{
			for (int index = rowToIndexList[i]; index < rowToIndexList[i + 1]; ++index)
			{
				tmpindexI.push_back(indexList[index]);
				tmpindexJ.push_back(i);
				tmpindexVal.push_back(valueList[index]);
			}
		}

		clear();

		blockInsert(tmpindexI, tmpindexJ, tmpindexVal);

		int numRowsTmp = numRows;
		numRows = numCols;
		numCols = numRowsTmp;

		
		/*flexVector<flexVector<int>> indexListTmp = indexList;
		flexVector<flexVector<T>> valueListTmp = valueList;

		indexList = indexListT;
		valueList = valueListT;

		indexListT = indexListTmp;
		valueListTmp = valueListTmp;


		return;

		flexVector<flexVector<int>> indexListOld;
		flexVector<flexVector<T>> valueListOld;
		indexListOld = indexList;
		valueListOld = valueList;

		//clear old content
		clear();

		int numRowsTmp = numRows;
		numRows = numCols;
		numCols = numRowsTmp;

		indexList.resize(numRows);
		valueList.resize(numRows);


		//iterate through old rows
		for (int i = 0; i < indexListOld.size(); i++)
		{
			for (int j = 0; j < indexListOld[i].size(); ++j)
			{
				//insert element at [j,i]
				insertElement(indexListOld[i][j], i, valueListOld[i][j]);
				//cout << indList->size() << endl;

				//indexList[indList->at(j)].push_back(i);
				//valueList[indList->at(j)].push_back(valList->at(j));
			}
		}

		//cout << "Fin";*/
	}
};

#endif