#ifndef flexMatrix_H
#define flexMatrix_H

#include "flexVector.h"

template < typename T >
class flexMatrix
{
private:
	flexVector<flexVector<int>> indexList;
	flexVector<flexVector<T>> valueList;
	int numCols;
	int numRows;
	
public:
	flexMatrix(void) : indexList(), valueList()
	{
		numCols = 0;
		numRows = 0;
	}

	flexMatrix(const int  _numRows,const int  _numCols) : indexList(), valueList()
	{
		numCols = _numCols;
		numRows = _numRows;

		indexList.resize(numRows);
		valueList.resize(numRows);
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
			output[i] = rowMulti(input,indexList[i], valueList[i]);
		}
	}

	void timesPlus(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			output[i] += rowMulti(input, indexList[i], valueList[i]);
		}
	}

	void timesMinus(const flexVector<T> &input, flexVector<T> &output)
	{
		#pragma omp parallel for
		for (int i = 0; i < numRows; i++)
		{
			output[i] -= rowMulti(input, indexList[i], valueList[i]);
		}
	}

	T rowMulti(const int i, const flexVector<T> &input)
	{
		// initialize result
		T result = static_cast<T>(0);

		int elements = indexList[i].size();

		if (elements > 0)
		{
			for (int j = 0; j < elements; j++)
			{
				result += input[indexList[i][j]] * valueList[i][j];
			}
		}
		return result;
	}

	T rowMulti(const flexVector<T> &input,const flexVector<int> &inputI, const flexVector<T> &inputV)
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
	}

	//inserts new matrix element val at position [i][j]
	void insertElement(int i, int j, T val, bool checkIfExists)
	{
		if (!checkIfExists)
		{
			insertElement(i,j,val);
		}
		else
		{

		}
	}

	//inserts new matrix element val at position [i][j]
	void insertElement(int i, int j, T val)
	{
		indexList[i].push_back(j);
		valueList[i].push_back(val);
	}

	//clear current matrix
	void clear()
	{
		//go through all rows and empty index and value list of corresponding row
		for (int i = 0; i < numRows; i++)
		{
			indexList[i].clear();
			valueList[i].clear();
		}
	}

	T getMaxRowSumAbs()
	{
		T maxSum = static_cast<T>(0);

		for (int i = 0; i < numRows; ++i)
		{
			flexVector<T> tmpRowSum;
			tmpRowSum = valueList[i];
			tmpRowSum.abs();

			T tmpSum = tmpRowSum.sum();
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
		for (int j = 0; j < valueList[i].size(); ++j)
		{
			printf("(%d,%d,%f)|",i,indexList[i][j], valueList[i][j]);
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
		/*flexVector<flexVector<int>> indexListTmp = indexList;
		flexVector<flexVector<T>> valueListTmp = valueList;

		indexList = indexListT;
		valueList = valueListT;

		indexListT = indexListTmp;
		valueListTmp = valueListTmp;


		return;*/
		
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

		//cout << "Fin";
	}
};

#endif