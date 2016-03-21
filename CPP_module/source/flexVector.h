#ifndef flexVector_H
#define flexVector_H

template<class T>
struct IsPrimitiveType {
	enum { VALUE = 0 };
};

template<>
struct IsPrimitiveType<int> {
	enum { VALUE = 1 };
};
template<>
struct IsPrimitiveType<float> {
	enum { VALUE = 1 };
};

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <new>
#include <algorithm>
#include <iostream>


template < typename T >
class flexVector
{

public:
	typedef typename int size_type;
	int _size;
	int _capacity;
	T* _data;

	static T* allocate(const int size)
	{
		//return static_cast<T*>(malloc(sizeof(T) * size));
		return static_cast<T*>(malloc(sizeof(T)* size));

	}

	static void copyRange(T* begin, T* end, T* dest)
	{
		while (begin != end)
		{
			new((void*)dest) T(*begin);
			++begin;
			++dest;
		}
	}
	static void deleteRange(T* begin, T* end)
	{
		if (IsPrimitiveType<T>::VALUE == 1)
		{
			//printf("Primitive \n");
			return;
		}
			

		while (begin != end)
		{
			begin->~T();
			++begin;
		}
	}

	void reallocate(const int newCapacity)
	{
		T* newData = allocate(newCapacity);
		copyRange(_data, _data + _size, newData);
		deleteRange(_data, _data + _size);

		//if (IsPrimitiveType<T>::VALUE != 1)
		//{
			free(_data);
		//}
		
		_data = newData;
		_capacity = newCapacity;
	}

	static void constructRange(T* begin, T* end, const T& fillWith)
	{
		while (begin != end)
		{
			new((void*)begin) T(fillWith);
			++begin;
		}
	}


	static void constructRange(T* begin, T* end)
	{
		if (IsPrimitiveType<T>::VALUE == 1)
			return;

		while (begin != end)
		{
			new((void*)begin) T();
			++begin;
		}
	}


	void resize(const size_type newSize)
	{
		if (newSize <= _size)
		{
			deleteRange(_data + newSize, _data + _size);
			_size = newSize;
			return;
		}
		if (newSize <= _capacity)
		{
			constructRange(_data + _size, _data + newSize);
			_size = newSize;
			return;
		}
		size_type newCapacity = newSize;
		if (newCapacity < _size * 2)
		{
			newCapacity = _size * 2;
		}
		reallocate(newCapacity);
		constructRange(_data + _size, _data + newSize);
		_size = newSize;
	}

	void resize(const size_type newSize, const T& fillWith)
	{
		/*if (newSize <= _size)
		{
			deleteRange(_data + newSize, _data + _size);
			_size = newSize;
			return;
		}
		if (newSize <= _capacity)
		{
			constructRange(_data + _size, _data + newSize, fillWith);
			_size = newSize;
			return;
		}*/
		size_type newCapacity = newSize;
		if (newCapacity < _size * 2)
		{
			newCapacity = _size * 2;
		}
		if (newCapacity < 2)
		{
			newCapacity = 2;
		}

		reallocate(newCapacity);
		constructRange(_data + _size, _data + newSize, fillWith);
		for (size_type i = _size; i < newSize; i++)
		{
			::new((void*)(_data + i)) T(fillWith);
		}
		_size = newSize;
	}

	flexVector(void) : _size(), _capacity(), _data()
	{
		_size = 0;
		_capacity = 1;
		_data = allocate(_capacity);

		//printf("Void constructor \n");
	}

	flexVector(const flexVector &classToCopy)
	{
		//printf("Copy constructor \n");
		_capacity = classToCopy._size;
		_data = allocate(_capacity);
		_size = classToCopy._size;

		if (IsPrimitiveType<T>::VALUE == 1)
		{
			// if data type is primitive, just copy the array
			std::copy(classToCopy._data, classToCopy._data + _size, _data);
		}
		else
		{
			copyRange(classToCopy._data, classToCopy._data + _size, _data);
		}
	}

	flexVector(const size_type size, const T& fillWith)
	{
		//printf("Standard constructor \n");

		if (size < 2)
		{
			_capacity = 2;
		}
		else
		{
			_capacity = size;
		}
		
		_size = size;
		_data = allocate(_capacity);

		constructRange(_data, _data + _size, fillWith);
	}

	~flexVector()
	{
		//printf("Destructor \n");
		//if (!IsPrimitiveType<T>::VALUE)
		//{
			deleteRange(_data, _data + _size);
		//}
		
		free(_data);
	}

	size_type size()
	{
		return _size;
	}

	T sum()
	{
		T result = static_cast<T>(0);

		for (size_type i = 0; i < _size; ++i)
		{
			result += _data[i];
		}

		return result;
	}

	T product()
	{
		T result = static_cast<T>(1);

		for (size_type i = 0; i < _size; ++i)
		{
			result *= _data[i];
		}

		return result;
	}

	void print()
	{
		printf("Vector contains: \n");
		for (size_type i = 0; i < _size; ++i)
		{
			printf("_data[%d]: %f\n", i, static_cast<float>(_data[i]));
		}
	}

	void abs()
	{
		T zero = static_cast<T>(0);

		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			if (_data[i] < zero)
			{
				_data[i] = -_data[i];
			}
		}
	}

	void pow2()
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] *= _data[i];
		}
	}

	void sqrt()
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] = std::sqrt(_data[i]);
		}
	}

	void clear()
	{
		resize(0);
	}

	void push_back(const T& value)
	{
		if (_size != _capacity)
		{
			new((void*)(_data + _size)) T(value);
			++_size;
			return;
		}
		size_type newCapacity = _capacity ? _capacity * 2 : 1;
		T* newData = allocate(newCapacity);
		copyRange(_data, _data + _size, newData);
		new((void*)(newData + _size)) T(value);
		deleteRange(_data, _data + _size);
		free(_data);
		_data = newData;
		_capacity = newCapacity;
		++_size;
	}

	T &operator[](const size_type index)
	{
		//printf("Brackets \n");
		return _data[index];
	}

	const T &operator[](const size_type index) const
	{
		//printf("Const brackets \n");
		return _data[index];
	}

	T* begin() const
	{
		return _data;
	}

	T* end() const
	{
		return _data + _size;
	}

	//projection:

	//L^\infty projection
	void project_minMax(T alpha)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			if (_data[i] < -alpha)
			{
				_data[i] = -alpha;
			}
			else if (_data[i] > alpha)
			{
				_data[i] = alpha;
			}
		}
	}

	//assignment operator
	flexVector<T> &operator=(const flexVector< T > &input)
	{
		if (input._size != _size)
		{
			resize(input._size);
		}

		if (IsPrimitiveType<T>::VALUE == 1)
		{
			// if data type is primitive, just copy the array
			std::copy(input._data, input._data + _size, _data);
		}
		else
		{
			copyRange(input._data, input._data + _size, _data);
		}
		

		return *this;
	}

	flexVector<T> &operator+=(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] += input[i];
		}

		return *this;
	}

	flexVector<T> &operator+(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] += input[i];
		}

		return *this;
	}

	flexVector<T> &operator-=(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] -= input[i];
		}

		return *this;
	}

	flexVector<T> &operator-(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] -= input[i];
		}

		return *this;
	}

	flexVector<T> &operator*(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] *= input[i];
		}

		return *this;
	}

	//sets all elements in vector to scalar input
	void scalarEquals(T input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] = input;
		}
	}

	//divides all elements in vector scalar input
	void scalarDivide(T input)
	{
		scalarMult(static_cast<T>(1) / input);
	}

	//calculates the maximum of every element and scalar input
	void scalarMult(T input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] *= input;
		}
	}

	//calculates the maximum of every element and scalar input
	void scalarMax(T input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] = max(_data[i],input);
		}
	}

	//divides all elements in vector by corresponding element in input
	void vectorDivide(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] = _data[i]/input[i];
		}
	}

	flexVector<T> &operator+=(const T &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] += input;
		}

		return *this;
	}

	flexVector<T> &operator-=(const T &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] -= input;
		}

		return *this;
	}



	flexVector<T> &operator*=(const flexVector< T > &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] *= input[i];
		}

		return *this;
	}

	flexVector<T> &operator*(const T &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] *= input;
		}

		return *this;
	}

	flexVector<T> &operator*=(const T &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] *= input;
		}

		return *this;
	}

	flexVector<T> &operator/=(const T &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] /= input;
		}

		return *this;
	}

	flexVector<T> &operator/(const T &input)
	{
		#pragma omp parallel for
		for (size_type i = 0; i < _size; ++i)
		{
			_data[i] /= input;
		}

		return *this;
	}

	T at(int i)
	{
		return _data[i];
	}

};

#endif