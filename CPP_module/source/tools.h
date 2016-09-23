#ifndef flexTools_H
#define flexTools_H

#include <vector>
#include <numeric>
#include <string>
#include <functional>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>

#define VERBOSE 0

#if __CUDACC__
	#include <thrust/device_vector.h> 
	#include <thrust/transform.h> 
	#include <thrust/sequence.h> 
	#include <thrust/copy.h> 
	#include <thrust/fill.h> 
	#include <thrust/replace.h> 
	#include <thrust/functional.h>
	#include <thrust/iterator/zip_iterator.h>
	#include <thrust/tuple.h> 
	#include <thrust/for_each.h>
	#include <thrust/transform_reduce.h>
	#include <thrust/extrema.h>
#endif


using std::cerr;
using std::cout;
using std::endl;
using std::string;


static const int SIGN_PLUS = 0;
static const int SIGN_MINUS = 1;
static const int SIGN_EQUALS = 2;

// Could be any number, but the whole array should fit into shared memory 
#define CONST_ARRAY_SIZE 512
#define BLOCK_SIZE (64)


enum mySign
{
	PLUS,
	MINUS,
	EQUALS
};

enum prox
{
	primalEmptyProx,
	dualL1AnisoProx,
	dualL1IsoProx,
	dualL2Prox,
	dualFrobeniusProx,
	dualHuberProx,
	dualL2DataProx,
	dualL1DataProx,
	dualKLDataProx
};

enum linOp
{
	linearOp,
	diagonalOp,
	gradientOp,
	identityOp,
	matrixOp,
	matrixGPUOp,
	zeroOp
};



template < typename T >
T myAbs(T x)
{
	return x > 0 ? x : -x;
}

template < typename T >
T myMin(T a, T b)
{
	return a > b ? b : a;
}

template < typename T >
T myMax(T a, T b)
{
	return a < b ? b : a;
}

template < typename T >
float myPow2(float x)
{
	return x * x;
}

float pow2(float x)
{
	return x * x;
}

template < typename T >
void calculateXYError(const std::vector<T> &x, const std::vector<T> &xOld, std::vector<T> &xError, T tau)
{
	std::transform(x.begin(), x.end(), xOld.begin(), xError.begin(), [tau](T x, T y) {return (x-y)/tau; });
}

template < typename T >
T vectorProduct(const std::vector<T> &v)
{
	return std::accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
}



template < typename T >
T vectorSum(const std::vector<T> &v)
{
	return std::accumulate(v.begin(), v.end(), (T)0);
}

float vectorMax(std::vector<float> &v)
{
	return *std::max_element(v.begin(), v.end());
}

template < typename T >
void vectorScalarProduct(std::vector<T> &v,T scalarValue)
{
	std::transform(v.begin(), v.end(), v.begin(), [scalarValue](T x) {return scalarValue*x;});
}

void vectorScalarSet(std::vector<float> &v, const float scalarValue)
{
	std::fill(v.begin(), v.end(), scalarValue);
}

void vectorPlus(std::vector<float> &v1, std::vector<float> &v2)
{
	std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::plus<float>());
}

void vectorMinus(std::vector<float> &v1, std::vector<float> &v2)
{
	std::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), std::minus<float>());
}

template < typename T >
void vectorAbs(std::vector<T> &v)
{
	std::transform(v.begin(), v.end(), v.begin(), [](T x) { return std::abs(x); });
}

template < typename T >
void vectorProjectL1Aniso(std::vector<T> &vTilde, std::vector<T> &v, T alpha)
{
	std::transform(vTilde.begin(), vTilde.end(), v.begin(), [alpha](T x) { return std::min(alpha, std::max(-alpha, x)); });
}

template < typename T >
void vectorProjectL1Iso2D(std::vector<T> &yTmp, std::vector<T> &y1Tilde, std::vector<T> &y2Tilde, std::vector<T> &y1, std::vector<T> &y2, T alpha)
{
	//square the first argument and add squared of second
	std::transform(y1Tilde.begin(), y1Tilde.end(), y2Tilde.begin(), yTmp.begin(), [](T x, T y) { return x*x + y*y; });

	//calculate norm projection
	std::transform(yTmp.begin(), yTmp.end(), yTmp.begin(), [alpha](T x) { return std::max((T)1, std::sqrt(x) / alpha); });

	//project onto norm
	std::transform(y1Tilde.begin(), y1Tilde.end(), yTmp.begin(), y1.begin(), [](T x, T y) { return x / y; });
	std::transform(y2Tilde.begin(), y2Tilde.end(), yTmp.begin(), y2.begin(), [](T x, T y) { return x / y; });
}

template < typename T >
struct HuberCalcTilder2D
{
	HuberCalcTilder2D(T sigma, T alpha) : factor((T)1.0 + sigma * (T)0.01 / alpha){};

	T operator()(T x, T y) const { return x/factor * x/factor + y/factor * y/factor; }
private:
	const T factor;
};

template < typename T >
void vectorProjectHuber2D(std::vector<T> &yTmp, std::vector<T> &y1Tilde, std::vector<T> &y2Tilde, std::vector<T> &y1, std::vector<T> &y2, T alpha, T sigma)
{
	//square the first argument and add squared of second
	std::transform(y1Tilde.begin(), y1Tilde.end(), y2Tilde.begin(), yTmp.begin(),HuberCalcTilder2D<T>(sigma,alpha) );

	//calculate norm projection
	std::transform(yTmp.begin(), yTmp.end(), yTmp.begin(), [alpha](T x) { return std::max((T)1, std::sqrt(x) / alpha); });

	//project onto norm
	std::transform(y1Tilde.begin(), y1Tilde.end(), yTmp.begin(), y1.begin(), [](T x, T y) { return x / y; });
	std::transform(y2Tilde.begin(), y2Tilde.end(), yTmp.begin(), y2.begin(), [](T x, T y) { return x / y; });
}

/*CPU prox for L2 regularizer*/
template < typename T >
void vectorProjectL2(std::vector<T> &yTilde, std::vector<T> &y, T alpha, T sigma)
{
	T factor = alpha / (sigma + alpha);
	std::transform(yTilde.begin(), yTilde.end(), y.begin(), [factor](T x) { return x * factor; });
}

template < typename T >
void vectorProjectL2data(std::vector<T> &y, std::vector<T> &yTilde, std::vector<T> &f, T alpha, T sigma)
{
	T factor = alpha / (sigma + alpha);
	std::transform(yTilde.begin(), yTilde.end(), f.begin(), y.begin(), [factor, sigma](T x, T y) { return factor*(x - sigma*y); });
}

/*CPU prox for KL data term (functor)*/
//main.y{dualNumber} = min(1,0.5*(1 + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}-1).^2 + 4*main.params.sigma{dualNumber}*obj.f(:) )));
template < typename T >
struct KLprojectionData
{
	KLprojectionData(T sigma) : sigma((T)4.0 * sigma){};

	T operator()(T x, T y) const { return (T)0.5 * ((T)1.0 + x - std::sqrt(myPow2<T>(x - (T)1.0) + sigma * y)); }
private:
	const T sigma;
};

/*GPU prox for KL data term */
template < typename T >
void vectorProjectKLdata(std::vector<T> &y, std::vector<T> &yTilde, std::vector<T> &f, T alpha, T sigma)
{
	std::transform(yTilde.begin(), yTilde.end(), f.begin(), y.begin(), KLprojectionData<T>(sigma));
}

/*CPU prox for L1 data term */
template < typename T >
void vectorProjectL1data(std::vector<T> &y, std::vector<T> &yTilde, std::vector<T> &f, T alpha, T sigma)
{
	T factor = alpha * sigma;
	std::transform(yTilde.begin(), yTilde.end(), f.begin(), y.begin(), [alpha,factor](T x, T y) { return std::min(alpha, std::max(-alpha, x - y * factor));});
}

/*CPU prox for Frobenius regularizer*/
template < typename T >
void vectorProjectFrobenius2D(std::vector<T> &yTmp, std::vector<T> &y1Tilde, std::vector<T> &y2Tilde, std::vector<T> &y1, std::vector<T> &y2, T alpha)
{
	//square the first argument and add squared of second
	std::transform(y1Tilde.begin(), y1Tilde.end(), y2Tilde.begin(), yTmp.begin(), [](T x,T y) { return x*x + y*y; });

	//sum up and take sqrt
	T norm = std::max((T)1, std::sqrt(std::accumulate(yTmp.begin(), yTmp.end(), (T)0)) / alpha);

	//project onto norm
	std::transform(y1Tilde.begin(), y1Tilde.end(), y1.begin(), [norm](T x) { return x / norm; });
	std::transform(y2Tilde.begin(), y2Tilde.end(), y2.begin(), [norm](T x) { return x / norm; });
}

template < typename T >
void doOverrelaxation(std::vector<T> &x, std::vector<T> &xOld, std::vector<T> &xBar)
{
	std::transform(x.begin(), x.end(), xOld.begin(), xBar.begin(), [](T x, T y) { return x + x - y; });
}

template < typename T >
void vectorPow2(std::vector<T> &v)
{
	std::transform(v.begin(), v.end(), v.begin(), [](T x) { return x *x ; });
}

template < typename T >
void vectorAddVectorTimesVector(std::vector<T> &result, const std::vector<T> &v1, const  std::vector<T> &v2, const int signRule)
{
	switch (signRule)
	{
		case SIGN_PLUS:
		{
			int numElements = result.size();
			for (int i = 0; i < numElements; ++i)
			{
				result[i] += v1[i] * v2[i];
			}
			break;
		}
		case SIGN_MINUS:
		{
			int numElements = result.size();
			for (int i = 0; i < numElements; ++i)
			{
				result[i] -= v1[i] * v2[i];
			}
			break;
		}
		case SIGN_EQUALS:
		{
			int numElements = result.size();
			for (int i = 0; i < numElements; ++i)
			{
				result[i] = v1[i] * v2[i];
			}
			break;
		}
	}
}

class Timer
{
public:
	Timer() : t_begin(omp_get_wtime())
	{
	};

	void reset()
	{
#if __CUDACC__
			cudaDeviceSynchronize();
		#endif
		t_begin = omp_get_wtime();
	}

	void end()
	{
#if __CUDACC__
			cudaDeviceSynchronize();
		#endif
		t_end = omp_get_wtime();
	}

	double elapsed() const
	{
		return t_end - t_begin;
	}

private:
	double t_begin;
	double t_end;
};





#if __CUDACC__

	template < typename T >
	void calculateXYError(thrust::device_vector<T> &x, thrust::device_vector<T> &xOld, thrust::device_vector<T> &xError, T tau)
	{
		thrust::transform(x.begin(), x.end(), xOld.begin(), xError.begin(), (thrust::placeholders::_1 - thrust::placeholders::_2) / tau);
	}

	template < typename T >
	__host__ __device__
	T myPow2GPU(T x)
	{
		return x * x;
	}

	template < typename T >
	struct myAbsGPU
	{
		__host__ __device__
		myAbsGPU() {}

		__host__ __device__
		T operator()(T x) const { return std::abs(x); }
	};

	template < typename T >
	void vectorAbs(thrust::device_vector<T> &v)
	{
		thrust::transform(v.begin(), v.end(), v.begin(), myAbsGPU<T>());
	}

	template < typename T >
	__host__ __device__
	T myMinGPU(T a, T b)
	{
		return a > b ? b : a;
	}

	__device__ float myMinGPUf(float a,float b)
	{
		return a > b ? b : a;
	}

	__device__ float myMaxGPUf(float a,float b)
	{
		return a > b ? a : b;
	}

	template < typename T >
	__host__ __device__
	T myMaxGPU(T a, T b)
	{
		return a < b ? b : a;
	}
	
	template < typename T >
	T vectorSum(thrust::device_vector<T> &v)
	{
		return thrust::reduce(v.begin(), v.end(), (T)0, thrust::plus<T>());
	}

	//binary functor for overrelax
	template < typename T >
	void doOverrelaxation(thrust::device_vector<T> &x, thrust::device_vector<T> &xOld, thrust::device_vector<T> &xBar)
	{
		thrust::transform(x.begin(), x.end(), xOld.begin(), xBar.begin(), (T)2 * thrust::placeholders::_1 - thrust::placeholders::_2 );
	}

	template < typename T >
	struct L1projectionGPU
	{
		__host__ __device__
		L1projectionGPU(T alpha) : alpha(alpha) 
		{}

		__host__ __device__
		T operator()(T x) const { return myMinGPU<T>(alpha, myMaxGPU<T>(-alpha, x)); }
	private:
		const T alpha;
	};

	template < typename T >
	void vectorProjectL1Aniso(thrust::device_vector<T> &vTilde, thrust::device_vector<T> &v, const T scalarValue)
	{
		thrust::transform(vTilde.begin(), vTilde.end(), v.begin(), L1projectionGPU<T>(scalarValue));
	}


	template<typename T>
	struct L1IsoProjectionGPU2D
	{
		__host__ __device__
		L1IsoProjectionGPU2D(T alpha) : alpha(alpha){}

		template <typename Tuple>
		__host__ __device__
			void operator()(Tuple t)
		{
			T norm = myMaxGPU<T>((T)1, std::sqrt(thrust::get<0>(t) * thrust::get<0>(t) +thrust::get<1>(t) * thrust::get<1>(t)) / this->alpha);

			thrust::get<2>(t) = thrust::get<0>(t) / norm;
			thrust::get<3>(t) = thrust::get<1>(t) / norm;
		}

		T alpha;
	};

	template < typename T >
	void vectorProjectL1Iso2D(thrust::device_vector<T> &yTmp, thrust::device_vector<T> &y1Tilde, thrust::device_vector<T> &y2Tilde, thrust::device_vector<T> &y1, thrust::device_vector<T> &y2, T alpha)
	{
		thrust::for_each(
			thrust::make_zip_iterator(thrust::make_tuple(y1Tilde.begin(), y2Tilde.begin(), y1.begin(), y2.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(y1Tilde.end(), y2Tilde.end(), y1.end(), y2.end())),
			L1IsoProjectionGPU2D<T>(alpha));
	}
	
	template < typename T >
	void vectorProjectHuber2D(thrust::device_vector<T> &yTmp, thrust::device_vector<T> &y1Tilde, thrust::device_vector<T> &y2Tilde, thrust::device_vector<T> &y1, thrust::device_vector<T> &y2, T alpha, T sigma)
	{
		/*
		//square the first argument and add squared of second
		std::transform(y1Tilde.begin(), y1Tilde.end(), y2Tilde.begin(), yTmp.begin(),HuberCalcTilder2D<T>(sigma,alpha) );

		//calculate norm projection
		std::transform(yTmp.begin(), yTmp.end(), yTmp.begin(), [alpha](T x) { return std::max((T)1, std::sqrt(x) / alpha); });

		//project onto norm
		std::transform(y1Tilde.begin(), y1Tilde.end(), yTmp.begin(), y1.begin(), [](T x, T y) { return x / y; });
		std::transform(y2Tilde.begin(), y2Tilde.end(), yTmp.begin(), y2.begin(), [](T x, T y) { return x / y; });*/
	}

	/*GPU prox for L2 regularizer*/
	template < typename T >
	void vectorProjectL2(thrust::device_vector<T> &yTilde, thrust::device_vector<T> &y, T alpha, T sigma)
	{
		T factor = alpha / (sigma + alpha);
		thrust::transform(yTilde.begin(), yTilde.end(), y.begin(), factor * thrust::placeholders::_1);
	}

	/*GPU prox for Frobenius regularizer*/
	template < typename T >
	void vectorProjectFrobenius2D(thrust::device_vector<T> &yTmp, thrust::device_vector<T> &y1Tilde, thrust::device_vector<T> &y2Tilde, thrust::device_vector<T> &y1, thrust::device_vector<T> &y2, T alpha)
	{
		//square the first argument and add squared of second
		thrust::transform(y1Tilde.begin(), y1Tilde.end(), y2Tilde.begin(), yTmp.begin(), thrust::placeholders::_1*thrust::placeholders::_1 + thrust::placeholders::_2 * thrust::placeholders::_2);

		//sum up and take sqrt
		T norm = std::max(static_cast<T>(1), std::sqrt(thrust::reduce(yTmp.begin(), yTmp.end(), (T)0, thrust::plus<T>())) / alpha);

		//project onto norm
		thrust::transform(y1Tilde.begin(), y1Tilde.end(), y1.begin(), thrust::placeholders::_1 / norm);
		thrust::transform(y2Tilde.begin(), y2Tilde.end(), y2.begin(), thrust::placeholders::_1 / norm);
	}

	/*GPU prox for L2 data term*/
	template < typename T >
	void vectorProjectL2data(thrust::device_vector<T> &y, thrust::device_vector<T> &yTilde, thrust::device_vector<T> &f, T alpha, T sigma)
	{
		T factor = alpha / (sigma + alpha);
		thrust::transform(yTilde.begin(), yTilde.end(), f.begin(), y.begin(), factor * (thrust::placeholders::_1 - sigma*thrust::placeholders::_2));
	}


	/*GPU prox for L1 data term (functor)*/
	template < typename T >
	struct L1projectionDataGPU
	{
		__host__ __device__
		L1projectionDataGPU(T alpha, T sigma) : alpha(alpha), factor(alpha * sigma){};

		__host__ __device__
			T operator()(T x, T y) const { return myMinGPU<T>(alpha, myMaxGPU<T>(-alpha, x - y * factor)); }
	private:
		const T factor;
		const T alpha;
	};

	/*GPU prox for L1 data term */
	template < typename T >
	void vectorProjectL1data(thrust::device_vector<T> &y, thrust::device_vector<T> &yTilde, thrust::device_vector<T> &f, T alpha, T sigma)
	{
		thrust::transform(yTilde.begin(), yTilde.end(), f.begin(), y.begin(), L1projectionDataGPU<T>(alpha,sigma));
	}


	/*GPU prox for KL data term (functor)*/
	//main.y{dualNumber} = min(1,0.5*(1 + main.yTilde{dualNumber} - sqrt( (main.yTilde{dualNumber}-1).^2 + 4*main.params.sigma{dualNumber}*obj.f(:) )));
	template < typename T >
	struct KLprojectionDataGPU
	{
		__host__ __device__ 
		KLprojectionDataGPU(T sigma) : sigma((T)4.0 * sigma){};

		__host__ __device__
		T operator()(T x, T y) const { return (T)0.5 * ((T)1.0 + x - std::sqrt( myPow2GPU(x - (T)1.0) + sigma * y)); }
	private:
		const T sigma;
	};

	/*GPU prox for KL data term */
	template < typename T >
	void vectorProjectKLdata(thrust::device_vector<T> &y, thrust::device_vector<T> &yTilde, thrust::device_vector<T> &f, T alpha, T sigma)
	{
		thrust::transform(yTilde.begin(), yTilde.end(), f.begin(), y.begin(), KLprojectionDataGPU<T>(sigma));
	}

	/*sets all elements in a vector to scalarValue*/
	template < typename T >
	void vectorScalarSet(thrust::device_vector<T> &v, T scalarValue)
	{
		thrust::fill(v.begin(), v.end(), scalarValue);
	}

	template < typename T >
	void vectorScalarProduct(thrust::device_vector<T> &v, const T scalarValue)
	{
		thrust::transform(v.begin(), v.end(), v.begin(), scalarValue * thrust::placeholders::_1);
	}

	template < typename T >
	void vectorMinus(thrust::device_vector<T> &v1, thrust::device_vector<T> &v2)
	{
		thrust::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), thrust::minus<T>());
	}

	template < typename T >
	void vectorAddSquared(thrust::device_vector<T> &v1, thrust::device_vector<T> &v2)
	{
		thrust::transform(v1.begin(), v1.end(), v2.begin(), v1.begin(), thrust::placeholders::_1 + thrust::placeholders::_2*thrust::placeholders::_2);
	}

	struct vectorAddVectorTimesVectorGPU
	{
		__host__ __device__ 
		vectorAddVectorTimesVectorGPU(const int signRule) : signRule(signRule){}

		template <typename Tuple>
		__host__ __device__
		void operator()(Tuple t)
		{
			switch (signRule)
			{
				case SIGN_PLUS:
				{
					thrust::get<0>(t) += thrust::get<1>(t) * thrust::get<2>(t);
					break;
				}
				case SIGN_MINUS:
				{
					thrust::get<0>(t) -= thrust::get<1>(t) * thrust::get<2>(t);
					break;
				}
				case SIGN_EQUALS:
				{
					thrust::get<0>(t) = thrust::get<1>(t) * thrust::get<2>(t);
					break;
				}
			}
		}

		const int signRule;
	};

	template < typename T >
	void vectorAddVectorTimesVector(thrust::device_vector<T> &result, const thrust::device_vector<T> &v1, const  thrust::device_vector<T> &v2, const int signRule)
	{
		thrust::for_each(
			thrust::make_zip_iterator(thrust::make_tuple(result.begin(), v1.begin(), v2.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(result.end(), v1.end(), v2.end())),
			vectorAddVectorTimesVectorGPU(signRule));
	}


	template < typename T >
	__host__ __device__
	T sqrtGPU(T x)
	{
		return std::sqrt(x);
	}

	template < typename T >
	void vectorSqrt(thrust::device_vector<T> &v1)
	{
		thrust::transform(v1.begin(), v1.end(), v1.begin(), sqrtGPU<T>());
	}

	template < typename T >
	float vectorMax(thrust::device_vector<T> &v)
	{
		return *thrust::max_element(v.begin(), v.end());
	}

#endif


	//kernels for CUDA operators

#if __CUDACC__

/*template<typename T>
__global__ void dxp2dCUDA(T* output, const T* input, int w, int h, const int signRule)
{
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;

	if (x >= w || y >= h)
		return;

	const int tmpIndex = x + y * w;

	if (x < w - 1)
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += input[tmpIndex + 1] - input[tmpIndex];
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= input[tmpIndex + 1] - input[tmpIndex];
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = input[tmpIndex + 1] - input[tmpIndex];
				break;
			}
		}
	}
}

template<typename T>
__global__ void dxp2dTransposedCUDA(T* output, const T* input, int w, int h, const int signRule)
{
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;

	if (x >= w || y >= h)
		return;

	const int tmpIndex = x + y * w;

	if (x < w - 1 && x > 0)
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += -(input[tmpIndex] - input[tmpIndex - 1]);
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= -(input[tmpIndex] - input[tmpIndex - 1]);
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = -(input[tmpIndex] - input[tmpIndex - 1]);
				break;
			}
		}
	}
	else if (x == 0)
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += -(input[tmpIndex]);
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= -(input[tmpIndex]);
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = -(input[tmpIndex]);
				break;
			}
		}
	}
	else
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += (input[tmpIndex - 1]);
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= (input[tmpIndex - 1]);
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = (input[tmpIndex - 1]);
				break;
			}
		}
	}
}

template<typename T>
__global__ void dyp2dCUDA(T* output, const T* input, int w, int h, const int signRule)
{
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;

	if (x >= w || y >= h)
		return;

	const int tmpIndex = x + y * w;

	if (y < h - 1)
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += input[tmpIndex + w] - input[tmpIndex];
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= input[tmpIndex + w] - input[tmpIndex];
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = input[tmpIndex + w] - input[tmpIndex];
				break;
			}
		}
	}
}

template<typename T>
__global__ void dyp2dTransposedCUDA(T* output, const T* input, int w, int h, const int signRule)
{
	int x = threadIdx.x + blockIdx.x * blockDim.x;
	int y = threadIdx.y + blockIdx.y * blockDim.y;

	if (x >= w || y >= h)
		return;

	const int tmpIndex = x + y * w;

	if (y < h - 1 && y > 0)
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += -(input[tmpIndex] - input[tmpIndex - w]);
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= -(input[tmpIndex] - input[tmpIndex - w]);
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = -(input[tmpIndex] - input[tmpIndex - w]);
				break;
			}
		}
	}
	else if (y == 0)
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += -(input[tmpIndex]);
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= -(input[tmpIndex]);
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = -(input[tmpIndex]);
				break;
			}
		}
	}
	else
	{
		switch (signRule)
		{
			case SIGN_PLUS:
			{
				output[tmpIndex] += (input[tmpIndex - w]);
				break;
			}
			case SIGN_MINUS:
			{
				output[tmpIndex] -= (input[tmpIndex - w]);
				break;
			}
			case SIGN_EQUALS:
			{
				output[tmpIndex] = (input[tmpIndex - w]);
				break;
			}
		}
	}
}*/


// cuda error checking
std::string prev_file = "";
int prev_line = 0;
void cuda_check(std::string file, int line)
{
	cudaError_t e = cudaGetLastError();
	if (e != cudaSuccess)
	{
		std::ofstream out("output.txt");
		out << endl << file << ", line " << line << ": " << cudaGetErrorString(e) << " (" << e << ")" << endl;
		if (prev_line>0) out << "Previous CUDA call:" << endl << prev_file << ", line " << prev_line << endl;
		out.close();

		cout << endl << file << ", line " << line << ": " << cudaGetErrorString(e) << " (" << e << ")" << endl;
		if (prev_line>0) cout << "Previous CUDA call:" << endl << prev_file << ", line " << prev_line << endl;
		system("pause");
		exit(1);
	}
	prev_file = file;
	prev_line = line;
}

void writeOutput(char* writeString)
{

	std::ofstream out("log.txt");
	out << endl << writeString << endl;
	out.close();
}

#if DO_CUDA_CHECK
	#define CUDA_CHECK cuda_check(__FILE__,__LINE__)
#else
	#define CUDA_CHECK 0
#endif

#endif

//




#endif