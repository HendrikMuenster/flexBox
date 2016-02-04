#include "math.h"
#include "tools.h"
#include "linearInterpolation.h"
#include "cubicInterpolation.h"

float cubicInterpolation(const float *image, float x, float y, const size_t *sizeMat)
{
	int x1 = (int)floorf(x)-1;
	int x2 = x1 + 1;
	int x3 = x2 + 1;
	int x4 = x3 + 1;
	
	int y1 = (int)floorf(y) -1;
	int y2 = y1 + 1;
	int y3 = y2 + 1;
	int y4 = y3 + 1;
	

	if (x1 < 0 || y1 < 0 || x4 >((int)sizeMat[0] - 1) || y4 >((int)sizeMat[1] - 1))
	{
		return linearInterpolation(image,x,y,sizeMat);
	}
	else
	{
		/*float p1[4] = {image[index2DtoLinear(sizeMat, x1, y1)],image[index2DtoLinear(sizeMat, x1, y2)],image[index2DtoLinear(sizeMat, x1, y3)],image[index2DtoLinear(sizeMat, x1, y4)]};
		float p2[4] = {image[index2DtoLinear(sizeMat, x2, y1)],image[index2DtoLinear(sizeMat, x2, y2)],image[index2DtoLinear(sizeMat, x2, y3)],image[index2DtoLinear(sizeMat, x2, y4)]};
		float p3[4] = {image[index2DtoLinear(sizeMat, x3, y1)],image[index2DtoLinear(sizeMat, x3, y2)],image[index2DtoLinear(sizeMat, x3, y3)],image[index2DtoLinear(sizeMat, x3, y4)]};
		float p4[4] = {image[index2DtoLinear(sizeMat, x4, y1)],image[index2DtoLinear(sizeMat, x4, y2)],image[index2DtoLinear(sizeMat, x4, y3)],image[index2DtoLinear(sizeMat, x4, y4)]};
		*/

		float p[4][4] = {
		{ image[index2DtoLinear(sizeMat, x1, y1)],image[index2DtoLinear(sizeMat, x1, y2)],image[index2DtoLinear(sizeMat, x1, y3)],image[index2DtoLinear(sizeMat, x1, y4)] },
		{ image[index2DtoLinear(sizeMat, x2, y1)],image[index2DtoLinear(sizeMat, x2, y2)],image[index2DtoLinear(sizeMat, x2, y3)],image[index2DtoLinear(sizeMat, x2, y4)] },
		{ image[index2DtoLinear(sizeMat, x3, y1)],image[index2DtoLinear(sizeMat, x3, y2)],image[index2DtoLinear(sizeMat, x3, y3)],image[index2DtoLinear(sizeMat, x3, y4)] },
		{ image[index2DtoLinear(sizeMat, x4, y1)],image[index2DtoLinear(sizeMat, x4, y2)],image[index2DtoLinear(sizeMat, x4, y3)],image[index2DtoLinear(sizeMat, x4, y4)] },
		};

		

		
		/*float p[4][4] = {
			{ image[index2DtoLinear(sizeMat, x1, y1)],image[index2DtoLinear(sizeMat, x2, y1)],image[index2DtoLinear(sizeMat, x3, y1)],image[index2DtoLinear(sizeMat, x4, y1)] },
			{ image[index2DtoLinear(sizeMat, x1, y2)],image[index2DtoLinear(sizeMat, x2, y2)],image[index2DtoLinear(sizeMat, x3, y2)],image[index2DtoLinear(sizeMat, x4, y2)] },
			{ image[index2DtoLinear(sizeMat, x1, y3)],image[index2DtoLinear(sizeMat, x2, y3)],image[index2DtoLinear(sizeMat, x3, y3)],image[index2DtoLinear(sizeMat, x4, y3)] },
			{ image[index2DtoLinear(sizeMat, x1, y4)],image[index2DtoLinear(sizeMat, x2, y4)],image[index2DtoLinear(sizeMat, x3, y4)],image[index2DtoLinear(sizeMat, x4, y4)] },
		};*/

		float val =  bicubicInterpolate (p, x - x2, y - y2);

		//mexPrintf("%f\n",val);

		return val;
	}
}


float cubicInterpolate (float p[4], float x)
{
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

float bicubicInterpolate (float p[4][4], float x, float y)
{
	float arr[4];
	arr[0] = cubicInterpolate(p[0], y);
	arr[1] = cubicInterpolate(p[1], y);
	arr[2] = cubicInterpolate(p[2], y);
	arr[3] = cubicInterpolate(p[3], y);
	return cubicInterpolate(arr, x);
}