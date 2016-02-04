#include "math.h"
#include "tools.h"

float linearInterpolation(const float *image, float x, float y, const size_t *sizeMat)
{
	int x1 = (int)floorf(x);
	int x2 = (int)ceilf(x);
	int y1 = (int)floorf(y);
	int y2 = (int)ceilf(y);

	if (x1 < 0)
	{
		x1 = 0;
		x2 = 0;
	}

	if (y1 < 0)
	{
		y1 = 0;
		y2 = 0;
	}

	if (x2 >((int)sizeMat[0] - 1))
	{
		x1 = (int)sizeMat[0] - 1;
		x2 = (int)sizeMat[0] - 1;
	}

	if (y2 > ((int)sizeMat[1] - 1))
	{
		y1 = (int)sizeMat[1] - 1;
		y2 = (int)sizeMat[1] - 1;
	}

	if (x1 == x2 && y1 == y2)
	{
		// we are on a regular grid point
		return image[index2DtoLinear(sizeMat, x1, y1)];
	}
	else if (x1 == x2)
	{
		// interpolate in y direction
		float frac = y - (float)y1;
		return (1.0f - frac) * image[index2DtoLinear(sizeMat, x1, y1)] + frac * image[index2DtoLinear(sizeMat, x1, y2)];
	}
	else if (y1 == y2)
	{
		// interpolate in x direction
		float frac = x - (float)x1;
		return (1.0f - frac) * image[index2DtoLinear(sizeMat, x1, y1)] + frac * image[index2DtoLinear(sizeMat, x2, y1)];
	}
	else
	{
		float frac1 = x - (float)x1;
		float frac2 = y - (float)y1;

		float inter1 = (1.0f - frac1) * image[index2DtoLinear(sizeMat, x1, y1)] + frac1 * image[index2DtoLinear(sizeMat, x2, y1)];
		float inter2 = (1.0f - frac1) * image[index2DtoLinear(sizeMat, x1, y2)] + frac1 * image[index2DtoLinear(sizeMat, x2, y2)];

		return (1.0f - frac2) * inter1 + frac2 * inter2;
	}
}