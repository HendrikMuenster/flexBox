int index2DtoLinear(const int* sizeMat, int i, int j)
{
	return (int)(i + j*sizeMat[0]);
}