#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>

#include "CImg.h"

//flexBox
#include "tools.h"
#include "flexBox.h"
#include "flexTermPrimal.h"
#include "flexLinearOperator.h"
#include "flexIdentityOperator.h"
#include "flexGradientOperator.h"
#include "flexProxDualDataL2.h"
#include "flexProxDualL1Iso.h"

#include "sampleFiles.h"

using namespace std;
using namespace cimg_library;

typedef float floatingType;
typedef thrust::device_vector<floatingType> vectorData;

int main()
{
	float weightDataTerm = 1.0;
	float weightRegularizer = 0.1;

	flexBox<floatingType, vectorData> mainObject;

	//read original image
	CImg<float> imageOriginal(sampleFiles.imgA);

	CImg<float>	imageOriginalGray(imageOriginal.width(), imageOriginal.height(), 1, 1, 0);
	//convert to gray
	for (int i = 0; i < imageOriginal.height(); ++i)
	{
		for (int j = 0; j < imageOriginal.width(); ++j)
		{
			int R = (int)imageOriginal(i, j, 0, 0);
			int G = (int)imageOriginal(i, j, 0, 1);
			int B = (int)imageOriginal(i, j, 0, 2);

			int grayValue = (int)(0.299*R + 0.587*G + 0.114*B);

			imageOriginalGray[j*imageOriginal.width() + i] = grayValue;
		}
	}

	CImg<float> imageNoise = imageOriginalGray;
	imageNoise.noise(30);

	//add primal variable
	std::vector<int> _dims;
	_dims.push_back(imageNoise.height()); _dims.push_back(imageNoise.width());
	mainObject.addPrimalVar(_dims);

	//number of elements in vector
	int nPx = imageNoise.height() * imageNoise.width();

	//add empty term for primal var, because we want to solve the fully dualized problem
	std::vector<int> _correspondingPrimals;
	_correspondingPrimals.push_back(0);

	mainObject.addPrimal(new flexTermPrimal<floatingType, vectorData>(1, 1, primalEmptyProx), _correspondingPrimals);

	//add dualized data term:
	std::vector<flexLinearOperator<floatingType, vectorData>*> operatorList;
	operatorList.push_back(new flexIdentityOperator<floatingType, vectorData>(nPx, nPx, false));

	//convert input image to grey value and add to vector data
	std::vector<floatingType> f(nPx, 0.0f);
	for (int i = 0; i < imageNoise.height(); ++i)
	{
		for (int j = 0; j < imageNoise.width(); ++j)
		{
			f[i*imageNoise.width() + j] = ((floatingType)imageOriginalGray[i*imageOriginal.width() + j]) / 255.0f;
		}
	}

	std::vector<std::vector<floatingType>> fList = std::vector<std::vector<floatingType>>();
	fList.push_back(f);

	//CImgDisplay main_disp(image, "Original"), draw_dispGr(gray, "Grey");

	flexProx<floatingType, vectorData>* myProx = new flexProxDualDataL2<floatingType, vectorData>();
	mainObject.addDual(new flexTermDual<floatingType, vectorData>(myProx, weightDataTerm, (int) _correspondingPrimals.size(), operatorList, fList), _correspondingPrimals);

	//add dualized regularizer
	operatorList.clear();
	//add gradient for x and y direction as operators
	operatorList.push_back(new flexGradientOperator<floatingType, vectorData>(mainObject.getDims(0), 0, 0));
	operatorList.push_back(new flexGradientOperator<floatingType, vectorData>(mainObject.getDims(0), 1, 0));

	flexProx<floatingType, vectorData>* myProx2 = new flexProxDualL1Iso<floatingType, vectorData>();
	mainObject.addDual(new flexTermDual<floatingType, vectorData>(myProx2, weightRegularizer, 1, operatorList), _correspondingPrimals);

	mainObject.runAlgorithm();

	std::vector<float> flexResult = mainObject.getPrimal(0);

	CImg<unsigned char>	imageResult(imageNoise.width(), imageNoise.height(), 1, 1, 0);
	for (int i = 0; i < imageNoise.height(); ++i)
	{
		for (int j = 0; j < imageNoise.width(); ++j)
		{
			imageResult(i, j, 0, 0) = (int)(flexResult[j*imageNoise.width() + i] * 255);
		}
	}

	CImgDisplay main_disp(imageOriginalGray, "Original"), draw_dispGr(imageNoise, "Noise"), draw_dispR(imageResult, "Result");
	while (!main_disp.is_closed() && !draw_dispGr.is_closed() &&
		!main_disp.is_keyESC() && !draw_dispGr.is_keyESC() && !main_disp.is_keyQ() && !draw_dispGr.is_keyQ()) {

		// Handle display window resizing (if any)
		if (main_disp.is_resized()) main_disp.resize().display(imageOriginalGray);
		if (draw_dispGr.is_resized()) draw_dispGr.resize().display(imageNoise);
		if (draw_dispR.is_resized()) draw_dispR.resize().display(imageResult);
	}

    return 0;
}
