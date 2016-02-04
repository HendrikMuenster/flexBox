// ConsoleApplication1.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//

#include "stdafx.h"

#include "opencv2/core/core.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <iostream>
#include <sstream>
#include <iomanip>
using namespace std;
using namespace cv;



#include "math.h"
#include <omp.h>
#include <stdio.h>
#include <sys/types.h>


//my includes
#include "flexTermPrimal.h"

#include "flexL2DataTerm.h"
#include "flexL2Gradient.h"

#include "helper.h"
#include "flexbox.h"
#include "flexMatrix.h"
#include "flexVector.h"
#include "flexBoxData.h"





#include <cstddef>
#include <ctime>

int main()
{
	int numImgs = 1;
	int startImg = 0;
	bool gray = false;

	flexVector<float> v1(2, 0.0f), v2(2, 0.0f), v3(2, 0.0f);

	v1[0] = 1;
	v2[1] = 1;

	v3 = v1;
	v3 += v2;

	v1.print();
	v2.print();
	v3.print();

	flexVector<flexVector<float>> test;
	test.push_back(v1);

	test[0].print();

	test[0] += v2;

	test[0].print();


	cv::Mat image, imageFloat, greyImage, greyImage2;

	//image = imread("D:/Dropbox/Uni/Projects/2015 - flexbox/CPP/a.png", 1);
	image = imread("C:/Users/Hendrik/Dropbox/Uni/Projects/2015 - flexbox/CPP/a.png",1);

	cv::cvtColor(image, greyImage, CV_BGR2GRAY);
	//greyImage.convertTo(imageFloat, CV_32F);

	imageFloat = greyImage;
	greyImage2 = greyImage;

	int sizeImage[2] = { greyImage.cols, greyImage.rows };

	int nPx = (int)(greyImage.rows * greyImage.cols);

	flexBox<float> mainObject;

	flexVector<int> _dims; _dims.push_back(sizeImage[0]); _dims.push_back(sizeImage[1]);

	float _alpha = 1.0;

	flexVector<float> _f(nPx, 0.0f);

	//define gradient matrix
	flexMatrix<float> Dx(nPx, nPx), DxT(nPx, nPx), Dy(nPx, nPx), DyT(nPx, nPx);
	for (int j = 0; j < sizeImage[1] - 1; ++j)
	{
		for (int i = 0; i < sizeImage[0] - 1; ++i)
		{
			_f[index2DtoLinear(sizeImage, i, j)] = imageFloat.data[index2DtoLinear(sizeImage, i, j)] / 255.0f;
			greyImage2.data[index2DtoLinear(sizeImage, i, j)] = 255.0f*_f[index2DtoLinear(sizeImage, i, j)];


			Dx.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i + 1, j), 1.0);
			Dx.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j), -1.0);

			Dy.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j + 1), 1.0);
			Dy.insertElement(index2DtoLinear(sizeImage, i, j), index2DtoLinear(sizeImage, i, j), -1.0);

		}
	}

	flexL2DataTerm<float> dataTerm(_dims, _alpha, _f);

	mainObject.addPrimal(&dataTerm);

	flexVector<flexMatrix<float>> _operatorList;
	_operatorList.push_back(Dx);
	_operatorList.push_back(Dy);

	flexVector<float> sigmaTau;
	sigmaTau.push_back(4);
	sigmaTau.push_back(4);

	flexVector<int> corrPrimals;
	corrPrimals.push_back(0);

	flexL2Gradient<float> regularizer(2, _operatorList, sigmaTau, sigmaTau);
	mainObject.addDual(&regularizer, corrPrimals);

	namedWindow("Before", CV_WINDOW_NORMAL);
	imshow("Before", greyImage);


	mainObject.runAlgorithm();

	flexVector<float> result = mainObject.getPrimal(0);



	for (int j = 0; j < sizeImage[1] - 1; ++j)
	{
		for (int i = 0; i < sizeImage[0] - 1; ++i)
		{
			//cout << result[index2DtoLinear(sizeImage, i, j)] << endl;

			greyImage2.data[index2DtoLinear(sizeImage, i, j)] = 255.0f*result[index2DtoLinear(sizeImage, i, j)];
		}
	}

	namedWindow("After", CV_WINDOW_NORMAL);
	imshow("After", greyImage2);


	waitKey(0);


	//imwrite("C:/Users/Hendrik/Dropbox/Uni/Projects/2015 - flexbox/CPP/result.png", greyImage);

    return 0;
}

