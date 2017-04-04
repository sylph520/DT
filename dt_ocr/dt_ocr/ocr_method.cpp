
#include "stdafx.h"
#include <opencv.hpp>
#include <ml.hpp>
#include "opencv_modules.hpp"
#define ATTRIBUTES 225  //每一个样本的像素总数.15*15
#define CLASSES    
#define TRAINING_SAMPLES 460  
#define TEST_SAMPLES 200
using namespace std;
using namespace cv;
using namespace cv::ml;


void opencv_builtin(char* img_path)
{
	//mlp
	
}

void test_opencv_mlp()
{
	Mat_<float> data(100, 100);
	randn(data, Mat::zeros(1, 1, data.type()), Mat::ones(1, 1, data.type()));

	//half of the samples for each class
	Mat_<float> responses(data.rows, 2);
	for (int i = 0; i<data.rows; ++i)
	{
		if (i < data.rows / 2)
		{
			responses(i, 0) = 1;
			responses(i, 1) = 0;
		}
		else
		{
			responses(i, 0) = 0;
			responses(i, 1) = 1;
		}
	}
	Mat_<int> layerSizes(1, 3);
	layerSizes(0, 0) = data.cols;
	layerSizes(0, 1) = 20;
	layerSizes(0, 2) = responses.cols;
	Ptr<ANN_MLP> network = ANN_MLP::create();
	network->setLayerSizes(layerSizes);
	network->setActivationFunction(ANN_MLP::SIGMOID_SYM, 0.1, 0.1);
	network->setTrainMethod(ANN_MLP::BACKPROP, 0.1, 0.1);
	Ptr<TrainData> TrainData = ml::TrainData::create(data, ROW_SAMPLE, responses);

	network->train(TrainData);
		if (network->isTrained())
		{
			printf("Predict one-vector:\n");
			Mat result;
			network->predict(Mat::ones(1, data.cols, data.type()), result);
			cout << result << endl;

			printf("Predict training data:\n");
			for (int i = 0; i<data.rows; ++i)
			{
				network->predict(data.row(i), result);
				cout << result << endl;
			}
		}
	
}
void opencv_tesseract(char* img_path)
{

}
void external_tesseract(char* img_path)
{
	
}
void gocr_test(char* img_path)
{
	
}

void opencv_cnn(char* img_path)
{
	
}