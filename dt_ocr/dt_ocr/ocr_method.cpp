
#include "stdafx.h"
#include <opencv.hpp>
#include "opencv_modules.hpp"
#include "tiny_dnn/tiny_dnn.h"
#include "boost/foreach.hpp"
#include "boost/filesystem.hpp"
using namespace boost::filesystem;
using namespace tiny_dnn;
using namespace tiny_dnn::activation;
using namespace tiny_dnn::layers;
#define ATTRIBUTES 225  //每一个样本的像素总数.15*15
#define CLASSES    
#define TRAINING_SAMPLES 460  
#define TEST_SAMPLES 200
using namespace std;
using namespace cv;

// convert image to vec_t
void convert_image(const std::string& imagefilename, double scale, int w, int h, std::vector<vec_t>& data)
{
	auto img = cv::imread(imagefilename, cv::IMREAD_GRAYSCALE);
	if (img.data == nullptr) return; // cannot open, or it's not an image

	cv::Mat_<uint8_t> resized;
	cv::resize(img, resized, cv::Size(w, h));
	vec_t d;
	std::transform(resized.begin(), resized.end(), std::back_inserter(d), [=](uint8_t c) { return c * scale; });
	data.push_back(d);
}

// convert all images found in directory to vec_t
void convert_images(const std::string& directory, double scale, int w, int h, std::vector<vec_t>& data)
{
	path dpath(directory);

	BOOST_FOREACH(const path& p, std::make_pair(directory_iterator(dpath), directory_iterator()))
	{
		if (is_directory(p)) continue;
		cout << dpath.leaf() << endl;
		convert_image(p.string(), scale, w, h, data);
	}
}
void construct_cnn()
{
	network<sequential> net;
	//add layers
	net << conv<tan_h>(32, 32, 5, 1, 6)  // in:32x32x1, 5x5conv, 6fmaps
		<< ave_pool<tan_h>(28, 28, 6, 2) // in:28x28x6, 2x2pooling
		<< fc<tan_h>(14 * 14 * 6, 120)   // in:14x14x6, out:120
		<< fc<activation::identity>(120, 8);        // in:120,     out:10

	assert(net.in_data_size() == 32 * 32);
	assert(net.out_data_size() == 8);

	// load MNIST dataset
	std::vector<label_t> train_labels;
	std::vector<vec_t> train_images;

	convert_images("D:\\data\\graph-DB\\nt4\\charImgs\\tiff\\charImgsep\\new2", 1, 32, 32, train_images);

//	parse_mnist_labels("train-labels.idx1-ubyte", &train_labels);
//	parse_mnist_images("train-images.idx3-ubyte", &train_images, -1.0, 1.0, 2, 2);


	// declare optimization algorithm
	adagrad optimizer;

	// train (50-epoch, 30-minibatch)
	net.train<mse>(optimizer, train_images, train_labels, 30, 50);

	// save
	net.save("net");

	// load
	// network<sequential> net2;
	// net2.load("net");
	
}

void opencv_builtin(char* img_path)
{
	
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