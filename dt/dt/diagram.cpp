# include "stdafx.h"
# include "diagram.h"
Mat image_binarizing(Mat input_image)
{
	/* This function is used for transform image to its binarized image */
	int block_size; double c;
	block_size = 13; c = 20;
	Mat binarized_image;
	//binarizing 
	adaptiveThreshold(input_image, binarized_image, 255, CV_ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY_INV, block_size, c);
	//imwrite("test_geo0.png", binarized_image);
	//namedWindow("binarized image");
	//imshow("binarized image", binarized_image);//test binarized image
	return binarized_image;
}
void image_labelling(Mat binarized_image,int &labeln, Mat &diagram_segment, vector<Mat> &label_segment)
{
	// this function is used to label image with connectedcomponnent analysis, and store the diagram 
	//and label segment
	Mat labeled(binarized_image.size(), CV_8UC3);
	Mat statsMat, centroidMat; Mat labeled_image; vector<Mat> segments;
	labeln = connectedComponentsWithStats(binarized_image, labeled_image, statsMat, centroidMat, 8, 4);
	//cout << statsMat << endl;
	// show the labelling image
	vector<Vec3b> colors(labeln);
	colors[0] = Vec3b(0, 0, 0);
	int dia_idx = 1;
	// random labelling color generation and segmented image logical matrix
	for (int label = 1; label < labeln; ++label)
	{
		//randomly generate the color the current label number
		
		colors[label] = Vec3b((rand() & 255), (rand() & 255), (rand() & 255));
		Mat seg_mat = ((labeled_image == label));
		int seg_area = statsMat.at<int>(label, 4);
		if (seg_area > 5)
		{
			segments.push_back(seg_mat);
			if (seg_area >= statsMat.at<int>(dia_idx, 4))
				dia_idx = label;
		}	
	}
	for (int r = 0; r < binarized_image.rows; ++r) {
		for (int c = 0; c < binarized_image.cols; ++c) {
			int label = labeled_image.at<int>(r, c);
			Vec3b &pixel = labeled.at<Vec3b>(r, c);
			pixel = colors[label];
		}
	}
	//namedWindow("labeled");
	//imshow("labeled", labeled);
	diagram_segment = (labeled_image == dia_idx);
	//namedWindow("diagram segment");
	//imshow("diagram segment", diagram_segment);
}

void primitive_parse(Mat diagram_segment, vector<pointX> &points, vector<lineX> &lines, vector<circleX> &circles)
{
	// primitive about points, lines, and circles
	//first detect the circle
	vector<Vec3f> circle_candidates; Mat color_img; cvtColor(diagram_segment, color_img, CV_GRAY2RGB);
	/*hough go*/
	//HoughCircles(diagram_segment, circle_candidates, HOUGH_GRADIENT, 2, diagram_segment.rows / 6, 60, 60);
	//for (size_t i = 0; i < circle_candidates.size(); i++)
	//{
	//	Point center(cvRound(circle_candidates[i][0]), cvRound(circle_candidates[i][1]));
	//	int radius = cvRound(circle_candidates[i][2]);
	//	// draw the circle center
	//	circle(color_img, center, 3, Scalar(0, 255, 0), -1, 8, 0);
	//	// draw the circle outline
	//	circle(color_img, center, radius, Scalar(0, 0, 255), 3, 8, 0);
	//}
	//namedWindow("circles", 1);
	//imshow("circles", color_img);
	/*mser blob*/
	Ptr<MSER> ms = MSER::create();
}

int test_diagram()
{
	//first load a image
	Mat image = imread("000.jpg", 0);
	//namedWindow("original image");
	//imshow("original image", image);
	// then binarize it
	Mat binarized_image = image_binarizing(image);
	// then go on a process of connectivity componnent analysis
	int labeln; Mat diagram_segment; vector<Mat> label_segment;
	image_labelling(binarized_image, labeln, diagram_segment, label_segment);
	vector<pointX> points; vector<lineX> lines; vector<circleX> circles;
	primitive_parse(diagram_segment, points, lines, circles);
	return 0;
}