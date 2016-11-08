# include "stdafx.h"
# include "diagram.h"


#define N 500
vector<Point2i> edgePoints;
int Sgn(double d)
{
	if (d<0)
		return -1;
	else
		return 1;
}

/*******************************load and segment phase***********************/
Mat image_binarizing(Mat input_image, bool showFlag=false)
{
	/* This function is used for transform image to its binarized image */
	
	int block_size; double c;
	block_size = 13; c = 20;
	Mat binarized_image;
	//binarizing 
	adaptiveThreshold(input_image, binarized_image, 255, CV_ADAPTIVE_THRESH_GAUSSIAN_C, THRESH_BINARY_INV, block_size, c);
	
	//imwrite("test_geo0.png", binarized_image);
	if (showFlag)
	{
		namedWindow("1.binarized image");
		imshow("1.binarized image", binarized_image);//test binarized image
	}
	return binarized_image;
}

void image_labelling(Mat binarized_image,int &labeln, Mat &diagram_segment, vector<Mat> &label_segment, bool showFlag = false)
{
	// this function is used to label image with connectedcomponnent analysis, and store the diagram 
	//and label segment
	Mat labeled(binarized_image.size(), CV_8UC3);
	Mat statsMat, centroidMat; Mat labeled_image; vector<Mat> segments;
	labeln = connectedComponentsWithStats(binarized_image, labeled_image, statsMat, centroidMat, 8, 4);
	
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
	diagram_segment = (labeled_image == dia_idx);
	int left = statsMat.at<int>(dia_idx, 0);
	int top = statsMat.at<int>(dia_idx, 1);
	int h = statsMat.at<int>(dia_idx, 3);
	//for (int label = 1; label < labeln; ++label)
	//{
	//	Mat seg_mat = ((labeled_image == label));
	//	int templeft = statsMat.at<int>(label, 0);
	//	int temptop = statsMat.at<int>(label, 1);
	//	if (templeft >= left && temptop >= top && temptop <= top + h)
	//		diagram_segment += (labeled_image == label);
	//	

	//}
	if (showFlag)
	{
		namedWindow("2.labeled");
		imshow("2.labeled", labeled);
		namedWindow("3.diagram segment");
		imshow("3.diagram segment", diagram_segment); 
	}

}


/********** *********************detection phase *************************/

/**********point part************/

vector<Point2i> getPointPositions(Mat bw)
{
	vector<Point2i> pointPositions;
	for (unsigned int y = 0; y < bw.rows; ++y)
	{
		for (unsigned int x = 0; x < bw.cols; ++x)
		{
			if (bw.at<unsigned char>(y, x)>0)
				pointPositions.push_back(Point2i(x, y));
		}
	}
	return pointPositions;
}

double p2pdistance(Vec2i pt1, Vec2i pt2)
{
	double distance;
	distance = sqrt(powf((pt1[0] - pt2[0]), 2) + powf((pt1[1] - pt2[1]), 2));
	return distance;
}

bool same_pt(Vec2i pt1, Vec2i pt2)
{
	double eps = 8;//to be set parameter
	if (p2pdistance(pt1, pt2) <= eps)
		return true;
	else
		return false;
}

Vec2i avg_pt(vector<Vec2i> pts)
{
	Vec2i new_pt; double xsum = 0.0; double ysum = 0.0;
	int ptsSize = pts.size();
	for (size_t i = 0; i < ptsSize; ++i)
	{
		xsum += pts[i][0];
		ysum += pts[i][1];
	}
	new_pt = { cvRound(xsum / ptsSize), cvRound(ysum / ptsSize) };
	return new_pt;
}

void computeNewPoints(Vec2i pt1, Vec2i pt2, Vec2i pt3, Vec2i pt4, Vec4i &newpts)
{
	/*choose the leftmost and rightmost points in 4 points coz the line is is other */
	int x[4] = { pt1[0], pt2[0], pt3[0], pt4[0] };
	int y[4] = { pt1[1], pt2[1], pt3[1], pt4[1] };
	int x1, y1, x2, y2; int idx1, idx2;
	idx1 = idx2 = 0;
	x1 = x[0]; x2 = x[0];
	for (int i = 1; i < 4; ++i)
	{
		if (x[i] <= x1)
		{
			x1 = x[i];
			idx1 = i;
		}
		if (x[i]>=x2)
		{
			x2 = x[i];
			idx2 = i;
		}
	}
	newpts = {x[idx1],y[idx1], x[idx2], y[idx2]};
}

/*********circle part********/
inline void getCircle(Point2i &p1, Point2i &p2, Point2i &p3, Point2f &center, float &radius)
{
	float x1 = p1.x;float x2 = p2.x;float x3 = p3.x;
	float y1 = p1.y;float y2 = p2.y;float y3 = p3.y;

	center.x = (x1*x1 + y1*y1)*(y2 - y3) + (x2*x2 + y2*y2)*(y3 - y1) + (x3*x3 + y3*y3)*(y1 - y2);
	center.x /= (2 * (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2));

	center.y = (x1*x1 + y1*y1)*(x3 - x2) + (x2*x2 + y2*y2)*(x1 - x3) + (x3*x3 + y3*y3)*(x2 - x1);
	center.y /= (2 * (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2));

	radius = sqrt((center.x - x1)*(center.x - x1) + (center.y - y1)*(center.y - y1));
}

float evaluateCircle(Mat dt, Point2f center, float radius)
{

	float completeDistance = 0.0f;
	int counter = 0;

	float maxDist = 1.0f;   //TODO: this might depend on the size of the circle!

	float minStep = 0.001f;
	// choose samples along the circle and count inlier percentage

	//HERE IS THE TRICK that no minimum/maximum circle is used, the number of generated points along the circle depends on the radius.
	// if this is too slow for you (e.g. too many points created for each circle), increase the step parameter, but only by factor so that it still depends on the radius

	// the parameter step depends on the circle size, otherwise small circles will create more inlier on the circle
	float step = 2 * 3.14159265359f / (6.0f * radius);
	if (step < minStep) step = minStep; // TODO: find a good value here.

	//for(float t =0; t<2*3.14159265359f; t+= 0.05f) // this one which doesnt depend on the radius, is much worse!
	for (float t = 0; t < 2 * 3.14159265359f; t += step)
	{
		float cX = radius*cos(t) + center.x;
		float cY = radius*sin(t) + center.y;

		if (cX < dt.cols)
			if (cX >= 0)
				if (cY < dt.rows)
					if (cY >= 0)
						if (dt.at<float>(cY, cX) <= maxDist)
						{
							completeDistance += dt.at<float>(cY, cX);
							counter++;
						}

	}

	return counter;
}

void detect_circle(Mat diagram_segment, Mat &color_img,Mat &diagram_segwithoutcircle, vector<Vec3f> &circle_candidates, vector<Point2i> &edgePositions, bool showFlag)
{
	unsigned int circleN_todetect = 2;
	diagram_segwithoutcircle = diagram_segment;
	for (unsigned int i = 0; i < circleN_todetect; ++i)
	{
		edgePositions = getPointPositions(diagram_segment);
		Mat dt;
		distanceTransform(255 - diagram_segment, dt, CV_DIST_L1, 3);
		unsigned int nIter = 0;
		Point2f bestCircleCenter; float bestCircleRadius;
		float bestCVal = -1;
		float minCircleRadius = 0.0f;
		for (unsigned int i = 0; i < 2000; ++i)
		{
			int edgepSize = edgePositions.size();
			unsigned int idx1 = rand() % edgepSize;
			unsigned int idx2 = rand() % edgepSize;
			unsigned int idx3 = rand() % edgepSize;

			if ((idx1 == idx2) || (idx1 == idx3) || (idx2 == idx3))
				continue;
			// create a circle from 3 points
			Point2f center; float radius;
			getCircle(edgePositions[idx1], edgePositions[idx2], edgePositions[idx3], center, radius);
			if (radius < minCircleRadius)
				continue;
			float cVal = evaluateCircle(dt, center, radius);
			if (cVal > bestCVal)
			{
				bestCVal = cVal;
				bestCircleRadius = radius;
				bestCircleCenter = center;
			}
			++nIter;
		}
		if (bestCVal<4 * bestCircleRadius) 
			break;
		if (bestCVal<(int)(edgePositions.size() / 8)) 
			break;
		if (bestCVal > 125)
		{
			//std::cout << "current best circle: " << bestCircleCenter << " with radius: " << bestCircleRadius << " and nInlier " << bestCVal << endl;
			Vec3f circle = { bestCircleCenter.x, bestCircleCenter.y, bestCircleRadius };
			circle_candidates.push_back(circle);
			//draw the cicle detected in red within the colorgeo blob image
			cv::circle(color_img, bestCircleCenter, bestCircleRadius, Scalar(0, 0, 255));
			//TODO: hold and save the detected circle.

			//TODO: instead of overwriting the mask with a drawn circle it might be better to hold and ignore detected circles and dont count new circles which are too close to the old one.
			// in this current version the chosen radius to overwrite the mask is fixed and might remove parts of other circles too!

			// update mask: remove the detected circle!
			cv::circle(diagram_segwithoutcircle, bestCircleCenter, bestCircleRadius, 0, 5); // here the radius is fixed which isnt so nice.
			//cv::circle(color_img, bestCircleCenter, bestCircleRadius, Scalar(255, 0, 255), 3);
		}
	}
	
	if (showFlag)
	{
		namedWindow("4.graygeo blob without circles"); cv::imshow("4.graygeo blob without circles", diagram_segwithoutcircle);
		namedWindow("5.colorgeo"); cv::imshow("5.colorgeo", color_img);
	}

}

bool on_circle(Vec2i pt, Vec3f circle)
{
	// check if the point pt is on one of the circles,or joints of multiple circles, or nothing to do with circles
	// joint_flag parameter: 0 means not on any circle, 1 means on a circle, 2 means on two circles and so on
	int dis = 5;// this parameter is to be set to check on the distance tolerants within the distance between radius and distance of pt and circle center point
	int count = 0;

	Vec2f center = { circle[0], circle[1] };
	double radius = circle[2];
	double distance = p2pdistance(center, pt);
	if (abs(distance - radius) <= dis)
		return true;
	else
		return false;

}

/**************line part****************/

bool on_line(Vec4i line, Vec2i pt)
{
	double linedis_eps = 3; double pointdis_eps = 5;
	Vec2i ref_point = { line[0], line[1] };
	Vec2i pt2 = { -pt[1], pt[0] };
	Vec2i ref_point_t = { line[1], -line[0] };
	Vec2i vec = { line[2] - line[0], line[3] - line[1] };
	double point2line = abs((pt2.dot(vec) + ref_point_t.dot(vec))) / sqrt(vec.dot(vec));
	if (point2line < linedis_eps)
		return true;
	else
		return false;
}

bool in_line(Vec4i line, Vec2i pt)
{
	if (on_line(line, pt))
	{
		Vec2i pt1, pt2; pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
		if (p2pdistance(pt1, pt) + p2pdistance(pt2, pt) - p2pdistance(pt1, pt2) < 20)
			return true;
		else
			return false;

	}
	else
		return false;
}
bool on_nin_line(Vec4i line, Vec2i pt)
{
	if (on_line(line, pt))
	{
		Vec2i pt1, pt2; pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
		if (p2pdistance(pt1, pt) + p2pdistance(pt2, pt) - p2pdistance(pt1, pt2) > 20)
			return true;
		else
			return false;
	}
	else
		return false;
}
bool on_other_noncollinearlines(vector<Vec4i> plainLines, Vec4i temp_line, Vec2i pt)
{
	for (size_t i = 0; i < plainLines.size(); i++)
	{
		Vec2i pt1 = { plainLines[i][0], plainLines[i][1] }; Vec2i pt2 = { plainLines[i][2], plainLines[i][3] };
		if ((pt == pt1) || (pt == pt2))
		{
			if ((!on_line(temp_line, pt1)) || (!on_line(temp_line, pt2)))
			{
				return true;
			}
		}
	}
	return false;
}

void findLineEnds(vector<Vec2i> colpoints, vector<Vec2i>& lineEnds, vector<Vec2i>& plainPoints, vector<Vec4i> plainLines)
{
	vector<distance_info> distance_infos;
	for (vector<Vec2i>::iterator iter1 = colpoints.begin(); iter1 != colpoints.end(); iter1++)
	{
		Vec2i pt1 = *iter1;
		for (vector<Vec2i>::iterator iter2 = iter1 + 1; iter2 != colpoints.end(); iter2++)
		{
			Vec2i pt2 = *iter2;
			double tempDistance = p2pdistance(pt1, pt2);
			distance_info dis_info = { pt1, pt2, tempDistance };
			distance_infos.push_back(dis_info);
		}
	}
	std::sort(distance_infos.begin(), distance_infos.end(),
		[](distance_info a, distance_info b){
		return (a.distance > b.distance);
	});
	Vec2i pt3 = distance_infos[0].pt1; Vec2i pt4 = distance_infos[0].pt2;
	lineEnds.push_back(pt3); lineEnds.push_back(pt4);
	int count0 = 0;
	//cout << "colpoints size " << colpoints.size();
	for (vector<Vec2i>::iterator iter3 = colpoints.begin(); iter3 != colpoints.end(); iter3++)
	{
		Vec2i tempPt1 = *iter3;
		bool flag = on_other_noncollinearlines(plainLines, { pt3[0], pt3[1], pt4[0], pt4[1] }, tempPt1);
		for (vector<Vec2i>::iterator iter4 = plainPoints.begin(); iter4 != plainPoints.end();)
		{
			Vec2i tempPt2 = *iter4;
			if (tempPt2 == tempPt1 && (tempPt2 != pt3) && (tempPt2 != pt4) && (!flag))
			{
				iter4 = plainPoints.erase(iter4);
				count0++;
				//cout << "iter4 value " << ++count0 <<" "<<tempPt2<< endl;
				break;
			}
			else
			{
				iter4++;
				continue;
			}
		}
	}
	//cout << " count0 " << count0 << endl;
}


void erase_pointfrom_set(vector<Vec2i>& set, Vec2i pt)
{
	vector<Vec2i>::iterator iter = find(set.begin(), set.end(), pt);
	if (iter != set.end())
		set.erase(iter);
}

Vec2i to_be_removed_similar_endpoint(Vec2i pt1, Vec2i pt2)
{
	Vec2i pt;
	if (pt1[0] < pt2[0])
		pt = pt2;
	else if (pt1[0]>pt2[0])
		pt = pt1;
	else if (pt1[1] < pt2[1])
		pt = pt2;
	else if (pt1[1] > pt2[1])
		pt = pt1;
	else
		cout << "error" << endl;
	return pt;

}

/***************************************ref code*****************************/
//void detect_line(Mat diagram_segwithoutcircle, Mat &color_img, vector<Vec4i> &plainLines, vector<Vec2i>& basicEndpoints)
//{
//	HoughLinesP(diagram_segwithoutcircle, plainLines, 1, CV_PI / 180, 30, 5, 10);
//	vector<Vec2i> plainPoints;
//	// store all the initial points to a vector for handling
//	for (size_t i = 0; i < plainLines.size(); ++i)
//	{
//		Vec2i pt1 = Vec2i(plainLines[i][0], plainLines[i][1]); Vec2i pt2 = Vec2i(plainLines[i][2], plainLines[i][3]);
//		plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//	}
//	// sort the points with a clear rule
//	std::sort(plainPoints.begin(), plainPoints.end(), [](Vec2i a, Vec2i b)
//	{
//		if (a[0] != b[0])
//			return (a[0] < b[0]);
//		else
//			return (a[1] < b[1]);
//	});
//	// define the vector to store the line endpoints along with the new lines
//	vector<Vec2i> lineEnds; vector<Vec4i> newlines;
//	int count = 0;
//	// a loop through the line candidates, first extend the line with all the points on this line
//	// and then find the line endpoints with the longest distance and generate the new lines
//	for (vector<Vec4i>::iterator iter3 = plainLines.begin(); iter3 != plainLines.end(); iter3++)
//	{
//		Vec4i temp_line = *iter3;
//		Vec2i pt1 = { temp_line[0], temp_line[1] }; Vec2i pt2 = { temp_line[2], temp_line[3] };
//
//		vector<Vec2i>::iterator iter5, iter6;
//		iter5 = find(plainPoints.begin(), plainPoints.end(), pt1);
//		iter6 = find(plainPoints.begin(), plainPoints.end(), pt2);
//
//		if (iter5 == plainPoints.end() || iter6 == plainPoints.end())
//		{
//			count++;
//			continue;
//		}
//		else
//		{
//			vector<Vec2i> tempCollinearPoints;
//			tempCollinearPoints.push_back(pt1); tempCollinearPoints.push_back(pt2);
//			for (vector<Vec2i>::iterator iter4 = plainPoints.begin(); iter4 != plainPoints.end(); iter4++)
//			{
//				Vec2i temp_pt = *iter4;
//				if (on_line(temp_line, temp_pt))
//				{
//					tempCollinearPoints.push_back(temp_pt);
//				}
//			}
//			std::sort(tempCollinearPoints.begin(), tempCollinearPoints.end(), [](Vec2i a, Vec2i b)
//			{
//				if (a[0] != b[0])
//					return (a[0] < b[0]);
//				else
//					return (a[1] < b[1]);
//			});
//			tempCollinearPoints.erase(unique(tempCollinearPoints.begin(), tempCollinearPoints.end()), tempCollinearPoints.end());
//			findLineEnds(tempCollinearPoints, lineEnds, plainPoints, plainLines);
//			size_t tempSize = lineEnds.size();
//			Vec4i newLine = { lineEnds[tempSize - 2][0], lineEnds[tempSize - 2][1], lineEnds[tempSize - 1][0], lineEnds[tempSize - 1][1] };
//			newlines.push_back(newLine);
//		}
//	}	
//	cout << "new line size: "<< newlines.size() << endl;
//	std::sort(lineEnds.begin(), lineEnds.end(), [](Vec2i a, Vec2i b)
//	{
//		if (a[0] != b[0])
//			return (a[0] < b[0]);
//		else
//			return (a[1] < b[1]);
//	});
//	lineEnds.erase(unique(lineEnds.begin(), lineEnds.end()), lineEnds.end());
//	cout << "line end size now : " << lineEnds.size() << endl;
//	for (size_t i = 0; i < newlines.size(); i++)
//	{
//		Vec2i pt_1 = { newlines[i][0], newlines[i][1] };
//		Vec2i pt_2 = { newlines[i][2], newlines[i][3] };
//		for (size_t j = i + 1; j < newlines.size(); j++)
//		{
//			Vec2i pt_3 = { newlines[j][0], newlines[j][1] };
//			Vec2i pt_4 = { newlines[j][2], newlines[j][3] };
//			if (same_pt(pt_1, pt_3) && (pt_3 != pt_1))
//			{
//				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_3, pt_1);
//				erase_pointfrom_set(lineEnds, tbrendpt);
//			}
//			else if (same_pt(pt_3, pt_2) && (pt_3 != pt_2))
//			{
//				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_3, pt_2);
//				erase_pointfrom_set(lineEnds, tbrendpt);
//			}
//			else if (same_pt(pt_4, pt_1) && (pt_4 != pt_1))
//			{
//				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_4, pt_1);
//				erase_pointfrom_set(lineEnds, tbrendpt);
//			}
//			else if (same_pt(pt_4, pt_2) && (pt_4 != pt_2))
//			{
//				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_4, pt_2);
//				erase_pointfrom_set(lineEnds, tbrendpt);
//			}
//		}
//	}
//
//	for (size_t i = 0; i < plainPoints.size(); ++i)
//	{
//		bool flag = true;
//		for (size_t j = 0; j < lineEnds.size(); ++j)
//		{
//			if (same_pt(plainPoints[i], lineEnds[j]))
//			{
//				flag = false;
//				break;
//			}
//		}
//		if (flag)
//		{
//			lineEnds.push_back(plainPoints[i]);
//		}
//	}
//	for (size_t i = 0; i < lineEnds.size(); i++)
//	{
//		Vec2i pt1 = lineEnds[i];
//		Scalar temp_color = Scalar((rand() % 255), (rand() % 255), (rand() % 255));
//		cv::circle(color_img, Point(pt1[0], pt1[1]), 1, temp_color, 10, 8, 0);
//	}
//	std::cout << "Now, the size of single points(not circle center) is: " << lineEnds.size() << endl;
//	for (size_t i = 0; i < newlines.size(); i++)
//	{
//		Vec2i pt1 = { newlines[i][0], newlines[i][1] }; Vec2i pt2 = { newlines[i][2], newlines[i][3] };
//		line(color_img, Point(pt1[0], pt1[1]), Point(pt2[0], pt2[1]), Scalar(0, 255, 0), 2, 8, 0);
//	}
//	//namedWindow("points test");
//	//cv::imshow("points test", color_img);
//	lines.clear();
//	lines.assign(newlines.begin(), newlines.end());
//	basicEndpoints.clear();
//	basicEndpoints.assign(lineEnds.begin(), lineEnds.end());
//}

//void detect_line2(Mat diagram_segwithoutcircle, Mat &color_img, vector<Vec4i> &plainLines, vector<Vec2i>& plainPoints)
//{
//	HoughLinesP(diagram_segwithoutcircle, plainLines, 1, CV_PI / 180, 30, 30, 10);
//	//vector<Vec2i> plainPoints;
//	for (size_t i = 0; i < plainLines.size(); ++i)
//	{
//		Vec2i pt1 = { plainLines[i][0], plainLines[i][1] }; Vec2i pt2 = { plainLines[i][2], plainLines[i][3] };
//		if (pt1[0] > 0 && pt1[1] > 0 && pt2[0] > 0 && pt2[1] > 0)
//		{
//			plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//		}
//
//	}
//	// sort the points with a clear rule
//	std::sort(plainPoints.begin(), plainPoints.end(), [](Vec2i a, Vec2i b)
//	{
//		if (a[0] != b[0])
//			return (a[0] < b[0]);
//		else
//			return (a[1] < b[1]);
//	});
//	cout << "test point" << endl;
//	// handle the similar points
//	for (vector<Vec2i>::iterator iter1 = plainPoints.begin(); iter1 != plainPoints.end(); iter1++)
//	{
//		for (vector<Vec2i>::iterator iter2 = iter1 + 1; iter2 != plainPoints.end();)
//		{
//			if (same_pt(*iter1, *iter2))
//			{
//
//				//erase iter2 element
//				if (iter2 != plainPoints.end())
//					iter2 = plainPoints.erase(iter2);
//				else
//				{
//					plainPoints.pop_back();
//				}
//					
//			}
//			else
//			{
//				iter2++;
//			}
//		}
//	}
//	//handle the intersection and add to the set
//	vector<Vec4i> newlines;
//	for (int i = 0; i < plainPoints.size(); ++i)
//	{
//		Vec2i pt1 = plainPoints[i];
//		for (int j = i + 1; j < plainPoints.size(); ++j)
//		{
//			Vec2i pt2 = plainPoints[j];
//			bool flag = false;
//			for (int k = 0; k < plainLines.size(); ++k)
//			{
//				Vec4i l = plainLines[k]; Vec2i pt3 = { l[0], l[1] }; Vec2i pt4 = { l[2], l[3] };
//				if ((same_pt(pt1, pt3) && same_pt(pt2, pt4)) || (same_pt(pt1, pt4) && same_pt(pt2, pt3)))
//				{
//					flag = true;
//					break;
//				}
//			}
//			if (flag)
//			{
//				Vec4i newline = { pt1[0], pt1[1], pt2[0], pt2[1] };
//				newlines.push_back(newline);
//			}
//		}
//	}
//	cout << newlines.size() << endl;
//	plainLines.clear();
//	plainLines.assign(newlines.begin(), newlines.end());
//	cout << "line size after removing the line with the removed points: "<<plainLines.size() << endl;
//	// this is meant to detect the crosses which is not the end point of line
//	for (size_t i = 0; i < plainLines.size(); i++)
//	{
//		Vec4i line1 = plainLines[i]; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
//
//		for (size_t j = i + 1; j < plainLines.size(); j++)
//		{
//			Vec4i line2 = plainLines[j]; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
//			Vec2i cross_cal;
//			bool have_line_cross = intersection(pt1, pt2, pt3, pt4, cross_cal);
//			if (have_line_cross)
//			{
//				bool similar_flag = false;
//				for (size_t k = 0; k < plainPoints.size(); ++k)
//				{
//					if (same_pt(cross_cal, plainPoints[k]))
//					{
//						similar_flag = true;
//						break;
//					}
//				}
//				if (!similar_flag)
//					plainPoints.push_back(cross_cal);
//			}
//		}
//	}
//	for (size_t i = 0; i < plainPoints.size(); ++i)
//	{
//		cout << plainPoints[i][0] << ", " << plainPoints[i][1] << endl;
//	}
//	cout << plainPoints.size() << endl;
//	
//	for (size_t i = 0; i < plainPoints.size(); i++)
//	{
//		Vec2i pt1 = plainPoints[i];
//		Scalar temp_color = Scalar((rand() % 255), (rand() % 255), (rand() % 255));
//		cv::circle(color_img, Point(pt1[0], pt1[1]), 1, temp_color, 5, 8, 0);
//	}
//	for (size_t i = 0; i < plainLines.size(); i++)
//	{
//		Vec2i pt1 = { plainLines[i][0], plainLines[i][1] }; Vec2i pt2 = { plainLines[i][2], plainLines[i][3] };
//		line(color_img, Point(pt1[0], pt1[1]), Point(pt2[0], pt2[1]), Scalar(0, 255, 0), 1, 8, 0);
//	}
//	//namedWindow("points test");
//	//cv::imshow("points test", color_img);
//}

/************concrete function*********/

void point_on_circle_line_check(vector<Vec2i> basicEndpoints, vector<Vec3f> circle_candidates, vector<circleX> &circles,
	vector<Vec4i> line_candidates, vector<lineX> &lines, vector<pointX> &points)
{
	bool flag = true;
	//cout << basicEndpoints.size() << endl;
	for (int i = 0; i < basicEndpoints.size(); ++i)
	{
		Vec2i bpoint = basicEndpoints[i];
		pointX point;
		point.p_idx = i; point.px = bpoint[0]; point.py = bpoint[1];
		for (int j = 0; j < circle_candidates.size(); ++j)
		{
			Vec3f bcircle = circle_candidates[j];
			circleX circle; circle.c_idx = j; circle.cx = bcircle[0]; circle.cy = bcircle[1]; circle.radius = bcircle[2];
			
			//Vec2i center = { cvRound(bcircle[0]),cvRound(bcircle[1])};
			if (flag)
			{
				circles.push_back(circle);
				pointX centerPoint; centerPoint.cflag = true; centerPoint.p_idx = -1; centerPoint.px = cvRound(bcircle[0]); centerPoint.py = cvRound(bcircle[1]);
				points.push_back(centerPoint);
			}
			if (on_circle(basicEndpoints[i], bcircle))
			{
				point.c_idx.push_back(j);
			}
			
		}
		flag = false;
		for (size_t k = 0; k < line_candidates.size(); ++k)
		{
			Vec4i bline = line_candidates[k];
			lineX line; Vec2i bpoint1, bpoint2; bpoint1 = { bline[0], bline[1] }; bpoint2 = { bline[2], bline[3] };
			line.l_idx = k; line.px1 = bline[0]; line.py1 = bline[1]; line.px2 = bline[2]; line.py2 = bline[3];
			line.length = p2pdistance(bpoint1, bpoint2);
			lines.push_back(line);
			point.cflag = false;
			if (on_line(line_candidates[k], bpoint))
			{
				point.l_idx.push_back(k);
			}
		}
		points.push_back(point);
	}
	//cout << points.size() << endl;
}

bool intersection(Vec2i o1, Vec2i p1, Vec2i o2, Vec2i p2,
	Vec2i &r)
{
	Vec2i x = o2 - o1;
	Vec2i d1 = p1 - o1;
	Vec2i d2 = p2 - o2;

	float cross = d1[0] * d2[1] - d1[1] * d2[0];
	if (abs(cross) < /*EPS*/1e-8)
		return false;

	double t1 = (x[0] * d2[1] - x[1] * d2[0]) / cross;
	r = o1 + d1 * t1;
	return true;
}
float evaluateLine(Mat diagram_segwithoutcircle, Vec4i rawLine)
{
	// compute the points number which is on the line
	float maxDist = 1.0f; int count = 0;
	vector<Point2i> edgePositions;
	edgePositions = getPointPositions(diagram_segwithoutcircle);
	Vec2i pt1 = { rawLine[0], rawLine[1] }; Vec2i pt2 = { rawLine[2], rawLine[3] };
	int smaller_x = (pt1[0] < pt2[0]) ? pt1[0] : pt2[0]; 
	int bigger_x = (pt1[0] > pt2[0]) ? pt1[0] : pt2[0]; 
	int smaller_y = (pt1[1] < pt2[1]) ? pt1[1] : pt2[1];
	int bigger_y = (pt1[1] > pt2[1]) ? pt1[1] : pt2[1];
	for (int i = smaller_x; i <= bigger_x; ++i)
	{
		for (int j = smaller_y; j <= bigger_y; ++j)
		{
			Point2i testPoint = CvPoint(i, j);
			Vec2i tp = { i, j };
			if (find(edgePositions.begin(), edgePositions.end(), testPoint) != edgePositions.end())
			{
				float d1 = norm(pt1 - tp);
				float d2 = norm(pt2 - tp);
				float d = norm(pt1 - pt2);
				float dist = d1 + d2 - d;
				if (dist < maxDist)
					count++;
			}
		}
	}
	return count;
}

//void dashlineRecovery(vector<Point2i> edgePositions, Vec4i line, Vec2i pt3, vector<Vec4i>::iterator iter1, vector<Vec2i>::iterator iter2)
//{
//	Vec2i pt1 = { line[0], line[1] }; Vec2i pt2 = { line[2], line[3] };
//	Vec2i pt4; double len1, len2; int gap;
//	len1 = p2pdistance(pt1, pt3); len2 = p2pdistance(pt2, pt3);
//	if (len1 > len2)
//	{
//		pt4 = pt1;
//		gap = abs(pt3[0] - pt1[0]);
//	}
//	else
//	{
//		pt4 = pt2;
//		gap = abs(pt3[0] - pt2[0]);
//	}
//	int bitmap[3000] = { 0 }; int count[10] = { 0 }; bool flag = true;
//	int step = int(gap / 10.0);
//	for (auto i = 0; i < edgePositions.size(); i++)
//	{
//		Vec2i tmp = { edgePositions[i].x, edgePositions[i].y };
//		if (on_line(line, tmp))
//		{
//			bitmap[tmp[0]] = 1;
//		}
//	}
//	for (int i = 0; i < 10; i++)
//	{
//		count[i] = count_if(bitmap + i*step, bitmap + (i + 1)*step, [](int a){return a == 1; });
//	}
//	for (int j = 0; j < 10; j++)
//	{
//		if (count[j] == 0)
//		{
//			flag = false;
//			break;
//		}
//	}
//	//iter1 = plainLines.erase(iter1);
//	//copy(count, count + 10, ostream_iterator<char>(cout, " "));
//
//	
//}
bool in_rect(Vec2i pt,int leftx, int rightx, int lowy, int highy)
{
	int eps = 5;
	if (pt[0] >= leftx - eps && pt[0] <= rightx + eps &&
		pt[1] >= lowy - eps && pt[1] <= highy + eps)
		return true;
	else
		return false;
}
bool dashLineRecovery(vector<Point2i> &edgePositions, Vec2i pt1, Vec2i pt2)
{
	Vec4i line = { pt1[0], pt1[1], pt2[0], pt2[1] }; int c = 0;
	float len = p2pdistance(pt1, pt2); int x1, x2,y1,y2;
	sort(edgePositions.begin(), edgePositions.end(), [](Vec2i a, Vec2i b){return a[0] < b[0]; });
	
	
	if (pt1[0] > pt2[0])
	{
		x1 = pt2[0];
		x2 = pt1[0];
	}
	else
	{
		x1 = pt1[0];
		x2 = pt2[0];
	}
	if (pt1[1] > pt2[1])
	{
		y1 = pt2[1];
		y2 = pt1[1];
	}
	else
	{
		y1 = pt1[1];
		y2 = pt2[1];
	}
	int bitmap[N] = { 0 }; int tmpc = 0;
	float alpha = 0.5;
	if (x2 - x1 > y2 - y1)
	{
		for (auto i = 0; i < edgePositions.size(); i++)
		{
			Vec2i tmpPt = edgePositions[i];
			int x = tmpPt[0];
			if (in_rect(tmpPt, x1, x2, y1, y2) && on_line(line, tmpPt))
			{
				bitmap[x] = 1; tmpc++;
			}
		}
		float low_threshold = alpha*(x2 - x1);
		int ptsCount = count(bitmap, bitmap + N, 1);
		if (ptsCount >= low_threshold)
			return true;
		else
			return false;
	}
	else
	{
		for (auto i = 0; i < edgePositions.size(); i++)
		{
			Vec2i tmpPt = edgePositions[i];
			int y = tmpPt[1];
			if (in_rect(tmpPt, x1, x2, y1, y2) && on_line(line, tmpPt))
			{
				bitmap[y] = 1; tmpc++;
			}
		}
		float low_threshold = alpha*(y2 - y1);
		int ptsCount = count(bitmap, bitmap + N, 1);
		if (ptsCount >= low_threshold)
			return true;
		else
			return false;
	}
		

}



bool isParallel(Vec4i line1, Vec4i line2)
{
	if (line1 == line2)
		return false;
	Vec2i line1V = { line1[2] - line1[0], line1[3] - line1[1] }; Vec2i line2V = { line2[2] - line2[0], line2[3] - line2[1] };
	double theta1 = (line1V[0] <= 3) ? CV_PI / 2.0 : atan2(line1V[1], line1V[0]);
	double theta2 = (line2V[0] <= 3) ? CV_PI / 2.0 : atan2(line2V[1], line2V[0]);
	double angle = abs(theta1 - theta2) / CV_PI * 180;
	return (angle < 5);
}
bool ptSortPred(Vec2i pt1, Vec2i pt2)
{
	return (pt1[0] < pt2[0]);
}
bool with_same_line(vector<Vec4i> &plainLines, Vec2i pt1, Vec2i pt2)
{
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i line = plainLines[i];
		Vec2i tmpPt1, tmpPt2; tmpPt1 = { line[0], line[1] }; tmpPt2 = { line[2], line[3] };
		if ((same_pt(pt1, tmpPt1) && same_pt(pt2, tmpPt2)) || (same_pt(pt1, tmpPt2) && (same_pt(pt2, tmpPt1))))
			return true;
		else if ((same_pt(pt1, tmpPt1) && on_line(line, pt2)) || (same_pt(pt2, tmpPt2) && on_line(line, pt1)) 
				||(same_pt(pt1, tmpPt2) && on_line(line, pt2))|| (same_pt(pt2, tmpPt1) && on_line(line, pt1)))
			return true;
		
	}
	return false;
}
void PointLineRevision(vector<Point2i> &edgePositions, vector<Vec4i> &plainLines, vector<Vec2i> &plainPoints)
{
	////obtain original line and points
	//vector<Vec2i> tmpPlainPoints;
	//for (size_t j = 0; j < plainLines.size(); ++j)
	//{
	//	Vec4i l = plainLines[j];
	//	tmpPlainPoints.push_back({ plainLines[j][0], plainLines[j][1] });
	//	tmpPlainPoints.push_back({ plainLines[j][2], plainLines[j][3] });
	//	//line(color_img, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(255, 255, 0), 1, 8);
	//}
	//// erase the points which is too close to each other
	//sort(tmpPlainPoints.begin(), tmpPlainPoints.end(), [](Vec2i a, Vec2i b){return a[0] < b[0]; });
	//tmpPlainPoints.erase(unique(tmpPlainPoints.begin(), tmpPlainPoints.end(), [](Vec2i a, Vec2i b){return same_pt(a, b); }));

	// dash line recovery
	//while(iter1 != plainLines.end())
	//{
	//	Vec4i line = *iter1; Vec2i pt1 = { line[0], line[1] }; Vec2i pt2 = { line[2], line[3] };
	//	for (auto iter2 = plainPoints.begin(); iter2 != plainPoints.end();)
	//	{
	//		Vec2i pt3 = *iter2;
	//		if (!same_pt(pt3, pt1) && !same_pt(pt3, pt2))
	//		{
	//			// dash line recovery
	//			Vec2i pt1 = { line[0], line[1] }; Vec2i pt2 = { line[2], line[3] };
	//			Vec2i pt4; double len1, len2; int gap;
	//			len1 = p2pdistance(pt1, pt3); len2 = p2pdistance(pt2, pt3); 
	//			bool zflag = true;//assume d(pt1, pt3) is bigger
	//			if (len1 > len2)
	//			{
	//				pt4 = pt1;
	//				gap = abs(pt3[0] - pt1[0]);
	//			}
	//			else
	//			{
	//				pt4 = pt2;
	//				gap = abs(pt3[0] - pt2[0]);
	//				zflag = false;
	//			}
	//			int bitmap[3000] = { 0 }; int count[10] = { 0 }; bool flag = true;
	//			int step = int(gap / 10.0);
	//			for (auto i = 0; i < edgePositions.size(); i++)
	//			{
	//				Vec2i tmp = { edgePositions[i].x, edgePositions[i].y };
	//				if (on_line(line, tmp))
	//				{
	//					bitmap[tmp[0]] = 1;
	//				}
	//			}
	//			for (int i = 0; i < 10; i++)
	//			{
	//				count[i] = count_if(bitmap + i*step, bitmap + (i + 1)*step, [](int a){return a == 1; });
	//			}
	//			for (int j = 0; j < 10; j++)
	//			{
	//				if (count[j] < 10)
	//				{
	//					flag = false;
	//					break;
	//				}
	//			}
	//			if (flag)
	//			{
	//				cout << "test" << endl;
	//				/*iter1 = plainLines.erase(iter1);
	//				Vec4i tmpline = { pt3[0], pt3[1], pt4[0], pt4[1] }; plainLines.push_back(tmpline);*/
	//				if (zflag)
	//				{
	//					line[2] = pt3[0]; line[3] = pt3[1];

	//				}
	//			}
	//			else
	//			{
	//				iter1++;
	//				plainPoints.push_back(pt1);
	//				plainPoints.push_back(pt2);
	//				iter2++;
	//			}
	//		
	//			//copy(count, count + 10, ostream_iterator<char>(cout, " "));	
	//		}
	//	}
	//}
	
	
	for (auto iter1 = plainLines.begin(); iter1 != plainLines.end();iter1++)
	{
		Vec4i line1 = *iter1; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
		for (auto iter2 = plainLines.begin(); iter2 != plainLines.end();)
		{
			Vec4i line2 = *iter2; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
			if (isParallel(line1, line2))
			{
				map<int, int> pMap; // pMap is a Map storing the x-axis and y-axis and from x-axis we can infer the y axis
				pMap[pt1[0]] = pt1[1]; pMap[pt2[0]] = pt2[1]; 
				pMap[pt3[0]] = pt3[1]; pMap[pt4[0]] = pt4[1];
				auto iter_begin = pMap.begin(); 
				auto iter_end = pMap.end(); iter_end--;
				//if two line is parallel, find the two closest point in 4 points, check if there're dash line
				Vec2i p1 = { iter_begin->first, iter_begin->second }; Vec2i p2 = { iter_end->first, iter_end->second };
				if (dashLineRecovery(edgePositions, p1, p2))
				{
					// check if there's a dash line to be recovered and erase the old line and points and add the new points and lines
					// the second and the third points in the map is to be erased
					plainPoints.push_back(p1); plainPoints.push_back(p2);
					// assign new line to the current iter1 line and erase the line which iter2 points to 
					Vec4i newline = { p1[0], p1[1], p2[0], p2[1] };
					*iter1 = newline;
					iter2 = plainLines.erase(iter2);
				}
				else
				{
					plainPoints.push_back(pt1); plainPoints.push_back(pt2);
					plainPoints.push_back(pt3); plainPoints.push_back(pt4);
					iter2++;
				}

				
			}
			else
			{
				//not parallel, choose the two end points pt3, pt4 to be tested for the dash line recovery
				if (same_pt(pt3, pt1)||same_pt(pt3, pt2)||same_pt(pt4,pt1)||same_pt(pt4,pt2))
				{
					plainPoints.push_back(pt1); plainPoints.push_back(pt2);
				}
				else if (on_nin_line(line1, pt3))
				{
					if (abs(pt1[0] - pt3[0]) > abs(pt2[0] - pt3[0]))
					{
						// pt1-pt3 longer
						if (dashLineRecovery(edgePositions, pt2, pt3))
						{
							Vec4i newline = { pt1[0], pt1[1], pt3[0], pt3[1] };
							*iter1 = newline;
							plainPoints.push_back(pt1); plainPoints.push_back(pt3);
						}
						else
						{
							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
						}
					}
					else
					{
						//pt2-pt3 longer
						if (dashLineRecovery(edgePositions, pt1, pt3))
						{
							Vec4i newline = { pt2[0], pt2[1], pt3[0], pt3[1] };
							*iter1 = newline;
							plainPoints.push_back(pt2); plainPoints.push_back(pt3);
						}
						else
						{
							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
						}
					}
				}
				else if (on_nin_line(line1, pt4))
				{
					if (abs(pt1[0] - pt4[0]) > abs(pt2[1] - pt4[0]))
					{
						if (dashLineRecovery(edgePositions, pt2, pt4))
						{
							Vec4i newline = { pt1[0], pt1[1], pt4[0], pt4[1] };
							*iter1 = newline;
							plainPoints.push_back(pt1); plainPoints.push_back(pt4);
						}
						else
						{
							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
						}
					}
					else
					{
						if (dashLineRecovery(edgePositions, pt1, pt4))
						{
							Vec4i newline = { pt2[0], pt2[1], pt4[0], pt4[1] };
							*iter1 = newline;
							plainPoints.push_back(pt2); plainPoints.push_back(pt4);
						}
						else
						{
							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
						}

					}
				}
				else
				{
					cout << "test" << endl;
				}
				iter2++;
			}
		}
	}
	// erase the points which is too close to each other
	sort(plainPoints.begin(), plainPoints.end(), [](Vec2i a, Vec2i b){return a[0] < b[0]; });
	plainPoints.erase(unique(plainPoints.begin(), plainPoints.end(), [](Vec2i a, Vec2i b){return same_pt(a, b); }), plainPoints.end());
	for (auto i = 0; i < plainPoints.size(); i++)
	{
		Vec2i pt1 = plainPoints[i];
		if (pt1[1] == 77)
			cout << "li" << endl;
		for (auto j = i+1; j < plainPoints.size(); j++)
		{
			Vec2i pt2 = plainPoints[j];
			if (pt2[1] == 77)
				cout << "li" << endl;
			if (with_same_line(plainLines, pt1, pt2))
				continue;
			else
			{
				if (dashLineRecovery(edgePoints, pt1, pt2))
				{
					Vec4i newline = { pt1[0], pt1[1], pt2[0], pt2[1] };
					plainLines.push_back(newline);
				}
			}
			
		}
	}
}


void detect_line3(vector<Point2i> &edgePositions, Mat diagram_segwithoutcircle, Mat &color_img, vector<Vec4i> &plainLines, vector<Vec2i>& plainPoints, Mat &drawedImages, bool showFlag=true, string fileName="")
{
#pragma region 

	/*detect line candidates by probabilistic hough transform */
	vector<Vec4i> rawLines;
	HoughLinesP(diagram_segwithoutcircle, rawLines, 1, CV_PI / 180, 15, 15, 10);
	if (showFlag)
	{
		for (size_t j = 0; j < rawLines.size(); ++j)
		{
			Vec4i l = rawLines[j];
			line(color_img, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 255, 255), 2, 8);
		}
	}
	if (showFlag)
	{ 
		namedWindow("6.lines first opt version now 2"); 
		imshow("6.lines first opt version now 2", color_img);

	};
	for (size_t i = 0; i < rawLines.size(); ++i)
	{
		//eliminate the false detected lines with few points on it
		Vec4i rawLine = rawLines[i];
		float lineEV = evaluateLine(diagram_segwithoutcircle, rawLine);
		if (lineEV > 10)//this threshold should be set a litter lower to generate enough line candidates
		{
			plainLines.push_back(rawLine);
			line(diagram_segwithoutcircle, Point(rawLine[0], rawLine[1]), Point(rawLine[2], rawLine[3]), Scalar(255, 0, 0), 2, 8);
		}
	}
	//display zone
	if (showFlag)
	{	
		for (size_t j = 0; j < plainLines.size(); ++j)
		{
			Vec4i l = plainLines[j];
			line(color_img, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 255, 0), 2, 8);
		}
		namedWindow("7.lines first opt version now 1"); imshow("7.lines first opt version now 1", color_img);
	}
	else if (fileName != "")
	{
		// write into txt
		ofstream ofile;
		ofile.open(fileName,ios_base::app);
		for (size_t k = 0; k < plainLines.size(); ++k)
		{
			Vec4i l = plainLines[k];
			cout << l[0] << "," << l[1] << endl << l[2] << "," << l[3] << endl;
			ofile << l[0] << "," << l[1] << "\n" << l[2] << "," << l[3] << "\n";
		}
		ofile << "\n";
		cout << endl;
		ofile.close();
	}


	/* then we handle the lines specifically*/
	/* for lines should be combined 1. colinear 2. */
	vector<double> angs; vector<double> slopes;
	// store the infomation of line slope and slope angle
	for (vector<Vec4i>::iterator iter = plainLines.begin(); iter != plainLines.end(); iter++)
	{
		// for collinear
		Vec4i l = *iter; Vec2i ld = { l[2] - l[0], l[3] - l[1] };
		double slope = (ld[0] <= 3) ?  10000000 : ld[1] * 1.0 / ld[0];
		double ang = (ld[0] <= 3) ? CV_PI / 2 : atan2(ld[1], ld[0]);
		slopes.push_back(slope);
		angs.push_back(ang);
	}
#pragma endregion

#pragma region 

	// loop througth, find the collinear line and combine
	for (vector<Vec4i>::iterator iter1 = plainLines.begin(); iter1 != plainLines.end();iter1++)
	{
		Vec4i l1 = *iter1;
		int i = iter1 - plainLines.begin();// log the location of current line
		//cout << " current the first line: " << i << ";" << endl;
		double slope1 = slopes[i];
		
		Vec2i pt1 = { l1[0], l1[1] }; Vec2i pt2 = { l1[2], l1[3] };
		Vec2i ld1 = { pt2[0] - pt1[0], pt2[1] - pt1[1] };

		double length1 = norm(ld1); Vec2f mp1 = (pt1 + pt2) / 2.0;
		//cout << pt1[0] << "," << pt1[1] << "  " << pt2[0] << "," << pt2[1] << endl;
		for (vector<Vec4i>::iterator iter2 = iter1 + 1; iter2 != plainLines.end();)
		{
			int ct = 0;
			if (iter1 != iter2)
			{			
				int j = iter2 - plainLines.begin() + ct;
				//cout << " currrent the second line: " << j << ";" << endl;
				double slope2 = slopes[j];		

				Vec2i pt3 = { plainLines[j][0], plainLines[j][1] }; Vec2i pt4 = { plainLines[j][2], plainLines[j][3] };
				Vec2i ld2 = { pt4[0] - pt3[0], pt4[1] - pt3[1] };

				double length2 = norm(ld2); Vec2f mp2 = (pt3 + pt4) / 2.0;

				//cout << pt3[0] << "," << pt3[1] << "  " << pt4[0] << "," << pt4[1] << endl;
				double sloped = abs(slope1 - slope2);
				//cout << "the slope diff is: "<< sloped << endl;
				
				if (sloped < 0.2)// the slope is considered to be equal 
				{
					if (on_line(plainLines[j], pt1))// the point on line i is on line j, then collinear
					{
						//cout << abs(norm(mp1 - mp2) - (length1 + length2) / 2) << endl;
						int dis = (norm(mp1 - mp2) - (length1 + length2) / 2);
						//cout << "dis: " << dis << endl;
						if (dis< 5 && dis > -5) // the two line are close to be combined
						{
							Vec4i newpts;
							computeNewPoints(pt1, pt2, pt3, pt4, newpts);
							*iter1 = newpts;
							//plainPoints.push_back({ newpts[0], newpts[1] }); plainPoints.push_back({ newpts[2], newpts[3] });
							iter2 = plainLines.erase(iter2);	
							slopes.erase(slopes.begin() + j);
							ct++;
						}
						else if (dis <= -5)
						{
							//cout << "test test" << endl;
							if (length1 < length2)
							{
								*iter1 = *iter2;
							}
							iter2 = plainLines.erase(iter2);
							slopes.erase(slopes.begin() + j);

						}
						else
						{
							iter2++;
						}
					}
					else
					{
						iter2++;
					}
				}
				else
				{
					
					iter2++;
				}
			}
			else
			{
				iter2++;
			}

		}
	}
	
#pragma endregion


	PointLineRevision(edgePositions, plainLines,plainPoints);
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		line(color_img, pt1, pt2, Scalar(0, 0, 255), 2, 8, 0);
	}
	if (showFlag)
	{
		namedWindow("8.lines first opt version now"); imshow("8.lines first opt version now", color_img);
	}
	drawedImages = color_img;
}

Mat preprocessing(Mat diagram_segment)
{
	//namedWindow("origianl diagram seg"); imshow("original diagram seg", diagram_segment);
	//now we'll do some preprocessing to make for the following detection procedure
	Mat eroded, dilated, opened, closed;
	Mat eroEle, dilEle;
	eroEle = getStructuringElement(MORPH_RECT, Size(3, 3));
	morphologyEx(diagram_segment, dilated, MORPH_DILATE, eroEle);
	//erode(diagram_segment, eroded, eroEle);
	namedWindow("dilated image"); imshow("dilated image", dilated);
	//namedWindow("eroded diagram"); imshow("eroded diagram", eroded);
	morphologyEx(diagram_segment, closed, MORPH_CLOSE, eroEle);
	namedWindow("closed diagram"); imshow("closed diagram", closed);
	return closed;
}


void primitive_parse(Mat &image, Mat diagram_segment, vector<pointX> &points, vector<lineX> &lines, vector<circleX> &circles, Mat &drawedImages, bool showFlag=true, string fileName="")
{
	/* primitive about points, lines, and circles
	first detect the circle, we can get the matrix of diagram segments without circle for the next
	 processing step*/
	vector<Vec3f> circle_candidates; Mat color_img; cvtColor(diagram_segment, color_img, CV_GRAY2RGB);
	Mat diagram_segwithoutcircle;
	
	//diagram_segment = preprocessing(diagram_segment);

	/*ransac go*/
	vector<Point2i> edgePositions;
	detect_circle(diagram_segment, color_img, diagram_segwithoutcircle,circle_candidates,edgePositions,showFlag);
	// then the line detection
	vector<Vec4i> line_candidates; vector<Vec2i> basicEndpoints;
	detect_line3(edgePoints ,diagram_segwithoutcircle, color_img, line_candidates, basicEndpoints, drawedImages,showFlag, fileName);
	
	//detect_line2(diagram_segwithoutcircle, color_img, line_candidates, basicEndpoints);
	//cout << "basic endpoints num: "<<basicEndpoints.size() << endl;
	
	// for points check if they are in circles
	point_on_circle_line_check(basicEndpoints, circle_candidates, circles, line_candidates, lines, points);
	
	/*display point text info*/
	//for (int i = 0; i < points.size(); ++i)
	//{
	//	pointX point = points[i];
	//	//cout << "point " << point.p_idx << " ("<<point.px<<", "<< point.py<<") " << " on circle ";
	//	for (int j = 0; j < point.c_idx.size(); ++j)
	//	{
	//		//cout << point.c_idx[j] << ", ";
	//	}
	//	cout << "on line ";
	//	for (int k = 0; k < point.l_idx.size(); ++k)
	//	{
	//		//cout << point.l_idx[k] << ", ";
	//	}
	//	//cout << endl;
	//}

	}

int test_diagram()
{
	//first load a image
	Mat image = imread("test1.jpg", 0);
	//namedWindow("original image");
	//imshow("original image", image);
	// then binarize it

	Mat binarized_image = image_binarizing(image, true);
	//binarized_image = preprocessing(binarized_image);
	// then go on a process of connectivity componnent analysis
	int labeln; Mat diagram_segment; vector<Mat> label_segment;
	
	edgePoints = getPointPositions(binarized_image);
	Mat pointss = Mat::zeros(1000, 1000, CV_8UC3);
	for (auto i = 0; i < edgePoints.size(); i++)
	{
		Point2i pt = edgePoints[i];
		circle(pointss, pt, 1, Scalar(0, 0, 255));
	}
	namedWindow("points"); imshow("points", pointss);
	image_labelling(binarized_image, labeln, diagram_segment, label_segment,true);
	vector<pointX> points; vector<lineX> lines; vector<circleX> circles;		Mat drawedImages;
	primitive_parse(image, diagram_segment, points, lines, circles, drawedImages);
	return 0;
}

int diagram()
{
	//a series of image
	//vector<Mat> images;
	char abs_path[100] = "D:\\data\\graph-DB";
	char imageName[150], saveimgName[150];
	string outputFN = "D:\\data\\graph-DB\\output.txt";
	for (int i = 1; i < 136; i++)
	{
		sprintf_s(imageName, "%s\\Sg-%d.jpg", abs_path, i);
		sprintf_s(saveimgName, "%s\\saveImage\\sgs-%d.jpg", abs_path, i);
		//first load a image
		Mat image = imread(imageName, 0);
		//namedWindow("original image");
		//imshow("original image", image);
		// then binarize it
		Mat binarized_image = image_binarizing(image);
		// then go on a process of connectivity componnent analysis
		int labeln; Mat diagram_segment; vector<Mat> label_segment;
		image_labelling(binarized_image, labeln, diagram_segment, label_segment);
		vector<pointX> points; vector<lineX> lines; vector<circleX> circles;
		Mat drawedImages;
		primitive_parse(image, diagram_segment, points, lines, circles, drawedImages, false, outputFN);
		//imwrite(saveimgName, drawedImages);
	}
	return 0;
}