# include "stdafx.h"
# include "diagram.h"

#define N 500

bool dashLineRecovery(vector<Point2i>&, Vec2i, Vec2i, Vec2i, vector<circle_class>&, bool plflag, bool pcflag, bool ppflag);

bool dashLineRecovery(vector<Point2i>&, Vec2i, Vec2i, vector<circle_class>&, bool plflag, bool pcflag, bool ppflag);

/*******************************load and segment phase***********************/
Mat image_binarizing(Mat input_image, bool showFlag = false)
{
	/* This function is used for transform image to its binarized image */

	int block_size;
	double c;
	block_size = 13;
	c = 20;
	Mat binarized_image = Mat::zeros(input_image.size(),CV_8UC1);
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

void image_labelling(Mat binarized_image, Mat& diagram_segment, bool showFlag = false)
{
	// this function is used to label image with connectedcomponnent analysis, and store the diagram 
	//and label segment
	Mat labeled(binarized_image.size(), CV_8UC3);
	Mat statsMat, centroidMat;
	Mat labeled_image;
	vector<Mat> segments;
	int labeln = connectedComponentsWithStats(binarized_image, labeled_image, statsMat, centroidMat, 8, 4);

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
	for (int r = 0; r < binarized_image.rows; ++r)
	{
		for (int c = 0; c < binarized_image.cols; ++c)
		{
			int label = labeled_image.at<int>(r, c);
			Vec3b& pixel = labeled.at<Vec3b>(r, c);
			pixel = colors[label];
		}
	}
	diagram_segment = (labeled_image == dia_idx);
	//int left = statsMat.at<int>(dia_idx, 0);
	//int top = statsMat.at<int>(dia_idx, 1);
	//int h = statsMat.at<int>(dia_idx, 3);
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
			if (bw.at<unsigned char>(y, x) > 0)
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

point_class* get_line_pt_by_id(vector<point_class>& pointxs, int id)
{
	auto pos = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
	                   {
		                   if (a.getPid() == id)
			                   return true;
		                   else
			                   return false;
	                   });
	point_class* p = nullptr;
	if (pos == pointxs.end())
	{
		cout << "error" << endl;
	}
	else
	{
		p = &(*pos);
	}
	return p;
}

Vec2i get_line_ptVec_by_id(vector<point_class>& pointxs, int id)
{
	auto pos = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
	                   {
		                   if (a.getPid() == id)
			                   return true;
		                   else
			                   return false;
	                   });
	if (pos == pointxs.end())
	{
		cout << "error" << endl;
		Vec2i ret = {-1, -1};
		return ret;
	}
	else
	{
		return pos->getXY();
	}
}


bool same_pt(point_class pt1, point_class pt2)
{
	double eps = 8;//to be set parameter
	double high_eps = 15;
	if (abs(pt1.getX() - pt2.getX()) < 3 && p2pdistance(pt1.getXY(), pt2.getXY()) < high_eps)
		return true;
	if (p2pdistance(pt1.getXY(), pt2.getXY()) <= eps)
		return true;
	else
		return false;
}

bool same_pt(Vec2i pt1, Vec2i pt2, double theta, bool in_line)
{
	double eps = 8;//to be set parametere
	double sin_val = sin(2 * CV_PI * theta / 360.0);
	double adapt_eps = eps / sin_val;
	double plainDis = p2pdistance(pt1, pt2);
	//cout << "plain dis " << plainDis << " ,adapt eps " << adapt_eps << endl;
	if (plainDis <= adapt_eps)
		return true;
	else
		return false;
}

bool same_pt(Vec2i pt1, Vec2i pt2, int adapt_eps)
{
	double plainDis = p2pdistance(pt1, pt2);
//	cout << "plain dis " << plainDis << " ,adapt eps " << adapt_eps << endl;
	if (plainDis <= adapt_eps)
		return true;
	else
		return false;
}

Vec2i avg_pt(vector<Vec2i> pts)
{
	Vec2i new_pt;
	double xsum = 0.0;
	double ysum = 0.0;
	int ptsSize = pts.size();
	for (size_t i = 0; i < ptsSize; ++i)
	{
		xsum += pts[i][0];
		ysum += pts[i][1];
	}
	new_pt = {cvRound(xsum / ptsSize), cvRound(ysum / ptsSize)};
	return new_pt;
}

void computeNewPoints(Vec2i pt1, Vec2i pt2, Vec2i pt3, Vec2i pt4, Vec4i& newpts)
{
	/*choose the leftmost and rightmost points in 4 points coz the line is is other */
	int x[4] = {pt1[0], pt2[0], pt3[0], pt4[0]};
	int y[4] = {pt1[1], pt2[1], pt3[1], pt4[1]};
	int x1, y1, x2, y2;
	int idx1, idx2;
	idx1 = idx2 = 0;
	x1 = x[0];
	x2 = x[0];
	for (int i = 1; i < 4; ++i)
	{
		if (x[i] <= x1)
		{
			x1 = x[i];
			idx1 = i;
		}
		if (x[i] >= x2)
		{
			x2 = x[i];
			idx2 = i;
		}
	}
	newpts = {x[idx1],y[idx1], x[idx2], y[idx2]};
}

void getRangePts(Vec2i pt1, Vec2i pt2, Vec2i pt3, Vec2i pt4, Vec2i& firstPt, Vec2i& fourthPt)
{
	map<int, int> tmpMap;
	tmpMap[pt1[0]] = pt1[1];
	tmpMap[pt2[0]] = pt2[1];
	tmpMap[pt3[0]] = pt3[1];
	tmpMap[pt4[0]] = pt4[1];
	auto iterBegin = tmpMap.begin();
	auto iterEnd = tmpMap.end();
	--iterEnd;
	firstPt = {iterBegin->first, iterBegin->second};
	fourthPt = {iterEnd->first, iterEnd->second};
}

void rm_point_by_id(vector<point_class>& pointxs, int id)
{
	auto pos = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
	                   {
		                   if (a.getPid() == id)
			                   return true;
		                   else
			                   return false;
	                   });
	if (pos == pointxs.end())
		cout << "error" << endl;
	else
	{
		pointxs.erase(pos);
	}
}

/*********circle part********/
inline void getCircle(Point2i p1, Point2i p2, Point2i p3, Point2f& center, float& radius)
{
	float x1 = p1.x;
	float x2 = p2.x;
	float x3 = p3.x;
	float y1 = p1.y;
	float y2 = p2.y;
	float y3 = p3.y;

	center.x = (x1 * x1 + y1 * y1) * (y2 - y3) + (x2 * x2 + y2 * y2) * (y3 - y1) + (x3 * x3 + y3 * y3) * (y1 - y2);
	center.x /= (2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2));

	center.y = (x1 * x1 + y1 * y1) * (x3 - x2) + (x2 * x2 + y2 * y2) * (x1 - x3) + (x3 * x3 + y3 * y3) * (x2 - x1);
	center.y /= (2 * (x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - x3 * y2));

	radius = p2pdistance(Point(p1.x, p1.y), Point(center.x, center.y));
}

float evaluateCircle(Mat dt, Point2f center, float radius)
{
	float completeDistance = 0.0f;
	int counter = 0;

	float maxDist = 1.0f; //TODO: this might depend on the size of the circle!

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
		float cX = radius * cos(t) + center.x;
		float cY = radius * sin(t) + center.y;

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

Vec3i radiusThicknessRev(Vec3f circle, Mat bw, int& width)
{
	Vec3i ret;
	Vec2i center0 = {int(circle[0]), int(circle[1])};
	int radius0 = int(circle[2]);
	//Vec2i center = {}; int radius = 0; 
	//int ptnumCurrent = (getPointPositions(bw)).size();
	//bool breakFlag = false;
	int minPtnum = 10000;
	Vec3i cCandidate = {};
	/*for (int w = 1; w <= 3; w++)
	{
		if (breakFlag)
			break;
		for (int i = center0[0]; i < center0[0] + 2; i++)
		{
			for (int j = center0[1]; j < center0[1] + 2; j++)
			{
				for (int k = radius0; k < radius0 + 2; k++)
				{
					Mat tmp = bw.clone();
					Vec2i tmpCenter = { i, j }; int tmpRadius = k; int tmpW = w;
					cv::circle(tmp, tmpCenter, tmpRadius, 0, tmpW);
					int tmpPtnum = (getPointPositions(tmp)).size();
					if (tmpPtnum < minPtnum)
					{
						minPtnum = tmpPtnum;
						cCandidate = { i, j, k };
						width = w;
						if (ptnumCurrent - minPtnum < 10)
						{
							breakFlag = true;
						}
					}
				}
			}
		}
		ptnumCurrent = minPtnum;
		ret = cCandidate;
	}*/
	for (int i = center0[0]; i < center0[0] + 2; i++)
	{
		for (int j = center0[1]; j < center0[1] + 2; j++)
		{
			for (int k = radius0; k < radius0 + 2; k++)
			{
				Mat tmp = bw.clone();
				Vec2i tmpCenter = {i, j};
				int tmpRadius = k;
				cv::circle(tmp, tmpCenter, tmpRadius, 0, 1);
				int tmpPtnum = (getPointPositions(tmp)).size();
				if (tmpPtnum < minPtnum)
				{
					minPtnum = tmpPtnum;
					cCandidate = {i, j, k};
				}
			}
		}
	}
	int a[3] = {};
	for (int w = 1; w <= 3; w++)
	{
		Mat tmp = bw.clone();
		Vec2i center_loc = {cCandidate[0], cCandidate[1]};
		int radius_loc = cCandidate[2];
		cv::circle(tmp, center_loc, radius_loc, 0, w);
		int tmpPtnum = (getPointPositions(tmp)).size();
		a[w - 1] = tmpPtnum;
		width = w;
		if (w != 1 && (a[w - 1] - a[w]) < 50)
			break;
	}
	ret = cCandidate;
	//int test = (getPointPositions(bw)).size();
	return ret;
}

void detect_circle(const Mat diagram_segment, Mat& color_img, Mat& diagram_segwithoutcircle, Mat& withoutCirBw, vector<circle_class>& circles, bool showFlag)
{
	unsigned int circleN_todetect = 2;
	diagram_segwithoutcircle = diagram_segment.clone();
	//int c_count = 0;
	for (unsigned int i = 0; i < circleN_todetect; ++i)
	{
		vector<Point2i> edgePositions = getPointPositions(diagram_segwithoutcircle);
		Mat dt;
		distanceTransform(255 - diagram_segwithoutcircle, dt, CV_DIST_L1, 3);
		unsigned int nIter = 0;
		Point2f bestCircleCenter = {};
		float bestCircleRadius = -1.0f;
		float bestCVal = -1;
		float minCircleRadius = 0.0f;
		for (unsigned int j = 0; j < 2000; ++j)
		{
			int edgepSize = edgePositions.size();
			unsigned int idx1 = rand() % edgepSize;
			unsigned int idx2 = rand() % edgepSize;
			unsigned int idx3 = rand() % edgepSize;

			if ((idx1 == idx2) || (idx1 == idx3) || (idx2 == idx3))
				continue;
			// create a circle from 3 points
			Point2f center = {};
			float radius = -1.0f;
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
		if (bestCVal < 4 * bestCircleRadius)
			break;
		if (bestCVal < static_cast<int>(edgePositions.size() / 8))
			break;
		if (bestCVal > 125)
		{
			//std::cout << "current best circle: " << bestCircleCenter << " with radius: " << bestCircleRadius << " and nInlier " << bestCVal << endl;

			Vec3f circlef = {bestCircleCenter.x, bestCircleCenter.y, bestCircleRadius};
			int width = 0;
			Vec3i circle = radiusThicknessRev(circlef, diagram_segment, width);
			//int c_id = c_count++;
			circle_class class_circle(circle, i, width);
			circles.push_back(class_circle);
			//draw the cicle detected in red within the colorgeo blob image
			//cv::circle(color_img, bestCircleCenter, bestCircleRadius, Scalar(0, 0, 255));
			//TODO: hold and save the detected circle.

			//TODO: instead of overwriting the mask with a drawn circle it might be better to hold and ignore detected circles and dont count new circles which are too close to the old one.
			// in this current version the chosen radius to overwrite the mask is fixed and might remove parts of other circles too!

			// update mask: remove the detected circle!
			cv::circle(diagram_segwithoutcircle, {circle[0], circle[1]}, circle[2], 0, width); // here the radius is fixed which isnt so nice.
			//edgePointsWithoutCircle = getPointPositions(diagram_segwithoutcircle);
			cv::circle(color_img, {circle[0], circle[1]}, circle[2], Scalar(255, 0, 255), 1);
			cv::circle(withoutCirBw, {circle[0], circle[1]}, circle[2], 0, width);
		}
	}

	if (showFlag)
	{
		namedWindow("4.graygeo blob without circles");
		cv::imshow("4.graygeo blob without circles", diagram_segwithoutcircle);
		namedWindow("5.colorgeo");
		cv::imshow("5.colorgeo", color_img);
	}
}

bool on_circle(Vec2i pt, Vec3f circle)
{
	// check if the point pt is on one of the circles,or joints of multiple circles, or nothing to do with circles
	// joint_flag parameter: 0 means not on any circle, 1 means on a circle, 2 means on two circles and so on
	int dis = 5;// this parameter is to be set to check on the distance tolerants within the distance between radius and distance of pt and circle center point
	//int count = 0;

	Vec2f center = {circle[0], circle[1]};
	double radius = circle[2];
	double distance = p2pdistance(center, pt);
	if (abs(distance - radius) <= dis)
		return true;
	else
		return false;
}

bool on_circle(Vec2f pt, Vec3f circle)
{
	int dis = 5;// this parameter is to be set to check on the distance tolerants within the distance between radius and distance of pt and circle center point
	//int count = 0;

	Vec2f center = {circle[0], circle[1]};
	double radius = circle[2];
	double distance = norm(center, pt);
	if (abs(distance - radius) <= dis)
		return true;
	else
		return false;
}

/**************line part****************/
double cross_product(Vec2i a, Vec2i b)
{
	return a[0] * b[1] - a[1] * b[0];
}

float pt2lineDis(Vec4i line, Vec2i pt)
{
	Vec2i bottom = {line[2] - line[0], line[3] - line[1]};
	Vec2i slope = {pt[0] - line[0], pt[1] - line[1]};
	float p2l = abs(cross_product(bottom, slope) / norm(bottom));
	return p2l;
}

bool on_line(Vec4i line, Vec2i pt)
{
	double linedis_eps = 3;// double pointdis_eps = 5;
	/*Vec2i ref_point = { line[0], line[1] };
	Vec2i pt2 = { -pt[1], pt[0] };
	Vec2i ref_point_t = { line[1], -line[0] };
	Vec2i vec = { line[2] - line[0], line[3] - line[1] };
	double point2line0 = abs((pt2.dot(vec) + ref_point_t.dot(vec))) / sqrt(vec.dot(vec));*/
	//point2line0 is always equal to point2line
	Vec2i bottom = {line[2] - line[0], line[3] - line[1]};
	Vec2i slope = {pt[0] - line[0], pt[1] - line[1]};
	double point2line = abs(cross_product(bottom, slope) / norm(bottom));
	//	cout << "test  " << point2line0 << "," << point2line << "  " << (point2line == point2line0) << endl;
	if (point2line < linedis_eps)
		return true;
	else
		return false;
}

bool in_line(Vec4i line, Vec2i pt)
{
	//if (on_line(line, pt))
	{
		Vec2i pt1, pt2;
		pt1 = {line[0], line[1]};
		pt2 = {line[2], line[3]};
		if (p2pdistance(pt1, pt) + p2pdistance(pt2, pt) - p2pdistance(pt1, pt2) < 10)
			return true;
		else
			return false;
	}
	//else
	//	return false;
}

int on_in_line(Vec4i line, Vec2i pt)
{
	Vec2i pt1, pt2;
	pt1 = {line[0], line[1]};
	pt2 = {line[2], line[3]};
	if (pt[0] - pt1[0] < 5 || pt[0] - pt2[0] < 5)
		return 1;
	if ((pt[0] - pt1[0] + 5) * (pt[0] - pt2[0] + 5) > 0)
	{
		return 0;
	}
	else
		return 2;
}

int ptLine(Vec4i line, Vec2i pt)
{
	//check the relationship between pt and line
	//returns: 0- not on line,1 - on but not in line, 2- same with pt1, 3 - same with pt2, 4 general in line
	Vec2i pt1, pt2;
	pt1 = {line[0], line[1]};
	pt2 = {line[2], line[3]};
	if (same_pt(pt, pt1))
	{
		return 2;
	}
	else if (same_pt(pt, pt2))
	{
		return 3;
	}
	else if (in_line(line, pt))
	{
		//cout << "pt is in line" << endl;
		return 4;
	}
	else if (on_line(line, pt))
		return 1;
	else
		return 0;
}

bool onLinePtInLine(Vec4i line, Vec2i pt)
{
	Vec2i pt1, pt2;
	pt1 = {line[0], line[1]};
	pt2 = {line[2], line[3]};
	int xLeft, xRight;
	if (pt1[0] > pt2[0])
	{
		xLeft = pt2[0];
		xRight = pt1[0];
	}
	else
	{
		xLeft = pt1[0];
		xRight = pt2[0];
	}
	if (pt[0] >= xLeft && pt[0] <= xRight)
		return true;
	else
		return false;
}

bool on_nin_line(Vec4i line, Vec2i pt)
{
	if (on_line(line, pt))
	{
		Vec2i pt1, pt2;
		pt1 = {line[0], line[1]};
		pt2 = {line[2], line[3]};
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
		Vec2i pt1 = {plainLines[i][0], plainLines[i][1]};
		Vec2i pt2 = {plainLines[i][2], plainLines[i][3]};
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

bool crossPtWithinLines(Vec4i line1, Vec4i line2, Vec2i cross)
{
	if (in_line(line1, cross))
	{
		//cross in line1
		if (in_line(line2, cross))
		{
			//cross in line2
			return true;
		}
		else
		{
			cout << "cross not in line2" << endl;
			return false;
		}
	}
	else
	{
		cout << "cross not in line1" << endl;
		return false;
	}
}

inline void line2pt(Vec4i line, Vec2i& pt1, Vec2i& pt2)
{
	pt1 = {line[0], line[1]};
	pt2 = {line[2], line[3]};
}

Vec2i ptAttachToCircle(Vec2f& tmp, vector<circle_class>& circles)
{
	Vec2i cross;
	if (circles.size() != 0)
	{
		for (int j = 0; j < circles.size(); j++)
		{
			Vec3i c = circles[j].getCircleVec();
			Vec2f center = circles[j].getCenter();
			float radius = circles[j].getRadius();
			/*			if (on_circle(pt1, c) || on_circle(pt2, c) || on_circle(pt3, c) || on_circle(pt4, c))
			{
			continue;
			}
			else */
			//To Note: there're changes here
			if (on_circle(tmp, c))
			{
				cout << "cross approximately on circle" << endl;
				float minDiff = 100;
				int cenx = int(tmp[0]);
				int ceny = int(tmp[1]);
				int offset = int(p2pdistance(tmp, center) - radius);
				offset = (offset < 1) ? 1 : offset;
				for (auto m = cenx - offset; m <= cenx + offset; m++)
				{
					for (auto n = ceny - offset; n <= ceny + offset; n++)
					{
						Vec2i temp2 = {m, n};
						float tmpDiff = abs(p2pdistance(temp2, center) - radius);
						if (tmpDiff < minDiff)
						{
							tmp = temp2;
							minDiff = tmpDiff;
						}
					}
				}
			}
		}
		cross = {int(tmp[0]), int(tmp[1])};
		cout << "after attach to circle the cross is " << cross << endl;
	}
	else
	{
		cross = {int(tmp[0]), int(tmp[1])};
		cout << "not circle on image, no attachment" << endl;
	}
	return cross;
}

void getCrossPt(Vec4i line1, Vec4i line2, Vec2f& tmpCross)
{
	int x1 = line1[0];
	int x2 = line1[2];
	int x3 = line2[0];
	int x4 = line2[2];
	int y1 = line1[1];
	int y2 = line1[3];
	int y3 = line2[1];
	int y4 = line2[3];
	tmpCross[0] = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) * 1.0 / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
	tmpCross[1] = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) * 1.0 / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
}

void findLineEnds(vector<Vec2i> colpoints, vector<Vec2i>& lineEnds, vector<Vec2i>& plainPoints, vector<Vec4i> plainLines)
{
	vector<distance_info> distance_infos;
	for (vector<Vec2i>::iterator iter1 = colpoints.begin(); iter1 != colpoints.end(); ++iter1)
	{
		Vec2i pt1 = *iter1;
		for (vector<Vec2i>::iterator iter2 = iter1 + 1; iter2 != colpoints.end(); ++iter2)
		{
			Vec2i pt2 = *iter2;
			double tempDistance = p2pdistance(pt1, pt2);
			distance_info dis_info = {pt1, pt2, tempDistance};
			distance_infos.push_back(dis_info);
		}
	}
	std::sort(distance_infos.begin(), distance_infos.end(),
	          [](distance_info a, distance_info b)
	          {
		          return (a.distance > b.distance);
	          });
	Vec2i pt3 = distance_infos[0].pt1;
	Vec2i pt4 = distance_infos[0].pt2;
	lineEnds.push_back(pt3);
	lineEnds.push_back(pt4);
	int count0 = 0;
	//cout << "colpoints size " << colpoints.size();
	for (vector<Vec2i>::iterator iter3 = colpoints.begin(); iter3 != colpoints.end(); ++iter3)
	{
		Vec2i tempPt1 = *iter3;
		bool flag = on_other_noncollinearlines(plainLines, {pt3[0], pt3[1], pt4[0], pt4[1]}, tempPt1);
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
				++iter4;
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
	else if (pt1[0] > pt2[0])
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

/************concrete function*********/

bool intersection(Vec2i o1, Vec2i p1, Vec2i o2, Vec2i p2,
                  Vec2i& r)
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
	float maxDist = 1.0f;
	int count = 0;
	vector<Point2i> edgePositions;
	edgePositions = getPointPositions(diagram_segwithoutcircle);
	Vec2i pt1 = {rawLine[0], rawLine[1]};
	Vec2i pt2 = {rawLine[2], rawLine[3]};
	int smaller_x = (pt1[0] < pt2[0]) ? pt1[0] : pt2[0];
	int bigger_x = (pt1[0] > pt2[0]) ? pt1[0] : pt2[0];
	int smaller_y = (pt1[1] < pt2[1]) ? pt1[1] : pt2[1];
	int bigger_y = (pt1[1] > pt2[1]) ? pt1[1] : pt2[1];
	for (int i = smaller_x; i <= bigger_x; ++i)
	{
		for (int j = smaller_y; j <= bigger_y; ++j)
		{
			Point2i testPoint = CvPoint(i, j);
			Vec2i tp = {i, j};
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

bool in_rect(Vec2i pt, int leftx, int rightx, int lowy, int highy)
{
	//int eps = 5;
	if (pt[0] >= leftx && pt[0] <= rightx &&
		pt[1] >= lowy && pt[1] <= highy)
		return true;
	else
		return false;
}

Vec4i pt2line(Vec2i pt1, Vec2i pt2)
{
	Vec4i ret = {pt1[0], pt1[1], pt2[0], pt2[1]};
	return ret;
}

bool in_circle(Vec2i center, int radius, Vec2i pt)
{
	if (radius - p2pdistance(center, pt) >= 3)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int ptWithCircle(Vec2i center, int radius, Vec2i pt)
{
	double pt2center = p2pdistance(center, pt);
	if (radius - pt2center > 4)
		return 0;//inside circle
	else if (abs(radius - pt2center) <= 4)
		return 1; // on cirlce
	else if (pt2center - radius < 8)
		return 2;//outside circle
	else
		return 2;
}

bool basicRev(vector<Point2i>& edgePt, Vec2i p1, Vec2i p2, vector<circle_class>& circles)
{
	Vec4i line = pt2line(p1, p2);
	int flag[1000] = {0};
	bool vertical_flag = (abs(p1[0] - p2[0]) < abs(p1[1] - p2[1])) ? true : false;
	int ranges = vertical_flag ? abs(p2[1] - p1[1]) : abs(p2[0] - p1[0]);
	for (auto i = 0; i < edgePt.size(); ++i)
	{
		Vec2i pt = edgePt[i];
		if (in_line(line, pt))
		{
			if (vertical_flag)
				flag[pt[1]] = 1;
			else
				flag[pt[0]] = 1;
		}
	}
	double ratio;
	int nums = count(flag, flag + 1000, 1);
	ratio = nums / ranges;
	cout << ratio * 100 << "%" << endl;
	int threshold_ratio = 0.8;
	if (ratio < threshold_ratio)
		return false;
	else
		return true;
}



bool dashLineRecovery(vector<Point2i>& edgePt, Vec2i col_p1, Vec2i col_p2, vector<circle_class>& circles, bool plflag = false, bool pcflag = false, bool ppflag = false)
{//check if there's dash line between
	Vec4i line = pt2line(col_p1, col_p2);
	int flag[1000] = {0};
	bool vertical_flag = (abs(col_p1[0] - col_p2[0]) < abs(col_p1[1] - col_p2[1])) ? true : false;
	int ranges = vertical_flag ? abs(col_p2[1] - col_p1[1]) : abs(col_p2[0] - col_p1[0]);
	cout << "line check between " << col_p1 << " and " << col_p2 << endl;
	if (circles.size() == 0)
	{
		cout << "non-circle" << endl;
		for (auto i = 0; i < edgePt.size(); ++i)
		{
			Vec2i pt = edgePt[i];
			if (in_line(line, pt))
			{
				if (vertical_flag)
					flag[pt[1]] = 1;
				else
					flag[pt[0]] = 1;
			}
		}
		double ratio;
		int nums = count(flag, flag + 1000, 1);
		ratio = nums / ranges;
		cout << ratio * 100 << "%" << endl;
		int threshold_ratio = 0.8;
		if (ratio < threshold_ratio)
			return false;
		else
			return true;
	}
	else
	{
		for (auto i = 0; i < circles.size(); ++i)
		{
			circle_class* c = &(circles[i]);
			Vec2i center = c->getCenter();
			int radius = c->getRadius();
			double center2lineDis = pt2lineDis(line, center);
			int flag1 = ptWithCircle(center, radius, col_p1);
			int flag2 = ptWithCircle(center, radius, col_p2);
			if (abs(center2lineDis - radius) <= 3)
			{
				if (p2pdistance(col_p1, col_p2) < 20)
				{
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				for (auto j = 0; j < edgePt.size(); ++j)
				{
					Vec2i pt = edgePt[j];
					if (in_line(line, pt))
					{
						if (vertical_flag)
							flag[pt[1]] = 1;
						else
							flag[pt[0]] = 1;
					}
				}
				double ratio;
				int nums = count(flag, flag + 1000, 1);
				ratio = nums / ranges;
				cout << ratio * 100 << "%" << endl;
				double threshold_ratio = 0.8;
				if (ratio < threshold_ratio)
					return false;
				else
					return true;;
			}
		}
	}
}

bool dashLineRecovery2(vector<Point2i>& withoutOnL_ept, Vec2i p_closer, Vec2i p_farther, Vec2i p_cross, vector<circle_class>& circles, bool plflag = false, bool pcflag = false, bool ppflag = false)
{
	//check if there's dash line between
	Vec4i line = pt2line(p_closer, p_cross);
	int flag[1000] = {0};
	bool vertical_flag = (abs(p_closer[0] - p_cross[0]) < abs(p_closer[1] - p_cross[1])) ? true : false;
	int ranges = vertical_flag ? abs(p_cross[1] - p_closer[1]) : abs(p_cross[0] - p_closer[0]);
	cout << "line check between " << p_closer << " and " << p_cross << endl;
	if (circles.size() == 0)
	{
		cout << "image without circle" << endl;
		for (auto i = 0; i < withoutOnL_ept.size(); ++i)
		{
			Vec2i pt = withoutOnL_ept[i];
			if (in_line(line, pt))
			{
				if (vertical_flag)
					flag[pt[1]] = 1;
				else
					flag[pt[0]] = 1;
			}
		}
		double ratio;
		int nums = count(flag, flag + 1000, 1);
		ratio = nums / ranges;
		cout << ratio * 100 << "%" << endl;
		double threshold_ratio = 0.8;
		if (ratio < threshold_ratio)
			return false;
		else 
			return true;
	}
	else
	{
		cout << "image with circle" << endl;
		for (auto i = 0; i < circles.size(); ++i)
		{
			circle_class* c = &(circles[i]);
			Vec2i center = c->getCenter();
			int radius = c->getRadius();
			double center2lineDis = pt2lineDis(line, center);
			int closer_flag = ptWithCircle(center, radius, p_closer);
			int cross_flag = ptWithCircle(center, radius, p_cross);
			double closerp_diff = abs(p2pdistance(p_closer, center) - radius);
			double crossp_diff = abs(p2pdistance(p_cross, center) - radius);
			cout << "closerp diff " << closerp_diff << " crossp diff " << crossp_diff << endl;
			cout << "line and center dis" << abs(center2lineDis - radius) << endl;
			if (abs(center2lineDis - radius) <= 5)
			{
				cout << "the line is almost tangent to the circle" << endl;
				if (closer_flag == 1 && cross_flag == 1)
				{
					// both are on the circle
					cout << "both points are on the circle" << endl;
					if (abs(p2pdistance(p_closer, center) - radius) < abs(p2pdistance(p_cross, center) - radius))
						return false;
					else
						return true;
				}
				else if (closer_flag == 1 && cross_flag == 2)
				{
					// point1 is on the circle
					cout << "point 1 is on the circle and point 2 is outside the circle" << endl;
					cout << "to go" << endl;
					return false;
				}
				else if (cross_flag == 1 && closer_flag == 2)
				{
					cout << "point 2 is on the circle" << endl;
					return true;
				}
				else
				{
					cout << "neither points are on the circle" << endl;
					for (auto j = 0; j < withoutOnL_ept.size(); ++j)
					{
						Vec2i pt = withoutOnL_ept[j];
						if (in_line(line, pt))
						{
							if (vertical_flag)
								flag[pt[1]] = 1;
							else
								flag[pt[0]] = 1;
						}
					}
					double ratio;
					int nums = count(flag, flag + 1000, 1);
					ratio = 1.0 * nums / ranges;
					cout << ratio * 100 << "%" << endl;
					double threshold_ratio = 0.7;
					if (ratio < threshold_ratio)
						return false;
					else
						return true;
				}
			}
			else
			{
				cout << "not tangent to the cirlce" << endl;
				for (auto j = 0; j < withoutOnL_ept.size(); ++j)
				{
					Vec2i pt = withoutOnL_ept[j];
					if (in_line(line, pt))
					{
						if (vertical_flag)
							flag[pt[1]] = 1;
						else
							flag[pt[0]] = 1;
					}
				}
				double ratio;
				int nums = count(flag, flag + 1000, 1);
				ratio = 1.0 * nums / ranges;
				cout << ratio * 100 << "%" << endl;
				double threshold_ratio = 0.7;
				if (ratio < threshold_ratio)
					return false;
				else
					return true;
				//return false;
			}
		}
	}
}

double frac_compute(Vec2i pt1, Vec2i pt2, bool vertical_flag)
{
	double ret = vertical_flag ? 1.0*(pt2[0] - pt1[0]) / (pt2[1] - pt1[1]) : 1.0*(pt2[1] - pt1[1]) / (pt2[0] - pt1[0]);
	return ret;
}

bool dashLineRecovery(vector<Point2i> &withoutOnL_ept, vector<Point2i>& withoutO_ept, Vec2i p_closer, Vec2i p_farther, Vec2i p_cross, vector<circle_class> &circles, bool plflag = false, bool pcflag = false, bool ppflag = false)
{
	//check if there's dash line between
	Vec4i line = pt2line(p_closer, p_cross);
	bool vertical_flag = (abs(p_closer[0] - p_cross[0]) < abs(p_closer[1] - p_cross[1])) ? true : false;
	int ranges = vertical_flag ? abs(p_cross[1] - p_closer[1]) - 1 : abs(p_cross[0] - p_closer[0]) - 1;
	cout << "range: " << ranges << endl;
	bool flag = true;
	if (vertical_flag)
	{
		cout << "measuring in the vertical form" << endl;
		int low_bound, high_bound;
		if (p_closer[1] < p_cross[1])
		{
			low_bound = p_closer[1];
			high_bound = p_cross[1];
		}
		else
		{
			low_bound = p_cross[1];
			high_bound = p_cross[1];
		}
		int step = ranges > 5 ? ranges / 5 : 1;
		int count = 0;
		for (int i = low_bound + 1; i < high_bound; ++i)
		{
			int j = (i - p_closer[1])*frac_compute(p_closer, p_farther, vertical_flag) + p_closer[0];
			Vec2i tmp_pt = { j, i };
			auto iter = find_if(withoutO_ept.begin(), withoutO_ept.end(), [&](Point2i a)
			{
				if (same_pt(a, tmp_pt, 3))
					return true;
				else
					return false;
			});
			if (iter == withoutO_ept.end())
			{
				cout << tmp_pt;
				cout << "       no recovery" << endl;
				return false;
			}
		}
		if (dashLineRecovery2(withoutOnL_ept, p_closer, p_farther, p_cross, circles))
		{
			cout << "recovery" << endl;
			return true;
		}
		else
			return false;
	}
	else
	{
		cout << "measuring in the horizional form" << endl;
		int low_bound, high_bound;
		if (p_closer[0] < p_cross[0])
		{
			low_bound = p_closer[0];
			high_bound = p_cross[0];
		}
		else
		{
			low_bound = p_cross[0];
			high_bound = p_cross[0];
		}
		int step = (ranges - 2) > 5 ? ranges / 5 : 1;
		int count = 0;
		for (int i = low_bound + 1; i < high_bound; ++i)
		{
			int j = (i - p_closer[0])*frac_compute(p_closer, p_farther, vertical_flag) + p_closer[1];
			Vec2i tmp_pt = { i, j };
			auto iter = find_if(withoutO_ept.begin(), withoutO_ept.end(), [&](Point2i a)
			{
				if (same_pt(a, tmp_pt, 3))
					return true;
				else
					return false;
			});
			if (iter == withoutO_ept.end())
			{
				cout << "no recovery" << endl;
				flag = false;
				break;
			}
		}
		if (dashLineRecovery2(withoutOnL_ept, p_closer, p_farther, p_cross, circles))
			return true;
		else
			return false;
	}

}

double lSlope(Vec4i line)
{
	Vec2i pt1, pt2;
	line2pt(line, pt1, pt2);
	Vec2i lineV = {line[2] - line[0], line[3] - line[1]};
	double theta = atan2(lineV[1], lineV[0]);
	return theta;
}

double llAngle(Vec4i line1, Vec4i line2)
{
	if (line1 == line2)
		return 0;
	double theta1, theta2;
	theta1 = lSlope(line1);
	theta2 = lSlope(line2);
	double angle = int(abs(theta1 - theta2) / CV_PI * 180);
	return angle;
}

bool isParallel(Vec4i line1, Vec4i line2)
{
	Vec2i pt1, pt2, pt3, pt4;
	line2pt(line1, pt1, pt2);
	line2pt(line2, pt3, pt4);
	Vec2i line1V = {line1[2] - line1[0], line1[3] - line1[1]};
	Vec2i line2V = {line2[2] - line2[0], line2[3] - line2[1]};
	double theta1, theta2;
	/*if (abs(line1V[0]) < 3)
	{
		theta1 = CV_PI / 2.0;
	}
	else if (abs(line1V[1]) < 3)
	{
		theta1 = 0;
	}
	else
	{
		theta1 = (atan2(line1V[1], line1V[0]));
	}
	if (abs(line2V[0]) < 3)
	{
		theta2 = CV_PI / 2.0;
	}
	else if (abs(line2V[1]) < 3)
	{
		theta2 = 0;
	}
	else
	{
		theta2 = (atan2(line2V[1], line2V[0]));
	}*/
	if (on_line(line1, pt3) && on_line(line1, pt4))
	{
		return true;
	}
	if ((abs(line1V[0]) < 3 && abs(line2V[0]) < 3) | (abs(line1V[1]) < 3 && abs(line2V[1]) < 3))
	{
		return true;
	}
	else
	{
		theta1 = (atan2(line1V[1], line1V[0]));
		theta2 = (atan2(line2V[1], line2V[0]));
	}

	/*double theta1 = abs((abs(line1V[0]) <= 3) ? CV_PI / 2.0 : atan2(line1V[1], line1V[0]));
	double theta2 = abs((abs(line2V[0]) <= 3) ? CV_PI / 2.0 : atan2(line2V[1], line2V[0]));*/
	double angle = llAngle(line1, line2);
	return (angle <= 5) || (angle >= 175);
}

bool ptSortPred(Vec2i pt1, Vec2i pt2)
{
	return (pt1[0] < pt2[0]);
}

bool with_same_line(vector<Vec4i>& plainLines, Vec2i pt1, Vec2i pt2)
{
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i line = plainLines[i];
		Vec2i tmpPt1, tmpPt2;
		tmpPt1 = {line[0], line[1]};
		tmpPt2 = {line[2], line[3]};
		if ((same_pt(pt1, tmpPt1) && same_pt(pt2, tmpPt2)) || (same_pt(pt1, tmpPt2) && (same_pt(pt2, tmpPt1))))
			return true;
		else if ((same_pt(pt1, tmpPt1) && on_line(line, pt2)) || (same_pt(pt2, tmpPt2) && on_line(line, pt1))
			|| (same_pt(pt1, tmpPt2) && on_line(line, pt2)) || (same_pt(pt2, tmpPt1) && on_line(line, pt1)))
			return true;
	}
	return false;
}

bool withinPtCRegion(Vec2i center, Vec2i pt)
{
	if (p2pdistance(center, pt) < 8)
	{
		cout << "within center area" << endl;
		return true;
	}
	else
	{
		cout << "not within center area" << endl;
		return false;
	}
}

bool isInImage(int xMax, int yMax, Vec2i pt)
{
	//the axis is supposed to be within[0,x]*[0,y]
	int x = pt[0];
	int y = pt[1];
	if (x <= 0 || x > xMax)
		return false;
	if (y <= 0 || y > yMax)
		return false;
	return true;
}

int isEndPoint(Vec4i line, Vec2i pt)
{
	Vec2i pt1 = {line[0], line[1]};
	Vec2i pt2 = {line[2], line[3]};
	if (same_pt(pt1, pt))
		return 1;
	else if (same_pt(pt2, pt))
		return 2;
	else
		return 0;
}

void chooseNearCircle(vector<Vec3f>& circle_candidates, Vec2i& cross, Vec2i& pt)
{
	for (auto i = 0; i < circle_candidates.size(); i++)
	{
		Vec3f circle = circle_candidates[i];
		Vec2i center = {int(circle[0]), int(circle[1])};
		if (p2pdistance(pt, center) - circle[2] < 5)
		{
			if (p2pdistance(pt, center) > p2pdistance(cross, center))
			{
				cross = pt;
			}
		}
	}
}

bool isPartOfLongerLine(Vec4i line1, Vec4i line2)
{
	//check if line1 is part of line2;
	Vec2i pt1, pt2;
	pt1 = {line1[0], line1[1]};
	pt2 = {line1[2], line1[3]};
	if (in_line(line2, pt1) && in_line(line2, pt2))
		return true;
	else
		return false;
}

bool existRealLineWithinPtxs(vector<line_class>& lineXs, vector<point_class> pointxs, point_class ptx1, point_class ptx2)
{
	//Vec4i tmpLxy = { ptx1.px, ptx1.py, ptx2.px, ptx2.py };
	//auto iter = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){return a.lxy == tmpLxy; });
	//if (iter != lineXs.end())
	//	return true;
	//else
	//{
	//	//not find exact line
	//	auto iter2 = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){return isPartOfLongerLine(tmpLxy, a.lxy); });
	//	if (iter2 != lineXs.end())
	//		return true;
	//	else
	//	{
	//		//not part of line
	//	}
	//}
	Vec4i line1, line2;
	line1 = {ptx1.getX(), ptx1.getY(), ptx2.getX(), ptx2.getY()};
	line2 = {ptx2.getX(), ptx2.getY(), ptx1.getX(), ptx1.getY()};
	//auto iter1 = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){return in_line(a.lxy, ptx1.pxy); });
	//auto iter2 = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){return in_line(a.lxy, ptx2.pxy); });
	auto iter = find_if(lineXs.begin(), lineXs.end(), [&](line_class a)
	                    {
		                    if (a.getLineVec(pointxs) == line1 || a.getLineVec(pointxs) == line2 || (on_line(a.getLineVec(pointxs), ptx1.getXY()) && on_line(a.getLineVec(pointxs), ptx2.getXY())))
			                    return true;
		                    else
			                    return false;
	                    });
	if (iter != lineXs.end())
		return true;
	else
		return false;
}

bool sameLine(Vec4i line1, Vec4i line2)
{
	Vec2i pt1 = {line1[0], line1[1]};
	Vec2i pt2 = {line1[2], line1[3]};
	Vec2i pt3 = {line2[0], line2[1]};
	Vec2i pt4 = {line2[0], line2[1]};
	bool flag1 = (same_pt(pt1, pt3) && same_pt(pt2, pt4));
	bool flag2 = (same_pt(pt1, pt4) && same_pt(pt2, pt3));
	if (flag1 || flag2)
		return true;
	else
		return false;
}

double slopeEval(Vec4i line)
{
	double diff1, diff2, diff;
	diff1 = lSlope(line) - 0;
	diff2 = lSlope(line) - 90;
	diff = (diff1 < diff2) ? diff1 : diff2;
	return diff;
}

void vecSwapValue(vector<Vec4i>::iterator iter1, vector<Vec4i>::iterator iter2)
{
	auto tmp = *iter1;
	*iter1 = *iter2;
	*iter2 = tmp;
}

bool nearToAcross(Vec2i pt, vector<Vec2i>& crossPts)
{
	auto iter = find_if(crossPts.begin(), crossPts.end(), [&](Vec2i a)
	                    {
		                    if (same_pt(pt, a))
			                    return true;
		                    else
			                    return false;
	                    });
	if (iter != crossPts.end())
		return true;
	else
		return false;
}

int line_recovery_process(line_class* linex, Vec2i p_cross, vector<Point2i>& withoutOnL_ept, vector<Point2i> &withoutO_ept, vector<point_class>& pointxs, vector<circle_class>& circlexs, bool id_change = false)
{
	Vec2i pt1 = linex->getpt1vec(pointxs);
	Vec2i pt2 = linex->getpt2vec(pointxs);
	double dis1 = p2pdistance(pt1, p_cross);
	double dis2 = p2pdistance(pt2, p_cross);
	int pos = 0;
	cout << "dis1 " << dis1 << " " << "dis2 " << dis2 << endl;
	if (dis1 < dis2)
	{
		if (dis1 < 8)
		{
			cout << linex->getpt1vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt1_vec(pointxs, p_cross);
			pos = 1;
		}
		else if (dashLineRecovery(withoutOnL_ept,withoutO_ept, pt1, pt2, p_cross, circlexs))
		{
			cout << linex->getpt1vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt1_vec(pointxs, p_cross);
			pos = 1;
		}
	}
	else
	{
		if (dis2 < 8)
		{
			cout << linex->getpt2vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt2_vec(pointxs, p_cross);
			pos = 2;
		}
		else if (dashLineRecovery(withoutOnL_ept,withoutO_ept, pt2, pt1, p_cross, circlexs))
		{
			cout << linex->getpt2vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt2_vec(pointxs, p_cross);
			pos = 2;
		}
	}
	return pos;
}

void cross_refinement(Vec2f& raw_cross, line_class* lx1, line_class* lx2, vector<circle_class>& circlexs, vector<point_class>& pointxs, vector<Point2i>& withoutOnL_ept, vector<Point2i> &withoutO_ept)
{
	Vec4i linex1_vec = lx1->getLineVec(pointxs);
	Vec4i linex2_vec = lx2->getLineVec(pointxs);
	Vec2i pt1, pt2, pt3, pt4;
	line2pt(linex1_vec, pt1, pt2);
	line2pt(linex2_vec, pt3, pt4);
	int id1, id2, id3, id4;
	id1 = lx1->getPt1Id();
	id2 = lx1->getPt2Id();
	id3 = lx2->getPt1Id();
	id4 = lx2->getPt2Id();
	cout << "pt1 id " << id1 << " pt2 id " << id2 << endl
		<< "pt3 id " << id3 << " pt4 id " << id4 << endl;
	// check according to the relationship between cross and the two lines

	bool in_line1, in_line2;
	in_line1 = in_line(linex1_vec, raw_cross);
	in_line2 = in_line(linex2_vec, raw_cross);
	double angle = llAngle(linex1_vec, linex2_vec);
	//	if (!in_line1 && !in_line2)
	//	{
	//		//if the cross points is not in either lines
	//		cout << "the cross points " << raw_cross << "  is not in either lines" << endl;
	//		int pos1 = line_recovery_process(lx1, raw_cross, edgePt, pointxs, circlexs);
	//		int pos2 = line_recovery_process(lx2, raw_cross, edgePt, pointxs, circlexs);
	//		if (pos1 != -1 && pos2 != -1)
	//		{
	//			if (pos2 == 1 && pos1 == 1)
	//			{
	//				lx2->setpt1Id(lx1->getPt1Id());
	//			}
	//			else if (pos2 == 1 && pos1 == 2)
	//			{
	//				lx2->setpt1Id(lx1->getPt2Id());
	//			}
	//			else if (pos2 == 2 && pos1 == 1)
	//			{
	//				lx2->setpt2Id(lx1->getPt1Id());
	//			}
	//			else if (pos2 == 2 && pos1 == 2)
	//			{
	//				lx2->setpt2Id(lx1->getPt2Id());
	//			}
	//		}
	//	}
	//	else if (in_line1 && in_line2)
	//	{
	//		cout << "the cross of two line" << endl;
	//	}
	//	else if (in_line1 && !in_line2)
	//	{
	//		cout << "the cross " << raw_cross << " is in line1 but not line2" << endl;
	//		line_recovery_process(lx2, raw_cross, edgePt, pointxs, circlexs ,true);
	//		if (same_pt(raw_cross, pt1))
	//		{
	//			
	//		}
	//	}
	//	else if (!in_line1 && in_line2)
	//	{
	//		cout << "the cross "<< raw_cross<< " is in line2 but not line1" << endl;
	//		line_recovery_process(lx1, raw_cross, edgePt, pointxs, circlexs ,true);
	//	}
	if (id1 == id3 || id1 == id4 || id2 == id3 || id2 == id4)
		cout << "already have same id point" << endl;
	else if (same_pt(raw_cross, pt1, angle, in_line1))
	{
		cout << "raw_cross =  pt1" << endl;
		//		cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
		if (same_pt(raw_cross, pt3, angle, in_line2))
		{
			cout << "also, raw_cross = pt3" << endl;
			//			if (same_pt(pt1, pt3))
			//				raw_cross = pt1;
			cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			double cpd1, cpd2, ld1, ld2;
			ld1 = p2pdistance(pt1, pt2);
			ld2 = p2pdistance(pt3, pt4);
			cpd1 = p2pdistance(raw_cross, pt2);
			cpd2 = p2pdistance(raw_cross, pt4);
			//			if (cpd1 >= ld1 && cpd2 >= p2pdistance(pt3,pt4))
			{
				lx2->setPt1_vec(pointxs, raw_cross);
				lx1->setPt1_vec(pointxs, raw_cross);
			}
			//			else
			//			{
			//				lx2->setPt1_vec(pointxs, lx1->getpt1vec(pointxs));
			//			}
			int tmp_id1, tmp_id2;
			tmp_id1 = lx1->getPt1Id();
			tmp_id2 = lx2->getPt1Id();
			if (tmp_id1 < tmp_id2)
			{
				cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
				cout << "point with id " << tmp_id2 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id2);
				lx2->setpt1Id(tmp_id1);
			}
			else
			{
				cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
				cout << "point with id " << tmp_id1 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id1);
				lx1->setpt1Id(tmp_id2);
			}
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else if (same_pt(raw_cross, pt4, angle, in_line2))
		{
			cout << "also, raw_cross = pt4" << endl;
			//			cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			//			if (p2pdistance(raw_cross, pt2) >= p2pdistance(pt1, pt2) && p2pdistance(raw_cross, pt3) >= p2pdistance(pt3, pt4))
			{
				lx2->setPt2_vec(pointxs, raw_cross);
				lx1->setPt1_vec(pointxs, raw_cross);
			}
			//			else
			//			{
			//				lx2->setPt2_vec(pointxs, lx1->getpt1vec(pointxs));
			//			}
			int tmp_id1, tmp_id2;
			tmp_id1 = lx1->getPt1Id();
			tmp_id2 = lx2->getPt2Id();
			if (tmp_id1 < tmp_id2)
			{
				cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
				cout << "point with id " << tmp_id2 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id2);
				lx2->setpt2Id(tmp_id1);
			}
			else
			{
				cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
				cout << "point with id " << tmp_id1 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id1);
				lx1->setpt1Id(tmp_id2);
			}
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else
		{
			if (in_line2)
			{
				cout << "cross in line2 and end of line1" << endl;
			}
			else
			{
				int pos = line_recovery_process(lx2, lx1->getpt1vec(pointxs), withoutOnL_ept,withoutO_ept, pointxs, circlexs);
				if (pos == 0)
				{
					cout << "no recovery" << endl;
				}
				else if (pos == 1)
				{
					cout << "recovery cross link to pt3" << endl;
					int tmp_id1, tmp_id2;
					tmp_id1 = lx1->getPt1Id();
					tmp_id2 = lx2->getPt1Id();
					if (tmp_id1 < tmp_id2)
					{
						cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
						lx2->setpt1Id(tmp_id1);
					}
					else
					{
						cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
						lx1->setpt1Id(tmp_id2);
					}
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
					cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
				}
				else if (pos == 2)
				{
					cout << "recovery cross link to pt4" << endl;
					int tmp_id1, tmp_id2;
					tmp_id1 = lx1->getPt1Id();
					tmp_id2 = lx2->getPt1Id();
					if (tmp_id1 < tmp_id2)
					{
						cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
						lx2->setpt2Id(tmp_id1);
					}
					else
					{
						cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
						lx1->setpt1Id(tmp_id2);
					}
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
					cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
				}
			}
		}
	}
	else if (same_pt(raw_cross, pt2, angle, in_line1))
	{
		cout << "raw_cross =  pt2" << endl;
		//		cout << lx1->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
		//		lx1->setPt2_vec(pointxs, raw_cross);
		if (same_pt(raw_cross, pt3, angle, in_line2))
		{
			cout << "also, raw_cross = pt3" << endl;
			//			cout << lx1->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			//			if (p2pdistance(raw_cross, pt1) >= p2pdistance(pt1, pt2) && p2pdistance(raw_cross, pt4) >= p2pdistance(pt3, pt4))
			{
				lx2->setPt1_vec(pointxs, raw_cross);
				lx1->setPt2_vec(pointxs, raw_cross);
			}
			//			else
			//			{
			//				lx2->setPt1_vec(pointxs, lx1->getpt2vec(pointxs));
			//			}
			int tmp_id1, tmp_id2;
			tmp_id1 = lx1->getPt2Id();
			tmp_id2 = lx2->getPt1Id();
			if (tmp_id1 < tmp_id2)
			{
				cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
				cout << "point with id " << tmp_id2 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id2);
				lx2->setpt1Id(lx1->getPt2Id());
			}
			else
			{
				cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
				cout << "point with id " << tmp_id1 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id1);
				lx1->setpt2Id(tmp_id2);
			}
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else if (same_pt(raw_cross, pt4, angle, in_line2))
		{
			cout << "also, raw_cross = pt4" << endl;
			//			cout << lx1->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			//			if (p2pdistance(raw_cross, pt1) >= p2pdistance(pt1, pt2) && p2pdistance(raw_cross, pt3) >= p2pdistance(pt3, pt4))
			{
				lx2->setPt2_vec(pointxs, raw_cross);
				lx1->setPt2_vec(pointxs, raw_cross);
			}
			//			else
			//			{
			//				lx2->setPt2_vec(pointxs, lx1->getpt2vec(pointxs));
			//				
			//			}
			int tmp_id1, tmp_id2;
			tmp_id1 = lx1->getPt2Id();
			tmp_id2 = lx2->getPt2Id();
			if (tmp_id1 < tmp_id2)
			{
				cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
				cout << "point with id " << tmp_id2 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id2);
				lx2->setpt2Id(tmp_id1);
			}
			else
			{
				cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
				cout << "point with id " << tmp_id1 << " is removed" << endl;
				rm_point_by_id(pointxs, tmp_id1);
				lx1->setpt2Id(tmp_id2);
			}
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else
		{
			if (in_line2)
			{
				cout << "sp" << endl;
			}
			else
			{
				int pos = line_recovery_process(lx2, lx1->getpt2vec(pointxs), withoutOnL_ept,withoutO_ept, pointxs, circlexs);
				if (pos == 0)
				{
					cout << "no recovery" << endl;
				}
				else if (pos == 1)
				{
					cout << "recovery cross link to pt3" << endl;
					int tmp_id1, tmp_id2;
					tmp_id1 = lx1->getPt2Id();
					tmp_id2 = lx2->getPt1Id();
					if (tmp_id1 < tmp_id2)
					{
						cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
						lx2->setpt1Id(tmp_id1);
					}
					else
					{
						cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
						lx1->setpt2Id(tmp_id2);
					}
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
					cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
				}
				else if (pos == 2)
				{
					cout << "recovery cross link to pt4" << endl;
					int tmp_id1, tmp_id2;
					tmp_id1 = lx1->getPt2Id();
					tmp_id2 = lx2->getPt2Id();
					if (tmp_id1 < tmp_id2)
					{
						cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
						lx2->setpt2Id(tmp_id1);
					}
					else
					{
						cout << "id" << tmp_id1 << " id change to id" << tmp_id2 << endl;
						lx1->setpt2Id(tmp_id2);
					}
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
					cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
				}
			}
		}
	}
	else if (same_pt(raw_cross, pt3, angle, in_line1))
	{
		cout << "raw_cross =  pt3" << endl;
		cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
		lx2->setPt1_vec(pointxs, raw_cross);
		if (in_line1)
		{
			cout << "half inner cross" << endl;
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else
		{
			int pos = line_recovery_process(lx1, lx2->getpt1vec(pointxs), withoutOnL_ept,withoutO_ept, pointxs, circlexs);
			if (pos == 0)
			{
				cout << "no recovery" << endl;
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else if (pos == 1)
			{
				//line rev between cross and first point in line	
				cout << "recovery cross link to pt1" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx2->getPt1Id();
				tmp_id2 = lx1->getPt1Id();
				if (tmp_id1 < tmp_id2)
				{
					cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
					lx1->setpt1Id(tmp_id1);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt1Id(tmp_id2);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else if (pos == 2)
			{
				//line rev between cross and second point in line	
				cout << "recovery cross link to pt2" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx2->getPt1Id();
				tmp_id2 = lx1->getPt2Id();
				if (tmp_id1 < tmp_id2)
				{
					cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
					lx1->setpt2Id(tmp_id1);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt2Id(tmp_id2);
				}

				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
		}
	}
	else if (same_pt(raw_cross, pt4, angle, in_line2))
	{
		cout << "raw_cross =  pt4" << endl;
		cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
		lx2->setPt2_vec(pointxs, raw_cross);
		if (in_line1)
		{
			cout << "half inner cross" << endl;
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else
		{
			int pos = line_recovery_process(lx1, lx2->getpt2vec(pointxs), withoutOnL_ept, withoutO_ept, pointxs, circlexs);
			if (pos == 0)
			{
				cout << "no recovery" << endl;
			}
			else if (pos == 1)
			{
				cout << "recovery cross link to pt1" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx2->getPt2Id();
				tmp_id2 = lx1->getPt1Id();
				if (tmp_id1 < tmp_id2)
				{
					cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
					lx1->setpt1Id(tmp_id1);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt2Id(tmp_id2);
				}

				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else if (pos == 2)
			{
				cout << "recovery cross link to pt2" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx2->getPt2Id();
				tmp_id2 = lx1->getPt2Id();
				if (tmp_id1 < tmp_id2)
				{
					cout << "id " << tmp_id2 << " id change to id " << tmp_id1 << endl;
					lx1->setpt2Id(tmp_id1);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt2Id(tmp_id2);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
		}
	}
	else
	{
		//cout << "two cross disjoint line" << endl;
		if (in_line1 && !in_line2)
		{
			cout << "cross in line1 but not in line2" << endl;
			line_recovery_process(lx2, raw_cross, withoutOnL_ept,withoutO_ept, pointxs, circlexs);
		}
		else if (in_line2 && !in_line1)
		{
			cout << "cross in line2 but not in line1" << endl;
			line_recovery_process(lx1, raw_cross, withoutOnL_ept,withoutO_ept, pointxs, circlexs);
		}
		else if (in_line1 && in_line2)
		{
			cout << "inner cross" << endl;
		}
		else
		{
			cout << "outer cross" << endl;
			int pos1 = line_recovery_process(lx1, raw_cross, withoutOnL_ept,withoutO_ept, pointxs, circlexs);
			int pos2 = line_recovery_process(lx2, raw_cross, withoutOnL_ept, withoutO_ept,pointxs, circlexs);
			if (pos1 == 0 || pos2 == 0)
			{
				cout << "no id set" << endl;
			}
			else if (pos1 == 1 && pos2 == 1)
			{
				cout << " point pt1 and point pt3 are to be recovered to the same" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx1->getPt1Id();
				tmp_id2 = lx2->getPt1Id();
				if (tmp_id1 < tmp_id2)
				{
					lx2->setpt1Id(tmp_id1);
				}
				else
				{
					lx1->setpt1Id(tmp_id2);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else if (pos1 == 1 && pos2 == 2)
			{
				cout << " point pt1 and point pt4 are to be recovered to the same" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx1->getPt1Id();
				tmp_id2 = lx2->getPt2Id();
				if (tmp_id1 < tmp_id2)
				{
					lx2->setpt2Id(tmp_id1);
				}
				else
				{
					lx1->setpt1Id(tmp_id2);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else if (pos1 == 2 && pos2 == 1)
			{
				cout << " point pt2 and point pt3 are to be recovered to the same" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx1->getPt2Id();
				tmp_id2 = lx1->getPt1Id();
				if (tmp_id1 < tmp_id2)
				{
					lx2->setpt1Id(tmp_id1);
				}
				else
				{
					lx1->setpt2Id(tmp_id2);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else if (pos1 == 2 && pos2 == 2)
			{
				cout << " point pt2 and point pt4 are to be recovered to the same" << endl;
				int tmp_id1, tmp_id2;
				tmp_id1 = lx1->getPt2Id();
				tmp_id2 = lx2->getPt2Id();
				if (tmp_id1 < tmp_id2)
				{
					lx2->setpt2Id(tmp_id1);
				}
				else
				{
					lx1->setpt2Id(tmp_id2);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
			else
			{
				cout << "pos error" << endl;
			}
		}
	}
	cout << "refinement stop" << endl;
}

void detect_line3(Mat diagram_segwithoutcircle, Mat& withoutCirBw, vector<point_class> pointXs, vector<circle_class>& circles, Mat& color_img, vector<line_class> lineXs, vector<Vec2i>& plainPoints, Mat& drawedImages, bool showFlag = true, string fileName = "")
{
	vector<Point2i> withoutO_ept = getPointPositions(withoutCirBw);
	vector<Vec4i> plainLines = {};

#pragma region raw detection

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
	/*display or write into file*/
	if (showFlag)
	{
		for (size_t j = 0; j < plainLines.size(); ++j)
		{
			Vec4i l = plainLines[j];
			line(color_img, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0, 255, 0), 2, 8);
		}
		namedWindow("7.lines first opt version now 1");
		imshow("7.lines first opt version now 1", color_img);
	}
	else if (fileName != "")
	{
		// write into txt
		ofstream ofile, ofile2;
		ofile.open(fileName, ios_base::app);
		//ofile2.open("t22222.txt", ios_base::app);
		for (size_t k = 0; k < plainLines.size(); ++k)
		{
			Vec4i l = plainLines[k];
			cout << l[0] << "," << l[1] << endl << l[2] << "," << l[3] << endl;
			ofile << l[0] << "," << l[1] << "\n" << l[2] << "," << l[3] << "\n";
			//ofile2 << l[0] << "," << l[1] << "," << l[2] << "," << l[3] << "\n";
		}
		ofile << "\n";
		cout << endl;
		ofile.close();
		//ofile2.close();
	}

#pragma endregion raw detection

	/* then we handle the lines specifically*/
	/* for lines should be combined 1. colinear 2. */

	/*display*/
	//	for (auto i = 0; i < plainLines.size(); i++)
	//	{
	//		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
	//		cout << "*********************" << pt1 << " " << pt2 << endl;
	//		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
	//		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
	//		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
	//		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	//	}
	//namedWindow("5.lines first opt version now", 0); imshow("5.lines first opt version now", color_img);

	cout << "stop and test point" << endl;

#pragma region rmParallel

	for (auto i = 0; i < plainLines.size(); i++)
	{
		// output the rough line infoes here
		cout << plainLines[i][0] << "," << plainLines[i][1] << "," << plainLines[i][2] << "," << plainLines[i][3] << endl;
	}
	// remove parallel lines
	for (auto iter1 = plainLines.begin(); iter1 != plainLines.end(); ++iter1)
	{
		// first combine collinear line
		//int maxXDiff = 0;
		for (auto iter2 = iter1 + 1; iter2 != plainLines.end();)
		{
			Vec4i line1 = *iter1;
			Vec2i pt1 = {line1[0], line1[1]};
			Vec2i pt2 = {line1[2], line1[3]};
			cout << pt1 << "pt" << pt2 << endl;

			Vec4i line2 = *iter2;
			Vec2i pt3 = {line2[0], line2[1]};
			Vec2i pt4 = {line2[2], line2[3]};
			cout << pt3 << "pt" << pt4 << endl;

			if (isParallel(line1, line2))
			{
				//if it's parallel, followed by checking if collinear
				cout << "line1 and line2 is parallel" << endl;
				vector<Vec2i> tmpPtVec;
				tmpPtVec.push_back(pt1);
				tmpPtVec.push_back(pt2);
				tmpPtVec.push_back(pt3);
				tmpPtVec.push_back(pt4);
				bool flag0 = false;

				//flag0 indicates whether line 1 is vertical, if vertical sort the tmpvec by y axis info, otherwise sort it  by x axis
				if ((flag0 = (abs(pt1[0] - pt2[0]) < 3)))
				{
					sort(tmpPtVec.begin(), tmpPtVec.end(), [](Vec2i a, Vec2i b) { return a[1] < b[1]; });
					cout << "line1 is vertical" << endl;
				}
				else
				{
					sort(tmpPtVec.begin(), tmpPtVec.end(), ptSortPred);
					cout << "line1 is not vertical" << endl;
				}
				Vec2i firstPt = tmpPtVec[0];
				Vec2i fourthPt = tmpPtVec[3];
				Vec2i secondPt = tmpPtVec[1];
				Vec2i thirdPt = tmpPtVec[2];
				//check;

				if (ptLine(line1, pt3)) //pt3 of line2 is on line 1
				//||sameLine(line1,line2)) 
				//TO Note: change here
				{
					// if it's collinear
					cout << "two line is collinear" << endl;
					int eps = 5;
					bool flag3, flag4;
					if (!flag0)
					{
						// line1 not vertical
						flag3 = ((pt3[0] - pt1[0]) * (pt3[0] - pt2[0]) > 0); //appoximate if pt3 is in line1 or on line1 but not in it
						flag4 = ((pt4[0] - pt1[0]) * (pt4[0] - pt2[0]) > 0);//approximate if pt4 is in  line1
					}
					else
					{
						flag3 = ((pt3[1] - pt1[1]) * (pt3[1] - pt2[1]) > 0); //appoximate if pt3 is in line1
						flag4 = ((pt4[1] - pt1[1]) * (pt4[1] - pt2[1]) > 0);//approximate if pt4 is in  line1					
					}
					// four points are all different
					if (flag3 && flag4 && !same_pt(pt1, pt3) && !same_pt(pt1, pt4) && !same_pt(pt2, pt3) && !same_pt(pt2, pt4))
					{
						// two disjoint collinear line
						cout << "four point are all different,";
						cout << firstPt << "range" << fourthPt << endl;
						if (dashLineRecovery(withoutO_ept, secondPt, thirdPt, circles, true, false, false))
						{
							*iter1 = {firstPt[0], firstPt[1], fourthPt[0], fourthPt[1]};
							cout << "collinear line recovery, now the line1 is " << *iter1 << "and erase line2" << endl << endl;
							iter2 = plainLines.erase(iter2);
						}
						else
						{
							cout << "collinear but two differnent line" << endl << endl;
							++iter2;
						}
						//iter2++;
					}
					else
					{
						cout << "colliner line cross over" << endl;
						cout << firstPt << "range" << fourthPt << endl;
						*iter1 = {firstPt[0], firstPt[1], fourthPt[0], fourthPt[1]};
						cout << "now the line1 is " << *iter1 << endl << endl;
						iter2 = plainLines.erase(iter2);
					}
				}
				else
				{
					cout << "parallel but not collinear" << endl << endl;
					++iter2;
				}
			}
			else
			{
				cout << "line1 and line2 is not parallel" << endl << endl;
				++iter2;
			}
		}
		cout << endl;
	}


#pragma endregion rmParallel

	Mat rmParaI = color_img.clone();
	/*display*/
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i];
		Vec2i pt1 = {l[0], l[1]};
		Vec2i pt2 = {l[2], l[3]};
		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(rmParaI, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(rmParaI, Point{pt1[0], pt1[1]}, 10, tmp);
		circle(rmParaI, Point{pt2[0], pt2[1]}, 10, tmp);
	}
	namedWindow("rm parallel", 0);
	imshow("rm parallel", rmParaI);

	// remove the line detected and store it as withoutCLOriBw

	/*************cover the detected line with white pixel*/
	Mat withoutCLOriBw = withoutCirBw.clone();
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i pl = plainLines[i];
		Vec2i pt1 = {pl[0], pl[1]};
		Vec2i pt2 = {pl[2], pl[3]};
		cv::line(withoutCLOriBw, pt1, pt2, 0, 3, 8, 0);
	}
	vector<Point2i> withoutOnL_ept = getPointPositions(withoutCLOriBw);

	bool testflag = false;
	/*write test block*/
	//	ofstream test; test.open("test/test.txt");
	//	for (auto i = 0; i < plainLines.size(); i++)
	//	{
	//		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
	//		cout << pt1 << " " << pt2 << endl;
	//		if (testflag)
	//			test << pt1[0] <<" "<<pt1[1]<<" "<< pt2[0] <<" "<<pt2[1] << endl;
	//	}
	//	test.close();

	cout << endl << "***********block stop" << endl << endl;

	/****************sort the line and mod the horionzal and vertical line*/
	vector<Vec4i> tmpLines;
	auto ordered_iter = plainLines.begin();
	for (auto iter = plainLines.begin() + 1; iter != plainLines.end(); ++iter)
	{
		Vec4i line = *iter;
		Vec2i pt1, pt2;
		line2pt(line, pt1, pt2);
		if (abs(pt1[0] - pt2[0]) < 3)
		{
			//the line is vertical
			(*iter)[0] = (*iter)[2];
			vecSwapValue(iter, ordered_iter);
			++ordered_iter;
		}
		else if (abs(pt1[1] - pt2[1]) < 3)
		{
			//the line is horizonal
			(*iter)[1] = (*iter)[3];
			vecSwapValue(iter, ordered_iter);
			++ordered_iter;
		}
	}

	int px_count = 0;
	int lx_count = 0;
	vector<line_class> linexs;
	vector<point_class> pointxs;

	/***********initialize lines and points**********/
	for (auto i = 0; i < plainLines.size(); i++)
	{
		auto p_l = plainLines[i];
		Vec2i pt1, pt2;
		line2pt(p_l, pt1, pt2);
		auto _pidx1 = px_count++;
		int _pidx2 = px_count++;
		int _lidx = lx_count++;
		point_class* ptx1 = new point_class(pt1, _pidx1);
		point_class* ptx2 = new point_class(pt2, _pidx2);
		//cout << &ptx1 << endl<<&ptx2<<endl;
		ptx1->pushLid(_lidx);
		ptx2->pushLid(_lidx);
		pointxs.push_back(*ptx1);
		pointxs.push_back(*ptx2);

		line_class* lx = new line_class(ptx1->getPid(), ptx2->getPid(), _lidx);
		//lx.p1 = &pointxs[_pidx1]; lx.p2 = &pointxs[_pidx2];
		linexs.push_back(*lx);
		delete ptx1 , ptx2 , lx;
	}

	/*display */
	//	for (auto i = 0; i < linexs.size(); i++)
	//	{
	//		Vec4i l = linexs[i].getLineVec(pointxs);
	//		Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
	//		//cout << "*********************" << pt1 << " " << pt2 << endl;
	//		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
	//		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
	//		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
	//		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	//	}
	//namedWindow("7.lines first opt version now", 0); imshow("7.lines first opt version now", color_img);

	/*print line axis and id infos*/
	for (auto i = 0; i < linexs.size(); i++)
	{
		Vec4i l = linexs[i].getLineVec(pointxs);
		Vec2i pt1 = {l[0], l[1]};
		Vec2i pt2 = {l[2], l[3]};
		cout << pt1 << " " << pt2 << endl;
		cout << linexs[i].getPt1Id() << " " << linexs[i].getPt2Id() << endl;
	}
	cout << "****************sep" << endl << endl;

	// with the raw detection result, we may get extra false line or omit some line(eg. dash line falsely removed previously or real line not detected)
	// in order to detect all lines, we suppose to first put all candidates into calculation then remove the false line at the last step.
	/*******************raw line candidates refinement*/
	map<int, int> change000;

	for (auto i = 0; i < linexs.size(); i++)
	{
		line_class* linex1 = &linexs[i];
		cout << endl << linex1->getLineVec(pointxs) << endl;
		for (auto j = i + 1; j < linexs.size(); j++)
		{
			line_class* linex2 = &linexs[j];
			Vec4i linex1_vec = linex1->getLineVec(pointxs);
			Vec4i linex2_vec = linex2->getLineVec(pointxs);
			//check whether the two line are parallel, if so, jump ahead
			if (isParallel(linex1_vec, linex2_vec))
			{
				//cout << linex1_vec << endl << linex2_vec << endl;
				cout << "the two line is parallel but not collinear" << endl;
				continue;
			}
			else
			{
				// not parallel, then calculate the rough cross first
				Vec2i pt1, pt2, pt3, pt4;
				Vec2f raw_cross;
				line2pt(linex1_vec, pt1, pt2);
				line2pt(linex2_vec, pt3, pt4);
				getCrossPt(linex1_vec, linex2_vec, raw_cross);
				cout << endl << "raw cross" << raw_cross << endl;

				if (!isInImage(diagram_segwithoutcircle.cols, diagram_segwithoutcircle.rows, raw_cross))
				{
					cout << raw_cross << endl;
					cout << "cross out of scope, then take it as no cross" << endl;
					continue;
				}
				else
				{
					//cross in scope
					cout << linex1_vec << endl << linex2_vec << endl;
					cross_refinement(raw_cross, linex1, linex2, circles, pointxs, withoutOnL_ept, withoutO_ept);
					//					cout << linex1->getLineVec(pointxs)<< linex1->getPt1Id() << "," << linex1->getPt2Id() << endl;
					//					cout << linex2->getLineVec(pointxs) << linex2->getPt1Id() << "," << linex2->getPt2Id() << endl;
				}
			}
		}
	}
	cout << endl;
	for (auto i = 0; i < linexs.size(); i++)
	{
		cout << linexs[i].getLineVec(pointxs) << linexs[i].getPt1Id() << ", " << linexs[i].getPt2Id() << endl;
	}
	for (auto j = 0; j < pointxs.size(); j++)
	{
		cout << pointxs[j].getPid() << "  " << pointxs[j].getXY() << endl;
	}
	//rm isolated short lines due to the non-complete circle removal
	// 1. two points almost on circles 2. the length is short 3. the end points is ont on other lines
	int rm_pt_num = 0;
	int rm_line_num = 0;
	map<int, int> changeMap0;
	for (auto iter = linexs.begin(); iter != linexs.end();)
	{
		Vec2i pt1, pt2;
		line2pt(iter->getLineVec(pointxs), pt1, pt2);
		bool rm_flag = false;
		for (auto i = 0; i < circles.size(); ++i)
		{
			Vec3i cir = circles[i].getCircleVec();
			if (on_circle(pt1, cir) && on_circle(pt2, cir))
			{
				double len = p2pdistance(pt1, pt2);
				if (len < 20)
				{
					auto iter1 = find_if(linexs.begin(), linexs.end(), [&](line_class a)
					                     {
						                     if (a.getLineVec(pointxs) != iter->getLineVec(pointxs))
						                     {
							                     if (a.getpt1vec(pointxs) == pt1 || a.getpt2vec(pointxs) == pt1)
								                     return true;
							                     else
								                     return false;
						                     }
						                     else
							                     return false;
					                     });
					auto iter2 = find_if(linexs.begin(), linexs.end(), [&](line_class a)
					                     {
						                     if (a.getLineVec(pointxs) != iter->getLineVec(pointxs))
						                     {
							                     if (a.getpt1vec(pointxs) == pt2 || a.getpt2vec(pointxs) == pt2)
								                     return true;
							                     else
								                     return false;
						                     }
						                     else
							                     return false;
					                     });
					bool not_find_pt1_flag = (iter1 == linexs.end()) ? true : false;
					bool not_find_pt2_flag = (iter2 == linexs.end()) ? true : false;
					if (not_find_pt1_flag && !not_find_pt2_flag)
					{
						iter = linexs.erase(iter);
						rm_flag = true;
						//						pointxs.erase(pointxs.begin() + iter->getPt1Id() - rm_pt_num);
						//						rm_pt_num++;
						rm_line_num++;
					}
					else if (not_find_pt2_flag && !not_find_pt1_flag)
					{
						iter = linexs.erase(iter);
						rm_flag = true;
						//						pointxs.erase(pointxs.begin() + iter->getPt2Id() - rm_pt_num);
						//						rm_pt_num++;
						rm_line_num++;
					}
					else if (not_find_pt1_flag && not_find_pt2_flag)
					{
						iter = linexs.erase(iter);
						rm_flag = true;
						//						pointxs.erase(pointxs.begin() + iter->getPt1Id() - rm_pt_num);
						//						rm_pt_num++;
						pointxs.erase(pointxs.begin() + iter->getPt2Id() - rm_pt_num);
						//						rm_pt_num++; rm_line_num++;
					}
				}
			}
		}
		if (!rm_flag)
		{
			++iter;
			changeMap0[iter - linexs.begin() + rm_line_num] = iter - linexs.begin();
		}
	}
	for (auto i = 0; i < linexs.size(); i++)
	{
		cout << linexs[i].getLineVec(pointxs) << linexs[i].getPt1Id() << ", " << linexs[i].getPt2Id() << endl;
	}
	for (auto j = 0; j < pointxs.size(); j++)
	{
		cout << pointxs[j].getPid() << "  " << pointxs[j].getXY() << endl;
	}

	/*set similar points identical*/


	/*remove indentical vec points and reorder the indices*/
	// point index indicate the line id in which it locate
	map<int, int> changeMap;
	for (auto i = 0; i < pointxs.size(); ++i)
	{
		//the location of the first pointx in the loop
		Vec2i px1 = pointxs[i].getXY();
		int erase_offset = 0; // the var log the erased points num and offset the indexing digit when indexing
		for (auto j = i + 1; j < pointxs.size(); ++j)
		{
			Vec2i px2 = pointxs[j].getXY();
			cout << pointxs[i].getXY() << "   " << pointxs[j].getXY() << endl;
			if (same_pt(px1, px2))
			{
				//the two points vector is the same, then erase the second point from the pointxs set
				// and ajust the respective line with previous erased point
				//				cout << "find the identical vec points " << px1 << " <- "<<pointxs[i].getPid() <<", "<< pointxs[j].getPid()<< endl;
				changeMap[pointxs[j].getPid()] = pointxs[i].getPid();
				int line_id = pointxs[j].getPid() / 2;
				int pos = pointxs[j].getPid() % 2;// if 0, means the first point in the line, otherwise 1, means the second point in the line
				pointxs.erase(pointxs.begin() + j);
				if (pos)
				{
					linexs[changeMap0[line_id]].setpt2Id(pointxs[i].getPid());
				}
				else
				{
					linexs[changeMap0[line_id]].setpt1Id(pointxs[i].getPid());
				}
				j--;
			}
		}
	}
	for (auto i = 0; i < 2 * linexs.size(); ++i)
	{
		if (changeMap[i] == 0)
		{
			changeMap[i] = i;
		}
	}

	//reordering index
	map<int, int> changeMap2;
	for (auto i = 0; i < pointxs.size(); i++)
	{
		//		cout << "point id " << pointxs[i].getPid() << pointxs[i].getXY() << endl;
		changeMap2[pointxs[i].getPid()] = i;
		pointxs[i].setPid(i);
	}
	for (auto j = 0; j < linexs.size(); j++)
	{
		cout << linexs[j].getPt1Id() << "->" << changeMap2[changeMap[linexs[j].getPt1Id()]] << endl;
		cout << linexs[j].getPt2Id() << "->" << changeMap2[changeMap[linexs[j].getPt2Id()]] << endl;
		linexs[j].setpt1Id(changeMap2[changeMap[linexs[j].getPt1Id()]]);
		linexs[j].setpt2Id(changeMap2[changeMap[linexs[j].getPt2Id()]]);
	}
	for (auto i = 0; i < linexs.size(); i++)
	{
		cout << linexs[i].getLineVec(pointxs) << linexs[i].getPt1Id() << ", " << linexs[i].getPt2Id() << endl;
	}
	for (auto j = 0; j < pointxs.size(); j++)
	{
		cout << pointxs[j].getPid() << "  " << pointxs[j].getXY() << endl;
	}

	/*display*/
	for (auto i = 0; i < linexs.size(); i++)
	{
		Vec4i l = linexs[i].getLineVec(pointxs);
		Vec2i pt1 = {l[0], l[1]};
		Vec2i pt2 = {l[2], l[3]};
		//		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{pt1[0], pt1[1]}, 10, tmp);
		circle(color_img, Point{pt2[0], pt2[1]}, 10, tmp);
	}
	cout << endl;


	//if (showFlag)
	{
		namedWindow("8.lines first opt version now", 0);
		imshow("8.lines first opt version now", color_img);
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
	namedWindow("dilated image");
	imshow("dilated image", dilated);
	//namedWindow("eroded diagram"); imshow("eroded diagram", eroded);
	morphologyEx(diagram_segment, closed, MORPH_CLOSE, eroEle);
	namedWindow("closed diagram");
	imshow("closed diagram", closed);
	return closed;
}


void primitive_parse(const Mat binarized_image, const Mat diagram_segment, vector<Point2i>& edgePoints, vector<point_class>& points, vector<line_class>& lines, vector<circle_class>& circles, Mat& drawedImages, bool showFlag = true, string fileName = "")
{
	/* primitive about points, lines, and circles
	first detect the circle, we can get the matrix of diagram segments without circle for the next
	 processing step*/
	//vector<Vec3i> circle_candidates = {}; 
	Mat color_img = Mat::zeros(diagram_segment.size(),CV_8UC3);
	cvtColor(diagram_segment, color_img, CV_GRAY2RGB);
	Mat diagram_segwithoutcircle = Mat::zeros(diagram_segment.size(), CV_8UC1);

	//diagram_segment = preprocessing(diagram_segment);

	/*ransac go*/
	Mat withoutCirBw = binarized_image.clone();
	// detect the circle and get the img without circle and bw img without cirlce and circle candidates
	detect_circle(diagram_segment, color_img, diagram_segwithoutcircle, withoutCirBw, circles, showFlag);
	// then the line detection
	vector<Vec2i> basicEndpoints = {};


	detect_line3(diagram_segwithoutcircle, withoutCirBw, points, circles, color_img, lines, basicEndpoints, drawedImages, showFlag, fileName);


	/*display point text info*/
	//for (int i = 0; i < points.size(); ++i)
	//{
	//	point_class point = points[i];
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
	int labeln;
	Mat diagram_segment = Mat::zeros(image.size(), CV_8UC1);
	vector<Mat> label_segment = {};
	vector<Point2i> oriEdgePoints = getPointPositions(binarized_image);

	image_labelling(binarized_image, diagram_segment, true);
	vector<point_class> points = {};
	vector<line_class> lines = {};
	vector<circle_class> circles = {};
	Mat drawedImages(image.size(), CV_8UC3);

	primitive_parse(binarized_image, diagram_segment, oriEdgePoints, points, lines, circles, drawedImages, false);
	return 0;
}

int diagram()
{
	//a series of image
	//vector<Mat> images;
	char abs_path[100] = "D:\\data\\graph-DB\\newtest33";
	char imageName[150], saveimgName[150];
	//string outputFN = "D:\\data\\graph-DB\\newtest6\\output.txt";
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
		vector<Mat> label_segment = {};
		Mat diagram_segment = Mat::zeros(binarized_image.size(), CV_8UC1);
		vector<Point2i> edgePoints = getPointPositions(binarized_image);
		/*Mat pointss = Mat::zeros(1000, 1000, CV_8UC3);
		for (auto i = 0; i < edgePoints.size(); i++)
		{
			Point2i pt = edgePoints[i];
			circle(pointss, pt, 1, Scalar(0, 0, 255));
		}
		namedWindow("points"); imshow("points", pointss);*/
		image_labelling(binarized_image, diagram_segment);
		vector<point_class> points = {};
		vector<line_class> lines = {};
		vector<circle_class> circles = {};
		Mat drawedImages = Mat::zeros(diagram_segment.size(), CV_8UC3);
		primitive_parse(binarized_image, diagram_segment, edgePoints, points, lines, circles, drawedImages, false);
		imwrite(saveimgName, drawedImages);
	}
	return 0;
}
