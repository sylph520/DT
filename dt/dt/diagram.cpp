# include "stdafx.h"
# include "diagram.h"

#define N 500

//vector<Point2i> edgePointsWithoutCircle;
int Sgn(double d)
{
	if (d<0)
		return -1;
	else
		return 1;
}
//point_class pt2x(Vec2i pt, int point_idx = -1, vector<int> circle_idxs = {}, vector<int> line_idxs = {}, string label = "")
//{
//	point_class ptx;
//	ptx.p_idx = point_idx; ptx.pxy = pt; ptx.px = pt[0]; ptx.py = pt[1];
//	ptx.c_idxs = circle_idxs; ptx.l_idxs = line_idxs;
//	ptx.label = label;
//	return ptx;
//}
//line_class line2x(point_class* p1, point_class* p2, int _l_idx = -1, string _label = "")
//{
//	line_class lx;
//	lx.l_idx = _l_idx; 
//	lx.p1 = p1; lx.p2 = p2;
//	lx.label = _label;
//	return lx;
//}
//circleX circle2x(Vec3i c, int c_idx = -1, int center_pid = -1, string label = "", int contour_width = 0, vector<int> p_idxs = {})
//{
//	circleX crx; crx.Circle = c;
//	crx.center = { c[0], c[1] }; crx.radius = c[2];
//	crx.center_pid = center_pid; crx.label = label; crx.contour_width = contour_width;
//	crx.cx = c[0]; crx.cy = c[1];
//	return crx;
//}

bool dashLineRecovery(vector<Point2i> &, Vec2i, Vec2i, Vec2i ,vector<circle_class> &, bool plflag, bool pcflag, bool ppflag);
bool dashLineRecovery(vector<Point2i> &, Vec2i, Vec2i, vector<circle_class> &, bool plflag, bool pcflag, bool ppflag);


/*******************************load and segment phase***********************/
Mat image_binarizing(Mat input_image, bool showFlag=false)
{
	/* This function is used for transform image to its binarized image */
	
	int block_size; double c;
	block_size = 13; c = 20;
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

void image_labelling(Mat binarized_image,Mat &diagram_segment,bool showFlag = false)
{
	// this function is used to label image with connectedcomponnent analysis, and store the diagram 
	//and label segment
	Mat labeled(binarized_image.size(), CV_8UC3);
	Mat statsMat, centroidMat; Mat labeled_image; vector<Mat> segments;
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
	for (int r = 0; r < binarized_image.rows; ++r) {
		for (int c = 0; c < binarized_image.cols; ++c) {
			int label = labeled_image.at<int>(r, c);
			Vec3b &pixel = labeled.at<Vec3b>(r, c);
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
	double high_eps = 15;
	double low_eps = 5;
	if (abs(pt1[0] - pt2[0]) < 3 && p2pdistance(pt1, pt2) < low_eps)
		return true;
//	if (abs(pt1[1] - pt2[1]) < 3 && p2pdistance(pt1, pt2) < high_eps)
//		return true;
	if (p2pdistance(pt1, pt2) <= eps)
		return true;
	else
		return false;
}

bool same_pt(point_class pt1, point_class pt2)
{
	double eps = 8;//to be set parameter
	double high_eps = 15;
	if (abs(pt1.getX() - pt2.getX()) < 3 && p2pdistance(pt1.getXY(),pt2.getXY()) < high_eps)
		return true;
	if (p2pdistance(pt1.getXY(), pt2.getXY()) <= eps)
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

void getRangePts(Vec2i pt1, Vec2i pt2, Vec2i pt3, Vec2i pt4, Vec2i &firstPt, Vec2i &fourthPt)
{
	map<int, int> tmpMap;
	tmpMap[pt1[0]] = pt1[1];tmpMap[pt2[0]] = pt2[1];
	tmpMap[pt3[0]] = pt3[1];tmpMap[pt4[0]] = pt4[1];
	auto iterBegin = tmpMap.begin(); auto iterEnd = tmpMap.end();
	--iterEnd;
	firstPt = { iterBegin->first, iterBegin->second }; fourthPt = { iterEnd->first, iterEnd->second };
}


/*********circle part********/
inline void getCircle(Point2i p1, Point2i p2, Point2i p3, Point2f &center, float &radius)
{
	float x1 = p1.x;float x2 = p2.x;float x3 = p3.x;
	float y1 = p1.y;float y2 = p2.y;float y3 = p3.y;

	center.x = (x1*x1 + y1*y1)*(y2 - y3) + (x2*x2 + y2*y2)*(y3 - y1) + (x3*x3 + y3*y3)*(y1 - y2);
	center.x /= (2 * (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2));

	center.y = (x1*x1 + y1*y1)*(x3 - x2) + (x2*x2 + y2*y2)*(x1 - x3) + (x3*x3 + y3*y3)*(x2 - x1);
	center.y /= (2 * (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2));

	radius = p2pdistance(Point(p1.x, p1.y), Point(center.x, center.y));
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
Vec3i radiusThicknessRev(Vec3f circle, Mat bw, int &width)
{
	Vec3i ret;
	Vec2i center0 = { int(circle[0]), int(circle[1]) };
	int radius0 = int(circle[2]);
	//Vec2i center = {}; int radius = 0; 
	//int ptnumCurrent = (getPointPositions(bw)).size();
	//bool breakFlag = false;
	int minPtnum = 10000; Vec3i cCandidate = {}; 
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
				Vec2i tmpCenter = { i, j }; int tmpRadius = k; 
				cv::circle(tmp, tmpCenter, tmpRadius, 0, 1);
				int tmpPtnum = (getPointPositions(tmp)).size();
				if (tmpPtnum < minPtnum)
				{
					minPtnum = tmpPtnum;
					cCandidate = { i, j, k };
				}
			}
		}
	}
	int a[3] = {};
	for (int w = 1; w <= 3; w++)
	{
		Mat tmp = bw.clone();
		Vec2i center_loc = { cCandidate[0], cCandidate[1] }; int radius_loc = cCandidate[2]; 
		cv::circle(tmp, center_loc, radius_loc, 0, w);
		int tmpPtnum = (getPointPositions(tmp)).size();
		a[w-1] = tmpPtnum;
		width = w;
		if (w != 1 && (a[w - 1] - a[w]) < 50)
			break;
		
	}
	ret = cCandidate;
	//int test = (getPointPositions(bw)).size();
	return ret;
}
void detect_circle(const Mat diagram_segment, Mat &color_img,Mat &diagram_segwithoutcircle,Mat &withoutCirBw, vector<circle_class> &circles,  bool showFlag)
{
	unsigned int circleN_todetect = 2;
	diagram_segwithoutcircle = diagram_segment.clone();
	//int c_count = 0;
	for (unsigned int i = 0; i < circleN_todetect; ++i)
	{
		vector<Point2i> edgePositions = getPointPositions(diagram_segwithoutcircle);
		Mat dt; distanceTransform(255 - diagram_segwithoutcircle , dt, CV_DIST_L1, 3);
		unsigned int nIter = 0;
		Point2f bestCircleCenter = {}; float bestCircleRadius = -1.0f;
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
			Point2f center = {}; float radius = -1.0f;
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
		if (bestCVal<static_cast<int>(edgePositions.size() / 8)) 
			break;
		if (bestCVal > 125)
		{
			//std::cout << "current best circle: " << bestCircleCenter << " with radius: " << bestCircleRadius << " and nInlier " << bestCVal << endl;
			
			Vec3f circlef = { bestCircleCenter.x, bestCircleCenter.y, bestCircleRadius };
			int width = 0;
			Vec3i circle = radiusThicknessRev(circlef, diagram_segment,width);
			//int c_id = c_count++;
			circle_class class_circle(circle,i , width);
			circles.push_back(class_circle);
			//draw the cicle detected in red within the colorgeo blob image
			//cv::circle(color_img, bestCircleCenter, bestCircleRadius, Scalar(0, 0, 255));
			//TODO: hold and save the detected circle.

			//TODO: instead of overwriting the mask with a drawn circle it might be better to hold and ignore detected circles and dont count new circles which are too close to the old one.
			// in this current version the chosen radius to overwrite the mask is fixed and might remove parts of other circles too!

			// update mask: remove the detected circle!
			cv::circle(diagram_segwithoutcircle, { circle[0], circle[1] }, circle[2], 0, width); // here the radius is fixed which isnt so nice.
			//edgePointsWithoutCircle = getPointPositions(diagram_segwithoutcircle);
			cv::circle(color_img, { circle[0], circle[1] }, circle[2], Scalar(255, 0, 255), 1);
			cv::circle(withoutCirBw, { circle[0], circle[1] }, circle[2], 0, width);
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
	//int count = 0;

	Vec2f center = { circle[0], circle[1] };
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

	Vec2f center = { circle[0], circle[1] };
	double radius = circle[2];
	double distance = norm(center, pt);
	if (abs(distance - radius) <= dis)
		return true;
	else
		return false;
}

/**************line part****************/
double cross_product(Vec2i a, Vec2i b){
	return a[0] * b[1] - a[1] * b[0];
}

float pt2lineDis(Vec4i line, Vec2i pt)
{
	Vec2i bottom = { line[2] - line[0], line[3] - line[1] };
	Vec2i slope = { pt[0] - line[0], pt[1] - line[1] };
	float p2l= abs(cross_product(bottom, slope) / norm(bottom));
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
	Vec2i bottom = { line[2] - line[0], line[3] - line[1] };
	Vec2i slope = { pt[0] - line[0], pt[1] - line[1]};
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
		Vec2i pt1, pt2; pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
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
	Vec2i pt1, pt2; pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
	if (pt[0] - pt1[0] < 5 || pt[0] - pt2[0] < 5)
		return 1;
	if ((pt[0] - pt1[0] + 5)*(pt[0] - pt2[0] + 5) > 0)
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
	Vec2i pt1, pt2; pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
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
	Vec2i pt1, pt2; pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
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
inline void line2pt(Vec4i line, Vec2i &pt1, Vec2i &pt2)
{
	pt1 = { line[0], line[1] }; pt2 = { line[2], line[3] };
}
Vec2i ptAttachToCircle(Vec2f &tmp ,vector<circle_class> &circles)
{
	Vec2i cross;
	if (circles.size() != 0)
	{
		for (int j = 0; j < circles.size(); j++)
		{
			Vec3i c = circles[j].getCircleVec();
			Vec2f center = circles[j].getCenter(); float radius = circles[j].getRadius();
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
				int cenx = int(tmp[0]); int ceny = int(tmp[1]);
				int offset = int(p2pdistance(tmp, center) - radius);
				offset = (offset < 1) ? 1 : offset;
				for (auto m = cenx - offset; m <= cenx + offset; m++)
				{
					for (auto n = ceny - offset; n <= ceny + offset; n++)
					{
						Vec2i temp2 = { m, n };
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
		cross = { int(tmp[0]), int(tmp[1]) };
		cout << "after attach to circle the cross is " << cross << endl;
	}
	else
	{
		cross = { int(tmp[0]), int(tmp[1]) };
		cout << "not circle on image, no attachment" << endl;
	}
	return cross;
}

void getCrossPt(Vec4i line1, Vec4i line2, Vec2f &tmpCross)
{
	int x1 = line1[0]; int x2 = line1[2]; int x3 = line2[0]; int x4 = line2[2];
	int y1 = line1[1]; int y2 = line1[3]; int y3 = line2[1]; int y4 = line2[3];
	tmpCross[0] = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) * 1.0 / ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
	tmpCross[1] = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) * 1.0 / ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
}
//void getCrossPtRev(vector<Point2i> edgePoints, line_class *linec1, line_class *linec2, point_class &cross_point, vector<circle_class> &circles, vector<line_class> &linexs, vector<point_class> &pointxs)
//{
//	Vec2f tmp; Vec2i pt1, pt2, pt3, pt4;
//	line2pt(linec1->getLineVec(pointxs), pt1, pt2); line2pt(linec2->getLineVec(pointxs), pt3, pt4);
//	int sharePts = 0; int cross_in_one_line = 0;
//	
//	getCrossPt(linec1->getLineVec(pointxs), linec2->getLineVec(pointxs), tmp);
//	cout << "initial cross point is " << tmp << endl;
//	if (same_pt(pt1, pt3))
//	{
//		tmp = (pt1 + pt3) / 2.0;
//		sharePts = 1;
//		cout << "the pt1 and pt3 are approximately the same" << endl;
//	}
//	else if (same_pt(pt1, pt4))
//	{
//		tmp = (pt1 + pt4) / 2.0;
//		sharePts = 2;
//		cout << "the pt1 and pt4 are approximately the same" << endl;
//	}
//	else if (same_pt(pt2, pt3))
//	{
//		tmp = (pt2 + pt3) / 2.0;
//		sharePts = 3;
//		cout << "the pt2 and pt3 are approximately the same" << endl;
//	}
//	else if (same_pt(pt2, pt4))
//	{
//		tmp = (pt2 + pt4) / 2.0;
//		sharePts = 4;
//		cout << "the pt2 and pt4 are approximately the same" << endl;
//	}
//	else if (same_pt(tmp, pt1))
//	{
//		if (dashLineRecovery(edgePoints, tmp, pt1, circles, false, false, false))
//		{
//			tmp = pt1;
//		}
//		cross_in_one_line = 1;
//		cout << "the cross and pt1 are approximately the same" << endl;
//	}
//	else if (same_pt(tmp, pt2))
//	{
//		cross_in_one_line = 2;
//		if (dashLineRecovery(edgePoints, tmp, pt1, circles, false, false, false))
//		{
//			tmp = pt2;
//		}
//		cout << "the cross and pt2 are approximately the same" << endl;
//	}
//	else if (same_pt(tmp, pt3))
//	{
//		cross_in_one_line = 3;
//		if (dashLineRecovery(edgePoints, tmp, pt1, circles, false, false, false))
//		{
//			tmp = pt3;
//		}
//		cout << "the cross and pt3 are approximately the same" << endl;
//	}
//	else if (same_pt(tmp, pt4))
//	{
//		cross_in_one_line = 4;
//		if (dashLineRecovery(edgePoints, tmp, pt1, circles, false,false,false))
//		{
//			tmp = pt4;
//		}
//		cout << "the cross and pt4 are approximately the same" << endl;
//	}
//	else
//	{
//		cout << "no special relationships" << endl;
//	}
//	/*getCrossPt(linec1.getLineVec(pointxs), linec2.getLineVec(pointxs), tmp);
//	if (!same_pt(tmp, tmp1))
//		tmp = tmp1;*/
//	Vec2i cross_vec = ptAttachToCircle(tmp, circles);
//	cross_point.setXY(cross_vec);
//	cout << "now the cross point is " << cross_point.getXY() << endl;
//	if (sharePts)
//	{
//		// if two line share a common point, the point's id should be asigned equal.
//		switch (sharePts)
//		{
//		default: break;
//		case 0:
//			break;
//		case 1:
//			{
//				cout << "set pt1 to cross" << " and set pt3 to cross" << endl;
//				cout << "change pt3 id " << linec2->getPt1Id() << " to pt1 id " << linec1->getPt1Id() << endl;
//				cross_point.setPid(linec1->getPt1Id());
//				linec1->setPt1_vec(pointxs, cross_vec); linec2->setPt1_vec(pointxs, cross_vec);
//				linec2->setpt1Id(linec1->getPt1Id());
//				break;
//			}
//		case 2:
//			{
//				cout << "set pt1 to cross" << " and set pt4 to cross" << endl;
//				cout << "change pt4 id " << linec2->getPt2Id() << " to pt1 id " << linec1->getPt1Id() << endl;
//				cross_point.setPid(linec1->getPt1Id());
//				linec1->setPt1_vec(pointxs, cross_vec); linec2->setPt2_vec(pointxs, cross_vec);
//				linec2->setpt2Id(linec1->getPt1Id());
//				break;
//			}
//		case 3:
//			{
//				cout << "set pt2 to cross" << " and set pt3 to cross" << endl;
//				cout << "change pt3 id " << linec2->getPt1Id() << " to pt2 id " << linec1->getPt2Id() << endl;
//				cross_point.setPid(linec1->getPt2Id());
//				linec1->setPt2_vec(pointxs, cross_vec); linec2->setPt1_vec(pointxs, cross_vec);
//				linec2->setpt1Id(linec1->getPt2Id());
//				break;
//			}
//		case 4:
//			{
//				cout << "set pt2 to cross" << " and set pt4 to cross" << endl;
//				cout << "change pt4 id " << linec2->getPt2Id() << " to pt2 id" << linec1->getPt2Id() << endl;
//				cross_point.setPid(linec1->getPt2Id());
//				linec1->setPt2_vec(pointxs, cross_vec); linec2->setPt2_vec(pointxs, cross_vec);
//				linec2->setpt2Id(linec1->getPt2Id());
//				break;
//			}
//		}
//	}
//	else if (cross_in_one_line)
//	{
//		// the cross points is only in one line, we'll only change the coordinate infos here, no id change issues here.
//		switch (cross_in_one_line)
//		{
//		default:break;
//		case 0:
//			break;
//		case 1:
//			{
//				linec1->setPt1_vec(pointxs, cross_vec);
//				cross_point.setPid(linec1->getPt1Id());
//				cout << "set pt1 to cross" << endl;
//				break;
//			}
//		case 2:
//			{
//				linec1->setPt2_vec(pointxs, cross_vec);
//				cross_point.setPid(linec1->getPt2Id());
//				cout << "set pt2 to cross" << endl;
//				break;
//			}
//		case 3:
//			{
//				linec2->setPt1_vec(pointxs, cross_vec);
//				cross_point.setPid(linec2->getPt1Id());
//				cout << "set pt3 to cross" << endl;
//				break;
//			}
//		case 4:
//			{
//				linec2->setPt2_vec(pointxs, cross_vec);
//				cross_point.setPid(linec2->getPt2Id());
//				cout << "set pt4 to cross" << endl;
//				break;
//			}
//		}
//	}
//	
//
////	if (circles.size() != 0)
////	{
////		for (int j = 0; j < circles.size(); j++)
////		{
////			Vec3i c = circles[j].getCircleVec();
////			Vec2f center = circles[j].getCenter(); float radius = circles[j].getRadius();
/////*			if (on_circle(pt1, c) || on_circle(pt2, c) || on_circle(pt3, c) || on_circle(pt4, c))
////			{
////				continue;
////			}
////			else */
////			//To Note: there're changes here
////			if (on_circle(tmp, c))
////			{
////				float minDiff = 100;
////				int cenx = int(tmp[0]); int ceny = int(tmp[1]);
////				int offset = int(p2pdistance(tmp, center) - radius);
////				offset = (offset < 1) ? 1 : offset;
////				for (auto m = cenx - offset ; m <= cenx + offset; m++)
////				{
////					for (auto n = ceny - offset ; n <= ceny + offset; n++)
////					{
////						Vec2i temp2 = { m, n };
////						float tmpDiff = abs(p2pdistance(temp2, center) - radius);
////						if (tmpDiff < minDiff)
////						{
////							tmp = temp2;
////							minDiff = tmpDiff;
////						}
////
////					}
////				}
////			}
////		}
////		cross = { int(tmp[0]), int(tmp[1]) };
////		linec1->setPt1_vec(cross); linec2->setPt2_vec(cross);
////	}
////	else if (!sharePts)
////	{
////		//if cross in one line
////		if (cross_in_line1)
////		{
////			cross = { int(tmp[0]), int(tmp[1]) };
////			linec1->setPt1_vec(cross);
////		}
////		else if (cross_in_line)
////		{
////			cross = { int(tmp[0]), int(tmp[1]) };
////			linec1->setPt2_vec(cross);
////		}
////		else if (cross_in_line)
////		{
////			cross = { int(tmp[0]), int(tmp[1]) };
////			linec2->setPt1_vec(cross);
////		}
////		else if (cross_in_line2_2)
////		{
////			cross = { int(tmp[0]), int(tmp[1]) };
////			linec2->setPt2_vec(cross);
////		}
////		else
////		{
////			/*for (int i = 0; i < linexs.size(); i++)
////			{
////				Vec4i line = linexs[i].getLineVec(pointxs);
////				if (in_line(line, tmp))
////				{
////					Vec2f temp1, temp2;
////					getCrossPt(line, linec1->getLineVec(pointxs), temp1); getCrossPt(line, linec2->getLineVec(pointxs), temp2);
////					Vec2f temp3 = (temp1 + temp2) / 2.0;
////					tmp = { temp3[0], temp3[1] };
////					break;
////				}
////
////			}*/
////			cross = { int(tmp[0]), int(tmp[1]) };
////		}
////	}
//
//}



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
	for (vector<Vec2i>::iterator iter3 = colpoints.begin(); iter3 != colpoints.end(); ++iter3)
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

//void point_on_circle_line_check(vector<Vec2i> basicEndpoints, vector<Vec3f> circle_candidates, vector<circleX> &circles,
//	vector<Vec4i> line_candidates, vector<line_class> &lines, vector<point_class> &points)
//{
//	bool flag = true;
//	//cout << basicEndpoints.size() << endl;
//	for (int i = 0; i < basicEndpoints.size(); ++i)
//	{
//		Vec2i bpoint = basicEndpoints[i];
//		point_class point;
//		point.p_idx = i; point.px = bpoint[0]; point.py = bpoint[1];
//		for (int j = 0; j < circle_candidates.size(); ++j)
//		{
//			Vec3f bcircle = circle_candidates[j];
//			circleX circle; circle.c_idx = j; circle.cx = bcircle[0]; circle.cy = bcircle[1]; circle.radius = bcircle[2];
//			
//			//Vec2i center = { cvRound(bcircle[0]),cvRound(bcircle[1])};
//			if (flag)
//			{
//				circles.push_back(circle);
//				point_class centerPoint; centerPoint.cflag = true; centerPoint.p_idx = -1; centerPoint.px = cvRound(bcircle[0]); centerPoint.py = cvRound(bcircle[1]);
//				points.push_back(centerPoint);
//			}
//			if (on_circle(basicEndpoints[i], bcircle))
//			{
//				point.c_idx.push_back(j);
//			}
//			
//		}
//		flag = false;
//		for (size_t k = 0; k < line_candidates.size(); ++k)
//		{
//			Vec4i bline = line_candidates[k];
//			line_class line; Vec2i bpoint1, bpoint2; bpoint1 = { bline[0], bline[1] }; bpoint2 = { bline[2], bline[3] };
//			line.l_idx = k; line.px1 = bline[0]; line.py1 = bline[1]; line.px2 = bline[2]; line.py2 = bline[3];
//			line.length = p2pdistance(bpoint1, bpoint2);
//			lines.push_back(line);
//			point.cflag = false;
//			if (on_line(line_candidates[k], bpoint))
//			{
//				point.l_idx.push_back(k);
//			}
//		}
//		points.push_back(point);
//	}
//	//cout << points.size() << endl;
//}

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
	//int eps = 5;
	if (pt[0] >= leftx && pt[0] <= rightx  &&
		pt[1] >= lowy  && pt[1] <= highy )
		return true;
	else
		return false;
}

Vec4i pt2line(Vec2i pt1, Vec2i pt2)
{
	Vec4i ret = { pt1[0], pt1[1], pt2[0], pt2[1] };
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
int ptWithCircle(Vec2i center, int radius ,Vec2i pt)
{
	double pt2center = p2pdistance(center, pt);
	if (radius - pt2center > 4)
		return 0;//inside circle
	else if (abs(radius - pt2center) <= 4)
		return 1; // on cirlce
	else if (pt2center - radius < 8)
		return 2;//outside circle
}
bool basicRev(vector<Point2i> &edgePt, Vec2i p1, Vec2i p2, vector<circle_class> &circles)
{
	Vec4i line = pt2line(p1, p2);
	int flag[1000] = { 0 };
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
	double ratio; int nums = count(flag, flag + 1000, 1);
	ratio = nums / ranges; cout << ratio * 100 << "%" << endl;
	int threshold_ratio = 0.8;
	if (ratio < threshold_ratio)
		return false;
	else
		return true;
}
bool dashLineRecovery(vector<Point2i> &edgePt, Vec2i col_p1,  Vec2i col_p2, vector<circle_class> &circles, bool plflag = false, bool pcflag = false, bool ppflag = false)
{//check if there's dash line between
	Vec4i line = pt2line(col_p1, col_p2);
	int flag[1000] = { 0 };
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
		double ratio; int nums = count(flag, flag + 1000, 1);
		ratio = nums / ranges; cout << ratio * 100 << "%" << endl;
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
			circle_class *c = &(circles[i]);
			Vec2i center = c->getCenter(); int radius = c->getRadius();
			double center2lineDis = pt2lineDis(line, center);
			int flag1 = ptWithCircle(center, radius, col_p1); int flag2 = ptWithCircle(center, radius, col_p2);
			if (abs(center2lineDis - radius) <= 3)
			{
				if (p2pdistance(col_p1,col_p2) < 20)
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
				double ratio; int nums = count(flag, flag + 1000, 1);
				ratio = nums / ranges; cout << ratio * 100 << "%" << endl;
				double threshold_ratio = 0.8;
				if (ratio < threshold_ratio)
					return false;
				else
					return true;;
			}
		}
	}
}
bool dashLineRecovery(vector<Point2i> &edgePt, Vec2i p_closer, Vec2i p_farther, Vec2i p_cross, vector<circle_class> &circles, bool plflag=false, bool pcflag=false, bool ppflag=false)
{		
	//check if there's dash line between
	Vec4i line = pt2line(p_closer, p_cross);
	int flag[1000] = { 0 };
	bool vertical_flag = (abs(p_closer[0] - p_cross[0]) < abs(p_closer[1] - p_cross[1])) ? true : false;
	int ranges = vertical_flag ? abs(p_cross[1] - p_closer[1]) : abs(p_cross[0] - p_closer[0]);
	cout << "line check between " << p_closer << " and " << p_cross << endl;
	if (circles.size()==0)
	{
		cout << "image without circle" << endl;
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
		double ratio; int nums =  count(flag, flag + 1000, 1);
		ratio = nums / ranges; cout << ratio * 100 << "%" << endl;
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
			circle_class *c = &(circles[i]);
			Vec2i center = c->getCenter(); int radius = c->getRadius();
			double center2lineDis = pt2lineDis(line, center);
			int closer_flag = ptWithCircle(center, radius, p_closer); int cross_flag = ptWithCircle(center, radius, p_cross);
			double closerp_diff = abs(p2pdistance(p_closer, center) - radius); double crossp_diff = abs(p2pdistance(p_cross, center) - radius);
			cout << "closerp diff " << closerp_diff << " crossp diff " << crossp_diff << endl;
			cout << abs(center2lineDis - radius) << endl;
			if (abs(center2lineDis - radius) <= 5)
			{
				// the line is almost tangent to the circle				
				if (closer_flag==1 && cross_flag == 1)
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
					return  true;
				}
				else
				{
					cout << "neither points are on the circle" << endl;
					return false;
				}
			}
			else
			{
				cout << "not tangent to the cirlce" << endl;
				return false;
			}
		}
	}
}


double lSlope(Vec4i line)
{
	Vec2i pt1, pt2;	line2pt(line, pt1, pt2);
	Vec2i lineV = { line[2] - line[0], line[3] - line[1] };
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
	line2pt(line1, pt1, pt2); line2pt(line2, pt3, pt4);
	Vec2i line1V = { line1[2] - line1[0], line1[3] - line1[1] }; Vec2i line2V = { line2[2] - line2[0], line2[3] - line2[1] };
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
	if ((abs(line1V[0]) < 3 && abs(line2V[0])<3) | (abs(line1V[1]) < 3 && abs(line2V[1])< 3) )
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
	return (angle <= 5)||(angle>=175);
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
//void PointLineRevision(vector<Point2i> &edgePositions, vector<vector<Vec4i> &plainLines, vector<Vec2i> &plainPoints)
//{
//	////obtain original line and points
//	//vector<Vec2i> tmpPlainPoints;
//	//for (size_t j = 0; j < plainLines.size(); ++j)
//	//{
//	//	Vec4i l = plainLines[j];
//	//	tmpPlainPoints.push_back({ plainLines[j][0], plainLines[j][1] });
//	//	tmpPlainPoints.push_back({ plainLines[j][2], plainLines[j][3] });
//	//	//line(color_img, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(255, 255, 0), 1, 8);
//	//}
//	//// erase the points which is too close to each other
//	//sort(tmpPlainPoints.begin(), tmpPlainPoints.end(), [](Vec2i a, Vec2i b){return a[0] < b[0]; });
//	//tmpPlainPoints.erase(unique(tmpPlainPoints.begin(), tmpPlainPoints.end(), [](Vec2i a, Vec2i b){return same_pt(a, b); }));
//
//	// dash line recovery
//	//while(iter1 != plainLines.end())
//	//{
//	//	Vec4i line = *iter1; Vec2i pt1 = { line[0], line[1] }; Vec2i pt2 = { line[2], line[3] };
//	//	for (auto iter2 = plainPoints.begin(); iter2 != plainPoints.end();)
//	//	{
//	//		Vec2i pt3 = *iter2;
//	//		if (!same_pt(pt3, pt1) && !same_pt(pt3, pt2))
//	//		{
//	//			// dash line recovery
//	//			Vec2i pt1 = { line[0], line[1] }; Vec2i pt2 = { line[2], line[3] };
//	//			Vec2i pt4; double len1, len2; int gap;
//	//			len1 = p2pdistance(pt1, pt3); len2 = p2pdistance(pt2, pt3); 
//	//			bool zflag = true;//assume d(pt1, pt3) is bigger
//	//			if (len1 > len2)
//	//			{
//	//				pt4 = pt1;
//	//				gap = abs(pt3[0] - pt1[0]);
//	//			}
//	//			else
//	//			{
//	//				pt4 = pt2;
//	//				gap = abs(pt3[0] - pt2[0]);
//	//				zflag = false;
//	//			}
//	//			int bitmap[3000] = { 0 }; int count[10] = { 0 }; bool flag = true;
//	//			int step = int(gap / 10.0);
//	//			for (auto i = 0; i < edgePositions.size(); i++)
//	//			{
//	//				Vec2i tmp = { edgePositions[i].x, edgePositions[i].y };
//	//				if (on_line(line, tmp))
//	//				{
//	//					bitmap[tmp[0]] = 1;
//	//				}
//	//			}
//	//			for (int i = 0; i < 10; i++)
//	//			{
//	//				count[i] = count_if(bitmap + i*step, bitmap + (i + 1)*step, [](int a){return a == 1; });
//	//			}
//	//			for (int j = 0; j < 10; j++)
//	//			{
//	//				if (count[j] < 10)
//	//				{
//	//					flag = false;
//	//					break;
//	//				}
//	//			}
//	//			if (flag)
//	//			{
//	//				cout << "test" << endl;
//	//				/*iter1 = plainLines.erase(iter1);
//	//				Vec4i tmpline = { pt3[0], pt3[1], pt4[0], pt4[1] }; plainLines.push_back(tmpline);*/
//	//				if (zflag)
//	//				{
//	//					line[2] = pt3[0]; line[3] = pt3[1];
//
//	//				}
//	//			}
//	//			else
//	//			{
//	//				iter1++;
//	//				plainPoints.push_back(pt1);
//	//				plainPoints.push_back(pt2);
//	//				iter2++;
//	//			}
//	//		
//	//			//copy(count, count + 10, ostream_iterator<char>(cout, " "));	
//	//		}
//	//	}
//	//}
//	
//	
//	for (auto iter1 = plainLines.begin(); iter1 != plainLines.end();iter1++)
//	{
//		Vec4i line1 = *iter1; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
//		for (auto iter2 = plainLines.begin(); iter2 != plainLines.end();)
//		{
//			Vec4i line2 = *iter2; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
//			if (isParallel(line1, line2))
//			{
//				map<int, int> pMap; // pMap is a Map storing the x-axis and y-axis and from x-axis we can infer the y axis
//				pMap[pt1[0]] = pt1[1]; pMap[pt2[0]] = pt2[1]; 
//				pMap[pt3[0]] = pt3[1]; pMap[pt4[0]] = pt4[1];
//				auto iter_begin = pMap.begin(); 
//				auto iter_end = pMap.end(); iter_end--;
//				//if two line is parallel, find the two closest point in 4 points, check if there're dash line
//				Vec2i p1 = { iter_begin->first, iter_begin->second }; Vec2i p2 = { iter_end->first, iter_end->second };
//				if (dashLineRecovery(edgePositions, p1, p2))
//				{
//					// check if there's a dash line to be recovered and erase the old line and points and add the new points and lines
//					// the second and the third points in the map is to be erased
//					plainPoints.push_back(p1); plainPoints.push_back(p2);
//					// assign new line to the current iter1 line and erase the line which iter2 points to 
//					Vec4i newline = { p1[0], p1[1], p2[0], p2[1] };
//					*iter1 = newline;
//					iter2 = plainLines.erase(iter2);
//				}
//				else
//				{
//					plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//					plainPoints.push_back(pt3); plainPoints.push_back(pt4);
//					iter2++;
//				}
//
//				
//			}
//			else
//			{
//				//not parallel, choose the two end points pt3, pt4 to be tested for the dash line recovery
//				if (same_pt(pt3, pt1)||same_pt(pt3, pt2)||same_pt(pt4,pt1)||same_pt(pt4,pt2))
//				{
//					plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//				}
//				else if (on_nin_line(line1, pt3))
//				{
//					if (abs(pt1[0] - pt3[0]) > abs(pt2[0] - pt3[0]))
//					{
//						// pt1-pt3 longer
//						if (dashLineRecovery(edgePositions, pt2, pt3))
//						{
//							Vec4i newline = { pt1[0], pt1[1], pt3[0], pt3[1] };
//							*iter1 = newline;
//							plainPoints.push_back(pt1); plainPoints.push_back(pt3);
//						}
//						else
//						{
//							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//						}
//					}
//					else
//					{
//						//pt2-pt3 longer
//						if (dashLineRecovery(edgePositions, pt1, pt3))
//						{
//							Vec4i newline = { pt2[0], pt2[1], pt3[0], pt3[1] };
//							*iter1 = newline;
//							plainPoints.push_back(pt2); plainPoints.push_back(pt3);
//						}
//						else
//						{
//							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//						}
//					}
//				}
//				else if (on_nin_line(line1, pt4))
//				{
//					if (abs(pt1[0] - pt4[0]) > abs(pt2[1] - pt4[0]))
//					{
//						if (dashLineRecovery(edgePositions, pt2, pt4))
//						{
//							Vec4i newline = { pt1[0], pt1[1], pt4[0], pt4[1] };
//							*iter1 = newline;
//							plainPoints.push_back(pt1); plainPoints.push_back(pt4);
//						}
//						else
//						{
//							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//						}
//					}
//					else
//					{
//						if (dashLineRecovery(edgePositions, pt1, pt4))
//						{
//							Vec4i newline = { pt2[0], pt2[1], pt4[0], pt4[1] };
//							*iter1 = newline;
//							plainPoints.push_back(pt2); plainPoints.push_back(pt4);
//						}
//						else
//						{
//							plainPoints.push_back(pt1); plainPoints.push_back(pt2);
//						}
//
//					}
//				}
//				else
//				{
//					cout << "test" << endl;
//				}
//				iter2++;
//			}
//		}
//	}
//	// erase the points which is too close to each other
//	sort(plainPoints.begin(), plainPoints.end(), [](Vec2i a, Vec2i b){return a[0] < b[0]; });
//	plainPoints.erase(unique(plainPoints.begin(), plainPoints.end(), [](Vec2i a, Vec2i b){return same_pt(a, b); }), plainPoints.end());
//	for (auto i = 0; i < plainPoints.size(); i++)
//	{
//		Vec2i pt1 = plainPoints[i];
//		if (pt1[1] == 77)
//			cout << "li" << endl;
//		for (auto j = i+1; j < plainPoints.size(); j++)
//		{
//			Vec2i pt2 = plainPoints[j];
//			if (pt2[1] == 77)
//				cout << "li" << endl;
//			if (with_same_line(plainLines, pt1, pt2))
//				continue;
//			else
//			{
//				if (dashLineRecovery(ept, pt1, pt2))
//				{
//					Vec4i newline = { pt1[0], pt1[1], pt2[0], pt2[1] };
//					plainLines.push_back(newline);
//				}
//			}
//			
//		}
//	}
//}

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
	int x = pt[0]; int y = pt[1];
	if (x<=0 || x>xMax)
		return false;
	if (y<=0 || y>yMax)
		return false;
	return true;
}
int isEndPoint(Vec4i line, Vec2i pt)
{
	Vec2i pt1 = { line[0], line[1] }; Vec2i pt2 = { line[2], line[3] };
	if (same_pt(pt1, pt))
		return 1;
	else if (same_pt(pt2, pt))
		return 2;
	else
		return 0;
}

void chooseNearCircle(vector<Vec3f> &circle_candidates, Vec2i &cross, Vec2i &pt)
{
	for (auto i = 0; i < circle_candidates.size(); i++)
	{
		Vec3f circle = circle_candidates[i];
		Vec2i center = { int(circle[0]), int(circle[1]) };
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
	pt1 = { line1[0], line1[1] }; pt2 = { line1[2], line1[3] };
	if (in_line(line2, pt1) && in_line(line2, pt2))
		return true;
	else
		return false;
}
bool existRealLineWithinPtxs(vector<line_class> &lineXs, vector<point_class> pointxs,point_class ptx1, point_class ptx2)
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
	line1 = { ptx1.getX(), ptx1.getY(), ptx2.getX(), ptx2.getY() }; line2 = { ptx2.getX(), ptx2.getY(), ptx1.getX(), ptx1.getY() };
	//auto iter1 = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){return in_line(a.lxy, ptx1.pxy); });
	//auto iter2 = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){return in_line(a.lxy, ptx2.pxy); });
	auto iter = find_if(lineXs.begin(), lineXs.end(), [&](line_class a){
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
	Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
	Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[0], line2[1] };
	bool flag1 = (same_pt(pt1, pt3) && same_pt(pt2, pt4)); bool flag2 = (same_pt(pt1, pt4) && same_pt(pt2, pt3));
	if (flag1 || flag2)
		return true;
	else
		return false;
}


double slopeEval(Vec4i line)
{
	double diff1, diff2, diff;
	diff1 = lSlope(line) - 0; diff2 = lSlope(line) - 90;
	diff = (diff1 < diff2) ? diff1 : diff2;
	return diff;
}
void vecSwapValue(vector<Vec4i>::iterator iter1, vector<Vec4i>::iterator iter2)
{
	auto tmp = *iter1;
	*iter1 = *iter2;
	*iter2 = tmp;
}

bool nearToAcross(Vec2i pt, vector<Vec2i> &crossPts)
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

int line_recovery_process(line_class *linex, Vec2i p_cross, vector<Point2i> &edgePt,vector<point_class> &pointxs,vector<circle_class> &circlexs, bool id_change = false)
{
	Vec2i pt1 = linex->getpt1vec(pointxs); Vec2i pt2 = linex->getpt2vec(pointxs);
	double dis1 = p2pdistance(pt1, p_cross); double dis2 = p2pdistance(pt2, p_cross);
	int pos = 0;
	if ( dis1 < dis2)
	{
		if (dis1 < 5)
		{
			cout << linex->getpt1vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt1_vec(pointxs, p_cross);
			pos = 1;
		}
		else if(dashLineRecovery(edgePt, pt1, pt2, p_cross, circlexs))
		{
			cout << linex->getpt1vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt1_vec(pointxs, p_cross);
			pos = 1;
		}
	}
	else
	{
		if (dis2 < 5 )
		{
			cout << linex->getpt2vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt2_vec(pointxs, p_cross); 
			pos = 2;
			
		}
		else if(dashLineRecovery(edgePt, pt2, pt1, p_cross, circlexs))
		{
			cout << linex->getpt2vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt2_vec(pointxs, p_cross);
			pos = 2;
		}
	}
	return pos;
}

void cross_refinement(Vec2f &raw_cross, line_class *lx1, line_class *lx2, vector<circle_class> &circlexs, vector<point_class> &pointxs, vector<Point2i> &edgePt)
{
	Vec4i linex1_vec = lx1->getLineVec(pointxs); Vec4i linex2_vec = lx2->getLineVec(pointxs);
	Vec2i pt1, pt2, pt3, pt4;  line2pt(linex1_vec, pt1, pt2); line2pt(linex2_vec, pt3, pt4);
	// check according to the relationship between cross and the two lines
	
	bool in_line1, in_line2; 
	in_line1 = in_line(linex1_vec, raw_cross);
	in_line2 = in_line(linex2_vec, raw_cross);
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

	if (same_pt(raw_cross, pt1))
	{
		cout << "raw_cross =  pt1" << endl;
		cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
		lx1->setPt1_vec(pointxs, raw_cross);
		if (same_pt(raw_cross, pt3))
		{
			cout << "also, raw_cross = pt3" << endl;
//			cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			lx2->setPt1_vec(pointxs, raw_cross);
			lx1->setPt1_vec(pointxs, raw_cross);
			lx2->setpt1Id(lx1->getPt1Id()); 
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
		}
		else if (same_pt(raw_cross, pt4))
		{
			cout << "also, raw_cross = pt4" << endl;
//			cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			lx1->setPt1_vec(pointxs, raw_cross);
			lx2->setPt2_vec(pointxs, raw_cross);
			lx2->setpt2Id(lx1->getPt1Id());
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
		}
		else
		{
			if (in_line2)
			{
				cout << "sp" << endl;
			}
			else
			{
				int pos = line_recovery_process(lx2, lx1->getpt1vec(pointxs), edgePt, pointxs, circlexs);
				if (pos == 0)
				{
					cout << "no recovery" << endl;
				}
				else if (pos == 1)
				{
					lx2->setpt1Id(lx1->getPt1Id());
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
				}
				else if (pos == 2)
				{
					lx2->setpt2Id(lx1->getPt1Id());
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
				}
			}
		}
	}
	else if (same_pt(raw_cross, pt2))
	{
		cout << "raw_cross =  pt2" << endl;
		cout << lx1->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
		lx1->setPt2_vec(pointxs, raw_cross);
		if (same_pt(raw_cross, pt3))
		{
			cout << "also, raw_cross = pt3" << endl;
//			cout << lx1->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			lx1->setPt2_vec(pointxs, raw_cross);
			lx2->setPt1_vec(pointxs, raw_cross);
			lx2->setpt1Id(lx1->getPt2Id());
			cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;
		}
		else if (same_pt(raw_cross, pt4))
		{
			cout << "also, raw_cross = pt4" << endl;
//			cout << lx1->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
			lx1->setPt2_vec(pointxs, raw_cross);
			lx2->setPt2_vec(pointxs, raw_cross);
			lx2->setpt2Id(lx1->getPt2Id());
			cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
		}
		else
		{
			if (in_line2)
			{
				cout << "sp" << endl;
			}
			else
			{
				int pos = line_recovery_process(lx2, lx1->getpt2vec(pointxs), edgePt, pointxs, circlexs);
				if (pos == 0)
				{
					cout << "no recovery" << endl;
				}
				else if (pos == 1)
				{
					lx2->setpt1Id(lx1->getPt2Id());
					cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
				}
				else if (pos == 2)
				{
					lx2->setpt2Id(lx1->getPt2Id());
					cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
				}
			}
		}
	}
	else if (same_pt(raw_cross, pt3))
	{
		cout << "raw_cross =  pt3" << endl;
		cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
		lx2->setPt1_vec(pointxs, raw_cross);
		if (in_line1)
		{
			cout << "half inner cross" << endl;
		}
		else
		{
			int pos = line_recovery_process(lx1, lx2->getpt1vec(pointxs), edgePt, pointxs, circlexs);
			if (pos == 0)
			{
				cout << "no recovery" << endl;
			}
			else if (pos == 1)
			{
				//line rev between cross and first point in line	
				lx1->setpt1Id(lx2->getPt1Id());
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
			}
			else if (pos == 2)
			{
				//line rev between cross and second point in line	
				lx1->setpt2Id(lx2->getPt1Id());
				cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
			}
		}
	}
	else if (same_pt(raw_cross, pt4))
	{
		cout << "raw_cross =  pt4" << endl;
		cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
		lx2->setPt2_vec(pointxs, raw_cross);
		if (in_line1)
		{
			cout << "half inner cross" << endl;
		}
		else
		{
			int pos = line_recovery_process(lx1, lx2->getpt2vec(pointxs), edgePt, pointxs, circlexs);
			if (pos == 0)
			{
				cout << "no recovery" << endl;
			}
			else if (pos == 1)
			{
				lx1->setpt1Id(lx2->getPt2Id());
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
			}
			else if (pos == 2)
			{
				lx1->setpt2Id(lx2->getPt2Id());
				cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
			}
		}
	}
	else
	{
		cout << "two cross disjoint line" << endl;
		if (in_line1 && !in_line2)
		{
			line_recovery_process(lx2, raw_cross, edgePt, pointxs, circlexs);
		}
		else if (in_line2 && !in_line1)
		{
			line_recovery_process(lx1, raw_cross, edgePt, pointxs, circlexs);
		}
		else if (in_line1 && in_line2)
		{
			cout << "inner cross" << endl;
		}
		else
		{
			cout << "outer cross" << endl;
			int pos1 = line_recovery_process(lx1, raw_cross, edgePt, pointxs, circlexs);
			int pos2 = line_recovery_process(lx2, raw_cross, edgePt, pointxs, circlexs);
			if (pos1==0 || pos2 == 0)
			{
				cout << "no id set" << endl;
			}
			else if (pos1 == 1 && pos2 == 1)
			{
				lx2->setpt1Id(lx1->getPt1Id());
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
			}
			else if (pos1 == 1 && pos2 == 2)
			{
				lx2->setpt2Id(lx1->getPt1Id());
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
			}
			else if (pos1 == 2 && pos2 == 1)
			{
				lx2->setpt1Id(lx1->getPt1Id());
				cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt1vec(pointxs) << " " << lx2->getPt1Id() << endl;;
			}
			else if (pos1 == 2 && pos2 == 2)
			{
				lx2->setpt2Id(lx1->getPt2Id());
				cout << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << ", " << lx2->getpt2vec(pointxs) << " " << lx2->getPt2Id() << endl;;
			}
			else
			{
				cout << "pos error" << endl;
			}
		}
	}
	cout << "test stop" << endl;
}

void detect_line3(Mat diagram_segwithoutcircle, Mat &withoutCirBw, vector<Point2i> &edgePoints,vector<point_class> pointXs, vector<circle_class> &circles, Mat &color_img, vector<line_class> lineXs, vector<Vec2i>& plainPoints, Mat &drawedImages, bool showFlag = true, string fileName = "")
{
	
	vector<Point2i> ept = getPointPositions(withoutCirBw);
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
		if (lineEV > 10 )//this threshold should be set a litter lower to generate enough line candidates
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
		namedWindow("7.lines first opt version now 1"); imshow("7.lines first opt version now 1", color_img);
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
	//for (auto iter1 = plainLines.begin(); iter1 != plainLines.end(); iter1++)
	//{
	//	Vec4i line1 = *iter1; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
	//	cout << pt1 << "pt" << pt2 << endl;
	//	int maxXDiff = 0;
	//	
	//	for (auto iter2 = iter1 + 1; iter2 != plainLines.end(); )
	//	{
	//		Vec4i line2 = *iter2; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
	//		cout << pt3 << "pt" << pt4 << endl;
	//		
	//		/*if (pt4[0] == 164)
	//			cout << "stop" << endl;*/
	//		if (isParallel(line1, line2))
	//		{
	//			//if it's parallel
	//			vector<Vec2i> tmpPtVec;
	//			tmpPtVec.push_back(pt1); tmpPtVec.push_back(pt2); 
	//			tmpPtVec.push_back(pt3); tmpPtVec.push_back(pt4);
	//			sort(tmpPtVec.begin(), tmpPtVec.end(), ptSortPred);	
	//			Vec2i firstPt = tmpPtVec[0]; Vec2i fourthPt = tmpPtVec[3];
	//			Vec2i secondPt = tmpPtVec[1]; Vec2i thirdPt = tmpPtVec[2];
	//			
	//			
	//			if (ptLine(line1, pt3))// pt3 of line2 is on line 1
	//			{
	//				// if it's collinear
	//				int eps = 5;
	//				bool flag3 = ((pt3[0] - pt1[0])*(pt3[0] - pt2[0]) > 0); //appoximate if pt3 is in line1
	//				bool flag4 = ((pt4[0] - pt1[0])*(pt4[0] - pt2[0]) > 0);//approximate if pt4 is in  line1
	//				if (same_pt(pt1, pt3))
	//				{
	//					if (abs(pt4[0] - pt1[0]) > abs(pt2[0] - pt1[0]))
	//					{
	//						*iter1 = { pt1[0], pt1[1], pt4[0], pt4[1] };
	//					}
	//					iter2 = plainLines.erase(iter2);
	//				}
	//				else if(same_pt(pt1, pt4))
	//				{
	//					if (abs(pt3[0] - pt1[0]) > abs(pt2[0] - pt1[0]))
	//					{
	//						*iter1 = { pt1[0], pt1[1], pt3[0], pt3[1] };
	//					}
	//					iter2 = plainLines.erase(iter2);
	//				}
	//				else if (same_pt(pt2, pt3))
	//				{
	//					if (abs(pt4[0] - pt2[0]) > abs(pt1[0] - pt2[0]))
	//					{
	//						*iter1 = { pt4[0], pt4[1], pt2[0], pt2[1] };
	//					}
	//					iter2 = plainLines.erase(iter2);
	//				}
	//				else if (same_pt(pt2, pt4))
	//				{
	//					if (abs(pt3[0] - pt2[0]) > abs(pt1[0] - pt2[0]))
	//					{
	//						*iter1 = { pt3[0], pt3[1], pt2[0], pt2[1] };
	//					}
	//					iter2 = plainLines.erase(iter2);
	//				}
	//				// four points are all different
	//				else if (flag3&&flag4)
	//				{
	//					// two disjoint collinear line
	//					cout << firstPt << "range" << fourthPt << endl;
	//					if (dashLineRecovery(ept,secondPt, thirdPt))
	//					{
	//						*iter1 = { firstPt[0], firstPt[1], fourthPt[0], fourthPt[1] };
	//						cout << "now the line1 is " << *iter1 << endl;
	//					}
	//					iter2++;
	//				}
	//				else
	//				{
	//					cout << firstPt << "range" << fourthPt << endl;
	//					*iter1 = { firstPt[0], firstPt[1], fourthPt[0], fourthPt[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//					iter2++;
	//				}
	//			}
	//			else
	//			{
	//				//parallel but not collinear;
	//				iter2++;
	//			}
	//		}
	//		else
	//		{
	//			//if it's not parallel then check the end points should be refined.
	//			Vec2i cross;getCrossPt(line1, line2, cross);
	//			if (cross[0] < 0 || cross[1] < 0)
	//			{
	//				iter2++;
	//				continue;
	//			}
	//			else
	//				cout << cross << endl;
	//			if (crossPtWithinLines(line1, line2, cross))
	//			{
	//				//cross within two lines
	//				plainPoints.push_back(cross);
	//				cout << "crossPtWithinLines " << endl;
	//				if (same_pt(cross, pt1))
	//				{
	//					*iter1 = { cross[0], cross[1], pt2[0], pt2[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//				}
	//				else if (same_pt(cross, pt2))
	//				{
	//					*iter1 = { pt1[0], pt1[1], cross[0], cross[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//				}
	//				else if (same_pt(cross, pt3))
	//				{
	//					*iter2 = { cross[0], cross[1], pt4[0], pt4[1] };
	//					cout << "now the line2 is " << *iter2 << endl;
	//				}
	//				else if (same_pt(cross, pt4))
	//				{
	//					*iter2 = { pt3[0], pt3[1], cross[0], cross[1] };
	//					cout << "now the line2 is " << *iter2 << endl;
	//				}
	//				else
	//				{
	//					cout << "not changed now" << endl;
	//				}
	//				iter2++;
	//				
	//				continue;
	//			}
	//			else
	//			{
	//				//otherwise check if we need to refine the line1's endpoints
	//				// if cross in a one line
	//				if (in_line(line1, cross))
	//				{
	//					// cross within line1
	//					if (abs(cross[0] - pt3[0]) < abs(cross[0] - pt4[0]))
	//					{
	//						// pt3 is closer to cross
	//						if (dashLineRecovery(ept, cross, pt3))
	//						{
	//							*iter2 = { cross[0], cross[1], pt4[0], pt4[1] };
	//							cout << "now the line2 is " << *iter2 << endl;
	//						}
	//					}
	//					else
	//					{
	//						// pt4 is closer to the cross
	//						if (dashLineRecovery(ept, cross, pt4))
	//						{
	//							*iter2 = { pt3[0], pt3[1], cross[0], cross[1] };
	//							cout << "now the line2 is " << *iter2 << endl;
	//						}
	//					}
	//				}
	//				else if (in_line(line2, cross))
	//				{
	//					// cross within line2
	//					if (abs(cross[0] - pt1[0]) < abs(cross[0] - pt2[0]))
	//					{
	//						// pt1 is closer to the cross
	//						if (dashLineRecovery(ept, cross, pt1))
	//						{
	//							*iter1 = { cross[0], cross[1], pt2[0], pt2[1] };
	//							cout << "now the line1 is " << *iter1 << endl;
	//						}
	//					}
	//					else
	//					{
	//						// pt2 is  closer to the cross
	//						if (dashLineRecovery(ept, cross, pt2))
	//						{
	//							*iter1 = { pt1[0], pt1[1], cross[0], cross[1] };
	//							cout << "now the line1 is " << *iter1 << endl;
	//						}
	//					}
	//				}
	//				// if cross is not within both the line here, then check if it's wihtin the small region around two endpoints from two lines
	//				else if (withinPtCRegion(cross,pt1)&& withinPtCRegion(cross,pt3))
	//				{
	//					//1,3
	//					*iter1 = { cross[0], cross[1], pt2[0], pt2[1] };
	//					*iter2 = { cross[0], cross[1], pt4[0], pt4[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//					cout << "now the line2 is " << *iter2 << endl;
	//				}
	//				else if ( withinPtCRegion(cross,pt2) && withinPtCRegion(cross,pt3))
	//				{
	//					//2, 3
	//					*iter1 = { pt1[0], pt1[1], cross[0], cross[1] };
	//					*iter2 = { cross[0], cross[1], pt4[0], pt4[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//					cout << "now the line2 is " << *iter2 << endl;
	//				}
	//				else if (  withinPtCRegion(cross,pt2) &&  withinPtCRegion(cross,pt4))
	//				{
	//					//2,4
	//					*iter1 = { pt1[0], pt1[1], cross[0], cross[1] };
	//					*iter2 = { pt3[0], pt3[1], cross[0], cross[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//					cout << "now the line2 is " << *iter2 << endl;
	//				}
	//				else if ( withinPtCRegion(cross,pt1) && withinPtCRegion(cross,pt4))
	//				{
	//					//1,4
	//					*iter1 = { cross[0], cross[1], pt2[0], pt2[1] };	
	//					*iter2 = { pt3[0], pt3[1], cross[0], cross[1] };
	//					cout << "now the line1 is " << *iter1 << endl;
	//					cout << "now the line2 is " << *iter2 << endl;
	//				}
	//				iter2++;
	//			}
	//		}
	//	}
	//	cout << endl;
	//}



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
			Vec4i line1 = *iter1; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
			cout << pt1 << "pt" << pt2 << endl;
		
			Vec4i line2 = *iter2; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
			cout << pt3 << "pt" << pt4 << endl;
			
			if (isParallel(line1, line2))
			{
				//if it's parallel, followed by checking if collinear
				cout << "line1 and line2 is parallel" << endl;
				vector<Vec2i> tmpPtVec;
				tmpPtVec.push_back(pt1); tmpPtVec.push_back(pt2);
				tmpPtVec.push_back(pt3); tmpPtVec.push_back(pt4);
				bool flag0 = false;
				
				//flag0 indicates whether line 1 is vertical, if vertical sort the tmpvec by y axis info, otherwise sort it  by x axis
				if ( (flag0 =(abs(pt1[0] - pt2[0]) < 3)))
				{
					sort(tmpPtVec.begin(), tmpPtVec.end(), [](Vec2i a, Vec2i b){return a[1] < b[1]; });
					cout << "line1 is vertical" << endl;
				}
				else
				{
					sort(tmpPtVec.begin(), tmpPtVec.end(), ptSortPred);
					cout << "line1 is not vertical" << endl;
					
				}
				Vec2i firstPt = tmpPtVec[0]; Vec2i fourthPt = tmpPtVec[3];
				Vec2i secondPt = tmpPtVec[1]; Vec2i thirdPt = tmpPtVec[2];
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
						flag3 = ((pt3[0] - pt1[0])*(pt3[0] - pt2[0]) > 0); //appoximate if pt3 is in line1 or on line1 but not in it
						flag4 = ((pt4[0] - pt1[0])*(pt4[0] - pt2[0]) > 0);//approximate if pt4 is in  line1
						
					}
					else
					{
						flag3 = ((pt3[1] - pt1[1])*(pt3[1] - pt2[1]) > 0); //appoximate if pt3 is in line1
						flag4 = ((pt4[1] - pt1[1])*(pt4[1] - pt2[1]) > 0);//approximate if pt4 is in  line1					
					}
					// four points are all different
					if (flag3&&flag4&&!same_pt(pt1,pt3)&&!same_pt(pt1,pt4)&&!same_pt(pt2,pt3)&&!same_pt(pt2,pt4))
					{
						// two disjoint collinear line
						cout << "four point are all different,";
						cout << firstPt << "range" << fourthPt << endl;
						if (dashLineRecovery(ept, secondPt, thirdPt,  circles, true,false,false))
						{
							*iter1 = { firstPt[0], firstPt[1], fourthPt[0], fourthPt[1] };
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
						*iter1 = { firstPt[0], firstPt[1], fourthPt[0], fourthPt[1] };
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
//	namedWindow("rm parallel", 0); imshow("rm parallel", color_img);
	
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
	vector<Point2i> ept0 = getPointPositions(withoutCLOriBw);

	bool testflag = false;
	/*write test block*/
	ofstream test; test.open("test/test.txt");
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		cout << pt1 << " " << pt2 << endl;
		if (testflag)
			test << pt1[0] <<" "<<pt1[1]<<" "<< pt2[0] <<" "<<pt2[1] << endl;
	}
	test.close();

	cout << endl << "***********block stop" << endl << endl;
	/****************sort the line and mod the horionzal and vertical line*/
	vector<Vec4i> tmpLines;
	auto ordered_iter = plainLines.begin();
	for (auto iter = plainLines.begin() + 1; iter != plainLines.end(); ++iter)
	{
		Vec4i line = *iter; Vec2i pt1, pt2; line2pt(line, pt1, pt2);
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

	int px_count = 0; int lx_count = 0;
	vector<line_class> linexs; vector<point_class> pointxs;
	
	/***********initialize lines and points**********/
	for (auto i = 0; i < plainLines.size(); i++)
	{
		auto p_l = plainLines[i];
		Vec2i pt1, pt2; line2pt(p_l, pt1, pt2);
		auto _pidx1 = px_count++; int _pidx2 = px_count++; int _lidx = lx_count++;
		point_class *ptx1 = new point_class(pt1, _pidx1); point_class *ptx2 = new point_class(pt2, _pidx2);
		//cout << &ptx1 << endl<<&ptx2<<endl;
		ptx1->pushLid(_lidx); ptx2->pushLid(_lidx);
		pointxs.push_back(*ptx1); pointxs.push_back(*ptx2);
		
		line_class *lx = new line_class(ptx1->getPid(), ptx2->getPid(), _lidx);
		//lx.p1 = &pointxs[_pidx1]; lx.p2 = &pointxs[_pidx2];
		linexs.push_back(*lx);
		delete ptx1, ptx2, lx;
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
		Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		cout << pt1 << " " << pt2 << endl;
		cout << linexs[i].getPt1Id() << " " << linexs[i].getPt2Id() << endl;
	}
	cout << "****************sep" << endl<<endl;

	// with the raw detection result, we may get extra false line or omit some line(eg. dash line falsely removed previously or real line not detected)
	// in order to detect all lines, we suppose to first put all candidates into calculation then remove the false line at the last step.
	/*******************raw line candidates refinement*/
	for (auto i = 0; i < linexs.size(); i++)
	{
		line_class *linex1 = &linexs[i];
		cout << linex1->getLineVec(pointxs) << endl;
		for (auto j = i + 1; j < linexs.size(); j++)
		{
			line_class *linex2 = &linexs[j];
			Vec4i linex1_vec = linex1->getLineVec(pointxs); Vec4i linex2_vec = linex2->getLineVec(pointxs);
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
				Vec2i pt1, pt2, pt3, pt4; Vec2f raw_cross;
				line2pt(linex1_vec, pt1, pt2); line2pt(linex2_vec, pt3, pt4);
				getCrossPt(linex1_vec, linex2_vec, raw_cross);
				cout << endl<<"raw cross" << raw_cross<< endl;

				if (!isInImage(diagram_segwithoutcircle.cols, diagram_segwithoutcircle.rows,raw_cross))
				{
					cout << raw_cross << endl;
					cout << "cross out of scope, then take it as no cross" << endl;
					continue;
				}
				else
				{
					//cross in scope
					cout << linex1_vec << endl << linex2_vec << endl;
					cross_refinement(raw_cross, linex1, linex2, circles, pointxs,ept0);
				}
			}
		}
	}
	cout << endl;
	/*set similar points identical*/


	/*remove indentical vec points and reorder the indices*/
	// point index indicate the line id in which it locate
	map<int, int> changeMap;
	for (auto i = 0; i < pointxs.size(); ++i)
	{
		//the location of the first pointx in the loop
		Vec2i px1 = pointxs[i].getXY();
		int erase_offset = 0; // the var log the erased points num and offset the indexing digit when indexing
		for (auto j = i + 1; j  < pointxs.size(); ++j)
		{
			Vec2i px2 = pointxs[j].getXY();
			cout << pointxs[i].getXY() << "   " << pointxs[j].getXY() << endl;
			if (px1 == px2)
			{
				//the two points vector is the same, then erase the second point from the pointxs set
				// and ajust the respective line with previous erased point
				cout << "find the identical vec points " << px1 << " <- "<<pointxs[i].getPid() <<", "<< pointxs[j].getPid()<< endl;
				changeMap[pointxs[j].getPid()] = pointxs[i].getPid();
				int line_id = pointxs[j].getPid() / 2;
				int pos = pointxs[j].getPid()  % 2;// if 0, means the first point in the line, otherwise 1, means the second point in the line
				pointxs.erase(pointxs.begin() + j);
				if (pos)
				{
					linexs[line_id].setpt2Id(pointxs[i].getPid());
				}
				else
				{
					linexs[line_id].setpt1Id(pointxs[i].getPid());
				}
				j--;
			}
		}
	}
	for (auto i = 0; i<2*linexs.size();++i)
	{
		if (changeMap[i]==0)
		{
			changeMap[i] = i;
		}
	}

	//reordering index
	map<int, int> changeMap2;
	for (auto i = 0; i < pointxs.size(); i++)
	{
		cout << "point id " << pointxs[i].getPid() << pointxs[i].getXY() << endl;
		changeMap2[pointxs[i].getPid()] = i;
		pointxs[i].setPid(i);
	}
	for (auto j = 0; j < linexs.size(); j++)
	{
		linexs[j].setpt1Id(changeMap2[changeMap[linexs[j].getPt1Id()]]);
		linexs[j].setpt2Id(changeMap2[changeMap[linexs[j].getPt2Id()]]);
	}

	/*display*/
	for (auto i = 0; i < linexs.size(); i++)
	{
		Vec4i l = linexs[i].getLineVec(pointxs);
		Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	}
	cout << endl;
	for (auto i = 0; i < pointxs.size();i++)
	{
		cout << pointxs[i].getXY() << endl;
	}
	

	//if (showFlag)
	{
		namedWindow("8.lines first opt version now",0); imshow("8.lines first opt version now", color_img);
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


void primitive_parse(const Mat binarized_image, const Mat diagram_segment, vector<Point2i> &edgePoints,vector<point_class> &points, vector<line_class> &lines, vector<circle_class> &circles, Mat &drawedImages, bool showFlag=true, string fileName="")
{
	/* primitive about points, lines, and circles
	first detect the circle, we can get the matrix of diagram segments without circle for the next
	 processing step*/
	//vector<Vec3i> circle_candidates = {}; 
	Mat color_img = Mat::zeros(diagram_segment.size(),CV_8UC3); cvtColor(diagram_segment, color_img, CV_GRAY2RGB);
	Mat diagram_segwithoutcircle = Mat::zeros(diagram_segment.size(), CV_8UC1);
	
	//diagram_segment = preprocessing(diagram_segment);

	/*ransac go*/
	Mat withoutCirBw = binarized_image.clone();
	// detect the circle and get the img without circle and bw img without cirlce and circle candidates
	detect_circle(diagram_segment, color_img, diagram_segwithoutcircle, withoutCirBw, circles, showFlag);
	// then the line detection
	vector<Vec2i> basicEndpoints = {};
	
	
	detect_line3(diagram_segwithoutcircle, withoutCirBw, edgePoints, points, circles, color_img, lines, basicEndpoints, drawedImages, showFlag, fileName);
	
	//detect_line2(diagram_segwithoutcircle, color_img, line_candidates, basicEndpoints);
	//cout << "basic endpoints num: "<<basicEndpoints.size() << endl;
	
	// for points check if they are in circles
	//point_on_circle_line_check(basicEndpoints, circles, circles, line_candidates, lines, points);
	
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
	Mat image = imread("Sg-18.jpg", 0);
	//namedWindow("original image");
	//imshow("original image", image);
	// then binarize it

	Mat binarized_image = image_binarizing(image, true);
	//binarized_image = preprocessing(binarized_image);
	// then go on a process of connectivity componnent analysis
	int labeln; Mat diagram_segment = Mat::zeros(image.size(), CV_8UC1);
	vector<Mat> label_segment = {};
	vector<Point2i> oriEdgePoints=getPointPositions(binarized_image);
	
	image_labelling(binarized_image, diagram_segment,true);
	vector<point_class> points = {}; vector<line_class> lines = {}; vector<circle_class> circles = {};		
	Mat drawedImages(image.size(), CV_8UC3);

	primitive_parse(binarized_image,diagram_segment, oriEdgePoints, points, lines, circles, drawedImages,false);
   	return 0;
}

int diagram()
{
	//a series of image
	//vector<Mat> images;
	char abs_path[100] = "D:\\data\\graph-DB\\newtest29";
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
		vector<point_class> points = {}; vector<line_class> lines = {}; vector<circle_class> circles = {};
		Mat drawedImages = Mat::zeros(diagram_segment.size(), CV_8UC3);
		primitive_parse(binarized_image,diagram_segment, edgePoints, points, lines, circles, drawedImages, false);
		imwrite(saveimgName, drawedImages);
	}
	return 0;
}

