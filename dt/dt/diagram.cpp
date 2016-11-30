# include "stdafx.h"
# include "diagram.h"

#include "time.h"
#define N 500

//vector<Point2i> edgePointsWithoutCircle;
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
	double distance = 0;
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

bool same_pt(pointX pt1, pointX pt2)
{
	double eps = 8;//to be set parameter
	double high_eps = 15;
	if (abs(pt1.px - pt2.px) < 3 && p2pdistance(pt1.pxy,pt2.pxy) < high_eps)
		return true;
	if (p2pdistance(pt1.pxy, pt2.pxy) <= eps)
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
	auto iterBegin = tmpMap.begin(); auto iterEnd = tmpMap.end(); iterEnd--;
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
	Vec3i ret = {};
	Vec2i center0 = { int(circle[0]), int(circle[1]) };
	int radius0 = int(circle[2]);
	Vec2i center = {}; int radius = 0; 
	int ptnumCurrent = (getPointPositions(bw)).size();
	bool breakFlag = false;
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
		Vec2i center = { cCandidate[0], cCandidate[1] }; int radius = cCandidate[2]; 
		cv::circle(tmp, center, radius, 0, w);
		int tmpPtnum = (getPointPositions(tmp)).size();
		a[w-1] = tmpPtnum;
		width = w;
		if (w != 1 && (a[w - 1] - a[w]) < 50)
			break;
		
	}
	ret = cCandidate;
	int test = (getPointPositions(bw)).size();
	return ret;
}
void detect_circle(const Mat diagram_segment, Mat &color_img,Mat &diagram_segwithoutcircle,Mat &withoutCirBw, vector<Vec3f> &circle_candidates,  bool showFlag)
{
	unsigned int circleN_todetect = 2;
	diagram_segwithoutcircle = diagram_segment;
	for (unsigned int i = 0; i < circleN_todetect; ++i)
	{
		vector<Point2i> edgePositions = getPointPositions(diagram_segwithoutcircle);
		Mat dt; distanceTransform(255 - diagram_segment, dt, CV_DIST_L1, 3);
		unsigned int nIter = 0;
		Point2f bestCircleCenter = {}; float bestCircleRadius = -1.0f;
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
		if (bestCVal<(int)(edgePositions.size() / 8)) 
			break;
		if (bestCVal > 125)
		{
			//std::cout << "current best circle: " << bestCircleCenter << " with radius: " << bestCircleRadius << " and nInlier " << bestCVal << endl;
			
			Vec3f circlef = { bestCircleCenter.x, bestCircleCenter.y, bestCircleRadius };
			int width = 0;
			Vec3i circle = radiusThicknessRev(circlef, diagram_segment,width);
			circle_candidates.push_back(circle);
			//draw the cicle detected in red within the colorgeo blob image
			//cv::circle(color_img, bestCircleCenter, bestCircleRadius, Scalar(0, 0, 255));
			//TODO: hold and save the detected circle.

			//TODO: instead of overwriting the mask with a drawn circle it might be better to hold and ignore detected circles and dont count new circles which are too close to the old one.
			// in this current version the chosen radius to overwrite the mask is fixed and might remove parts of other circles too!

			// update mask: remove the detected circle!
			cv::circle(diagram_segwithoutcircle, { circle[0], circle[1] }, circle[2], 0, width); // here the radius is fixed which isnt so nice.
			//edgePointsWithoutCircle = getPointPositions(diagram_segwithoutcircle);
			cv::circle(color_img, { circle[0], circle[1] }, circle[2], Scalar(255, 0, 255), width);
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
	int count = 0;

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
	int count = 0;

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

float point2Line(Vec4i line, Vec2i pt)
{
	Vec2i bottom = { line[2] - line[0], line[3] - line[1] };
	Vec2i slope = { pt[0] - line[0], pt[1] - line[1] };
	float p2l= abs(cross_product(bottom, slope) / norm(bottom));
	return p2l;
}
bool on_line(Vec4i line, Vec2i pt)
{
	double linedis_eps = 3; double pointdis_eps = 5;
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

void getCrossPt(Vec4i line1, Vec4i line2, Vec2f &tmpCross)
{
	int x1 = line1[0]; int x2 = line1[2]; int x3 = line2[0]; int x4 = line2[2];
	int y1 = line1[1]; int y2 = line1[3]; int y3 = line2[1]; int y4 = line2[3];
	tmpCross[0] = ((x1*y2 - y1*x2)*(x3 - x4) - (x1 - x2)*(x3*y4 - y3*x4)) * 1.0 / ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
	tmpCross[1] = ((x1*y2 - y1*x2)*(y3 - y4) - (y1 - y2)*(x3*y4 - y3*x4)) * 1.0 / ((x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4));
}
void getCrossPtRev(Vec4i line1, Vec4i line2, Vec2i &cross, vector<Vec3f> &circle_candidates, vector<Vec4i> &plainLines)
{
	Vec2f tmp;
	getCrossPt(line1, line2, tmp);
	if (circle_candidates.size() != 0)
	{
		for (int j = 0; j < circle_candidates.size(); j++)
		{
			Vec3f c = circle_candidates[j];
			Vec2f center = { c[0], c[1] }; float radius = c[2];
			if (on_circle(tmp, c))
			{
				float minDiff = 100;
				int cenx = int(tmp[0]); int ceny = int(tmp[1]);
				int offset = int(p2pdistance(tmp, center) - radius);
				for (auto m = cenx - offset; m <= cenx + offset; m++)
				{
					for (auto n = ceny - offset; n <= ceny + offset; n++)
					{
						Vec2i tmp2 = { m, n };
						float tmpDiff = abs(p2pdistance(tmp2, center) - radius);
						if (tmpDiff < minDiff)
						{
							tmp = tmp2;
							minDiff = tmpDiff;
						}

					}
				}
			}
		}
		cross = { int(tmp[0]), int(tmp[1]) };
	}
	else
	{
		if (on_line(line1, tmp) || on_line(line2, tmp))
		{
			cross = { int(tmp[0]), int(tmp[1]) };
		}
		else
		{
			for (int i = 0; i < plainLines.size(); i++)
			{
				Vec4i line = plainLines[i];
				if (in_line(line, tmp))
				{
					Vec2f tmp1, tmp2;
					getCrossPt(line, line1, tmp1); getCrossPt(line, line2, tmp2);
					Vec2f tmp3 = (tmp1 + tmp2) / 2.0;
					tmp = { tmp3[0], tmp3[1] };
					break;
				}

			}
			cross = { int(tmp[0]), int(tmp[1]) };
		}
	}
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
	if (pt[0] >= leftx && pt[0] <= rightx  &&
		pt[1] >= lowy  && pt[1] <= highy )
		return true;
	else
		return false;
}
bool dashLineRecovery(vector<Point2i> &edgePositions, Vec2i pt1, Vec2i pt2, vector<Vec3f> &circle_candidates, bool plflag = false, bool pcflag = false, bool ppflag = false)
{
	Vec4i line = { pt1[0], pt1[1], pt2[0], pt2[1] }; int c = 0;
	float len = p2pdistance(pt1, pt2); int x1, x2,y1,y2;
	Vec2i mid = (pt1 + pt2) / 2;
	for (auto i = 0; i < circle_candidates.size(); i++)
	{
		Vec3f c = circle_candidates[i];
		Vec2i center = { int(c[0]), int(c[1]) };
		float radius = c[2];
		if (on_circle(pt1, c) && on_circle(pt2, c) && len < 25 && !ppflag)
			return true;
		if (on_circle(mid, c) && len < 20)
		{
			return true;
		}
		if (p2pdistance(center, pt1) - radius < 10 && p2pdistance(center, pt2)  - radius < 10 && abs(point2Line(line, center) - radius) < 5 && plflag)
			return true;
		if (pcflag&&on_circle(pt1, c) && len < 30)
			return true;
	}
	
	
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
	float alpha = 0.45;
	if (x2 - x1 > y2 - y1)
	{
		sort(edgePositions.begin(), edgePositions.end(), [](Vec2i a, Vec2i b){return a[0] < b[0]; });
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
		if (low_threshold < 5)
			low_threshold = x2 - x1;
		int ptsCount = count(bitmap, bitmap + N, 1);
		if (ptsCount > low_threshold)
		{
			cout << "exist dashline between " << pt1<<" "<<pt2<<endl;
			return true;
		}
		else
		{
			cout << "not exist dashline between " << pt1<<" "<<pt2<<endl;
			return false;
		}
	}
	else
	{
		sort(edgePositions.begin(), edgePositions.end(), [](Vec2i a, Vec2i b){return a[1] < b[1]; });
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
		if (low_threshold < 8)
			low_threshold = y2 - y1;
		int ptsCount = count(bitmap, bitmap + N, 1);
		if (ptsCount > low_threshold)
		{
			cout << "exist dashline between " << pt1<<" "<<pt2<<endl;
			return true;
		}
		else
		{
			cout << "not exist dashline between " << pt1<<" "<<pt2<<endl;
			return false;
		}
	}
		

}



bool isParallel(Vec4i line1, Vec4i line2)
{
	if (line1 == line2)
		return false;
	Vec2i line1V = { line1[2] - line1[0], line1[3] - line1[1] }; Vec2i line2V = { line2[2] - line2[0], line2[3] - line2[1] };
	double theta1, theta2;
	if (abs(line1V[0]) <= 3)
	{
		theta1 = CV_PI / 2.0;
	}
	else if (abs(line1V[1]) <= 3)
	{
		theta1 = 0;
	}
	else
	{
		theta1 = (atan2(line1V[1], line1V[0]));
	}
	if (abs(line2V[0]) <= 3)
	{
		theta2 = CV_PI / 2.0;
	}
	else if (abs(line2V[1]) <= 3)
	{
		theta2 = 0;
	}
	else
	{
		theta2 = (atan2(line2V[1], line2V[0]));
	}

	/*double theta1 = abs((abs(line1V[0]) <= 3) ? CV_PI / 2.0 : atan2(line1V[1], line1V[0]));
	double theta2 = abs((abs(line2V[0]) <= 3) ? CV_PI / 2.0 : atan2(line2V[1], line2V[0]));*/
	double angle = abs(theta1 - theta2) / CV_PI * 180;
	return (angle < 8);
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
bool existRealLineWithinPtxs(vector<lineX> &lineXs, pointX ptx1, pointX ptx2)
{
	//Vec4i tmpLxy = { ptx1.px, ptx1.py, ptx2.px, ptx2.py };
	//auto iter = find_if(lineXs.begin(), lineXs.end(), [&](lineX a){return a.lxy == tmpLxy; });
	//if (iter != lineXs.end())
	//	return true;
	//else
	//{
	//	//not find exact line
	//	auto iter2 = find_if(lineXs.begin(), lineXs.end(), [&](lineX a){return isPartOfLongerLine(tmpLxy, a.lxy); });
	//	if (iter2 != lineXs.end())
	//		return true;
	//	else
	//	{
	//		//not part of line
	//	}
	//}
	Vec4i line1, line2;
	line1 = { ptx1.px, ptx1.py, ptx2.px, ptx2.py }; line2 = { ptx2.px, ptx2.py, ptx1.px, ptx1.py };
	//auto iter1 = find_if(lineXs.begin(), lineXs.end(), [&](lineX a){return in_line(a.lxy, ptx1.pxy); });
	//auto iter2 = find_if(lineXs.begin(), lineXs.end(), [&](lineX a){return in_line(a.lxy, ptx2.pxy); });
	auto iter = find_if(lineXs.begin(), lineXs.end(), [&](lineX a){
		if (a.lxy == line1 || a.lxy == line2 || (on_line(a.lxy, ptx1.pxy) && on_line(a.lxy, ptx2.pxy)))
			return true;
		else
			return false;
	});
	if (iter != lineXs.end())
		return true;
	else
		return false;
}


void detect_line3(Mat diagram_segwithoutcircle, Mat &withoutCirBw, vector<Point2i> &edgePoints,vector<Vec3f> &circle_candidates, Mat &color_img, vector<Vec4i> &plainLines, vector<Vec2i>& plainPoints, Mat &drawedImages, bool showFlag = true, string fileName = "")
{
	vector<Point2i> ept = getPointPositions(withoutCirBw);
#pragma region region1
	
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

#pragma endregion region1

	/* then we handle the lines specifically*/
	/* for lines should be combined 1. colinear 2. */
	
	/*for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	}*/
	//namedWindow("5.lines first opt version now", 0); imshow("5.lines first opt version now", color_img);

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
		cout << plainLines[i][0] << "," << plainLines[i][1] << "," << plainLines[i][2] << "," << plainLines[i][3] << endl;
	}
	for (auto iter1 = plainLines.begin(); iter1 != plainLines.end(); iter1++)
	{
		// first combine collinear line
		int maxXDiff = 0;
		for (auto iter2 = iter1 + 1; iter2 != plainLines.end();)
		{
			Vec4i line1 = *iter1; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };
			cout << pt1 << "pt" << pt2 << endl;
		
			Vec4i line2 = *iter2; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
			cout << pt3 << "pt" << pt4 << endl;

			if (isParallel(line1, line2))
			{
				//if it's parallel
				cout << "line1 and line2 is parallel" << endl;
				vector<Vec2i> tmpPtVec;
				tmpPtVec.push_back(pt1); tmpPtVec.push_back(pt2);
				tmpPtVec.push_back(pt3); tmpPtVec.push_back(pt4);
				bool flag0 = false;
				if ( flag0 =(abs(pt1[0] - pt2[0]) < 3))
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

				if (ptLine(line1, pt3))// pt3 of line2 is on line 1
				{
					// if it's collinear
					cout << "two line is collinear" << endl;
					int eps = 5;
					bool flag3, flag4, flag5;
					if (!flag0)
					{
						flag3 = ((pt3[0] - pt1[0])*(pt3[0] - pt2[0]) > 0); //appoximate if pt3 is in line1
						flag4 = ((pt4[0] - pt1[0])*(pt4[0] - pt2[0]) > 0);//approximate if pt4 is in  line1
						
					}
					else
					{
						flag3 = ((pt3[1] - pt1[1])*(pt3[1] - pt2[1]) > 0); //appoximate if pt3 is in line1
						flag4 = ((pt4[1] - pt1[1])*(pt4[1] - pt2[1]) > 0);//approximate if pt4 is in  line1					
					}

					//if (same_pt(pt1, pt3))
					//{
					//	cout << "pt1 and pt3 seems the same" << endl;
					//	if (!flag0)
					//	{
					//		if (abs(pt4[0] - pt1[0]) > abs(pt2[0] - pt1[0]))
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt4[0], pt4[1] };
					//			cout << "pt4 is farther to pt1 than pt2, and the line1 changed to " << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt2 is farther to pt1 than pt4, and the line1 changed to " << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	else
					//	{
					//		if (abs(pt4[1] - pt1[1]) > abs(pt2[1] - pt1[1]))
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt4[0], pt4[1] };
					//			cout << "pt4 is farther to pt1 than pt2, and the line1 changed to " << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt2 is farther to pt1 than pt4, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	cout << "erase line2" << endl << endl;
					//	iter2 = plainLines.erase(iter2);
					//}
					//else if (same_pt(pt1, pt4))
					//{
					//	if (!flag0)
					//	{
					//		if (abs(pt3[0] - pt1[0]) > abs(pt2[0] - pt1[0]))
					//		{
					//			cout << "pt3 is farther to pt1 than pt2, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//			*iter1 = { pt1[0], pt1[1], pt3[0], pt3[1] };
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt2 is farther to pt1 than pt2, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	else
					//	{
					//		if (abs(pt3[1] - pt1[1]) > abs(pt2[1] - pt1[1]))
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt3[0], pt3[1] };
					//			cout << "pt3 is farther to pt1 than pt2, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt2 is farther to pt1 than pt2, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	cout << "erase line2" << endl << endl;
					//	iter2 = plainLines.erase(iter2);
					//}
					//else if (same_pt(pt2, pt3))
					//{
					//	if (!flag0)
					//	{
					//		if (abs(pt4[0] - pt2[0]) > abs(pt1[0] - pt2[0]))
					//		{
					//			*iter1 = { pt4[0], pt4[1], pt2[0], pt2[1] };
					//			cout << "pt4 is farther to pt2 than pt1, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt1 is farther to pt2 than pt4, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	else
					//	{
					//		if (abs(pt4[1] - pt2[1]) > abs(pt1[1] - pt2[1]))
					//		{
					//			*iter1 = { pt4[0], pt4[1], pt2[0], pt2[1] };
					//			cout << "pt1 is farther to pt2 than pt4, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;

					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt1 is farther to pt2 than pt4, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	cout << "erase line2" << endl << endl;
					//	iter2 = plainLines.erase(iter2);
					//}
					//else if (same_pt(pt2, pt4))
					//{
					//	if (!flag0)

					//	{
					//		if (abs(pt3[0] - pt2[0]) > abs(pt1[0] - pt2[0]))
					//		{
					//			*iter1 = { pt3[0], pt3[1], pt2[0], pt2[1] };
					//			cout << "pt3 is farther to pt2 than pt1, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt1 is farther to pt2 than pt3, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	else
					//	{
					//		if (abs(pt3[1] - pt2[1]) > abs(pt1[1] - pt2[1]))
					//		{
					//			*iter1 = { pt3[0], pt3[1], pt2[0], pt2[1] };
					//			cout << "pt3 is farther to pt2 than pt1, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//		else
					//		{
					//			*iter1 = { pt1[0], pt1[1], pt2[0], pt2[1] };
					//			cout << "pt1 is farther to pt2 than pt3, and the line1 changed to" << (*iter1)[0] << "," << (*iter1)[1]
					//				<< "," << (*iter1)[2] << "," << (*iter1)[3] << endl;
					//		}
					//	}
					//	cout << "erase line2" << endl << endl;
					//	iter2 = plainLines.erase(iter2);
					//}


					// four points are all different
					if (flag3&&flag4&&!same_pt(pt1,pt3)&&!same_pt(pt1,pt4)&&!same_pt(pt2,pt3)&&!same_pt(pt2,pt4))
					{
						// two disjoint collinear line
						cout << "four point are all different,";
						cout << firstPt << "range" << fourthPt << endl;
						if (dashLineRecovery(ept, secondPt, thirdPt, circle_candidates, true,false,false))
						{
							*iter1 = { firstPt[0], firstPt[1], fourthPt[0], fourthPt[1] };
							cout << "collinear line recovery, now the line1 is " << *iter1 << "and erase line2" << endl << endl;
							iter2 = plainLines.erase(iter2);
						}
						else
						{
							cout << "collinear but two differnent line" << endl << endl;
							iter2++;
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
					iter2++;
				}
			}
			else
			{
				cout << "line1 and line2 is not parallel" << endl << endl;
				iter2++;
			}
		}
		cout << endl;
	}
	

#pragma endregion rmParallel
	
	//for (auto i = 0; i < plainLines.size(); i++)
	//{
	//	Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
	//	cout << "*********************" << pt1 << " " << pt2 << endl;
	//	line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
	//	Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
	//	circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
	//	circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	//}
	//namedWindow("6.lines first opt version now", 0); imshow("6.lines first opt version now", color_img);

	bool testflag = true;
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
#pragma region recover line-based dash line and point refinement
	
#pragma region cross point combination
	vector<Vec4i> tmpLines; vector<lineX> lineXs; vector<pointX> pointXs;
	for (int i = 0; i < plainLines.size(); i++)
	{
		Vec4i line1; lineX lineX1; Vec2i pt1, pt2;
		pointX p1, p2; p1.l_idx.push_back(i); p2.l_idx.push_back(i); 
		for (int j = i + 1; j < plainLines.size(); j++)
		{
			line1 = plainLines[i];
			pt1 = { line1[0], line1[1] }; pt2 = { line1[2], line1[3] }; 
			Vec4i line2 = plainLines[j]; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
			bool flag00 = (abs(pt1[0] - pt2[0]) < 3); bool flag01 = (abs(pt3[0] - pt4[0]) < 3);
			int maxD1, maxD2;
			maxD1 = (abs(pt1[0] - pt2[0]) < 3) ? abs(pt2[1] - pt1[1]) : abs(pt2[0] - pt1[0]);
			maxD2 = (abs(pt4[0] - pt3[0]) < 3) ? abs(pt4[1] - pt3[1]) : abs(pt4[0] - pt3[0]);
			/*if (!flag00)
			{
				maxD1 = abs(pt2[0] - pt1[0]);
			}
			else
			{
				maxD1 = abs(pt2[1] - pt1[1]);
			}
			if (!flag01)
			{
				maxD2 = abs(pt4[0] - pt3[0]);
			}
			else
			{
				maxD2 = abs(pt4[1] - pt3[1]);
			}*/

			cout << endl<<"line1 is " << line1 << endl << "line2 is " << line2 << endl;
			if (isParallel(line1, line2))
			{
				cout << "line1 and line2 is parallel" << endl;
				continue;
			}
			else
			{
				// not parallel
				cout << "line1 and line2 is not parallel" << endl;
				Vec2i cross; getCrossPtRev(line1, line2, cross,circle_candidates, plainLines);
				//Vec2f cross; getCrossPt(line1, line2, cross);
				cout << "cross point is now " << cross << endl;
				if (!isInImage(color_img.cols, color_img.rows, cross))
				{
					// take it as no cross
					cout << "the cross point is out of scope" << endl;
					continue;
				}
				else
				{
					//cross point in image
					//first check if the cross is same with 
					int flag1, flag2; flag1 = flag2 = -1;
					if (!flag00)
					{
						// use x;
						cout << "line1 is not vertical, use x" << endl;
						bool tmpFlag;
						if (abs(cross[0] - pt1[0]) < abs(cross[0] - pt2[0]))
						{

							// pt1 is closer to cross
							cout << "pt1 is closer to the cross" << endl;
							flag1 = 0;
							if ( (tmpFlag = withinPtCRegion(cross, pt1)) || dashLineRecovery(ept, cross, pt1, circle_candidates,false,true,false))
							{
								if (tmpFlag||abs(cross[0] - pt2[0]) >= maxD1)
								{
									plainLines[i][0] = cross[0]; plainLines[i][1] = cross[1];
									line1 = plainLines[i];
									cout << "the line1 is now " << line1 << endl;
									maxD1 = abs(cross[0] - pt2[0]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}


						}
						else
						{
							// pt2 is closer to the cross
							cout << "pt2 is closer to the cross" << endl;
							flag1 = 1;
							if ((tmpFlag = withinPtCRegion(cross, pt2)) || dashLineRecovery(ept, cross, pt2, circle_candidates, false, true, false))
							{
								if (tmpFlag || abs(cross[0] - pt1[0]) >= maxD1)
								{
									//chooseNearCircle(circle_candidates, cross, pt2);
									plainLines[i][2] = cross[0]; plainLines[i][3] = cross[1];
									line1 = plainLines[i];
									cout << "the line1 is now " << line1 << endl;
									maxD1 = abs(cross[0] - pt1[0]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}
						}
					}
					else
					{
						//use y
						cout << "line1 is vertical, use y" << endl;
						bool tmpFlag;
						if (abs(cross[1] - pt1[1]) < abs(cross[1] - pt2[1]))
						{

							// pt1 is closer to cross
							cout << "pt1 is closer to the cross" << endl;
							flag1 = 0;
							if ((tmpFlag = withinPtCRegion(cross, pt1)) || dashLineRecovery(ept, cross, pt1, circle_candidates, false, true, false))
							{
								if (tmpFlag || abs(cross[1] - pt2[1]) >= maxD1)
								{
									plainLines[i][0] = cross[0]; plainLines[i][1] = cross[1];
									line1 = plainLines[i];
									cout << "the line1 is now " << line1 << endl;
									maxD1 = abs(cross[1] - pt2[1]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}


						}
						else
						{
							// pt2 is closer to the cross
							cout << "pt2 is closer to the cross" << endl;
							flag1 = 1;
							if ((tmpFlag = withinPtCRegion(cross, pt2)) || dashLineRecovery(ept, cross, pt2, circle_candidates, false, true, false))
							{
								if (tmpFlag || abs(cross[1] - pt1[1]) >= maxD1)
								{
									//chooseNearCircle(circle_candidates, cross, pt2);
									plainLines[i][2] = cross[0]; plainLines[i][3] = cross[1];
									line1 = plainLines[i];
									cout << "the line1 is now " << line1 << endl;
									maxD1 = abs(cross[1] - pt1[1]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}
						}
					}
					if (!flag01)
					{
						cout << "line2 is not vertical, use x" << endl;
						bool tmpFlag;
						if (abs(cross[0] - pt3[0]) < abs(cross[0] - pt4[0]))
						{
							//pt3 is closer to the cross
							cout << "pt3 is closer to the cross" << endl;
							flag2 = 0;
							if ((tmpFlag = withinPtCRegion(cross, pt3)) || dashLineRecovery(ept, cross, pt3, circle_candidates, false, true, false))
							{
								if (tmpFlag || abs(cross[0] - pt4[0]) >= maxD2)
								{
									//chooseNearCircle(circle_candidates, cross, pt3);
									plainLines[j][0] = cross[0]; plainLines[j][1] = cross[1];
									line2 = plainLines[j];
									cout << "the line2 is now " << line2 << endl;
									maxD2 = abs(cross[0] - pt3[0]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}
						}
						else
						{
							//pt4 is closer to the cross
							cout << "pt4 is closer to the cross" << endl;
							flag2 = 1;
							if ((tmpFlag = withinPtCRegion(cross, pt4)) || dashLineRecovery(ept, cross, pt4, circle_candidates, false, true, false))
							{
								if (tmpFlag || abs(cross[0] - pt3[0]) >= maxD2)
								{
									//chooseNearCircle(circle_candidates, cross, pt4);
									plainLines[j][2] = cross[0]; plainLines[j][3] = cross[1];
									line2 = plainLines[j];
									cout << "the line2 is now " << line2 << endl;
									maxD2 = abs(cross[0] - pt3[0]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}
						}
					}
					else
					{
						cout << "line2 is vertical, use y" << endl;
						bool tmpFlag;
						if (abs(cross[1] - pt3[1]) < abs(cross[1] - pt4[1]))
						{
							//pt3 is closer to the cross
							cout << "pt3 is closer to the cross" << endl;
							flag2 = 0;
							if ((tmpFlag = withinPtCRegion(cross, pt3)) || dashLineRecovery(ept, cross, pt3, circle_candidates, false, true, false))
							{
								if (tmpFlag || abs(cross[1] - pt4[1]) >= maxD2)
								{
									//chooseNearCircle(circle_candidates, cross, pt3);
									plainLines[j][0] = cross[0]; plainLines[j][1] = cross[1];
									line2 = plainLines[j];
									cout << "the line2 is now " << line2 << endl;
									maxD2 = abs(cross[1] - pt3[1]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}
						}
						else
						{
							//pt4 is closer to the cross
							cout << "pt4 is closer to the cross" << endl;
							flag2 = 1;
							if ((tmpFlag = withinPtCRegion(cross, pt4)) || dashLineRecovery(ept, cross, pt4, circle_candidates, false, true, false))
							{
								if (tmpFlag||abs(cross[1] - pt3[1]) >= maxD2)
								{
									//chooseNearCircle(circle_candidates, cross, pt4);
									plainLines[j][2] = cross[0]; plainLines[j][3] = cross[1];
									line2 = plainLines[j];
									cout << "the line2 is now " << line2 << endl;
									maxD2 = abs(cross[1] - pt3[1]);
								}
								else
								{
									cout << "no change" << endl;
								}
							}
						}
					}
					
					//{
					//	if (abs(cross[0] - pt1[1]) < abs(cross[1] - pt2[1]))
					//	{
					//		// pt1 is closer to cross
					//		flag1 = 0;
					//		if (withinPtCRegion(cross, pt1) || dashLineRecovery(edgePoints, cross, pt1))
					//		{
					//			if (abs(cross[1] - pt2[1]) >= maxD1)
					//			{
					//				//chooseNearCircle(circle_candidates, cross, pt1);
					//				plainLines[i][0] = cross[0]; plainLines[i][1] = cross[1];
					//				maxD1 = abs(cross[1] - pt2[1]);
					//			}
					//		}
					//	}
					//	else
					//	{
					//		// pt2 is closer to the cross
					//		flag1 = 1;
					//		if (withinPtCRegion(cross, pt2) || dashLineRecovery(edgePoints, cross, pt2))
					//		{
					//			if (abs(cross[1] - pt1[1]) >= maxD1)
					//			{
					//				//chooseNearCircle(circle_candidates, cross, pt2);
					//				plainLines[i][2] = cross[0]; plainLines[i][3] = cross[1];
					//				maxD1 = abs(cross[1] - pt1[1]);
					//			}
					//		}
					//	}
					//	if (abs(cross[1] - pt3[1]) < abs(cross[1] - pt4[1]))
					//	{
					//		//pt3 is closer to the cross
					//		flag2 = 0;
					//		if (withinPtCRegion(cross, pt3) || dashLineRecovery(edgePoints, cross, pt3))
					//		{
					//			if (abs(cross[1] - pt4[1]) >= maxD2)
					//			{
					//				//chooseNearCircle(circle_candidates, cross, pt3);
					//				plainLines[j][0] = cross[0]; plainLines[j][1] = cross[1];
					//				maxD2 = abs(cross[1] - pt3[1]);
					//			}
					//		}
					//	}
					//	else
					//	{
					//		//pt4 is closer to the cross
					//		flag2 = 1;
					//		if (withinPtCRegion(cross, pt4) || dashLineRecovery(edgePoints, cross, pt4))
					//		{
					//			if (abs(cross[1] - pt3[1]) >= maxD2)
					//			{
					//				//chooseNearCircle(circle_candidates, cross, pt4);
					//				plainLines[j][2] = cross[0]; plainLines[j][3] = cross[1];
					//				maxD2 = abs(cross[1] - pt3[1]);
					//			}
					//		}
					//	}
					//}
				}


			}
		}
		line1 = plainLines[i]; pt1 = { line1[0], line1[1] }; pt2 = { line1[2], line1[3] };
		
		p1.px = pt1[0]; p1.py = pt1[1];p1.pxy =pt1 ; p1.p_idx = 2*i; pointXs.push_back(p1);
		p2.px = pt2[0]; p2.py = pt2[1];p2.pxy =pt2 ; p2.p_idx = 2*i + 1; pointXs.push_back(p2);

		
		
		lineX1.l_idx = i; lineX1.lxy = line1;
		
		lineX1.pt1 = pt1; lineX1.pt2 = pt2; lineX1.pidx1 = p1.p_idx; lineX1.pidx2 = p2.p_idx;
		
		lineX1.px1 = line1[0]; lineX1.py1 = line1[1]; lineX1.px2 = line1[2]; lineX1.py2 = line1[3]; 

		
		
		lineXs.push_back(lineX1);
	}
	/*for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	}*/
	int pxsize = pointXs.size();
	int erasenum = 0;
	for (auto i = 0; i < pointXs.size(); i++)
	{
		pointX p1 = pointXs[i];
		int erasenum1 = 0;
		//int div1 = i / 2; int mod1 = i % 2;
		for (int j = i + 1; j < pointXs.size(); j++)
		{
			pointX p2 = pointXs[j];
			if (same_pt(p1,p2))
			{
				//erase point2
				int div2 = p2.p_idx/ 2;
				int mod2 = p2.p_idx % 2;
				if (mod2)
				{
					//second point in line to be changed
					lineXs[div2].pt2 = p1.pxy;
					lineXs[div2].px2 = p1.pxy[0]; lineXs[div2].py2 = p1.pxy[1];
					lineXs[div2].lxy[2] = p1.pxy[0];  lineXs[div2].lxy[3] = p1.pxy[1];
					lineXs[div2].pidx2 = p1.p_idx;
					//pointXs[j - erasenum].pxy = p1.pxy;
					pointXs.erase(pointXs.begin() + j);
					j--;
					erasenum++; erasenum1++;
					
				}
				else
				{
					//first point in line to be changed
					lineXs[div2].pt1 = p1.pxy; 
					lineXs[div2].px1 = p1.pxy[0]; lineXs[div2].py1 = p1.pxy[1];
					lineXs[div2].lxy[0] = p1.pxy[0];  lineXs[div2].lxy[1] = p1.pxy[1];
					lineXs[div2].pidx1 = p1.p_idx;
					pointXs.erase(pointXs.begin() + j );
					j--;
					erasenum++; erasenum1++;
				}
			}
		}
	}
//#pragma region drag close to circle
//	for (auto i = 0; i < pointXs.size(); i++)
//	{
//		pointX ptx = pointXs[i];
//		for (auto j = 0; j < circle_candidates.size(); j++)
//		{
//			Vec3f c = circle_candidates[j]; 
//			Vec2f center = { c[0], c[1] }; float radius = c[2];
//			if (on_circle(ptx.pxy, c))
//			{
//				 c is thought to be on circle
//				float minDiff = 100;
//				for (auto m = ptx.px - 3; m <= ptx.px + 3; m++)
//				{
//					for (auto n = ptx.py - 3; n <= ptx.py + 3; n++)
//					{
//						Vec2i tmp = { m, n };
//						float tmpDiff = abs(p2pdistance(tmp, center) - radius);
//						if (tmpDiff< minDiff)
//						{
//							minDiff = tmpDiff;
//							ptx.pxy = tmp;
//							int div = ptx.p_idx / 2;
//							int mod = ptx.p_idx % 2; 
//							if (!mod)
//							{
//								mod == 0
//								lineXs[div].pt1 = tmp;
//							}
//							else
//							{
//								lineXs[div].pt2 = tmp;
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//#pragma endregion drag close to circle
	
#pragma endregion cross point combination
	
	//for (auto i = 0; i < lineXs.size(); i++)
	//{
	//	Vec4i l = lineXs[i].lxy;
	//	Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
	//	//cout << "*********************" << pt1 << " " << pt2 << endl;
	//	line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
	//	Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
	//	circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
	//	circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	//}
	//namedWindow("7.lines first opt version now", 0); imshow("7.lines first opt version now", color_img);

	for (auto i = 0; i < lineXs.size(); i++)
	{
		Vec4i l = lineXs[i].lxy;
		Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		cout << pt1 << " " << pt2 << endl;
	}
	cout << endl<<"block stop" << endl << endl;

	Mat withoutCLBw = withoutCirBw;
	for (auto i = 0; i < lineXs.size(); i++)
	{
		line(withoutCLBw, lineXs[i].pt1, lineXs[i].pt2, 0, 2);
	}
	vector<Point2i> ept2 = getPointPositions(withoutCLBw);
	for (auto i = 0; i < pointXs.size(); i++)
	{
		pointX ptx1 = pointXs[i];
		for (auto j = i + 1; j < pointXs.size(); j++)
		{
			pointX ptx2 = pointXs[j];
			if (!existRealLineWithinPtxs(lineXs, ptx1, ptx2))
			{
				if (dashLineRecovery(ept2, ptx1.pxy, ptx2.pxy, circle_candidates, false,false,true))
				{
					Vec4i tmpLine = { ptx1.px, ptx1.py, ptx2.px, ptx2.py };
					lineX tmpLx;
					tmpLx.lxy = tmpLine; tmpLx.pidx1 = ptx1.p_idx; tmpLx.pidx2 = ptx2.p_idx;
					tmpLx.l_idx = lineXs.size(); tmpLx.pt1 = ptx1.pxy; tmpLx.pt2 = ptx2.pxy;
					tmpLx.px1 = ptx1.px; tmpLx.py1 = ptx1.py; tmpLx.px2 = ptx2.px; tmpLx.py2 = ptx2.py;
					lineXs.push_back(tmpLx);
				}
			}
			else
			{
				continue;
			}
		}
	}
#pragma endregion recover line-based dash line and point refinement
	
	//for (auto i = 0; i < plainLines.size(); i++)
	//{
	//	Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
	//	cout << "*********************"<<pt1 << " " << pt2 << endl;
	//	line(color_img, pt1, pt2, Scalar(rand()%255, rand()%255, rand()%255), 2, 8, 0);
	//	Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
	//	circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
	//	circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	//}
	//


	for (auto i = 0; i < lineXs.size(); i++)
	{
		Vec4i l = lineXs[i].lxy;
		Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		//cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp);
		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp);
	}
	//
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


void primitive_parse(const Mat binarized_image, const Mat diagram_segment, vector<Point2i> &edgePoints,vector<pointX> &points, vector<lineX> &lines, vector<circleX> &circles, Mat &drawedImages, bool showFlag=true, string fileName="")
{
	/* primitive about points, lines, and circles
	first detect the circle, we can get the matrix of diagram segments without circle for the next
	 processing step*/
	vector<Vec3f> circle_candidates = {}; 
	Mat color_img = Mat::zeros(diagram_segment.size(),CV_8UC3); cvtColor(diagram_segment, color_img, CV_GRAY2RGB);
	Mat diagram_segwithoutcircle = Mat::zeros(diagram_segment.size(), CV_8UC1);
	
	//diagram_segment = preprocessing(diagram_segment);

	/*ransac go*/
	Mat withoutCirBw = binarized_image;
	detect_circle(diagram_segment, color_img, diagram_segwithoutcircle, withoutCirBw, circle_candidates, showFlag);
	// then the line detection
	vector<Vec4i> line_candidates = {}; vector<Vec2i> basicEndpoints = {};
	detect_line3(diagram_segwithoutcircle, withoutCirBw, edgePoints, circle_candidates, color_img, line_candidates, basicEndpoints, drawedImages, showFlag, fileName);
	
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
	Mat image = imread("Sg-90.jpg", 0);
	//namedWindow("original image");
	//imshow("original image", image);
	// then binarize it

	Mat binarized_image = image_binarizing(image, true);
	//binarized_image = preprocessing(binarized_image);
	// then go on a process of connectivity componnent analysis
	int labeln; Mat diagram_segment; vector<Mat> label_segment;
	vector<Point2i> edgePoints;
	edgePoints = getPointPositions(binarized_image);
	//Mat pointss = Mat::zeros(1000, 1000, CV_8UC3);
	//for (auto i = 0; i < edgePoints.size(); i++)
	//{
	//	Point2i pt = edgePoints[i];
	//	circle(pointss, pt, 1, Scalar(0, 0, 255));
	//}
	//namedWindow("points"); imshow("points", pointss);
	image_labelling(binarized_image, diagram_segment,true);
	vector<pointX> points; vector<lineX> lines; vector<circleX> circles;		Mat drawedImages;
	primitive_parse(binarized_image,diagram_segment, edgePoints, points, lines, circles, drawedImages,false);
	return 0;
}

int diagram()
{
	//a series of image
	//vector<Mat> images;
	char abs_path[100] = "D:\\data\\graph-DB\\newtest13";
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
		vector<pointX> points = {}; vector<lineX> lines = {}; vector<circleX> circles = {};
		Mat drawedImages = Mat::zeros(diagram_segment.size(), CV_8UC3);
		primitive_parse(binarized_image,diagram_segment, edgePoints, points, lines, circles, drawedImages, false);
		imwrite(saveimgName, drawedImages);
	}
	return 0;
}