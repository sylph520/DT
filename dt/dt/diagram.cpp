# include "stdafx.h"
# include "diagram.h"
#include <filesystem>

#define N 500

bool dashLineRecovery(vector<Point2i>&, Vec2i, Vec2i, Vec2i, vector<circle_class>&, bool plflag, bool pcflag, bool ppflag);

bool dashLineRecovery(vector<Point2i>&, Vec2i, Vec2i, vector<circle_class>&, bool plflag, bool pcflag, bool ppflag);

/*deserted now*/
//void chooseNearCircle(vector<Vec3f>& circle_candidates, Vec2i& cross, Vec2i& pt)
//{
//	for (auto i = 0; i < circle_candidates.size(); i++)
//	{
//		Vec3f circle = circle_candidates[i];
//		Vec2i center = {int(circle[0]), int(circle[1])};
//		if (p2pdistance(pt, center) - circle[2] < 5)
//		{
//			if (p2pdistance(pt, center) > p2pdistance(cross, center))
//			{
//				cross = pt;
//			}
//		}
//	}
//}
//bool isPartOfLongerLine(Vec4i line1, Vec4i line2)
//{
//	//check if line1 is part of line2;
//	Vec2i pt1, pt2;
//	pt1 = {line1[0], line1[1]};
//	pt2 = {line1[2], line1[3]};
//	if (in_line(line2, pt1) && in_line(line2, pt2))
//		return true;
//	else
//		return false;
//}
//bool nearToAcross(Vec2i pt, vector<Vec2i>& crossPts)
//{
//	auto iter = find_if(crossPts.begin(), crossPts.end(), [&](Vec2i a)
//	                    {
//		                    if (same_pt(pt, a))
//			                    return true;
//		                    else
//			                    return false;
//	                    });
//	if (iter != crossPts.end())
//		return true;
//	else
//		return false;
//}
//Mat preprocessing(Mat diagram_segment)
//{
//	//namedWindow("origianl diagram seg"); imshow("original diagram seg", diagram_segment);
//	//now we'll do some preprocessing to make for the following detection procedure
//	Mat eroded, dilated, opened, closed;
//	Mat eroEle, dilEle;
//	eroEle = getStructuringElement(MORPH_RECT, Size(3, 3));
//	morphologyEx(diagram_segment, dilated, MORPH_DILATE, eroEle);
//	//erode(diagram_segment, eroded, eroEle);
//	namedWindow("dilated image");
//	imshow("dilated image", dilated);
//	//namedWindow("eroded diagram"); imshow("eroded diagram", eroded);
//	morphologyEx(diagram_segment, closed, MORPH_CLOSE, eroEle);
//	namedWindow("closed diagram");
//	imshow("closed diagram", closed);
//	return closed;
//}
//bool sameLine(Vec4i line1, Vec4i line2)
//{
//	Vec2i pt1 = {line1[0], line1[1]};
//	Vec2i pt2 = {line1[2], line1[3]};
//	Vec2i pt3 = {line2[0], line2[1]};
//	Vec2i pt4 = {line2[0], line2[1]};
//	bool flag1 = (same_pt(pt1, pt3) && same_pt(pt2, pt4));
//	bool flag2 = (same_pt(pt1, pt4) && same_pt(pt2, pt3));
//	if (flag1 || flag2)
//		return true;
//	else
//		return false;
//}
//double slopeEval(Vec4i line)
//{
//	double diff1, diff2, diff;
//	diff1 = lSlope(line) - 0;
//	diff2 = lSlope(line) - 90;
//	diff = (diff1 < diff2) ? diff1 : diff2;
//	return diff;
//}
//bool on_other_noncollinearlines(vector<Vec4i> plainLines, Vec4i temp_line, Vec2i pt)
//{
//	for (size_t i = 0; i < plainLines.size(); i++)
//	{
//		Vec2i pt1 = { plainLines[i][0], plainLines[i][1] };
//		Vec2i pt2 = { plainLines[i][2], plainLines[i][3] };
//		if ((pt == pt1) || (pt == pt2))
//		{
//			if ((!on_line(temp_line, pt1)) || (!on_line(temp_line, pt2)))
//			{
//				return true;
//			}
//		}
//	}
//	return false;
//}
//void findLineEnds(vector<Vec2i> colpoints, vector<Vec2i>& lineEnds, vector<Vec2i>& plainPoints, vector<Vec4i> plainLines)
//{
//	vector<distance_info> distance_infos;
//	for (vector<Vec2i>::iterator iter1 = colpoints.begin(); iter1 != colpoints.end(); ++iter1)
//	{
//		Vec2i pt1 = *iter1;
//		for (vector<Vec2i>::iterator iter2 = iter1 + 1; iter2 != colpoints.end(); ++iter2)
//		{
//			Vec2i pt2 = *iter2;
//			double tempDistance = p2pdistance(pt1, pt2);
//			distance_info dis_info = {pt1, pt2, tempDistance};
//			distance_infos.push_back(dis_info);
//		}
//	}
//	std::sort(distance_infos.begin(), distance_infos.end(),
//	          [](distance_info a, distance_info b)
//	          {
//		          return (a.distance > b.distance);
//	          });
//	Vec2i pt3 = distance_infos[0].pt1;
//	Vec2i pt4 = distance_infos[0].pt2;
//	lineEnds.push_back(pt3);
//	lineEnds.push_back(pt4);
//	int count0 = 0;
//	//cout << "colpoints size " << colpoints.size();
//	for (vector<Vec2i>::iterator iter3 = colpoints.begin(); iter3 != colpoints.end(); ++iter3)
//	{
//		Vec2i tempPt1 = *iter3;
//		bool flag = on_other_noncollinearlines(plainLines, {pt3[0], pt3[1], pt4[0], pt4[1]}, tempPt1);
//		for (vector<Vec2i>::iterator iter4 = plainPoints.begin(); iter4 != plainPoints.end();)
//		{
//			Vec2i tempPt2 = *iter4;
//			if (tempPt2 == tempPt1 && (tempPt2 != pt3) && (tempPt2 != pt4) && (!flag))
//			{
//				iter4 = plainPoints.erase(iter4);
//				count0++;
//				//cout << "iter4 value " << ++count0 <<" "<<tempPt2<< endl;
//				break;
//			}
//			else
//			{
//				++iter4;
//				continue;
//			}
//		}
//	}
//	//cout << " count0 " << count0 << endl;
//bool isLSeg(Mat seg)
//{
////	morphologyEx(seg, seg, MORPH_OPEN,getStructuringElement(MORPH_RECT,Size(3,3)));
//	Rect tmpRect = boundingRect(seg); Rect centerRect = tmpRect + Point(tmpRect.width / 4, tmpRect.height / 4) + Size(-tmpRect.width / 2, -tmpRect.height / 2);
//	Rect leftRect = tmpRect + Size(-tmpRect.width / 2, 0); Rect rightRect = leftRect+Point(tmpRect.width/4, 0);
//	Rect topRect = tmpRect + Size(0, -tmpRect.height / 2);Rect bottomRect = topRect + Point(0, tmpRect.height / 4);
//	
//	Mat lRect, rRect, tRect, bRect, cRect,oRect;
//	int fac = 1;
////	resize(Mat(seg, leftRect), lRect, Size(),fac,fac,INTER_NEAREST);resize(Mat(seg, rightRect), rRect, Size(), fac, fac, INTER_NEAREST);
////	resize(Mat(seg, topRect), tRect, Size(), fac, fac, INTER_NEAREST);resize(Mat(seg, bottomRect), bRect, Size(), fac, fac, INTER_NEAREST);
////	resize(Mat(seg, centerRect), cRect, Size(), fac, fac, INTER_NEAREST);	resize(Mat(seg, tmpRect), oRect, Size(), fac, fac, INTER_NEAREST);
//	lRect= Mat(seg, leftRect);rRect= Mat(seg, rightRect);
//	tRect= Mat(seg, topRect);bRect= Mat(seg, bottomRect);
//	cRect= Mat(seg, centerRect);oRect=	Mat(seg, tmpRect);
//	int LRNum, RRNum, TRNum, BRNum,tNum,cNum;
//	tNum = countNonZero(oRect); cNum = countNonZero(cRect);
//	LRNum = countNonZero(lRect);RRNum = countNonZero(rRect);
//	TRNum = countNonZero(tRect);BRNum = countNonZero(bRect);
//	if (abs(LRNum - RRNum) < 5 && abs(TRNum - BRNum) < 5)
//		return true;
//	else
//		return false;
//}
//}
//bool isLSeg(Mat seg)
//{
////	morphologyEx(seg, seg, MORPH_OPEN,getStructuringElement(MORPH_RECT,Size(3,3)));
//	Rect tmpRect = boundingRect(seg); Rect centerRect = tmpRect + Point(tmpRect.width / 4, tmpRect.height / 4) + Size(-tmpRect.width / 2, -tmpRect.height / 2);
//	Rect leftRect = tmpRect + Size(-tmpRect.width / 2, 0); Rect rightRect = leftRect+Point(tmpRect.width/4, 0);
//	Rect topRect = tmpRect + Size(0, -tmpRect.height / 2);Rect bottomRect = topRect + Point(0, tmpRect.height / 4);
//	
//	Mat lRect, rRect, tRect, bRect, cRect,oRect;
//	int fac = 1;
////	resize(Mat(seg, leftRect), lRect, Size(),fac,fac,INTER_NEAREST);resize(Mat(seg, rightRect), rRect, Size(), fac, fac, INTER_NEAREST);
////	resize(Mat(seg, topRect), tRect, Size(), fac, fac, INTER_NEAREST);resize(Mat(seg, bottomRect), bRect, Size(), fac, fac, INTER_NEAREST);
////	resize(Mat(seg, centerRect), cRect, Size(), fac, fac, INTER_NEAREST);	resize(Mat(seg, tmpRect), oRect, Size(), fac, fac, INTER_NEAREST);
//	lRect= Mat(seg, leftRect);rRect= Mat(seg, rightRect);
//	tRect= Mat(seg, topRect);bRect= Mat(seg, bottomRect);
//	cRect= Mat(seg, centerRect);oRect=	Mat(seg, tmpRect);
//	int LRNum, RRNum, TRNum, BRNum,tNum,cNum;
//	tNum = countNonZero(oRect); cNum = countNonZero(cRect);
//	LRNum = countNonZero(lRect);RRNum = countNonZero(rRect);
//	TRNum = countNonZero(tRect);BRNum = countNonZero(bRect);
//	if (abs(LRNum - RRNum) < 5 && abs(TRNum - BRNum) < 5)
//		return true;
//	else
//		return false;
//}



/**************************************************************geometry funcs*******************************/
/********basic eff funcs*******/
Vec4i pt2line(Vec2i pt1, Vec2i pt2)
{
	Vec4i ret = { pt1[0], pt1[1], pt2[0], pt2[1] };
	return ret;
}
Vec4f pt2line(Vec2f pt1, Vec2f pt2)
{
	Vec4f ret = { pt1[0], pt1[1], pt2[0], pt2[1] };
	return ret;
}

inline void line2pt(Vec4i line, Vec2i& pt1, Vec2i& pt2)
{
	pt1 = { line[0], line[1] };
	pt2 = { line[2], line[3] };
}
inline void line2pt(Vec4f line, Vec2f& pt1, Vec2f& pt2)
{
	pt1 = { line[0], line[1] };
	pt2 = { line[2], line[3] };
}

Vec2i f2i(Vec2f fp)
{
	Vec2i retP = { int(fp[0]), int(fp[1]) }; 
	return retP;
}
Vec4i f2i(Vec4f fl)
{
	Vec4i retL= { int(fl[0]), int(fl[1]), int(fl[2]), int(fl[3]) };
	return retL;
}

/**p2p**/
double p2pdistance(Vec2i pt1, Vec2i pt2)
{
	double distance;
	distance = sqrt(powf((pt1[0] - pt2[0]), 2) + powf((pt1[1] - pt2[1]), 2));
	return distance;
}

bool same_pt(point_class pt1, point_class pt2)
{
	double eps = 8;//to be set parameter
	double high_eps = 10;
	if (abs(pt1.getX() - pt2.getX()) < 3 && p2pdistance(pt1.getXY(), pt2.getXY()) < high_eps)
		return true;
	if (p2pdistance(pt1.getXY(), pt2.getXY()) <= eps)
		return true;
	else
		return false;
}
bool same_pt(Vec2i pt1, Vec2i pt2, double theta, bool in_line)
{
	double eps = 6;//to be set parametere
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

bool close_pt(Vec2i pt1, Vec2i pt2)
{
	double eps = 10;//to be set parametere
	double plainDis = p2pdistance(pt1, pt2);
	//cout << "plain dis " << plainDis << " ,adapt eps " << adapt_eps << endl;
	if (plainDis <= eps)
		return true;
	else
		return false;
}

double cross_product(Vec2i a, Vec2i b)
{
	return a[0] * b[1] - a[1] * b[0];
}


/**  p2l **/
float pt2lineDis(Vec4i line, Vec2i pt)
{
	Vec2i bottom = { line[2] - line[0], line[3] - line[1] };
	Vec2i slope = { pt[0] - line[0], pt[1] - line[1] };
	float p2l = abs(cross_product(bottom, slope) / norm(bottom));
	return p2l;
}

float pt2lineDis(Vec4f line, Vec2f pt)
{
	Vec2f pt1, pt2;
	line2pt(line, pt1, pt2);
	Vec2f bottom = { line[2] - line[0], line[3] - line[1] };
	Vec2f slope = { pt[0] - line[0], pt[1] - line[1] };
	float p2l = abs(cross_product(bottom, slope) / p2pdistance(pt1, pt2));
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
	Vec2i slope = { pt[0] - line[0], pt[1] - line[1] };
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
		pt1 = { line[0], line[1] };
		pt2 = { line[2], line[3] };
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
	pt1 = { line[0], line[1] };
	pt2 = { line[2], line[3] };
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
	pt1 = { line[0], line[1] };
	pt2 = { line[2], line[3] };
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
	pt1 = { line[0], line[1] };
	pt2 = { line[2], line[3] };
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
		pt1 = { line[0], line[1] };
		pt2 = { line[2], line[3] };
		if (p2pdistance(pt1, pt) + p2pdistance(pt2, pt) - p2pdistance(pt1, pt2) > 20)
			return true;
		else
			return false;
	}
	else
		return false;
}

double frac_compute(Vec2i pt1, Vec2i pt2, bool vertical_flag)
{
	double ret = vertical_flag ? 1.0*(pt2[0] - pt1[0]) / (pt2[1] - pt1[1]) : 1.0*(pt2[1] - pt1[1]) / (pt2[0] - pt1[0]);
	return ret;
}

double lSlope_r(Vec4i line)
{
	Vec2i pt1, pt2;
	line2pt(line, pt1, pt2);
	Vec2i lineV = { line[2] - line[0], line[3] - line[1] };
	if ((abs(lineV[1]) < 3))
	{
		return 0;
	}
	else if (abs(lineV[0]) < 3)
	{
		return CV_PI / 2;
	}
	int ydelta = lineV[1];
	int xdelta = lineV[0];
	double theta_r = atan2(ydelta, xdelta);
	return theta_r;
}

double llAngle(Vec4i line1, Vec4i line2)
{
	if (line1 == line2)
		return 0;
	double theta1_r, theta2_r;
	theta1_r = lSlope_r(line1);
	theta2_r = lSlope_r(line2);
	if (theta1_r < 0 || theta2_r < 0)
		cout << "stop" << endl;
	double angle = int(abs(theta1_r - theta2_r) / CV_PI * 180);
	return angle;
}

bool isParallel(Vec4i line1, Vec4i line2)
{
	Vec2i pt1, pt2, pt3, pt4;
	line2pt(line1, pt1, pt2);
	line2pt(line2, pt3, pt4);
	Vec2i line1V = { line1[2] - line1[0], line1[3] - line1[1] };
	Vec2i line2V = { line2[2] - line2[0], line2[3] - line2[1] };
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
	if ((abs(line1V[0]) < 3 && abs(line2V[0]) < 3) || (abs(line1V[1]) < 3 && abs(line2V[1]) < 3))
	{
		return true;
	}
//	else
//	{
//		theta1 = (atan2(line1V[1], line1V[0]));
//		theta2 = (atan2(line2V[1], line2V[0]));
//	}

	/*double theta1 = abs((abs(line1V[0]) <= 3) ? CV_PI / 2.0 : atan2(line1V[1], line1V[0]));
	double theta2 = abs((abs(line2V[0]) <= 3) ? CV_PI / 2.0 : atan2(line2V[1], line2V[0]));*/
	double angle = llAngle(line1, line2);
	if (angle > 180)
		cout << "stop" << endl;
	double threshold_ori = 10;
	double threshold_adp = (norm(line1V) < 15 || norm(line2V) < 15) ? 20 : threshold_ori;
	return (angle <= threshold_adp) || (angle >= 180-threshold_adp);
}

bool isParallel2(Vec4i line1,Vec4i line2)
{
	Vec2i line1V = { line1[2] - line1[0], line1[3] - line1[1] };
	Vec2i line2V = { line2[2] - line2[0], line2[3] - line2[1] };
	double len1 = norm(line1V); double len2 = norm(line2V);
	double ret = (line1V[0] * line2V[0] + line1V[1] * line2V[1])/(len1*len2);
	double angle = acos(ret) / CV_PI * 180;
	return angle < 20;
}


/********************preprocessing funcs*************************************************/
Mat image_binarizing(Mat input_image, bool showFlag = false)
{
	/* This function is used for transform image to its binarized image */

	int block_size;
	double c;
	block_size = 11;
	c = 20;
	Mat binarized_image = Mat::zeros(input_image.size(), CV_8UC1);
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

void image_labelling(Mat binarized_image, Mat& diagram_segment, vector<Mat> &char_imgs, bool showFlag = false)
{
	// this function is used to label image with connectedcomponnent analysis, and store the diagram 
	//and label segment
	Mat labeled(binarized_image.size(), CV_8UC3);
	Mat boundImg = binarized_image.clone();
	Mat statsMat, centroidMat;
	Mat labeled_image;
	int labeln = connectedComponentsWithStats(binarized_image, labeled_image, statsMat, centroidMat, 8, 4);
	vector<Mat> charSegs;
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
		int wholeIdx = -1;
		if (seg_area > 5)
		{
			//segments.push_back(seg_mat);
			if (seg_area >= statsMat.at<int>(dia_idx, 4))
				dia_idx = label;
		}
	}

	vector<Rect> bRects;
	map<int, int> matLabelp;
	for (auto label = 1; label < labeln; ++label)
	{
		if (label != dia_idx)
		{
			Mat charImgCandicate = labeled_image == label;
			int seg_area = statsMat.at<int>(label, 4);
			if (seg_area > 10)
			{
				matLabelp[charSegs.size()] = label;
				charSegs.push_back(charImgCandicate);
				Rect oriRect = boundingRect(charImgCandicate);
				int xoffset, yoffset; xoffset = (oriRect.tl()).x >= 1 ? 1 : 0;  yoffset = (oriRect.tl()).y >= 1 ? 1 : 0;
				Rect adjustRect = oriRect + Point(-xoffset, -yoffset) + Size(2 * xoffset, 2 * yoffset);
				if (adjustRect.x + adjustRect.width > charImgCandicate.cols)
				{
					adjustRect.width = charImgCandicate.cols - adjustRect.x;
				}
				if (adjustRect.y + adjustRect.height > charImgCandicate.rows)
				{
					adjustRect.height = charImgCandicate.rows - adjustRect.y;
				}
				Mat pushM = Mat(charImgCandicate, adjustRect);
				char_imgs.push_back(pushM);
				//				Rect brect = boundingRect(charImgCandicate);
				//				bRects.push_back(brect);
			}
		}
	}
	/*test seg img check*/

	for (auto i = 0; i < charSegs.size(); ++i)
	{
		int label1 = matLabelp[i];
		Vec2i tmpCentroid1 = centroidMat.row(label1);
		Mat tmp1 = charSegs[i].clone();

		for (auto j = i + 1; j < charSegs.size(); ++j)
		{
			int label2 = matLabelp[j];
			Vec2i tmpCentroid2 = centroidMat.row(label2);
			double regiondis = norm(tmpCentroid1, tmpCentroid2);
			cout << regiondis << endl;
			if (regiondis < 15)
			{
				//merge
				Mat tmp2 = charSegs[j].clone();
				Mat tmp3 = tmp1 + tmp2;
				cout << label1 << "  " << label2 << "    " << tmpCentroid1 << "  " << tmpCentroid2 << endl;
				charSegs[i] = tmp3; charSegs[j] = tmp3;
				Mat tmp4 = charSegs[i]; Mat tmp5 = charSegs[j];
				cout << "test" << endl;
			}
			else
			{
				//nothing
			}
		}
	}
	sort(charSegs.begin(), charSegs.end(), [](Mat a, Mat b)
	{
		if (countNonZero(a) < countNonZero(b))
			return true;
		else
			return false;
	});
	charSegs.erase(unique(charSegs.begin(), charSegs.end(), [](Mat a, Mat b)
	{
		if (countNonZero(a != b) == 0)
			return true;
		else
			return false;
	}), charSegs.end());
	for (auto iter = charSegs.begin(); iter != charSegs.end();)
	{
		Rect tmp = boundingRect(*iter);
		double tmpArea = tmp.area();
		//		if (tmpArea < 30||isLSeg(*iter))
		//			iter = char_imgs.erase(iter);
		//		else
		iter++;
	}
	for (auto i = 0; i < charSegs.size(); ++i)
	{
		Rect brect = boundingRect(charSegs[i]);
		bRects.push_back(brect);
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
	for (auto i = 0; i < bRects.size(); ++i)
	{
		rectangle(boundImg, bRects[i], (0, 255, 255), 2);
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
		namedWindow("2.bounded");
		imshow("2.bounded", boundImg);
	}
}

/******************************detection funcs*******************************************/
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

//Vec2i avg_pt(vector<Vec2i> pts)
//{
//	Vec2i new_pt;
//	double xsum = 0.0;
//	double ysum = 0.0;
//	int ptsSize = pts.size();
//	for (size_t i = 0; i < ptsSize; ++i)
//	{
//		xsum += pts[i][0];
//		ysum += pts[i][1];
//	}
//	new_pt = {cvRound(xsum / ptsSize), cvRound(ysum / ptsSize)};
//	return new_pt;
//}
//void computeNewPoints(Vec2i pt1, Vec2i pt2, Vec2i pt3, Vec2i pt4, Vec4i& newpts)
//{
//	/*choose the leftmost and rightmost points in 4 points coz the line is is other */
//	int x[4] = {pt1[0], pt2[0], pt3[0], pt4[0]};
//	int y[4] = {pt1[1], pt2[1], pt3[1], pt4[1]};
//	int x1, y1, x2, y2;
//	int idx1, idx2;
//	idx1 = idx2 = 0;
//	x1 = x[0];
//	x2 = x[0];
//	for (int i = 1; i < 4; ++i)
//	{
//		if (x[i] <= x1)
//		{
//			x1 = x[i];
//			idx1 = i;
//		}
//		if (x[i] >= x2)
//		{
//			x2 = x[i];
//			idx2 = i;
//		}
//	}
//	newpts = {x[idx1],y[idx1], x[idx2], y[idx2]};
//}
//void getRangePts(Vec2i pt1, Vec2i pt2, Vec2i pt3, Vec2i pt4, Vec2i& firstPt, Vec2i& fourthPt)
//{
//	map<int, int> tmpMap;
//	tmpMap[pt1[0]] = pt1[1];
//	tmpMap[pt2[0]] = pt2[1];
//	tmpMap[pt3[0]] = pt3[1];
//	tmpMap[pt4[0]] = pt4[1];
//	auto iterBegin = tmpMap.begin();
//	auto iterEnd = tmpMap.end();
//	--iterEnd;
//	firstPt = {iterBegin->first, iterBegin->second};
//	fourthPt = {iterEnd->first, iterEnd->second};
//}

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
	float step = 2 * CV_PI / (6.0f * radius);
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

void detect_circle3(Mat diagram_segment, Mat& color_img, Mat& diagram_segwithoutcircle, Mat& withoutCirBw, vector<circle_class>& circles, bool showFlag)
{
	//Hough circle detection go
	vector<Vec3f> houghc;
	Mat diagram_segment2 = color_img.clone();

	HoughCircles(diagram_segment, houghc, HOUGH_GRADIENT,1, 20, 1,22, 50,diagram_segment.rows/2);
	for (size_t i = 0; i < houghc.size(); i++)
	{
		Point center(cvRound(houghc[i][0]), cvRound(houghc[i][1]));
		int radius = cvRound(houghc[i][2]);
		// circle center
		circle(diagram_segment2, center, 3, Scalar(0, 255, 0), -1, 8, 0);
		// circle outline
		circle(diagram_segment2, center, radius, Scalar(0, 0, 255), 3, 8, 0);
	}
}

void circleRANSAC(Mat &image, std::vector<Vec3f> &circles, double canny_threshold, double circle_threshold, int numIterations)
{
	CV_Assert(image.type() == CV_8UC1 || image.type() == CV_8UC3);
	circles.clear();

	// Edge Detection
	Mat edges;
	Canny(image, edges, MAX(canny_threshold / 2, 1), canny_threshold, 3);

	// Create point set from Canny Output
	std::vector<Point2d> points;
	for (int r = 0; r < edges.rows; r++)
	{
		for (int c = 0; c < edges.cols; c++)
		{
			if (edges.at<unsigned char>(r, c) == 255)
			{
				points.push_back(cv::Point2d(c, r));
			}
		}
	}

	// 4 point objects to hold the random samples
	Point2d pointA;
	Point2d pointB;
	Point2d pointC;
	Point2d pointD;

	// distances between points
	double AB;
	double BC;
	double CA;
	double DC;

	// varibales for line equations y = mx + b
	double m_AB;
	double b_AB;
	double m_BC;
	double b_BC;

	// varibles for line midpoints
	double XmidPoint_AB;
	double YmidPoint_AB;
	double XmidPoint_BC;
	double YmidPoint_BC;

	// variables for perpendicular bisectors
	double m2_AB;
	double m2_BC;
	double b2_AB;
	double b2_BC;

	// RANSAC
	cv::RNG rng;
	int min_point_separation = 10; // change to be relative to image size?
	int colinear_tolerance = 1; // make sure points are not on a line
	int radius_tolerance = 3; // change to be relative to image size?
	int points_threshold = 10; //should always be greater than 4
	//double min_circle_separation = 10; //reject a circle if it is too close to a previously found circle
	//double min_radius = 10.0; //minimum radius for a circle to not be rejected

	int x, y;
	Point2d center;
	double radius;

	// Iterate
	for (int iteration = 0; iteration < numIterations; iteration++)
	{
		//std::cout << "RANSAC iteration: " << iteration << std::endl;

		// get 4 random points
		pointA = points[rng.uniform((int)0, (int)points.size())];
		pointB = points[rng.uniform((int)0, (int)points.size())];
		pointC = points[rng.uniform((int)0, (int)points.size())];
		pointD = points[rng.uniform((int)0, (int)points.size())];

		// calc lines
		AB = norm(pointA - pointB);
		BC = norm(pointB - pointC);
		CA = norm(pointC - pointA);
		DC = norm(pointD - pointC);

		// one or more random points are too close together
		if (AB < min_point_separation || BC < min_point_separation || CA < min_point_separation || DC < min_point_separation) continue;

		//find line equations for AB and BC
		//AB
		m_AB = (pointB.y - pointA.y) / (pointB.x - pointA.x + 0.000000001); //avoid divide by 0
		b_AB = pointB.y - m_AB*pointB.x;

		//BC
		m_BC = (pointC.y - pointB.y) / (pointC.x - pointB.x + 0.000000001); //avoid divide by 0
		b_BC = pointC.y - m_BC*pointC.x;


		//test colinearity (ie the points are not all on the same line)
		if (abs(pointC.y - (m_AB*pointC.x + b_AB + colinear_tolerance)) < colinear_tolerance) continue;

		//find perpendicular bisector
		//AB
		//midpoint
		XmidPoint_AB = (pointB.x + pointA.x) / 2.0;
		YmidPoint_AB = m_AB * XmidPoint_AB + b_AB;
		//perpendicular slope
		m2_AB = -1.0 / m_AB;
		//find b2
		b2_AB = YmidPoint_AB - m2_AB*XmidPoint_AB;

		//BC
		//midpoint
		XmidPoint_BC = (pointC.x + pointB.x) / 2.0;
		YmidPoint_BC = m_BC * XmidPoint_BC + b_BC;
		//perpendicular slope
		m2_BC = -1.0 / m_BC;
		//find b2
		b2_BC = YmidPoint_BC - m2_BC*XmidPoint_BC;

		//find intersection = circle center
		x = (b2_AB - b2_BC) / (m2_BC - m2_AB);
		y = m2_AB * x + b2_AB;
		center = Point2d(x, y);
		radius = cv::norm(center - pointB);

		/// geometry debug image
		if (true)
		{
			Mat debug_image = edges.clone();
			cvtColor(debug_image, debug_image, CV_GRAY2RGB);

			Scalar pink(255, 0, 255);
			Scalar blue(255, 0, 0);
			Scalar green(0, 255, 0);
			Scalar yellow(0, 255, 255);
			Scalar red(0, 0, 255);

			// the 3 points from which the circle is calculated in pink
			circle(debug_image, pointA, 3, pink);
			circle(debug_image, pointB, 3, pink);
			circle(debug_image, pointC, 3, pink);

			// the 2 lines (blue) and the perpendicular bisectors (green)
			line(debug_image, pointA, pointB, blue);
			line(debug_image, pointB, pointC, blue);
			line(debug_image, Point(XmidPoint_AB, YmidPoint_AB), center, green);
			line(debug_image, Point(XmidPoint_BC, YmidPoint_BC), center, green);

			circle(debug_image, center, 3, yellow); // center
			circle(debug_image, center, radius, yellow);// circle

			// 4th point check
			circle(debug_image, pointD, 3, red);

			imshow("ransac debug", debug_image);
			waitKey(0);
		}

		//check if the 4 point is on the circle
		if (abs(cv::norm(pointD - center) - radius) > radius_tolerance) continue;

		// vote
		std::vector<int> votes;
		std::vector<int> no_votes;
		for (int i = 0; i < (int)points.size(); i++)
		{
			double vote_radius = norm(points[i] - center);

			if (abs(vote_radius - radius) < radius_tolerance)
			{
				votes.push_back(i);
			}
			else
			{
				no_votes.push_back(i);
			}
		}

		// check votes vs circle_threshold
		if ((float)votes.size() / (2.0*CV_PI*radius) >= circle_threshold)
		{
			circles.push_back(Vec3f(x, y, radius));

			// voting debug image
			if (false)
			{
				Mat debug_image2 = edges.clone();
				cvtColor(debug_image2, debug_image2, CV_GRAY2RGB);

				Scalar yellow(0, 255, 255);
				Scalar green(0, 255, 0);

				circle(debug_image2, center, 3, yellow); // center
				circle(debug_image2, center, radius, yellow);// circle

				// draw points that voted
				for (int i = 0; i < (int)votes.size(); i++)
				{
					circle(debug_image2, points[votes[i]], 1, green);
				}

				imshow("ransac debug", debug_image2);
				waitKey(0);
			}

			// remove points from the set so they can't vote on multiple circles
			std::vector<Point2d> new_points;
			for (int i = 0; i < (int)no_votes.size(); i++)
			{
				new_points.push_back(points[no_votes[i]]);
			}
			points.clear();
			points = new_points;
		}

		// stop RANSAC if there are few points left
		if ((int)points.size() < points_threshold)
			break;
	}

	return;
}

void ransac_circle(Mat diagram_segment, Mat& color_img, Mat& diagram_segwithoutcircle, Mat& withoutCirBw, vector<circle_class>& circles, bool showFlag)
{
	vector<Vec3f> ransac_cs;

	const clock_t start = clock();
	circleRANSAC(diagram_segment, ransac_cs, 100, 10, 2000);
	clock_t end = clock();

	cout << "Found " << (int)ransac_cs.size() << " Circles." << endl;

	double time = ((double)(end - start)) / (double)CLOCKS_PER_SEC;
	std::cout << "RANSAC runtime: " << time << " seconds" << std::endl;

	// Draw Circles
	cvtColor(diagram_segment, diagram_segment, CV_GRAY2RGB);
	for (int i = 0; i < (int)ransac_cs.size(); i++)
	{
		int x = ransac_cs[i][0];
		int y = ransac_cs[i][1];
		float rad = ransac_cs[i][2];

		circle(diagram_segment, Point(x, y), rad, Scalar(0, 255, 0));
	}

	imshow("circles", diagram_segment);
}

/**************line part****************/


//bool crossPtWithinLines(Vec4i line1, Vec4i line2, Vec2i cross)
//{
//	if (in_line(line1, cross))
//	{
//		//cross in line1
//		if (in_line(line2, cross))
//		{
//			//cross in line2
//			return true;
//		}
//		else
//		{
//			cout << "cross not in line2" << endl;
//			return false;
//		}
//	}
//	else
//	{
//		cout << "cross not in line1" << endl;
//		return false;
//	}
//}
//Vec2i ptAttachToCircle(Vec2f& tmp, vector<circle_class>& circles)
//{
//	Vec2i cross;
//	if (circles.size() != 0)
//	{
//		for (int j = 0; j < circles.size(); j++)
//		{
//			Vec3i c = circles[j].getCircleVec();
//			Vec2f center = circles[j].getCenter();
//			float radius = circles[j].getRadius();
//			/*			if (on_circle(pt1, c) || on_circle(pt2, c) || on_circle(pt3, c) || on_circle(pt4, c))
//			{
//			continue;
//			}
//			else */
//			//To Note: there're changes here
//			if (on_circle(tmp, c))
//			{
//				cout << "cross approximately on circle" << endl;
//				float minDiff = 100;
//				int cenx = int(tmp[0]);
//				int ceny = int(tmp[1]);
//				int offset = int(p2pdistance(tmp, center) - radius);
//				offset = (offset < 1) ? 1 : offset;
//				for (auto m = cenx - offset; m <= cenx + offset; m++)
//				{
//					for (auto n = ceny - offset; n <= ceny + offset; n++)
//					{
//						Vec2i temp2 = {m, n};
//						float tmpDiff = abs(p2pdistance(temp2, center) - radius);
//						if (tmpDiff < minDiff)
//						{
//							tmp = temp2;
//							minDiff = tmpDiff;
//						}
//					}
//				}
//			}
//		}
//		cross = {int(tmp[0]), int(tmp[1])};
//		cout << "after attach to circle the cross is " << cross << endl;
//	}
//	else
//	{
//		cross = {int(tmp[0]), int(tmp[1])};
//		cout << "not circle on image, no attachment" << endl;
//	}
//	return cross;
//}

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

//bool intersection(Vec2i o1, Vec2i p1, Vec2i o2, Vec2i p2,
//                  Vec2i& r)
//{
//	Vec2i x = o2 - o1;
//	Vec2i d1 = p1 - o1;
//	Vec2i d2 = p2 - o2;
//
//	float cross = d1[0] * d2[1] - d1[1] * d2[0];
//	if (abs(cross) < /*EPS*/1e-8)
//		return false;
//
//	double t1 = (x[0] * d2[1] - x[1] * d2[0]) / cross;
//	r = o1 + d1 * t1;
//	return true;
//}

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
	cout << "radius - pt2center: " << radius - pt2center << endl;
	if (radius - pt2center > 5)
		return 0;//inside circle
	else if (abs(radius - pt2center) <= 5)
		return 1; // on cirlce
	else if (pt2center - radius < 8)
		return 2;//outside circle with a close distance
	else
		return 3;//outside circle with large distance
}

bool dashLineRecovery(vector<Point2i>& edgePt, Vec2i col_p1, Vec2i col_p2, vector<circle_class>& circles, bool plflag = false, bool pcflag = false, bool ppflag = false)
{//check if there's dash line between
	Vec4i line = pt2line(col_p1, col_p2);
	int flag[1000] = {0};
	bool vertical_flag = (abs(col_p1[0] - col_p2[0]) < abs(col_p1[1] - col_p2[1])) ? true : false;
	int ranges = vertical_flag ? abs(col_p2[1] - col_p1[1]) : abs(col_p2[0] - col_p1[0]);
	cout << "line check between " << col_p1 << " and " << col_p2 << endl;
	int xleft, xright, ysmall, ylarge;
	if (col_p1[0]<col_p2[0])
	{
		xleft = col_p1[0];
		xright = col_p2[0];
	}
	else
	{
		xleft = col_p2[0];
		xright = col_p1[0];
	}
	if (col_p1[1]<col_p2[1])
	{
		ysmall = col_p1[1];
		ylarge = col_p2[1];
	}
	else
	{
		ysmall = col_p2[1];
		ylarge = col_p1[1];
	}
	int offset = 3;
	xleft = xleft - offset; xright = xright + offset;
	ysmall = ysmall - offset; ylarge = ylarge + offset;
	if (circles.size() == 0)
	{
		cout << "non-circle" << endl;
		for (auto i = 0; i < edgePt.size(); ++i)
		{
			Vec2i pt = edgePt[i];
			if ((pt[0]> xleft&& pt[0]<xright&&pt[1]>ysmall&&pt[1]<ylarge)&&on_line(line, pt))
			{
				if (vertical_flag)
					flag[pt[1]] = 1;
				else
					flag[pt[0]] = 1;
			}
		}
		double ratio;
		int nums = count(flag, flag + 1000, 1);
		ratio = 1.0* nums / ranges;
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
//				if (p2pdistance(col_p1, col_p2) < 20)
//				{
					return true;
//				}
//				else
//				{
//					return false;
//				}
			}
			else
			{
				for (auto j = 0; j < edgePt.size(); ++j)
				{
					Vec2i pt = edgePt[j];
					   if ((pt[0]> xleft&& pt[0]<xright&&pt[1]>ysmall&&pt[1]<ylarge)&&on_line(line, pt))
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
				double threshold_ratio = 0.8;
				if (ratio < threshold_ratio)
					return false;
				else
					return true;;
			}
		}
	}
}

bool basicDashLineRev(vector<Point2i>& withoutOnL_ept,bool vertical_flag, Vec4i line,int flag[], int ranges)
{
	Vec2i pt1, pt2; line2pt(line, pt1, pt2);
	int xleft, xright, ysmall, ylarge;
	if (pt1[0]<pt2[0])
	{
		xleft = pt1[0];
		xright = pt2[0];
	}
	else
	{
		xleft = pt2[0];
		xright = pt1[0];
	}
	if (pt1[1]<pt2[1])
	{
		ysmall = pt1[1];
		ylarge = pt2[1];
	}
	else
	{
		ysmall = pt2[1];
		ylarge = pt1[1];
	}
	int offset = 3;
	xleft = xleft - offset; xright = xright + offset;
	ysmall = ysmall - offset; ylarge = ylarge + offset;
	for (auto j = 0; j < withoutOnL_ept.size(); ++j)
	{
		Vec2i pt = withoutOnL_ept[j];
		   if ((pt[0]> xleft&& pt[0]<xright&&pt[1]>ysmall&&pt[1]<ylarge)&&on_line(line, pt))
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
	double threshold_ratio = 0.6;
	if (ratio < threshold_ratio)
		return false;
	else
		return true;
}

bool basicDashLineRev(vector<Point2i>& ept, Vec2i pt1, Vec2i pt2,double changeThres=0)
{
	int flag[1000] = { 0 };
	bool vertical_flag = (abs(pt1[0] - pt2[0]) < abs(pt1[1] - pt2[1])) ? true : false;
	int ranges = vertical_flag ? abs(pt1[1] - pt2[1]) : abs(pt1[0] - pt2[0]);
	Vec4i line =  pt2line(pt1, pt2);
	int xleft, xright, ysmall, ylarge;
	if (pt1[0]<pt2[0])
	{
		xleft = pt1[0];
		xright = pt2[0];
	}
	else
	{
		xleft = pt2[0];
		xright = pt1[0];
	}
	if (pt1[1]<pt2[1])
	{
		ysmall = pt1[1];
		ylarge = pt2[1];
	}
	else
	{
		ysmall = pt2[1];
		ylarge = pt1[1];
	}
	int offset = 3;
	xleft = xleft - offset; xright = xright + offset;
	ysmall = ysmall - offset; ylarge = ylarge + offset;
	for (auto j = 0; j < ept.size(); ++j)
	{
		Vec2i pt = ept[j];
		   if ((pt[0]> xleft&& pt[0]<xright&&pt[1]>ysmall&&pt[1]<ylarge)&&on_line(line, pt))
		{
			if (in_line(line, pt))
			{
				if (vertical_flag)
					flag[pt[1]] = 1;
				else
					flag[pt[0]] = 1;
			}
		}
	}
	double ratio;
	int nums = count(flag, flag + 1000, 1);
	ratio = 1.0 * nums / ranges;
	cout << ratio * 100 << "%" << endl;
	double threshold_ratio = (changeThres == 0) ? 0.6 : changeThres;
	if (ratio < threshold_ratio)
		return false;
	else
		return true;
}

bool dashLineRecovery2(vector<Point2i>& withoutOnL_ept, Vec2i p_closer, Vec2i p_farther, Vec2i p_cross, vector<circle_class>& circles, bool plflag = false, bool pcflag = false, bool ppflag = false)
{
	//check if there's dash line between
	Vec4i line = pt2line(p_closer, p_cross);
	int flag[1000] = {0};
	bool vertical_flag = (abs(p_closer[0] - p_cross[0]) < abs(p_closer[1] - p_cross[1])) ? true : false;
	int ranges = vertical_flag ? abs(p_cross[1] - p_closer[1]) : abs(p_cross[0] - p_closer[0]);
	cout << "line check between " << p_closer << " and " << p_cross << endl;
	int xleft, xright, ysmall, ylarge;
	if (p_closer[0]<p_cross[0])
	{
		xleft = p_closer[0];
		xright = p_cross[0];
	}
	else
	{
		xleft = p_cross[0];
		xright = p_closer[0];
	}
	if (p_closer[1]<p_cross[1])
	{
		ysmall = p_closer[1];
		ylarge= p_cross[1];
	}
	else
	{
		ysmall = p_cross[1];
		ylarge = p_closer[1];
	}
	int offset = 3;
	xleft = xleft - offset; xright = xright + offset;
	ysmall = ysmall - offset; ylarge = ylarge + offset;
	if (circles.size() == 0)
	{
		cout << "image without circle" << endl;
		for (auto i = 0; i < withoutOnL_ept.size(); ++i)
		{
			Vec2i pt = withoutOnL_ept[i];
			   if ((pt[0]> xleft&& pt[0]<xright&&pt[1]>ysmall&&pt[1]<ylarge)&&on_line(line, pt))
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
		bool tmpFlag[2] = {false, false};
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
//			if (abs(center2lineDis - radius) <= 5)
//			{
//				cout << "the line is almost tangent to the circle" << endl;
				if (closer_flag == 1 && cross_flag == 1)
				{
					// both are on the circle
					cout << "both points are on the circle" << endl;
					tmpFlag[i] = true;
				}
				else if (closer_flag == 1 && cross_flag == 2)
				{
					// point1 is on the circle
					cout << "closer point is on the circle and cross point is outside the circle with close distance" << endl;
					cout << "to go" << endl;
					tmpFlag[i] = true;
				}
				else if (closer_flag ==1 && cross_flag == 3)
				{
					cout << "closer point is on the circle and cross point is outside the circle with large distance" << endl;
					cout << "norm handle" << endl;
					tmpFlag[i] = basicDashLineRev(withoutOnL_ept, vertical_flag, line, flag, ranges);
				}
				else if (cross_flag == 1 && closer_flag == 2)
				{
					cout << "cross point is on the circle and closer point is outside the circle with a close distance" << endl;
					tmpFlag[i] = true;

				}
				else if (cross_flag == 1 && closer_flag == 0)
				{
					cout << "cross point is on the circle and closer point is inside the circle" <<"norm handle"<< endl;
					tmpFlag[i] = basicDashLineRev(withoutOnL_ept, vertical_flag, line, flag, ranges);
				}
				else if (cross_flag == 1 && closer_flag == 3)
				{
					
					cout << "cross point is on the circle and closer point is outside the circle" <<"norm handle"<< endl;
					tmpFlag[i] = basicDashLineRev(withoutOnL_ept, vertical_flag, line, flag, ranges);
				}
				else
				{
					cout << "neither points are on the circle" << endl;
//					for (auto j = 0; j < withoutOnL_ept.size(); ++j)
//					{
//						Vec2i pt = withoutOnL_ept[j];
//						   if ((pt[0]> xleft&& pt[0]<xright&&pt[1]>ysmall&&pt[1]<ylarge)&&on_line(line, pt))
//						{
//							if (vertical_flag)
//								flag[pt[1]] = 1;
//							else
//								flag[pt[0]] = 1;
//						}
//					}
//					double ratio;
//					int nums = count(flag, flag + 1000, 1);
//					ratio = 1.0 * nums / ranges;
//					cout << ratio * 100 << "%" << endl;
//					double threshold_ratio = 0.7;
//					if (ratio < threshold_ratio)
//						tmpFlag[i]= false;
//					else
//						tmpFlag[i] = true;
					tmpFlag[i] = basicDashLineRev(withoutOnL_ept, vertical_flag, line, flag, ranges);
				}
//			}
//			else
//			{
//				cout << "not tangent to the cirlce" << endl;
//				for (auto j = 0; j < withoutOnL_ept.size(); ++j)
//				{
//					Vec2i pt = withoutOnL_ept[j];
//					if (in_line(line, pt))
//					{
//						if (vertical_flag)
//							flag[pt[1]] = 1;
//						else
//							flag[pt[0]] = 1;
//					}
//				}
//				double ratio;
//				int nums = count(flag, flag + 1000, 1);
//				ratio = 1.0 * nums / ranges;
//				cout << ratio * 100 << "%" << endl;
//				double threshold_ratio = 0.7;
//				if (ratio < threshold_ratio)
//					return false;
//				else
//					return true;
//				//return false;
//			}
		}
		if (tmpFlag[0] == true || tmpFlag[1] == true)
			return true;
		else
			return false;
	}
}

bool dashLineRecovery(vector<Point2i> &withoutOnL_ept, vector<Point2i>& oriEdgePoints, Vec2i p_closer, Vec2i p_farther, Vec2i p_cross, vector<circle_class> &circles, bool plflag = false, bool pcflag = false, bool ppflag = false, bool secCheck = true)
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
			auto iter = find_if(oriEdgePoints.begin(), oriEdgePoints.end(), [&](Point2i a)
			{
				if (same_pt(a, tmp_pt, 3))
					return true;
				else
					return false;
			});
			if (iter == oriEdgePoints.end())
			{
				cout << tmp_pt;
				cout << "       no recovery" << endl;
				return false;
			}
		}
		if (secCheck)
		{
			cout << "through the first recovery check" << endl;
			if (dashLineRecovery2(withoutOnL_ept, p_closer, p_farther, p_cross, circles))
			{
				cout << "recovery" << endl;
				return true;
			}
			else
				return false;
		}
		else
			return true;
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
			auto iter = find_if(oriEdgePoints.begin(), oriEdgePoints.end(), [&](Point2i a)
			{
				if (same_pt(a, tmp_pt, 3))
					return true;
				else
					return false;
			});
			if (iter == oriEdgePoints.end())
			{
				cout << "no recovery" << endl;
				flag = false;
				return false;
			}
		}
		if (secCheck)
		{
			cout << "through the first recovery check" << endl;
			if (dashLineRecovery2(withoutOnL_ept, p_closer, p_farther, p_cross, circles))
				return true;
			else
				return false;
		}
		else
			return true;
	}
}

bool ptXiSortPred(Vec2i pt1, Vec2i pt2)
{
	return (pt1[0] < pt2[0]);
}
bool ptYiSortPred(Vec2i pt1, Vec2i pt2)
{
	return (pt1[1] < pt2[1]);
}
bool ptXfSortPred(Vec2f pt1, Vec2f pt2)
{
	return (pt1[0] < pt2[0]);
}
bool ptYfSortPred(Vec2f pt1, Vec2f pt2)
{
	return (pt1[1] < pt2[1]);
}

//bool with_same_line(vector<Vec4i>& plainLines, Vec2i pt1, Vec2i pt2)
//{
//	for (auto i = 0; i < plainLines.size(); i++)
//	{
//		Vec4i line = plainLines[i];
//		Vec2i tmpPt1, tmpPt2;
//		tmpPt1 = {line[0], line[1]};
//		tmpPt2 = {line[2], line[3]};
//		if ((same_pt(pt1, tmpPt1) && same_pt(pt2, tmpPt2)) || (same_pt(pt1, tmpPt2) && (same_pt(pt2, tmpPt1))))
//			return true;
//		else if ((same_pt(pt1, tmpPt1) && on_line(line, pt2)) || (same_pt(pt2, tmpPt2) && on_line(line, pt1))
//			|| (same_pt(pt1, tmpPt2) && on_line(line, pt2)) || (same_pt(pt2, tmpPt1) && on_line(line, pt1)))
//			return true;
//	}
//	return false;
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

void vecSwapValue(vector<Vec4i>::iterator iter1, vector<Vec4i>::iterator iter2)
{
	auto tmp = *iter1;
	*iter1 = *iter2;
	*iter2 = tmp;
}
void vecSwapValue(vector<Vec4f>::iterator iter1, vector<Vec4f>::iterator iter2)
{
	auto tmp = *iter1;
	*iter1 = *iter2;
	*iter2 = tmp;
}

int line_recovery_process(line_class* linex, Vec2i p_cross, vector<Point2i>& withoutOnL_ept, vector<Point2i> &oriEdgePoints, vector<point_class>& pointxs, vector<circle_class>& circlexs, bool id_change = false)
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
		else if (dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt1, pt2, p_cross, circlexs))
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
		else if (dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt2, pt1, p_cross, circlexs))
		{
			cout << linex->getpt2vec(pointxs) << "  ->  " << p_cross << endl;
			linex->setPt2_vec(pointxs, p_cross);
			pos = 2;
		}
	}
	return pos;
}

void cross_refinement(Vec2f& raw_cross, line_class* lx1, line_class* lx2, vector<circle_class>& circlexs, vector<point_class>& pointxs, vector<point_class>& TACpointxs ,vector<Point2i>& withoutOnL_ept, vector<Point2i> &oriEdgePoints)
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
			if (same_pt(pt1, pt3))
				raw_cross = (pt1+pt3)/2;
			cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			double cpd1, cpd2, ld1, ld2;
			ld1 = p2pdistance(pt1, pt2);
			ld2 = p2pdistance(pt3, pt4);
			cpd1 = p2pdistance(raw_cross, pt2);
			cpd2 = p2pdistance(raw_cross, pt4);
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
			if (same_pt(pt1, pt4))
				raw_cross = (pt1 + pt4) / 2;
			//			cout << lx1->getpt1vec(pointxs) << "  ->   " << raw_cross << endl;
			cout << lx2->getpt2vec(pointxs) << "  ->   " << raw_cross << endl;
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
				cout << "cross disjoint with line1" << endl;
				int pos = line_recovery_process(lx2, lx1->getpt1vec(pointxs), withoutOnL_ept,oriEdgePoints, pointxs, circlexs);
//				if (close_pt(raw_cross, pt3))
//				{
//					bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt3, pt4, raw_cross, circlexs, false, false, false, false);
//					if (FalseFlag)
//					{
//						cout << "id " << lx2->getPt1Id() << " id change to id " << lx1->getPt1Id() << endl;
//						lx2->setpt1Id(lx1->getPt1Id());
//					}
//				}
//				else if (close_pt(raw_cross,pt4))
//				{
//					bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt4, pt3, raw_cross, circlexs, false, false, false, false);
//					if (FalseFlag)
//					{
//						cout << "id " << lx2->getPt2Id() << " id change to id " << lx1->getPt1Id() << endl;
//						lx2->setpt2Id(lx1->getPt1Id());
//					}
//				}
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
						rm_point_by_id(pointxs, tmp_id2);
					}
					else
					{
						cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
						lx1->setpt1Id(tmp_id2);
						rm_point_by_id(pointxs, tmp_id1);
					}
					cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
					cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
				}
				else if (pos == 2)
				{
					cout << "recovery cross link to pt4" << endl;
					int tmp_id1, tmp_id2;
					tmp_id1 = lx1->getPt1Id();
					tmp_id2 = lx2->getPt2Id();
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
			if (same_pt(pt2, pt3))
				raw_cross = (pt2 + pt3) / 2;
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
			if (same_pt(pt2, pt4))
				raw_cross = (pt2 + pt4) / 2;
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
				cout <<  "cross in line2 and end of line1" << endl;
			}
			else
			{
				int pos = line_recovery_process(lx2, lx1->getpt2vec(pointxs), withoutOnL_ept,oriEdgePoints, pointxs, circlexs);
//				if (close_pt(raw_cross, pt3))
//				{
//					bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt3, pt4, raw_cross, circlexs, false, false, false, false);
//					if (FalseFlag)
//					{
//						cout << "id " << lx2->getPt1Id() << " id change to id " << lx1->getPt2Id() << endl;
//						lx2->setpt1Id(lx1->getPt2Id());
//					}
//				}
//				else if (close_pt(raw_cross, pt4))
//				{
//					bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt4, pt3, raw_cross, circlexs, false, false, false, false);
//					if (FalseFlag)
//					{
//						cout << "id " << lx2->getPt2Id() << " id change to id " << lx1->getPt2Id() << endl;
//						lx2->setpt2Id(lx1->getPt2Id());
//					}
//				}
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
						rm_point_by_id(pointxs, tmp_id2);
					}
					else
					{
						cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
						lx1->setpt2Id(tmp_id2);
						rm_point_by_id(pointxs, tmp_id1);
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
		raw_cross = { (raw_cross[0] + pt3[0]) / 2, (raw_cross[1] + pt3[1]) / 2 };
		lx2->setPt1_vec(pointxs, raw_cross);
		if (in_line1)
		{
			cout <<  "cross in line1 and end of line2" << endl;
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else
		{
			int pos = line_recovery_process(lx1, lx2->getpt1vec(pointxs), withoutOnL_ept,oriEdgePoints, pointxs, circlexs);
//			if (close_pt(raw_cross, pt1))
//			{
//				bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt1, pt2, raw_cross, circlexs, false, false, false, false);
//				if (FalseFlag)
//				{
//					cout << "id " << lx1->getPt1Id() << " id change to id " << lx2->getPt1Id() << endl;
//					lx2->setpt1Id(lx1->getPt1Id());
//				}
//			}
//			else if (close_pt(raw_cross, pt2))
//			{
//				bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt2, pt1, raw_cross, circlexs, false, false, false, false);
//				if (FalseFlag)
//				{
//					cout << "id " << lx1->getPt2Id() << " id change to id " << lx2->getPt1Id() << endl;
//					lx2->setpt2Id(lx1->getPt1Id());
//				}
//			}
			
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
					rm_point_by_id(pointxs, tmp_id2);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt1Id(tmp_id2);
					rm_point_by_id(pointxs, tmp_id1);
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
					lx2->setpt1Id(tmp_id2);
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
		raw_cross = { (raw_cross[0] + pt4[0]) / 2, (raw_cross[1] + pt4[1]) / 2 };
		lx2->setPt2_vec(pointxs, raw_cross);
		if (in_line1)
		{
			cout << "cross in line1 and end of line2" << endl;
			cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
			cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
		}
		else
		{
			int pos = line_recovery_process(lx1, lx2->getpt2vec(pointxs), withoutOnL_ept, oriEdgePoints, pointxs, circlexs);
//			if (close_pt(raw_cross, pt1))
//			{
//				bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt1, pt2, raw_cross, circlexs, false, false, false, false);
//				if (FalseFlag)
//				{
//					cout << "False Flag" << endl;
//					cout << "id " << lx1->getPt1Id() << " id change to id " << lx2->getPt2Id() << endl;
//					lx2->setpt1Id(lx1->getPt1Id());
//				}
//			}
//			else if (close_pt(raw_cross, pt2))
//			{
//				bool FalseFlag = !dashLineRecovery(withoutOnL_ept, oriEdgePoints, pt2, pt1, raw_cross, circlexs, false, false, false, false);
//				if (FalseFlag)
//				{
//					cout << "False Flag" << endl;
//					cout << "id " << lx1->getPt2Id() << " id change to id " << lx2->getPt2Id() << endl;
//					lx2->setpt2Id(lx1->getPt1Id());
//				}
//			}
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
					rm_point_by_id(pointxs,tmp_id2);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt2Id(tmp_id2);
					rm_point_by_id(pointxs, tmp_id1);
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
					rm_point_by_id(pointxs,tmp_id2);
				}
				else
				{
					cout << "id " << tmp_id1 << " id change to id " << tmp_id2 << endl;
					lx2->setpt2Id(tmp_id2);
					rm_point_by_id(pointxs,tmp_id1);
				}
				cout << lx1->getPt1Id() << " " << lx1->getpt1vec(pointxs) << ", " << lx1->getPt2Id() << " " << lx1->getpt2vec(pointxs) << endl;
				cout << lx2->getPt1Id() << " " << lx2->getpt1vec(pointxs) << ", " << lx2->getPt2Id() << " " << lx2->getpt2vec(pointxs) << endl;
			}
		}
	}
	else
	{
		cout << "two cross disjoint line" << endl;
		if (in_line1 && !in_line2)
		{
			cout << "cross in line1 but not in line2" << endl;
			line_recovery_process(lx2, raw_cross, withoutOnL_ept,oriEdgePoints, pointxs, circlexs);
		}
		else if (in_line2 && !in_line1)
		{
			cout << "cross in line2 but not in line1" << endl;
			line_recovery_process(lx1, raw_cross, withoutOnL_ept,oriEdgePoints, pointxs, circlexs);
		}
		else if (in_line1 && in_line2)
		{
			cout << "inner cross" << endl;
			point_class a; a.setPid(TACpointxs.size()); a.setXY(raw_cross); a.pushLid(lx1->getLid());
			TACpointxs.push_back(a);
		}
		else
		{
			cout << "outer cross" << endl;
			int pos1 = line_recovery_process(lx1, raw_cross, withoutOnL_ept,oriEdgePoints, pointxs, circlexs);
			int pos2 = line_recovery_process(lx2, raw_cross, withoutOnL_ept, oriEdgePoints,pointxs, circlexs);
			cout << "*********************pos1 " << pos1 << ", pos2 " << pos2 << endl;
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
				tmp_id2 = lx2->getPt1Id();////////////
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

void detect_line3(Mat diagram_segment, Mat diagram_segwithoutcircle, Mat& withoutCirBw, vector<point_class> pointXs, vector<circle_class>& circles, Mat& color_img, vector<line_class> lineXs, vector<Point2i>& oriEdgePoints, Mat& drawedImages, bool showFlag = true, string fileName = "")
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
	Mat oriL_img = color_img.clone();
		for (auto i = 0; i < plainLines.size(); i++)
		{
			Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
//			cout << "*********************" << pt1 << " " << pt2 << endl;
			line(oriL_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
			Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
//			circle(oriL_img, Point{ pt1[0], pt1[1] }, 10, tmp);
//			circle(oriL_img, Point{ pt2[0], pt2[1] }, 10, tmp);
		}
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
					sort(tmpPtVec.begin(), tmpPtVec.end(), ptXiSortPred);
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
	Mat withoutCnLBw = withoutCirBw.clone();
	Mat withoutLBw = diagram_segment.clone();
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i pl = plainLines[i];
		Vec2i pt1 = {pl[0], pl[1]};
		Vec2i pt2 = {pl[2], pl[3]};
		cv::line(withoutCnLBw, pt1, pt2, 0, 3, 8, 0);
		cv::line(withoutLBw, pt1, pt2, 0, 3, 8, 0);
	}
	vector<Point2i> withouCnlBw_ept = getPointPositions(withoutCnLBw);
	vector<Point2i> withoutL_ept = getPointPositions(withoutLBw);

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
		Vec2i tmpPt1, tmpPt2; line2pt(plainLines[i], tmpPt1, tmpPt2);
		if (norm(tmpPt1, tmpPt2) > 20)
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
			delete ptx1, ptx2, lx;
		}
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
	vector<point_class> toAddCrossPoints;
	for (auto i = 0; i < linexs.size(); i++)
	{
		line_class* linex1 = &linexs[i];
		cout << endl << "**********"<<linex1->getLineVec(pointxs) << endl;
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
				double angle = llAngle(linex1_vec, linex2_vec);
				if ((angle>15)&&(!isInImage(diagram_segwithoutcircle.cols, diagram_segwithoutcircle.rows, raw_cross)))
				{
					cout << raw_cross << endl;
					cout << "cross out of scope, then take it as no cross" << endl;
					continue;
				}
				else
				{
					//cross in scope
					cout << linex1_vec << endl << linex2_vec << endl;
					cross_refinement(raw_cross, linex1, linex2, circles, pointxs,toAddCrossPoints, withouCnlBw_ept, oriEdgePoints);
					//					cout << linex1->getLineVec(pointxs)<< linex1->getPt1Id() << "," << linex1->getPt2Id() << endl;
					//					cout << linex2->getLineVec(pointxs) << linex2->getPt1Id() << "," << linex2->getPt2Id() << endl;
					cout << "step separate" << endl;
				}
			}
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
	cout << "after cross points refinement"<<endl;
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
				if (len < cir[2]/2)
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
						auto rm_iter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt1Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;

					}
					else if (not_find_pt2_flag && !not_find_pt1_flag)
					{
						auto rm_iter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt2Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;

					}
					else if (not_find_pt1_flag && not_find_pt2_flag)
					{
						auto rm_iter1 = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt1Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter1);
						auto rm_iter2 = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt2Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter2);
						iter = linexs.erase(iter);
						rm_flag = true; rm_line_num++;
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
	cout << "after remove isolated short lines due to the non-complete circle removal " << endl;
	for (auto i = 0; i < linexs.size(); i++)
	{
		cout << linexs[i].getLineVec(pointxs) << linexs[i].getPt1Id() << ", " << linexs[i].getPt2Id() << endl;
	}
	for (auto j = 0; j < pointxs.size(); j++)
	{
		cout << pointxs[j].getPid() << "  " << pointxs[j].getXY() << endl;
	}

	/*set similar points identical*/


//	cout << "now remove indentical vec points and reorder the indices";
	// point index indicate the line id in which it locate
	map<int, int> changeMap;
//	for (auto i = 0; i < pointxs.size(); ++i)
//	{
//		//the location of the first pointx in the loop
//		Vec2i px1 = pointxs[i].getXY();
//		int erase_offset = 0; // the var log the erased points num and offset the indexing digit when indexing
//		for (auto j = i + 1; j < pointxs.size(); ++j)
//		{
//			Vec2i px2 = pointxs[j].getXY();
//			cout << pointxs[i].getXY() << "   " << pointxs[j].getXY() << endl;
//			if (same_pt(px1, px2))
//			{
//				//the two points vector is the same, then erase the second point from the pointxs set
//				// and ajust the respective line with previous erased point
//				//				cout << "find the identical vec points " << px1 << " <- "<<pointxs[i].getPid() <<", "<< pointxs[j].getPid()<< endl;
//				changeMap[pointxs[j].getPid()] = pointxs[i].getPid();
//				int line_id = pointxs[j].getPid() / 2;
//				int pos = pointxs[j].getPid() % 2;// if 0, means the first point in the line, otherwise 1, means the second point in the line
//				pointxs.erase(pointxs.begin() + j);
//				if (pos)
//				{
//					linexs[changeMap0[line_id]].setpt2Id(pointxs[i].getPid());
//				}
//				else
//				{
//					linexs[changeMap0[line_id]].setpt1Id(pointxs[i].getPid());
//				}
//				j--;
//			}
//		}
//	}
//	for (auto i = 0; i < 2 * linexs.size(); ++i)
//	{
//		if (changeMap[i] == 0)
//		{
//			changeMap[i] = i;
//		}
//	}

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
		cout << linexs[j].getPt1Id() << "->" << changeMap2[linexs[j].getPt1Id()] << endl;
		cout << linexs[j].getPt2Id() << "->" << changeMap2[linexs[j].getPt2Id()] << endl;
		linexs[j].setpt1Id(changeMap2[linexs[j].getPt1Id()]);
		linexs[j].setpt2Id(changeMap2[linexs[j].getPt2Id()]);
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
	
	cout << endl;
	for (auto j = 0; j < toAddCrossPoints.size(); j++)
	{
		auto tmpIter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
		{
			if (same_pt(a,toAddCrossPoints[j]))
				return true;
			else
				return false;
		});
		if (tmpIter == pointxs.end())
//		if (true)
		{
			Vec2i pt = toAddCrossPoints[j].getXY();
//			cout << pt << endl;
			circle(color_img, Point(pt[0], pt[1]), 7, (255, 255, 255), 2);
		}
		else
			toAddCrossPoints.erase(toAddCrossPoints.begin() + j--);
	}

//	for (auto m = 0; m < linexs.size(); ++m)
//	{
//		for (auto n = m + 1; n < linexs.size(); ++n)
//		{
//			Vec2i pt1, pt2, pt3, pt4;
//			line2pt(linexs[m].getLineVec(pointxs), pt1, pt2);
//			line2pt(linexs[n].getLineVec(pointxs), pt3, pt4);
//			int id1, id2, id3, id4;
//			id1 = linexs[m].getPt1Id(); id2 = linexs[m].getPt2Id(); id3 = linexs[n].getPt1Id(); id4 = linexs[n].getPt2Id();
//			if (pt1 == pt3)
//			{
//				cout << "pt1 and pt3 same point" << endl;
//				if (basicDashLineRev(oriEdgePoints, pt2, pt4))
//				{
//					int pid1 = (linexs[m].getPt2Id());
//					int pid2 = (linexs[n].getPt2Id());
//					int lid = int(linexs.size());
//					line_class tmpNewline(pid1, pid2, lid);
//					linexs.push_back(tmpNewline);
//				}
//			}
//			else if (pt1 == pt4)
//			{
//				cout << "pt1 and pt4 same point" << endl;
//				if (basicDashLineRev(oriEdgePoints, pt2, pt3))
//				{
//					int pid1 = (linexs[m].getPt2Id());
//					int pid2 = (linexs[n].getPt1Id());
//					int lid = int(linexs.size());
//					line_class tmpNewline(pid1, pid2, lid);
//					linexs.push_back(tmpNewline);
//				}
//			}
//			else if (pt2 == pt3)
//			{
//				cout << "pt2 and pt3 same point" << endl;
//				if (basicDashLineRev(oriEdgePoints, pt1, pt4))
//				{
//					int pid1 = (linexs[m].getPt1Id());
//					int pid2 = (linexs[n].getPt2Id());
//					int lid = int(linexs.size());
//					line_class tmpNewline(pid1, pid2, lid);
//					linexs.push_back(tmpNewline);
//				}
//			}
//			else if (pt2 == pt4)
//			{
//				cout << "pt2 and pt4 same point" << endl;
//				if (basicDashLineRev(oriEdgePoints, pt1, pt3))
//				{
//					int pid1 = (linexs[m].getPt1Id());
//					int pid2 = (linexs[n].getPt1Id());
//					int lid = int(linexs.size());
//					line_class tmpNewline(pid1, pid2, lid);
//					linexs.push_back(tmpNewline);
//				}
//			}
//			else
//			{
	
	for (auto iter1 = pointxs.begin(); iter1 != pointxs.end(); ++iter1)
	{ 
		for (auto iter2 = pointxs.begin(); iter2 != pointxs.end(); ++iter2)
		{
			Vec2i pt1 = iter1->getXY(); Vec2i pt2 = iter2->getXY();
			if (!existRealLineWithinPtxs(linexs, pointxs,pt1,pt2))
			{
				Vec2i tmpPt1 = iter1->getXY(); Vec2i tmpPt2 = iter2->getXY();
				int id1 = iter1->getPid(); int id2 = iter2->getPid();
				if (basicDashLineRev(withouCnlBw_ept, tmpPt1,tmpPt2, 0.4))
				{
					int pid1 = id1;
					int pid2 = id2;
					int lid = int(linexs.size());
					line_class tmpNewline(pid1, pid2, lid);
					linexs.push_back(tmpNewline);
				}
			}
		}
	}
	for (auto i = 0; i < linexs.size(); i++)
	{
		Vec4i l = linexs[i].getLineVec(pointxs);
		Vec2i pt1 = { l[0], l[1] };
		Vec2i pt2 = { l[2], l[3] };
		//		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp, 2);
		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp, 2);
	}
	cout << endl;
	//if (showFlag)
	{
		namedWindow("8.lines first opt version now", 0);
		imshow("8.lines first opt version now", color_img);
	}
	drawedImages = color_img;
}

void detect_line_lsd1(Mat diagram_segment, Mat diagram_segwithoutcircle, Mat& withoutCirBw, vector<point_class> pointXs, vector<circle_class>& circles, Mat& color_img, vector<line_class> lineXs, vector<Point2i>& oriEdgePoints, Mat& drawedImages, bool showFlag = true, string fileName = "")
{
	vector<Point2i> withoutO_ept = getPointPositions(withoutCirBw);

	vector<Vec4i> plainLines = {};

#pragma region raw detection

	Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_STD);
	vector<Vec4f> line_std;
	//	ofstream tmpLogFile;
	//	tmpLogFile.open("tmpLog.txt");
	ls->detect(diagram_segwithoutcircle, line_std);
	
	for (size_t i = 0; i < line_std.size(); ++i)
	{
		//eliminate the false detected lines with few points on it
		Vec4i rawLine = f2i(line_std[i]);

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

#pragma endregion raw detection

	/* then we handle the lines specifically*/
	/* for lines should be combined 1. colinear 2. */

	/*display*/
	Mat oriL_img = color_img.clone();
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i]; Vec2i pt1 = { l[0], l[1] }; Vec2i pt2 = { l[2], l[3] };
		//			cout << "*********************" << pt1 << " " << pt2 << endl;
		line(oriL_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
	}

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
		cout << "stop" << endl;
		for (auto iter2 = iter1 + 1; iter2 != plainLines.end();)
		{
			Vec4i line1 = *iter1;
			Vec2i pt1 = { line1[0], line1[1] };
			Vec2i pt2 = { line1[2], line1[3] };
			cout << pt1 << "pt" << pt2 << endl;

			Vec4i line2 = *iter2;
			Vec2i pt3 = { line2[0], line2[1] };
			Vec2i pt4 = { line2[2], line2[3] };
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
					sort(tmpPtVec.begin(), tmpPtVec.end(), ptXiSortPred);
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

	Mat rmParaI = color_img.clone();
	/*display*/
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i l = plainLines[i];
		Vec2i pt1 = { l[0], l[1] };
		Vec2i pt2 = { l[2], l[3] };
		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(rmParaI, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 2, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(rmParaI, Point{ pt1[0], pt1[1] }, 10, tmp);
		circle(rmParaI, Point{ pt2[0], pt2[1] }, 10, tmp);
	}
	namedWindow("rm parallel", 0);
	imshow("rm parallel", rmParaI);

	// remove the line detected and store it as withoutCLOriBw

	/*************cover the detected line with white pixel*/
	Mat withoutCnLBw = withoutCirBw.clone();
	Mat withoutLBw = diagram_segment.clone();
	for (auto i = 0; i < plainLines.size(); i++)
	{
		Vec4i pl = plainLines[i];
		Vec2i pt1 = { pl[0], pl[1] };
		Vec2i pt2 = { pl[2], pl[3] };
		cv::line(withoutCnLBw, pt1, pt2, 0, 3, 8, 0);
		cv::line(withoutLBw, pt1, pt2, 0, 3, 8, 0);
	}
	vector<Point2i> withouCnlBw_ept = getPointPositions(withoutCnLBw);
	vector<Point2i> withoutL_ept = getPointPositions(withoutLBw);

	bool testflag = false;

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
		Vec2i tmpPt1, tmpPt2; line2pt(plainLines[i], tmpPt1, tmpPt2);
		if (norm(tmpPt1, tmpPt2) > 20)
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
			delete ptx1, ptx2, lx;
		}
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
		Vec2i pt1 = { l[0], l[1] };
		Vec2i pt2 = { l[2], l[3] };
		cout << pt1 << " " << pt2 << endl;
		cout << linexs[i].getPt1Id() << " " << linexs[i].getPt2Id() << endl;
	}
	cout << "****************sep" << endl << endl;

	// with the raw detection result, we may get extra false line or omit some line(eg. dash line falsely removed previously or real line not detected)
	// in order to detect all lines, we suppose to first put all candidates into calculation then remove the false line at the last step.
	/*******************raw line candidates refinement*/
	map<int, int> change000;
	vector<point_class> toAddCrossPoints;
	for (auto i = 0; i < linexs.size(); i++)
	{
		line_class* linex1 = &linexs[i];
		cout << endl << "**********" << linex1->getLineVec(pointxs) << endl;
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
				double angle = llAngle(linex1_vec, linex2_vec);
				if ((angle>15) && (!isInImage(diagram_segwithoutcircle.cols, diagram_segwithoutcircle.rows, raw_cross)))
				{
					cout << raw_cross << endl;
					cout << "cross out of scope, then take it as no cross" << endl;
					continue;
				}
				else
				{
					//cross in scope
					cout << linex1_vec << endl << linex2_vec << endl;
					cross_refinement(raw_cross, linex1, linex2, circles, pointxs, toAddCrossPoints, withouCnlBw_ept, oriEdgePoints);
					cout << "step separate" << endl;
				}
			}
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
	cout << "after cross points refinement" << endl;
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
				if (len < cir[2] / 2)
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
						auto rm_iter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt1Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;

					}
					else if (not_find_pt2_flag && !not_find_pt1_flag)
					{
						auto rm_iter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt2Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;

					}
					else if (not_find_pt1_flag && not_find_pt2_flag)
					{
						auto rm_iter1 = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt1Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter1);
						auto rm_iter2 = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						{
							if (a.getPid() == iter->getPt2Id())
								return true;
							else
								return false;
						});
						pointxs.erase(rm_iter2);
						iter = linexs.erase(iter);
						rm_flag = true; rm_line_num++;
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
	cout << "after remove isolated short lines due to the non-complete circle removal " << endl;
	for (auto i = 0; i < linexs.size(); i++)
	{
		cout << linexs[i].getLineVec(pointxs) << linexs[i].getPt1Id() << ", " << linexs[i].getPt2Id() << endl;
	}
	for (auto j = 0; j < pointxs.size(); j++)
	{
		cout << pointxs[j].getPid() << "  " << pointxs[j].getXY() << endl;
	}

	/*set similar points identical*/


	//	cout << "now remove indentical vec points and reorder the indices";
	// point index indicate the line id in which it locate
	map<int, int> changeMap;

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
		cout << linexs[j].getPt1Id() << "->" << changeMap2[linexs[j].getPt1Id()] << endl;
		cout << linexs[j].getPt2Id() << "->" << changeMap2[linexs[j].getPt2Id()] << endl;
		linexs[j].setpt1Id(changeMap2[linexs[j].getPt1Id()]);
		linexs[j].setpt2Id(changeMap2[linexs[j].getPt2Id()]);
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

	cout << endl;
	for (auto j = 0; j < toAddCrossPoints.size(); j++)
	{
		auto tmpIter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
		{
			if (same_pt(a, toAddCrossPoints[j]))
				return true;
			else
				return false;
		});
		if (tmpIter == pointxs.end())
			//		if (true)
		{
			Vec2i pt = toAddCrossPoints[j].getXY();
			//			cout << pt << endl;
			circle(color_img, Point(pt[0], pt[1]), 7, (255, 255, 255), 2);
		}
		else
			toAddCrossPoints.erase(toAddCrossPoints.begin() + j--);
	}

	for (auto iter1 = pointxs.begin(); iter1 != pointxs.end(); ++iter1)
	{
		for (auto iter2 = pointxs.begin(); iter2 != pointxs.end(); ++iter2)
		{
			Vec2i pt1 = iter1->getXY(); Vec2i pt2 = iter2->getXY();
			if (!existRealLineWithinPtxs(linexs, pointxs, pt1, pt2))
			{
				Vec2i tmpPt1 = iter1->getXY(); Vec2i tmpPt2 = iter2->getXY();
				int id1 = iter1->getPid(); int id2 = iter2->getPid();
				if (basicDashLineRev(withouCnlBw_ept, tmpPt1, tmpPt2, 0.4))
				{
					int pid1 = id1;
					int pid2 = id2;
					int lid = int(linexs.size());
					line_class tmpNewline(pid1, pid2, lid);
					linexs.push_back(tmpNewline);
				}
			}
		}
	}
	for (auto i = 0; i < linexs.size(); i++)
	{
		Vec4i l = linexs[i].getLineVec(pointxs);
		Vec2i pt1 = { l[0], l[1] };
		Vec2i pt2 = { l[2], l[3] };
		//		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{ pt1[0], pt1[1] }, 10, tmp, 2);
		circle(color_img, Point{ pt2[0], pt2[1] }, 10, tmp, 2);
	}
	cout << endl;
	//if (showFlag)
	{
		namedWindow("8.lines first opt version now", 0);
		imshow("8.lines first opt version now", color_img);
	}
	drawedImages = color_img;
}


Vec4f plineRet(Vec4f a, Vec4f b)
{
	Vec2f pt1, pt2, pt3, pt4;
	pt1 = { a[0], a[1] };pt2 = { a[2], a[3] };
	pt3 = { b[0], b[1] };pt4 = { b[2], b[3] };
	float x1, x2, x3, x4, y1, y2, y3, y4,x0,y0;
	x1 = pt1[0];x2 = pt2[0];x3 = pt3[0];x4 = pt4[0];
	y1 = pt1[1];y2 = pt2[1];y3 = pt3[1];y4 = pt4[1];
	float x_4[4] = { x1, x2, x3, x4 };
//	sort(x_4, x_4 + 4);
	
	bool vertical_flag = (abs(pt1[0] - pt2[0]) < abs(pt1[1] - pt2[1])) ? true : false;
	Vec4f ret;

//	bool overlap = (x3 - x1)*(x4 - x2) >= 0 ? true : false;
//	double frac = frac_compute(pt1, pt2, vertical_flag);
//	ret[2] = (x4 - frac*y1 + frac*frac*x1) / (1 + frac*frac);
//	x0 = ret[2];
//	ret[3] = y1 + frac*(x0 - x1);
	if (!vertical_flag)
	{

	}
	else
	{
		
	}
	return ret;
}

double lldis(Vec4f l1,Vec4f l2)
{
	Vec2f pt1, pt2;
	line2pt(l1, pt1, pt2);
	Vec2f mid = (pt1 + pt2) / 2;
	return pt2lineDis(l2, mid);
}

bool get_vertical_flag(Vec4f lf)
{
	Vec2f pt1, pt2; line2pt(lf, pt1, pt2);
	if (abs(pt1[0] - pt2[0]) < 3)
		return true;
	else
		return false;
}
bool get_vertical_flag(Vec4i lf)
{
	Vec2i pt1, pt2; line2pt(lf, pt1, pt2);
	if (abs(pt1[0] - pt2[0]) < 3)
		return true;
	else
		return false;
}

void line_l(Mat img,Vec4i line,const Scalar& scalar,int thickness)
{
	Vec2i pt1, pt2; line2pt(line,pt1, pt2);
	cv::line(img, pt1, pt2, scalar, thickness);
}

void detect_line_lsd(Mat diagram_segment, Mat diagram_segwithoutcircle, Mat& withoutCirBw, vector<point_class> pointXs, vector<circle_class>& circles, Mat& color_img, vector<line_class> lineXs, vector<Point2i>& oriEdgePoints, Mat& drawedImages, bool showFlag = true, string fileName = "")
{
	Ptr<LineSegmentDetector> ls = createLineSegmentDetector(LSD_REFINE_STD);
	vector<Vec4f> line_std;
	//	ofstream tmpLogFile;
	//	tmpLogFile.open("tmpLog.txt");
	ls->detect(diagram_segwithoutcircle, line_std);
	Mat drawLines(diagram_segwithoutcircle);
	cvtColor(drawLines, drawLines, CV_GRAY2RGB);
	Mat drawLines0 = drawLines.clone();
	ls->drawSegments(drawLines0, line_std);
	imshow("standard refinement", drawLines);
	for (auto iter = line_std.begin(); iter != line_std.end(); )
	{
		//rm the too short lines
		Vec2f pt1, pt2;
		line2pt(*iter, pt1, pt2);
		if (p2pdistance(pt1, pt2) < 10)
			iter = line_std.erase(iter);
		else
		{
			if (get_vertical_flag(*iter))
			{
				cout << " the line is thought to be vertical" << endl;
				if (pt2[1] < pt1[1])
					*iter = pt2line(pt2, pt1);
			}
			else
			{
				cout << "not verical" << endl;
				if (pt2[0] < pt1[0])
					*iter = pt2line(pt2, pt1);
			}
			++iter;
		}

	}


	int corlinear_num = 0;
	Mat tmpImg;
	for (auto iter1 = line_std.begin(); iter1 != line_std.end(); ++iter1)
	{
		Vec4f line1 = *iter1;
		for (auto iter2 = iter1+1; iter2 != line_std.end(); )
		{
			//			Vec4f line1 = *iter1; Vec4f line2 = *iter2;
			Vec2f pt1, pt2;  line2pt(line1, pt1, pt2);
			tmpImg = drawLines.clone();
			line(tmpImg, f2i(pt1), f2i(pt2), Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 0), 1);
			circle(tmpImg, f2i(pt1), 10, Scalar(255, 255, 0), 1);
			circle(tmpImg, f2i(pt2), 10, Scalar(255, 255, 0), 1);
			cout << "sep" << endl;
			Vec4f line2 = *iter2;
			if (abs(line1[3] - line1[1]) < 3 && abs(line2[3] - line2[1])< 3)
				cout << "stop" << endl;
			int id1, id2; id1 = int(iter1 - line_std.begin()); id2 = int(iter2 - line_std.begin());
			Vec2f pt3, pt4;
			line2pt(line2, pt3, pt4);
			line(tmpImg, f2i(pt3), f2i(pt4), Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 0), 1);
			circle(tmpImg, f2i(pt3), 5, Scalar(255, 0, 255), 1);
			circle(tmpImg, f2i(pt4), 5, Scalar(255, 0, 255), 1);
			Vec4i test1 = { 57, 34, 170, 59 };
			Vec4i test2 = { 41, 54, 68, 68 };
			if (test1 == f2i(line1))
				cout << "test stop" << endl;
			if (test2 == f2i(line2))
				cout << "stop" << endl;
			if (isParallel(line1, line2))
			{
				cout << "the two line is parallel" << endl;
				double tmp_lldis = lldis(line1, line2);
				//				Vec2f pt1, pt2, pt3, pt4;
				//				line(tmpImg, f2i(pt3), f2i(pt4), Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 0), 1);
				//				circle(tmpImg, f2i(pt3), 5, Scalar(255, 0, 255), 1);
				//				circle(tmpImg, f2i(pt4), 5, Scalar(255, 0, 255), 1);
				cout << "the distance between the two line is " << tmp_lldis << endl;
				if (tmp_lldis < 5)
				{
					cout << "collinear,then combine" << endl;
					cout << "line " << id1 << ": " << *iter1 << endl << "line " << id2 << ": " << *iter2 << endl;
					//					Vec2f pt1, pt2, pt3, pt4;
					//					line2pt(line1, pt1, pt2); line2pt(line2, pt3, pt4);
					//					Mat tmpImg = drawLines.clone();
					//					line(tmpImg, f2i(pt1), f2i(pt2), Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 0), 1);
					//					line(tmpImg, f2i(pt3), f2i(pt4), Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 0), 1);
					//					circle(tmpImg, f2i(pt1), 10, Scalar(255, 255, 0), 1);
					//					circle(tmpImg, f2i(pt2), 10, Scalar(255, 255, 0), 1);
					//					circle(tmpImg, f2i(pt3), 5, Scalar(255, 0, 255), 1);
					//					circle(tmpImg, f2i(pt4), 5, Scalar(255, 0, 255), 1);

//					vector<Vec2f> tmp;
//					tmp.push_back(pt1);
//					tmp.push_back(pt2);
//					tmp.push_back(pt3);
//					tmp.push_back(pt4);
//					sort(tmp.begin(), tmp.end(), ptXfSortPred);

					Vec2f newPt1, newPt2;
//					newPt1 = (p2pdistance(tmp[0], tmp[3]) < p2pdistance(tmp[1], tmp[3])) ? tmp[1] : tmp[0];
//					newPt2 = (p2pdistance(tmp[3], tmp[0]) < p2pdistance(tmp[2], tmp[0])) ? tmp[2] : tmp[3];
					bool tmpVFlag = get_vertical_flag(line1);
					
					if (!tmpVFlag)
					{
						if (pt1[0] < pt3[0])
						{
							newPt1 = (p2pdistance(pt1, pt4) > p2pdistance(pt3, pt4)) ? pt1 : pt3;
						}
						else
						{
							if (p2pdistance(pt3, pt4) > p2pdistance(pt1, pt4))
								newPt1 = pt3;
							else
								newPt1 = pt1;
							newPt1 = (p2pdistance(pt3, pt4) > p2pdistance(pt1, pt4)) ? pt3 : pt1;
						}
						if (pt4[0] > pt2[0])
						{
							if (p2pdistance(pt4, pt1) > p2pdistance(pt2, pt1))
								newPt2 = pt4;
							else
								newPt2 = pt2;
						}
						else
						{
							if (p2pdistance(pt2, pt1) > p2pdistance(pt4, pt1))
								newPt2 = pt2;
							else
								newPt2 = pt4;
						}
					}
					else
					{
						if (pt1[1] < pt3[1])
						{
							if (p2pdistance(pt1, pt4) > p2pdistance(pt3, pt4))
								newPt1 = pt1;
							else
								newPt1 = pt3;
						}
						else
						{
							if (p2pdistance(pt3, pt4) > p2pdistance(pt1, pt4))
								newPt1 = pt3;
							else
								newPt1 = pt1;
						}
						if (pt4[1] > pt2[1])
						{
							if (p2pdistance(pt4, pt1) > p2pdistance(pt2, pt1))
								newPt2 = pt4;
							else
								newPt2 = pt2;
						}
						else
						{
							if (p2pdistance(pt2, pt1) > p2pdistance(pt4, pt1))
								newPt2 = pt2;
							else
								newPt2 = pt4;
						}
					}
					Vec4f newline = pt2line(newPt1, newPt2);

					circle(tmpImg, f2i(newPt1), 7, Scalar(0, 255, 255), 1);
					circle(tmpImg, f2i(newPt2), 7, Scalar(0, 255, 255), 1);
					line(tmpImg, f2i(newPt1), f2i(newPt2), Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 0), 2);
					cout << "both line -> " << newline << endl;
					*iter1 = newline; 
//					*iter2 = newline;
					
					cout << line_std[id1] << "  " << line_std[id2] << endl;
					iter2 = line_std.erase(iter2);
					cout << "sep" << endl;
				}
				else
				{
					++iter2;
				}
			}
			else
			{
				++iter2;
				cout << "not parallel and no change" << endl;
			}
		}
	}

	sort(line_std.begin(), line_std.end(), [](Vec4f a, Vec4f b){
		if (a[0] < b[0])
			return true;
		else
			return false;
	});
	line_std.erase(unique(line_std.begin(), line_std.end(), [](Vec4f a, Vec4f b)
	{
		Vec2i pt1, pt2, pt3, pt4;
		pt1 = { int(a[0]), int(a[1]) }; pt2 = { int(a[2]), int(a[3]) };
		pt3 = { int(b[0]), int(b[1]) }; pt4 = { int(b[2]), int(b[3]) };
		if ((a == b) || (same_pt(pt1, pt3) && same_pt(pt2, pt4)) || (same_pt(pt1, pt4)&same_pt(pt2, pt3)))
		{
			return true;
		}
		else
			return false;
	}), line_std.end());

	
//	for (auto iter = line_std.begin(); iter != line_std.end(); ++iter)
//	{
//		Vec2i p1, p2;
//		Scalar tmp = Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)));
//		p1 = { int((*iter)[0]), int((*iter)[1]) };
//		p2 = { int((*iter)[2]), int((*iter)[3]) };
//		circle(drawLines, p1, 10 * (rand() / double(RAND_MAX)) + 5, tmp, 1);
//		circle(drawLines, p2, 5, tmp, 1);
//		line(drawLines, p1, p2, Scalar(0, 255, 255), 1);
//	}
//	cout << "stop" << endl;

//	// initialize linex and pointx
//	int px_count = 0; int lx_count = 0;
//	vector<point_class> pointxs;
//	vector<line_class> linexs;
//	for (auto i = 0; i < line_std.size(); i++)
//	{
//		Vec2i tmpPt1, tmpPt2; line2pt(line_std[i], tmpPt1, tmpPt2);
//		//		if (norm(tmpPt1, tmpPt2) > 20)
//		{
//			auto p_l = line_std[i];
//			Vec2i pt1, pt2;
//			line2pt(p_l, pt1, pt2);
//			auto _pidx1 = px_count++;
//			int _pidx2 = px_count++;
//			int _lidx = lx_count++;
//			point_class* ptx1 = new point_class(pt1, _pidx1);
//			point_class* ptx2 = new point_class(pt2, _pidx2);
//			//cout << &ptx1 << endl<<&ptx2<<endl;
//			ptx1->pushLid(_lidx);
//			ptx2->pushLid(_lidx);
//			pointxs.push_back(*ptx1);
//			pointxs.push_back(*ptx2);
//
//			line_class* lx = new line_class(ptx1->getPid(), ptx2->getPid(), _lidx);
//			//lx.p1 = &pointxs[_pidx1]; lx.p2 = &pointxs[_pidx2];
//			linexs.push_back(*lx);
//			delete ptx1, ptx2, lx;
//		}
//	}
//	for (auto iter1 = linexs.begin(); iter1 != linexs.end(); ++iter1)
//	{
//		Vec4i lx1 = iter1->getLineVec(pointxs);
//		Vec2i pt1, pt2; line2pt(lx1, pt1, pt2);
//		for (auto iter2 = iter1 + 1; iter2 != linexs.end(); ++iter2)
//		{
//			Vec4i lx2 = iter2->getLineVec(pointxs);
//			Vec2i pt3, pt4; line2pt(lx2, pt3, pt4);
//			if (same_pt(pt1, pt3))
//			{
//				int old_id = iter2->getPt1Id();
//				int new_id = iter1->getPt1Id();
//				iter2->setpt1Id(new_id);
//				rm_point_by_id(pointxs, old_id);
//				cout << "id " << old_id << " -> " << "id " << new_id << endl;
//			}
//			else if (same_pt(pt1, pt4))
//			{
//				int old_id = iter2->getPt2Id();
//				int new_id = iter1->getPt1Id();
//				iter2->setpt2Id(new_id);
//				rm_point_by_id(pointxs, old_id);
//				cout << "id " << old_id << " -> " << "id " << new_id << endl;
//			}
//			else if (same_pt(pt2, pt3))
//			{
//				int old_id = iter2->getPt1Id();
//				int new_id = iter1->getPt2Id();
//				iter2->setpt1Id(new_id);
//				rm_point_by_id(pointxs, old_id);
//				cout << "id " << old_id << " -> " << "id " << new_id << endl;
//			}
//			else if (same_pt(pt2, pt4))
//			{
//				int old_id = iter2->getPt2Id();
//				int new_id = iter1->getPt2Id();
//				iter2->setpt2Id(new_id);
//				rm_point_by_id(pointxs, old_id);
//				cout << "id " << old_id << " -> " << "id " << new_id << endl;
//			}
//			else
//			{
//				cout << "no id change here" << endl;
//			}
//		}
//	}
//	//recovery
////	for (auto iter1 = linexs.begin(); iter1 != linexs.end(); ++iter1)
////	{
////		Vec4i lx1 = iter1->getLineVec(pointxs);
////		Vec2i pt1, pt2; line2pt(lx1, pt1, pt2);
////		for (auto iter2 = iter1 + 1; iter2 != linexs.end(); ++iter2)
////		{
////			Vec4i lx2 = iter2->getLineVec(pointxs);
////			Vec2i pt3, pt4; line2pt(lx2, pt3, pt4);
////			int id1, id2, id3, id4;
////			id1 = iter1->getPt1Id(); id2 = iter1->getPt2Id(); id3 = iter2->getPt1Id(); id4 = iter2->getPt2Id();
////			if (id1 == id3)
////			{
////				if (!existRealLineWithinPtxs(linexs, pointxs, pointxs[id2], pointxs[id4]) && dashLineRecovery(oriEdgePoints, pt2, pt4, circles))
////				{
////					line_class newlx(id2, id4, linexs.size());
////					linexs.push_back(newlx);
////				}
////			}
////			else if (id1 == id4)
////			{
////				if (!existRealLineWithinPtxs(linexs, pointxs, pointxs[id2], pointxs[id3]) && dashLineRecovery(oriEdgePoints, pt2, pt3, circles))
////				{
////					line_class newlx(id2, id3, linexs.size());
////					linexs.push_back(newlx);
////				}
////			}
////			else if (id2 == id3)
////			{
////				if (!existRealLineWithinPtxs(linexs, pointxs, pointxs[id1], pointxs[id4]) && dashLineRecovery(oriEdgePoints, pt1, pt4, circles))
////				{
////					line_class newlx(id1, id4, linexs.size());
////					linexs.push_back(newlx);
////				}
////			}
////			else if (id2 == id4)
////			{
////				if (!existRealLineWithinPtxs(linexs, pointxs, pointxs[id1], pointxs[id3]) && dashLineRecovery(oriEdgePoints, pt1, pt3, circles))
////				{
////					line_class newlx(id2, id3, linexs.size());
////					linexs.push_back(newlx);
////				}
////			}
////			else
////			{
////				Vec2f tmp_cross;
////				getCrossPt(lx1, lx2, tmp_cross);
////				Vec2i i_tmp_cross = f2i(tmp_cross);
////				//line1
////				if (same_pt(i_tmp_cross, pt1))
////				{
////					cout << "cross == pt1" << endl;
////					double dis1 = norm(pt3, i_tmp_cross);
////					double dis2 = norm(pt4, i_tmp_cross);
////					int pt1id = iter1->getPt1Id();
////					if (dis1 < dis2)
////					{
////						if (dashLineRecovery(oriEdgePoints, pt3, i_tmp_cross, circles))
////						{
////							cout << "assign pt1 to line2 pt3" << endl;
////							iter2->setpt1Id(pt1id);
////						}
////					}
////					else
////					{
////						if (dashLineRecovery(oriEdgePoints, pt4, i_tmp_cross, circles))
////						{
////							cout << "assign pt1 to line2 pt4" << endl;
////							iter2->setpt2Id(pt1id);
////						}
////					}
////				}
////				else if (same_pt(i_tmp_cross, pt2))
////				{
////					cout << "cross == pt2" << endl;
////					double dis1 = norm(pt3, i_tmp_cross);
////					double dis2 = norm(pt4, i_tmp_cross);
////					int pt2id = iter1->getPt2Id();
////					if (dis1 < dis2)
////					{
////						if (dashLineRecovery(oriEdgePoints, pt3, i_tmp_cross, circles))
////						{
////							cout << "assign pt2 to line2 pt3" << endl;
////							iter2->setpt1Id(pt2id);
////						}
////					}
////					else
////					{
////						if (dashLineRecovery(oriEdgePoints, pt4, i_tmp_cross, circles))
////						{
////							cout << "assign pt2 to line2 pt4" << endl;
////							iter2->setpt2Id(pt2id);
////						}
////					}
////				}
////				else if (same_pt(i_tmp_cross, pt3))
////				{
////					cout << "cross == pt3" << endl;
////					double dis1 = norm(pt1, i_tmp_cross);
////					double dis2 = norm(pt2, i_tmp_cross);
////					int pt3id = iter2->getPt1Id();
////					if (dis1 < dis2)
////					{
////						if (dashLineRecovery(oriEdgePoints, pt3, i_tmp_cross, circles))
////						{
////							cout << "assign pt3 to line1 pt1" << endl;
////							iter1->setpt1Id(pt3id);
////						}
////					}
////					else
////					{
////						if (dashLineRecovery(oriEdgePoints, pt4, i_tmp_cross, circles))
////						{
////							cout << "assign pt3 to line1 pt2" << endl;
////							iter1->setpt2Id(pt3id);
////						}
////					}
////				}
////				else if (same_pt(i_tmp_cross, pt4))
////				{
////					cout << "cross == pt4" << endl;
////					double dis1 = norm(pt1, i_tmp_cross);
////					double dis2 = norm(pt2, i_tmp_cross);
////					int pt4id = iter2->getPt2Id();
////					if (dis1 < dis2)
////					{
////						if (dashLineRecovery(oriEdgePoints, pt3, i_tmp_cross, circles))
////						{
////							cout << "assign pt4 to line1 pt1" << endl;
////							iter1->setpt1Id(pt4id);
////						}
////					}
////					else
////					{
////						if (dashLineRecovery(oriEdgePoints, pt4, i_tmp_cross, circles))
////						{
////							cout << "assign pt4 to line1 pt2" << endl;
////							iter1->setpt2Id(pt4id);
////						}
////					}
////				}
////				else if (in_line(lx1, i_tmp_cross))
////				{
////					cout << "cross in line1" << endl;
////					if (in_line(lx2, i_tmp_cross))
////					{
////						cout << "cross in line2, inner cross" << endl;
////
////					}
////					else
////					{
////						cout << "" << endl;
////					}
////				}
////
////			}
////		}
////	}

	Mat withoutCnLBw = withoutCirBw.clone();
	Mat withoutLBw = diagram_segment.clone();
	for (auto i = 0; i < line_std.size(); i++)
	{
		Vec4i pl = line_std[i];
		Vec2i pt1 = {pl[0], pl[1]};
		Vec2i pt2 = {pl[2], pl[3]};
		cv::line(withoutCnLBw, pt1, pt2, 0, 3, 8, 0);
		cv::line(withoutLBw, pt1, pt2, 0, 3, 8, 0);
	}
	vector<Point2i> withouCnlBw_ept = getPointPositions(withoutCnLBw);
	int px_count = 0;
	int lx_count = 0;
	vector<line_class> linexs;
	vector<point_class> pointxs;
	/***********initialize lines and points**********/
	for (auto i = 0; i < line_std.size(); i++)
	{
		Vec2i tmpPt1, tmpPt2;
		line2pt(line_std[i], tmpPt1, tmpPt2);
		if (norm(tmpPt1, tmpPt2) > 20)
		{
			auto p_l = line_std[i];
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
	}
	
	for (auto iter = linexs.begin(); iter != linexs.end(); ++iter)
	{
		Vec2i p1, p2;
		Scalar tmp = Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)));
		p1 = { int((iter->getLineVec(pointxs))[0]), int((iter->getLineVec(pointxs))[1]) };
		p2 = { int((iter->getLineVec(pointxs))[2]), int((iter->getLineVec(pointxs))[3]) };
		circle(drawLines, p1, 10 * (rand() / double(RAND_MAX)) + 5, tmp, 1);
		circle(drawLines, p2, 5, tmp, 1);
		line(drawLines, p1, p2, tmp, 1);
		cout << p1 << " " << p2 << endl;
	}
	cout << "stop" << endl;

	map<int, int> change000;
	vector<point_class> toAddCrossPoints;
	Mat tmpRef = color_img.clone();
	for (auto i = 0; i < linexs.size(); i++)
	{
		line_class* linex1 = &linexs[i];
		cout << endl << "**********" << linex1->getLineVec(pointxs) << endl;
		for (auto j = i + 1; j < linexs.size(); j++)
		{
			line_class* linex2 = &linexs[j];
			Vec4i linex1_vec = linex1->getLineVec(pointxs);
			Vec4i linex2_vec = linex2->getLineVec(pointxs);

			//check whether the two line are parallel, if so, jump ahead
			Mat tmp = color_img.clone();
			Vec2i pt1, pt2, pt3, pt4;
			line2pt(linex1_vec, pt1, pt2);
			line2pt(linex2_vec, pt3, pt4);
			line(tmp, pt1, pt2, Scalar(0, 255, 255), 1);
			line(tmp, pt3, pt4, Scalar(0, 255, 255), 1);
			if (isParallel(linex1_vec, linex2_vec))
			{
				//cout << linex1_vec << endl << linex2_vec << endl;
				cout << "the two line is parallel but not collinear" << endl;
				continue;
			}
			else
			{
				// not parallel, then calculate the rough cross first
				Vec2f raw_cross;
//				Vec2i pt1, pt2, pt3, pt4;
//				line2pt(linex1_vec, pt1, pt2);
//				line2pt(linex2_vec, pt3, pt4);
				getCrossPt(linex1_vec, linex2_vec, raw_cross);
				cout << linex1_vec << endl << linex2_vec << endl;
//				cout << endl << "raw cross" << raw_cross << endl;
				double angle = llAngle(linex1_vec, linex2_vec);
				
				if ((angle > 15) && (!isInImage(diagram_segwithoutcircle.cols, diagram_segwithoutcircle.rows, raw_cross)))
				{
					cout << raw_cross << endl;
					cout << "cross out of scope, then take it as no cross" << endl;
					continue;
				}
				else
				{
					//cross in scope
//					cout << linex1_vec << endl << linex2_vec << endl;
					circle(tmp, f2i(raw_cross), 5, Scalar(255, 255, 0), 1);
					cross_refinement(raw_cross, linex1, linex2, circles, pointxs, toAddCrossPoints, withouCnlBw_ept, oriEdgePoints);

					//					cout << linex1->getLineVec(pointxs)<< linex1->getPt1Id() << "," << linex1->getPt2Id() << endl;
					//					cout << linex2->getLineVec(pointxs) << linex2->getPt1Id() << "," << linex2->getPt2Id() << endl;
					cout << "step separate" << endl << endl;
				}
			}
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
	cout << "after cross points refinement" << endl;
	//rm isolated short lines due to the non-complete circle removal
	// 1. two points almost on circles 2. the length is short 3. the end points is ont on other lines
	int rm_pt_num = 0;
	int rm_line_num = 0;
	map<int, int> changeMap0;
	Mat test = color_img.clone();
	for (auto iter = linexs.begin(); iter != linexs.end(); ++iter)
	{
		Vec2i p1, p2;
		Scalar tmp = Scalar(255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)), 255 * (rand() / double(RAND_MAX)));
		p1 = { int((iter->getLineVec(pointxs))[0]), int((iter->getLineVec(pointxs))[1]) };
		p2 = { int((iter->getLineVec(pointxs))[2]), int((iter->getLineVec(pointxs))[3]) };
		circle(test, p1, 10 * (rand() / double(RAND_MAX)) + 5, tmp, 1);
		circle(test, p2, 5, tmp, 1);
		line(test, p1, p2, tmp, 1);
		cout << p1 << " " << p2 << endl;
	}
	cout << "stop" << endl;
	for (auto iter = linexs.begin(); iter != linexs.end();)
	{
		Vec2i pt1, pt2;
		line2pt(iter->getLineVec(pointxs), pt1, pt2);
		cout << "pt1 " << pt1 << " and pt2 " << pt2 << endl;
		Mat tmp = color_img.clone();
		line(tmp, pt1, pt2, Scalar(0, 255, 255), 1);
		bool rm_flag = false;
		for (auto i = 0; i < circles.size(); ++i)
		{
			Vec3i cir = circles[i].getCircleVec();
			bool flag1, flag2;
			flag1 = on_circle(pt1, cir); flag2 = on_circle(pt2, cir);
			if (flag1&&flag2)
			{
				double len = p2pdistance(pt1, pt2);
				if (len < cir[2] / 3)
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
						auto rm_iter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						                       {
							                       if (a.getPid() == iter->getPt1Id())
								                       return true;
							                       else
								                       return false;
						                       });
						pointxs.erase(rm_iter);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;
					}
					else if (not_find_pt2_flag && !not_find_pt1_flag)
					{
						auto rm_iter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						                       {
							                       if (a.getPid() == iter->getPt2Id())
								                       return true;
							                       else
								                       return false;
						                       });
						pointxs.erase(rm_iter);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;
					}
					else if (not_find_pt1_flag && not_find_pt2_flag)
					{
						auto rm_iter1 = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						                        {
							                        if (a.getPid() == iter->getPt1Id())
								                        return true;
							                        else
								                        return false;
						                        });
						pointxs.erase(rm_iter1);
						auto rm_iter2 = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
						                        {
							                        if (a.getPid() == iter->getPt2Id())
								                        return true;
							                        else
								                        return false;
						                        });
						pointxs.erase(rm_iter2);
						iter = linexs.erase(iter);
						rm_flag = true;
						rm_line_num++;
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
	cout << "after remove isolated short lines due to the non-complete circle removal " << endl;
	for (auto i = 0; i < linexs.size(); i++)
	{
		cout << linexs[i].getLineVec(pointxs) << linexs[i].getPt1Id() << ", " << linexs[i].getPt2Id() << endl;
	}
	for (auto j = 0; j < pointxs.size(); j++)
	{
		cout << pointxs[j].getPid() << "  " << pointxs[j].getXY() << endl;
	}

	/*set similar points identical*/


	//	cout << "now remove indentical vec points and reorder the indices";
	// point index indicate the line id in which it locate
	map<int, int> changeMap;


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
		cout << linexs[j].getPt1Id() << "->" << changeMap2[linexs[j].getPt1Id()] << endl;
		cout << linexs[j].getPt2Id() << "->" << changeMap2[linexs[j].getPt2Id()] << endl;
		linexs[j].setpt1Id(changeMap2[linexs[j].getPt1Id()]);
		linexs[j].setpt2Id(changeMap2[linexs[j].getPt2Id()]);
	}


	/*display*/

	cout << endl;
	for (auto j = 0; j < toAddCrossPoints.size(); j++)
	{
		auto tmpIter = find_if(pointxs.begin(), pointxs.end(), [&](point_class a)
		                       {
			                       if (same_pt(a, toAddCrossPoints[j]))
				                       return true;
			                       else
				                       return false;
		                       });
		if (tmpIter == pointxs.end())
		//		if (true)
		{
			Vec2i pt = toAddCrossPoints[j].getXY();
			//			cout << pt << endl;
			circle(color_img, Point(pt[0], pt[1]), 7, (255 , 255 , 255), 2);
		}
		else
			toAddCrossPoints.erase(toAddCrossPoints.begin() + j--);
	}


	for (auto iter1 = pointxs.begin(); iter1 != pointxs.end(); ++iter1)
	{
		for (auto iter2 = pointxs.begin(); iter2 != pointxs.end(); ++iter2)
		{
			Vec2i pt1 = iter1->getXY();
			Vec2i pt2 = iter2->getXY();
			if (!existRealLineWithinPtxs(linexs, pointxs, pt1, pt2))
			{
				Vec2i tmpPt1 = iter1->getXY();
				Vec2i tmpPt2 = iter2->getXY();
				int id1 = iter1->getPid();
				int id2 = iter2->getPid();
				if (basicDashLineRev(withouCnlBw_ept, tmpPt1, tmpPt2, 0.4))
				{
					int pid1 = id1;
					int pid2 = id2;
					int lid = int(linexs.size());
					line_class tmpNewline(pid1, pid2, lid);
					linexs.push_back(tmpNewline);
				}
			}
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
	for (auto i = 0; i < linexs.size(); i++)
	{
		Vec4i l = linexs[i].getLineVec(pointxs);
		Vec2i pt1 = {l[0], l[1]};
		Vec2i pt2 = {l[2], l[3]};
		//		cout << "*********************" << pt1 << " " << pt2 << endl;
		line(color_img, pt1, pt2, Scalar(rand() % 255, rand() % 255, rand() % 255), 1, 8, 0);
		Scalar tmp = Scalar(rand() % 255, rand() % 255, rand() % 255);
		circle(color_img, Point{pt1[0], pt1[1]}, 10, tmp, 2);
		circle(color_img, Point{pt2[0], pt2[1]}, 10, tmp, 2);
	}
	cout << endl;
	//if (showFlag)
	{
		namedWindow("8.lines first opt version now", 0);
		imshow("8.lines first opt version now", color_img);
	}
	drawedImages = color_img;

}

void primitive_parse(const Mat binarized_image, const Mat diagram_segment, vector<Point2i>& oriEdgePoints, vector<point_class>& points, vector<line_class>& lines, vector<circle_class>& circles, Mat& drawedImages, bool showFlag = true, string fileName = "")
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
//	detect_circle3(diagram_segment, color_img, diagram_segwithoutcircle, withoutCirBw, circles, showFlag);
//	ransac_circle(diagram_segment, color_img, diagram_segwithoutcircle, withoutCirBw, circles, showFlag);
	detect_circle(diagram_segment, color_img, diagram_segwithoutcircle, withoutCirBw, circles, showFlag);
	
	// then the line detection
	//vector<Vec2i> basicEndpoints = {};


	detect_line3(diagram_segment, diagram_segwithoutcircle, withoutCirBw, points, circles, color_img, lines, oriEdgePoints, drawedImages, showFlag, fileName);
//	detect_line_lsd(diagram_segment, diagram_segwithoutcircle, withoutCirBw, points, circles, color_img, lines, oriEdgePoints, drawedImages, showFlag, fileName);

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
	vector<Mat> char_imgs;
	image_labelling(binarized_image, diagram_segment,char_imgs, true);
	for (auto i = 0; i < char_imgs.size();++i)
	{
		char fullNameStr[20];
		sprintf_s(fullNameStr, "charImg-%d.png", i);
		imwrite(fullNameStr, char_imgs[i]);
	}
	vector<point_class> points = {};
	vector<line_class> lines = {};
	vector<circle_class> circles = {};
	Mat drawedImages(image.size(), CV_8UC3);

	primitive_parse(binarized_image, diagram_segment, oriEdgePoints, points, lines, circles, drawedImages, false);
	// then we turn to handle the characters

	return 0;
}

int diagram()
{
	//a series of image
	//vector<Mat> images;
	char abs_path[100] = "D:\\data\\graph-DB\\nt6";
	char imageName[150], saveimgName[150];
	//string outputFN = "D:\\data\\graph-DB\\newtest6\\output.txt";
	int charCount = 0;
	for (int i = 1; i < 87; i++)
	{
		sprintf_s(imageName, "%s\\graph-%d.jpg", abs_path, i);
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
		vector<Point2i> oriEdgePoints = getPointPositions(binarized_image);
		/*Mat pointss = Mat::zeros(1000, 1000, CV_8UC3);
		for (auto i = 0; i < edgePoints.size(); i++)
		{
			Point2i pt = edgePoints[i];
			circle(pointss, pt, 1, Scalar(0, 0, 255));
		}
		namedWindow("points"); imshow("points", pointss);*/
		vector<Mat> char_imgs;
		image_labelling(binarized_image, diagram_segment,char_imgs);
		char tmpOutFolder[100]; char tmpCmd[100];
		sprintf_s(tmpOutFolder, "%s\\charImgs\\%d", abs_path, i);
		sprintf_s(tmpCmd, "mkdir %s", tmpOutFolder);
		system(tmpCmd);
		for (auto j = 0; j < char_imgs.size(); ++j)
		{
			char fullNameStr[100];
			sprintf_s(fullNameStr, "%s\\charImgs\\charImg-%d.png",abs_path, charCount++);
			imwrite(fullNameStr, char_imgs[j]);
			char subNameStr[100];
			sprintf_s(subNameStr, "%s\\charImgs\\%d\\charImg-%d.png",abs_path, i,j);
			imwrite(subNameStr, char_imgs[j]);

		}

		vector<point_class> points = {};
		vector<line_class> lines = {};
		vector<circle_class> circles = {};
		Mat drawedImages(image.size(), CV_8UC3);
		primitive_parse(binarized_image, diagram_segment, oriEdgePoints, points, lines, circles, drawedImages, false);
		imwrite(saveimgName, drawedImages);
	}
	return 0;
}
