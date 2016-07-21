# include "stdafx.h"
# include "diagram.h"
#include "zhangsuen.h"
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

vector<Point2f> getPointPositions(Mat bw)
{
	vector<Point2f> pointPositions;
	for (unsigned int y = 0; y < bw.rows; ++y)
	{
		for (unsigned int x = 0; x < bw.cols; ++x)
		{
			if (bw.at<unsigned char>(y, x)>0)
				pointPositions.push_back(Point2f(x, y));
		}
	}
	return pointPositions;
}

inline void getCircle(Point2f &p1, Point2f &p2, Point2f &p3, Point2f &center, float &radius)
{
	float x1 = p1.x;float x2 = p2.x;float x3 = p3.x;
	float y1 = p1.y;float y2 = p2.y;float y3 = p3.y;

	center.x = (x1*x1 + y1*y1)*(y2 - y3) + (x2*x2 + y2*y2)*(y3 - y1) + (x3*x3 + y3*y3)*(y1 - y2);
	center.x /= (2 * (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2));

	center.y = (x1*x1 + y1*y1)*(x3 - x2) + (x2*x2 + y2*y2)*(x1 - x3) + (x3*x3 + y3*y3)*(x2 - x1);
	center.y /= (2 * (x1*(y2 - y3) - y1*(x2 - x3) + x2*y3 - x3*y2));

	radius = sqrt((center.x - x1)*(center.x - x1) + (center.y - y1)*(center.y - y1));
}

// specific fucntions
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

void detect_circle(Mat diagram_segment, Mat &color_img,Mat &diagram_segwithoutcircle, vector<Vec3f> &circle_candidates)
{
	unsigned int circleN_todetect = 2;
	diagram_segwithoutcircle = diagram_segment;
	for (unsigned int i = 0; i < circleN_todetect; ++i)
	{
		vector<Point2f> edgePositions;
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
		if (bestCVal > 300)
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
			cv::circle(diagram_segwithoutcircle, bestCircleCenter, bestCircleRadius, 0, 3); // here the radius is fixed which isnt so nice.
			//cv::circle(color_img, bestCircleCenter, bestCircleRadius, Scalar(255, 0, 255), 3);
		}
	}
	namedWindow("graygeo blob without circles"); cv::imshow("graygeo blob without circles", diagram_segwithoutcircle);
	namedWindow("colorgeo"); cv::imshow("colorgeo", color_img);
}

int Sgn(double d)
{
	if (d<0)
		return -1;
	else
		return 1;
}

double p2pdistance(Vec2i pt1, Vec2i pt2)
{
	double distance;
	distance = sqrt(powf((pt1[0] - pt2[0]), 2) + powf((pt1[1] - pt2[1]), 2));
	return distance;
}

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


bool same_pt(Vec2i pt1, Vec2i pt2)
{
	double eps = 8;//to be set parameter
	if (p2pdistance(pt1, pt2) <= eps)
		return true;
	else
		return false;
}

bool on_other_noncollinearlines(vector<Vec4i> lines, Vec4i temp_line, Vec2i pt)
{
	for (size_t i = 0; i < lines.size(); i++)
	{
		Vec2i pt1 = { lines[i][0], lines[i][1] }; Vec2i pt2 = { lines[i][2], lines[i][3] };
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

void findLineEnds(vector<Vec2i> colpoints, vector<Vec2i>& lineEnds, vector<Vec2i>& temp_points, vector<Vec4i> lines)
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
		bool flag = on_other_noncollinearlines(lines, { pt3[0], pt3[1], pt4[0], pt4[1] }, tempPt1);
		for (vector<Vec2i>::iterator iter4 = temp_points.begin(); iter4 != temp_points.end();)
		{
			Vec2i tempPt2 = *iter4;
			if (tempPt2 == tempPt1 && (tempPt2 != pt3) && (tempPt2 != pt4) && (!flag))
			{
				iter4 = temp_points.erase(iter4);
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
void detect_line(Mat diagram_segwithoutcircle, Mat &color_img, vector<Vec4i> &lines, vector<Vec2i>& basicEndpoints)
{
	HoughLinesP(diagram_segwithoutcircle, lines, 1, CV_PI / 180, 30, 5, 10);
	vector<Vec2i> temp_points;
	// store all the initial points to a vector for handling
	for (size_t i = 0; i < lines.size(); ++i)
	{
		Vec2i pt1 = Vec2i(lines[i][0], lines[i][1]); Vec2i pt2 = Vec2i(lines[i][2], lines[i][3]);
		temp_points.push_back(pt1); temp_points.push_back(pt2);
	}
	// sort the points with a clear rule
	std::sort(temp_points.begin(), temp_points.end(), [](Vec2i a, Vec2i b)
	{
		if (a[0] != b[0])
			return (a[0] < b[0]);
		else
			return (a[1] < b[1]);
	});
	// define the vector to store the line endpoints along with the new lines
	vector<Vec2i> lineEnds; vector<Vec4i> newlines;
	int count = 0;
	// a loop through the line candidates, first extend the line with all the points on this line
	// and then find the line endpoints with the longest distance and generate the new lines
	for (vector<Vec4i>::iterator iter3 = lines.begin(); iter3 != lines.end(); iter3++)
	{
		Vec4i temp_line = *iter3;
		Vec2i pt1 = { temp_line[0], temp_line[1] }; Vec2i pt2 = { temp_line[2], temp_line[3] };

		vector<Vec2i>::iterator iter5, iter6;
		iter5 = find(temp_points.begin(), temp_points.end(), pt1);
		iter6 = find(temp_points.begin(), temp_points.end(), pt2);

		if (iter5 == temp_points.end() || iter6 == temp_points.end())
		{
			count++;
			continue;
		}
		else
		{
			vector<Vec2i> tempCollinearPoints;
			tempCollinearPoints.push_back(pt1); tempCollinearPoints.push_back(pt2);
			for (vector<Vec2i>::iterator iter4 = temp_points.begin(); iter4 != temp_points.end(); iter4++)
			{
				Vec2i temp_pt = *iter4;
				if (on_line(temp_line, temp_pt))
				{
					tempCollinearPoints.push_back(temp_pt);
				}
			}
			std::sort(tempCollinearPoints.begin(), tempCollinearPoints.end(), [](Vec2i a, Vec2i b)
			{
				if (a[0] != b[0])
					return (a[0] < b[0]);
				else
					return (a[1] < b[1]);
			});
			tempCollinearPoints.erase(unique(tempCollinearPoints.begin(), tempCollinearPoints.end()), tempCollinearPoints.end());
			findLineEnds(tempCollinearPoints, lineEnds, temp_points, lines);
			size_t tempSize = lineEnds.size();
			Vec4i newLine = { lineEnds[tempSize - 2][0], lineEnds[tempSize - 2][1], lineEnds[tempSize - 1][0], lineEnds[tempSize - 1][1] };
			newlines.push_back(newLine);
		}
	}	
	cout << "new line size: "<< newlines.size() << endl;
	std::sort(lineEnds.begin(), lineEnds.end(), [](Vec2i a, Vec2i b)
	{
		if (a[0] != b[0])
			return (a[0] < b[0]);
		else
			return (a[1] < b[1]);
	});
	lineEnds.erase(unique(lineEnds.begin(), lineEnds.end()), lineEnds.end());
	cout << "line end size now : " << lineEnds.size() << endl;
	for (size_t i = 0; i < newlines.size(); i++)
	{
		Vec2i pt_1 = { newlines[i][0], newlines[i][1] };
		Vec2i pt_2 = { newlines[i][2], newlines[i][3] };
		for (size_t j = i + 1; j < newlines.size(); j++)
		{
			Vec2i pt_3 = { newlines[j][0], newlines[j][1] };
			Vec2i pt_4 = { newlines[j][2], newlines[j][3] };
			if (same_pt(pt_1, pt_3) && (pt_3 != pt_1))
			{
				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_3, pt_1);
				erase_pointfrom_set(lineEnds, tbrendpt);
			}
			else if (same_pt(pt_3, pt_2) && (pt_3 != pt_2))
			{
				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_3, pt_2);
				erase_pointfrom_set(lineEnds, tbrendpt);
			}
			else if (same_pt(pt_4, pt_1) && (pt_4 != pt_1))
			{
				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_4, pt_1);
				erase_pointfrom_set(lineEnds, tbrendpt);
			}
			else if (same_pt(pt_4, pt_2) && (pt_4 != pt_2))
			{
				Vec2i tbrendpt = to_be_removed_similar_endpoint(pt_4, pt_2);
				erase_pointfrom_set(lineEnds, tbrendpt);
			}
		}
	}

	for (size_t i = 0; i < temp_points.size(); ++i)
	{
		bool flag = true;
		for (size_t j = 0; j < lineEnds.size(); ++j)
		{
			if (same_pt(temp_points[i], lineEnds[j]))
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			lineEnds.push_back(temp_points[i]);
		}
	}
	for (size_t i = 0; i < lineEnds.size(); i++)
	{
		Vec2i pt1 = lineEnds[i];
		Scalar temp_color = Scalar((rand() % 255), (rand() % 255), (rand() % 255));
		cv::circle(color_img, Point(pt1[0], pt1[1]), 1, temp_color, 10, 8, 0);
	}
	std::cout << "Now, the size of single points(not circle center) is: " << lineEnds.size() << endl;
	for (size_t i = 0; i < newlines.size(); i++)
	{
		Vec2i pt1 = { newlines[i][0], newlines[i][1] }; Vec2i pt2 = { newlines[i][2], newlines[i][3] };
		line(color_img, Point(pt1[0], pt1[1]), Point(pt2[0], pt2[1]), Scalar(0, 255, 0), 2, 8, 0);
	}
	namedWindow("points test");
	cv::imshow("points test", color_img);
	lines.clear();
	lines.assign(newlines.begin(), newlines.end());
	basicEndpoints.clear();
	basicEndpoints.assign(lineEnds.begin(), lineEnds.end());
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
void point_on_circle_line_check(vector<Vec2i> basicEndpoints, vector<Vec3f> circle_candidates, vector<circleX> &circles,
	vector<Vec4i> line_candidates, vector<lineX> &lines, vector<pointX> &points)
{
	bool flag = true;
	cout << basicEndpoints.size() << endl;
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
	cout << points.size() << endl;
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
void detect_line2(Mat diagram_segwithoutcircle, Mat &color_img, vector<Vec4i> &lines, vector<Vec2i>& temp_points)
{
	HoughLinesP(diagram_segwithoutcircle, lines, 1, CV_PI / 180, 30, 30, 10);
	//vector<Vec2i> temp_points;
	for (size_t i = 0; i < lines.size(); ++i)
	{
		Vec2i pt1 = { lines[i][0], lines[i][1] }; Vec2i pt2 = { lines[i][2], lines[i][3] };
		if (pt1[0] > 0 && pt1[1] > 0 && pt2[0] > 0 && pt2[1] > 0)
		{
			temp_points.push_back(pt1); temp_points.push_back(pt2);
		}

	}
	// sort the points with a clear rule
	std::sort(temp_points.begin(), temp_points.end(), [](Vec2i a, Vec2i b)
	{
		if (a[0] != b[0])
			return (a[0] < b[0]);
		else
			return (a[1] < b[1]);
	});
	cout << "test point" << endl;
	// handle the similar points
	for (vector<Vec2i>::iterator iter1 = temp_points.begin(); iter1 != temp_points.end(); iter1++)
	{
		for (vector<Vec2i>::iterator iter2 = iter1 + 1; iter2 != temp_points.end();)
		{
			if (same_pt(*iter1, *iter2))
			{

				//erase iter2 element
				if (iter2 != temp_points.end())
					iter2 = temp_points.erase(iter2);
				else
				{
					temp_points.pop_back();
				}
					
			}
			else
			{
				iter2++;
			}
		}
	}
	//handle the intersection and add to the set
	vector<Vec4i> newlines;
	for (int i = 0; i < temp_points.size(); ++i)
	{
		Vec2i pt1 = temp_points[i];
		for (int j = i + 1; j < temp_points.size(); ++j)
		{
			Vec2i pt2 = temp_points[j];
			bool flag = false;
			for (int k = 0; k < lines.size(); ++k)
			{
				Vec4i l = lines[k]; Vec2i pt3 = { l[0], l[1] }; Vec2i pt4 = { l[2], l[3] };
				if ((same_pt(pt1, pt3) && same_pt(pt2, pt4)) || (same_pt(pt1, pt4) && same_pt(pt2, pt3)))
				{
					flag = true;
					break;
				}
			}
			if (flag)
			{
				Vec4i newline = { pt1[0], pt1[1], pt2[0], pt2[1] };
				newlines.push_back(newline);
			}
		}
	}
	cout << newlines.size() << endl;
	lines.clear();
	lines.assign(newlines.begin(), newlines.end());
	cout << "line size after removing the line with the removed points: "<<lines.size() << endl;
	// this is meant to detect the crosses which is not the end point of line
	for (size_t i = 0; i < lines.size(); i++)
	{
		Vec4i line1 = lines[i]; Vec2i pt1 = { line1[0], line1[1] }; Vec2i pt2 = { line1[2], line1[3] };

		for (size_t j = i + 1; j < lines.size(); j++)
		{
			Vec4i line2 = lines[j]; Vec2i pt3 = { line2[0], line2[1] }; Vec2i pt4 = { line2[2], line2[3] };
			Vec2i cross_cal;
			bool have_line_cross = intersection(pt1, pt2, pt3, pt4, cross_cal);
			if (have_line_cross)
			{
				bool similar_flag = false;
				for (size_t k = 0; k < temp_points.size(); ++k)
				{
					if (same_pt(cross_cal, temp_points[k]))
					{
						similar_flag = true;
						break;
					}
				}
				if (!similar_flag)
					temp_points.push_back(cross_cal);
			}
		}
	}
	for (size_t i = 0; i < temp_points.size(); ++i)
	{
		cout << temp_points[i][0] << ", " << temp_points[i][1] << endl;
	}
	cout << temp_points.size() << endl;
	
	for (size_t i = 0; i < temp_points.size(); i++)
	{
		Vec2i pt1 = temp_points[i];
		Scalar temp_color = Scalar((rand() % 255), (rand() % 255), (rand() % 255));
		cv::circle(color_img, Point(pt1[0], pt1[1]), 1, temp_color, 5, 8, 0);
	}
	for (size_t i = 0; i < lines.size(); i++)
	{
		Vec2i pt1 = { lines[i][0], lines[i][1] }; Vec2i pt2 = { lines[i][2], lines[i][3] };
		line(color_img, Point(pt1[0], pt1[1]), Point(pt2[0], pt2[1]), Scalar(0, 255, 0), 1, 8, 0);
	}
	namedWindow("points test");
	cv::imshow("points test", color_img);
}
float evaluateLine(Mat diagram_segwithoutcircle, Vec4i rawLine)
{
	// compute the points number which is on the line
	float maxDist = 1.0f; int count = 0;
	vector<Point2f> edgePositions;
	edgePositions = getPointPositions(diagram_segwithoutcircle);
	Vec2i pt1 = { rawLine[0], rawLine[1] }; Vec2i pt2 = { rawLine[2], rawLine[3] };
	int smaller_x = (pt1[0] < pt2[0]) ? pt1[0] : pt2[0]; 
	int bigger_x = (pt1[0] > pt2[0]) ? pt1[0] : pt2[0]; 
	int smaller_y = (pt1[1] < pt2[1]) ? pt1[1] : pt2[1];
	int bigger_y = (pt1[1] > pt2[1]) ? pt1[1] : pt2[1];
	for (int i = smaller_x; i < bigger_x; ++i)
	{
		for (int j = smaller_y; j < bigger_y; ++j)
		{
			Point2f testPoint = CvPoint(i, j);
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
void detect_line3(Mat diagram_segwithoutcircle, Mat &color_img, vector<Vec4i> &lines, vector<Vec2i>& temp_points)
{
	vector<Vec4i> rawLines;
	HoughLinesP(diagram_segwithoutcircle, rawLines, 1, CV_PI / 180, 15, 15, 10);
	for (size_t i = 0; i < rawLines.size(); ++i)
	{
		Vec4i rawLine = rawLines[i];
		float lineEV = evaluateLine(diagram_segwithoutcircle, rawLine);
		if (lineEV > 15)
		{
			lines.push_back(rawLine);
			line(diagram_segwithoutcircle, Point(rawLine[0], rawLine[1]), Point(rawLine[2], rawLine[3]), Scalar(0, 0, 0), 6, 8);
		}
	}
}
void primitive_parse(Mat diagram_segment, vector<pointX> &points, vector<lineX> &lines, vector<circleX> &circles)
{
	// primitive about points, lines, and circles
	//first detect the circle, we can get the matrix of diagram segments without circle for the next
	// processing step
	vector<Vec3f> circle_candidates; Mat color_img; cvtColor(diagram_segment, color_img, CV_GRAY2RGB);
	Mat diagram_segwithoutcircle;

	/*ransac go*/
	detect_circle(diagram_segment, color_img, diagram_segwithoutcircle,circle_candidates);

	// then the line detection
	vector<Vec4i> line_candidates; vector<Vec2i> basicEndpoints;
	detect_line3(diagram_segwithoutcircle, color_img, line_candidates, basicEndpoints);
	//detect_line2(diagram_segwithoutcircle, color_img, line_candidates, basicEndpoints);
	cout << "basic endpoints num: "<<basicEndpoints.size() << endl;
	
	// for points check if they are in circles
	point_on_circle_line_check(basicEndpoints, circle_candidates, circles, line_candidates, lines, points);
	for (int i = 0; i < points.size(); ++i)
	{
		pointX point = points[i];
		cout << "point " << point.p_idx << " ("<<point.px<<", "<< point.py<<") " << " on circle ";
		for (int j = 0; j < point.c_idx.size(); ++j)
		{
			cout << point.c_idx[j] << ", ";
		}
		cout << "on line ";
		for (int k = 0; k < point.l_idx.size(); ++k)
		{
			cout << point.l_idx[k] << ", ";
		}
		cout << endl;
	}
}

int test_diagram()
{
	//first load a image
	Mat image = imread("Sg-2.jpg", 0);
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