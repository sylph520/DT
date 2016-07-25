// class and data structure definition
struct pointX
{
	bool cflag;// wheather it's  circle center or not
	int p_idx; // the sequence, if -1 , it's a circle center.
	int px, py;
	vector<int> l_idx;
	vector<int> c_idx;
	string label = "";
};
struct lineX
{
	int l_idx; 
	int px1, py1, px2, py2;
	double length;
	string label = "";
};
struct circleX
{
	int c_idx; 
	float cx, cy;
	float radius;
	string label = "";
};
struct distance_info
{
	Vec2i pt1; Vec2i pt2;// the points to calculate on
	double distance;//computed distance results with pt1 and pt2
};
//functions
int image_parse(Mat image);
int test_diagram();
int diagram();