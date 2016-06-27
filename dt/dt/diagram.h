// class and data structure definition
struct pointX
{
	int p_idx; string label;
	int px, py;
	vector<int> l_idx;
	vector<int> c_idx;
};
struct lineX
{
	int l_idx; string label;
	int p_idx1, p_idx2;
	double length;
};
struct circleX
{
	int c_idx; string label;
	int center_idx;
	double radius;
};
//functions
int image_parse(Mat image);
int test_diagram();