// class and data structure definition
struct pointX
{
	bool cflag;// wheather it's  circle center or not
	int p_idx; // the sequence, if -1 , it's a circle center.
	int px, py; Vec2i pxy;
	vector<int> l_idx;
	vector<int> c_idx;
	string label = "";
};
struct lineX
{
	int l_idx;
	
	int px1, py1, px2, py2;
	Vec2i pt1, pt2; Vec4i lxy;
	int pidx1, pidx2;
	
	double length;
	string label = "";
};
struct circleX
{
	int c_idx; int center_pid;
	
	Vec3i Circle; Vec2i center; int radius; int contour_width;
	int cx, cy;

	vector<int> p_idxs;
	string label = "";
};
struct distance_info
{
	Vec2i pt1; Vec2i pt2;// the points to calculate on
	double distance;//computed distance results with pt1 and pt2
};
class point_class
{
public:
	point_class(Vec2i _pxy, int _p_idx = -1, string _label = "", vector<int> _l_idxs = {}, vector<int> _c_idxs = {})
	{
		pxy = _pxy; 
		px = _pxy[0]; py = _pxy[1];
		p_idx = _p_idx;
		label = _label; l_idxs = _l_idxs; c_idxs = _c_idxs;
	}
private:
	int p_idx; Vec2i pxy; string label; vector<int> l_idxs; vector<int> c_idxs;
	int px, py;
};
class line_class
{
public:
	line_class(Vec4i _lxy, int _l_idx = -1, int _p_idx1 = -1, int _p_idx2 = -1, double _length = 0.0, string _label = "")
	{
		lxy = _lxy; pt1 = { _lxy[0], _lxy[1] }; pt2 = { _lxy[2], _lxy[3] }; px1 = _lxy[0]; py1 = _lxy[1]; px2 = _lxy[2]; py2 = _lxy[3];
		l_idx = _l_idx; p_idx1 = _p_idx1; p_idx2 = _p_idx2;
		length = _length;
		label = _label;
	}
private:
	Vec4i lxy; Vec2i pt1, pt2; int px1, py1, px2, py2;
	int l_idx, p_idx1, p_idx2;
	double length;
	string label;
};
class circle_class
{
public:
	circle_class(Vec3i _c, int _c_idx = -1, int _center_pidx = -1, int _contour_width = 0, vector<int> _p_idxs = {})
	{
		Circle = _c; Center = { _c[0], _c[1] }; px = _c[0]; py = _c[1]; radius = _c[2];
		c_idx = _c_idx; center_pidx = _center_pidx; contour_width = _contour_width;
		p_idxs = _p_idxs;
	}
private:
	Vec3i Circle; Vec2i Center;
	int radius, px, py;
	int c_idx, center_pidx, contour_width;
	vector<int> p_idxs;
	
};
//functions
int image_parse(Mat image);
int test_diagram();
int diagram();