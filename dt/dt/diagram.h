// class and data structure definition
double p2pdistance(Vec2i pt1, Vec2i pt2);
struct pointX
{
	//bool cflag;// wheather it's  circle center or not
	int p_idx; // the sequence, if -1 , it's a circle center.
	int px, py; Vec2i pxy;
	vector<int> l_idxs;
	vector<int> c_idxs;
	string label = "";
};
struct lineX
{
	int l_idx;
	
//	int px1, py1, px2, py2;
//	Vec2i pt1, pt2; 
	//Vec4i lxy;
	//int p_idx1, p_idx2;
	pointX *p1, *p2;
	
	//double length;
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
/*class definition*/
class point_class
{
private:
	Vec2i pxy; int p_idx; string label;
	vector<int> l_idxs, c_idxs;
public:
	point_class(Vec2i _pxy, int _p_idx = -1, string _label = "", vector<int> _l_idxs = {}, vector<int> _c_idxs = {})
	{
		pxy = _pxy; p_idx = _p_idx; label = _label;
		l_idxs = _l_idxs; c_idxs = _c_idxs;
	}
	point_class(const point_class &b)
	{
		pxy = b.pxy; p_idx = b.p_idx; label = b.label;
		l_idxs = b.l_idxs; c_idxs = b.c_idxs;
	}
	int getPid() const
	{
		return p_idx;
	}
	int getX() const
	{
		return pxy[0];
	}
	int getY() const
	{
		return pxy[1];
	}
	Vec2i getXY() const
	{
		return pxy;
	}
	vector<int> getLidxs() const
	{
		return l_idxs;
	}
	vector<int> getCidxs() const
	{
		return c_idxs;
	}
	string get_label() const
	{
		return label;
	}
};


class line_class
{
private:
	int l_id; string label;
	point_class *p1, *p2;
public:
	line_class(point_class *_p1, point_class *_p2, int _l_id = -1, string _label = "")
	{
		l_id = _l_id; label = _label;
		p1 = _p1; p2 = _p2;
	}
	line_class(const line_class &b)
	{
		l_id = b.l_id; label = b.label;
		p1 = b.p1; p2 = b.p2;
	}
	Vec4i getPlainLineVec() const
	{
		Vec4i plainLineVec = { (*p1).getX(), (*p1).getY(), (*p2).getX(), (*p2).getY() };
		return plainLineVec;
	}
};

int image_parse(Mat image);
int test_diagram();
int diagram();