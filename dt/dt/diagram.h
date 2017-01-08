// class and data structure definition
double p2pdistance(Vec2i pt1, Vec2i pt2);
//struct pointX
//{
//	//bool cflag;// wheather it's  circle center or not
//	int p_idx; // the sequence, if -1 , it's a circle center.
//	int px, py; Vec2i pxy;
//	vector<int> l_idxs;
//	vector<int> c_idxs;
//	string label = "";
//};
//struct lineX
//{
//	int l_idx;
//	
////	int px1, py1, px2, py2;
////	Vec2i pt1, pt2; 
//	//Vec4i lxy;
//	//int p_idx1, p_idx2;
//	point_class *p1, *p2;
//	
//	//double length;
//	string label = "";
//};
//struct circleX
//{
//	int c_idx; int center_pid;
//	Vec3i Circle; Vec2i center; int radius; int contour_width;
//	int cx, cy;
//
//	vector<int> p_idxs;
//	string label = "";
//};

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
	vector<int> l_idxs, c_idxs; vector<int> con_pidxs;
public:
	point_class(Vec2i _pxy = { 0, 0}, int _p_idx = -1, vector<int> _l_idxs = {}, string _label = "", vector<int> _c_idxs = {}, vector<int> _con_pidxs = {})
	{
		pxy = _pxy; p_idx = _p_idx; label = _label;
		l_idxs = _l_idxs; c_idxs = _c_idxs; con_pidxs = _con_pidxs;
	}
	point_class(const point_class &b)
	{
		pxy = b.pxy; p_idx = b.p_idx; label = b.label;
		l_idxs = b.l_idxs; c_idxs = b.c_idxs; con_pidxs = b.con_pidxs;
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
	vector<int> getConidxs() const
	{
		return con_pidxs;
	}
	string get_label() const
	{
		return label;
	}
	void setPid(int a)
	{
		p_idx = a;
	}
	void setXY(Vec2i x)
	{
		pxy = x;
	}
	void pushLid(int a)
	{
		l_idxs.push_back(a);
	}
	void pushCid(int a)
	{
		c_idxs.push_back(a);
	}
	void pushConid(int a)
	{
		con_pidxs.push_back(a);
	}
	
	~point_class(){}
};

class line_class
{
private:
	int l_id; string label;
	int p1_id, p2_id;
	//point_class p1, p2;
	vector<point_class> ptcs;
public:
	line_class(int _p1_id, int _p2_id, int _l_id = -1, string _label = "", vector<point_class> _ptcs={})
	{
		l_id = _l_id; label = _label;
		p1_id = _p1_id; p2_id = _p2_id;
		ptcs = _ptcs;
	}
	line_class(const line_class &b)
	{
		l_id = b.l_id; label = b.label;
		p1_id = b.p1_id; p2_id = b.p2_id;
		ptcs = b.ptcs;
	}
	Vec4i getLineVec(vector<point_class> &pointxs) const
	{
		Vec4i lineVec = { pointxs[p1_id].getX(), pointxs[p1_id].getY(), pointxs[p2_id].getX(), pointxs[p2_id].getY() };
		return lineVec;
	}
	int getPt1Id() const
	{
		return p1_id;
	}
	int getPt2Id() const
	{
		return p2_id;
	}
	int getPtId(int pos) const
	{
		if (pos == 1)
			return p1_id;
		else
			return p2_id;
	}
	/*point_class getPt1(vector<point_class> &pointxs) const
	{
		return  pointxs[p1_id];
	}
	point_class getPt2(vector<point_class> &pointxs) const
	{
		return  pointxs[p2_id];
	}*/
	Vec2i getpt1vec(vector<point_class> &pointxs) const
	{
		return pointxs[p1_id].getXY();
	}
	Vec2i getpt2vec(vector<point_class> &pointxs) const
	{
		return pointxs[p2_id].getXY();
	}
	double getLen(vector<point_class> &pointxs) const
	{
		double length = p2pdistance(pointxs[p1_id].getXY(), pointxs[p2_id].getXY());
		return length;
	}
	void setPt1_vec(vector<point_class> &pointxs,Vec2i pv) const
	{
		pointxs[p1_id].setXY(pv);
	}
	void setPt2_vec(vector<point_class> &pointxs, Vec2i pv) const
	{
		pointxs[p2_id].setXY(pv);
	}
	void setPt_vec(vector<point_class> &pointxs, Vec2i pv, int ind) const
	{
		if (ind==1)
			pointxs[p1_id].setXY(pv);
		else
			pointxs[p2_id].setXY(pv);

	}
	/*void setPt1(vector<point_class> pointxs, int p) 
	{
		p1_id = p;
	}
	void setPt2(vector<point_class> pointxs, int p) 
	{
		p2_id = p;
	}*/
	void setpt1Id(int p)
	{
		p1_id = p;
	}
	void setpt2Id(int p)
	{
		p2_id = p;
	}
	void setptId(int p, int ind)
	{
		if (ind == 1)
			p1_id = p;
		else
			p2_id = p;
	}
	~line_class(){}
};

class circle_class
{
private:
	int c_idx; Vec2i center; int radius; 
	int contour_width;
	vector<point_class> p_idxs;
public:
	circle_class(Vec2i _center, int _radius, int _c_idx = -1, int _contour_width = 0, vector<point_class> _p_idxs = {})
	{
		c_idx = _c_idx; center = _center; radius = _radius;
		contour_width = _contour_width; 
		p_idxs = _p_idxs;
	}

	circle_class(Vec3i c, int _c_idx = -1, int _contour_width = 0, vector<point_class> _p_idxs = {})
	{
		c_idx = _c_idx; center = { c[0], c[1] }; radius = c[2];
		contour_width = _contour_width;
		p_idxs = _p_idxs;
		
	}

	circle_class(const circle_class &b)
	{
		c_idx = b.c_idx; center = b.center; radius = b.radius;
		contour_width = b.contour_width;
		p_idxs = b.p_idxs;
	}
	Vec3i getCircleVec()
	{
		Vec3i circleVec = { center[0], center[1], radius };
		return circleVec;
	}
	Vec2i getCenter() const
	{
		return center;
	}
	int getRadius() const
	{
		return radius;
	}
	~circle_class(){}
};
int image_parse(Mat image);
int test_diagram();
int diagram();