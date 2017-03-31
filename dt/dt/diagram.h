// class and data structure definition

#pragma once

double p2pdistance(Vec2i pt1, Vec2i pt2);


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
	
	~point_class()
	{
		l_idxs.clear();
		c_idxs.clear();
		con_pidxs.clear();
	}
};

Vec2i get_line_ptVec_by_id(vector<point_class>& pointxs, int id);
point_class* get_line_pt_by_id(vector<point_class> & pointxs, int id);
Vec4i pt2line(Vec2i pt1, Vec2i pt2);

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
		Vec2i pt1, pt2;
		pt1 = get_line_ptVec_by_id(pointxs, p1_id); pt2 = get_line_ptVec_by_id(pointxs, p2_id);
		Vec4i lineVec = pt2line(pt1, pt2);
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

	int getLid() const
	{
		return l_id;
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
		return get_line_ptVec_by_id(pointxs, p1_id);
	}

	Vec2i getpt2vec(vector<point_class> &pointxs) const
	{
		return get_line_ptVec_by_id(pointxs, p2_id);
	}

	double getLen(vector<point_class> &pointxs) const
	{
		double length = p2pdistance(pointxs[p1_id].getXY(), pointxs[p2_id].getXY());
		return length;
	}

	void setPt1_vec(vector<point_class> &pointxs,Vec2i pv) const
	{
		point_class *p = get_line_pt_by_id(pointxs, p1_id);
		p->setXY(pv);
	}

	void setPt2_vec(vector<point_class> &pointxs, Vec2i pv) const
	{
		point_class *p = get_line_pt_by_id(pointxs, p2_id);
		p->setXY(pv);
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

	~line_class()
	{
		ptcs.clear();
	}
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