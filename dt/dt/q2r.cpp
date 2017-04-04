#include "stdafx.h"
#include "diagram.h"

bool same_len_line(line_class a, line_class b, vector<point_class> &pcs)
{
	double threshold = 3;
	if (abs(a.getLen(pcs) - b.getLen(pcs)) < threshold)
		return true;
	else
		return false;
}
bool _same_len_line(Vec2i pt, Vec2i pt1, Vec2i pt2, vector<point_class> &pcs)
{
	double dis1, dis2; 
	dis1 = p2pdistance(pt, pt1);
	dis2 = p2pdistance(pt, pt2);
	if (abs(dis1 - dis2) < 3)
		return true;
	else
		return false;
}
bool same_len_line(point_class a,point_class b, point_class c,vector<point_class> &pcs)
{
	//this is for two line with a shared point a
	Vec2i axy, bxy, cxy;
	axy = a.getXY(); bxy = b.getXY(); cxy = c.getXY();
	return _same_len_line(axy, bxy, cxy, pcs);
	

}
bool same_ang_lines(line_class a, line_class b,line_class c, line_class d, vector<point_class> &pcs)
{
	double ang1, ang2;
	ang1 = llAngle(a.getLineVec(pcs), b.getLineVec(pcs));
	ang2 = llAngle(c.getLineVec(pcs), d.getLineVec(pcs));

	if (abs(ang1 - ang2) < 3)
		return true;
	else
		return false;
}

bool pt_line_mid(line_class l,point_class p,vector<point_class> &pcs)
{
	Vec2i pt1, pt2; line2pt(l.getLineVec(pcs), pt1, pt2);
	if (_same_len_line(p.getXY(), pt1, pt2, pcs))
	{
		//same dis  to the end point
		if (on_line(l.getLineVec(pcs), p.getXY()))
			return true;
		else
			return false;
	}
	else
		return false;
}

//bool line_tan_cirlce(line_class l, circle_class c,vector<point_class> pcs)
//{
//	
//}