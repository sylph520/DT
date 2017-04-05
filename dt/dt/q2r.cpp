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


string pOnCircle2string(point_class p, circle_class c, map<string,int> vars)
{
	string p_label = p.getLabel();
	Vec2i center = c.getCenter(); double radius = c.getRadius();
	int var_id = vars[p_label];
	string equation = "(x" + to_string(var_id) + " - " + to_string(center[0]) + ")^2 + " 
		+ "(y" + to_string(var_id) + " - " + to_string(center[1]) + ")^2 - " + to_string(c.getRadius()) + "^2 = 0";
	return equation;
}

string pInLine2string(point_class p,line_class l, map<string, int> vars, vector<point_class> &points)
{
	int pid1, pid2; pid1 = l.getPt1Id(); pid2 = l.getPt2Id();
	string p_label1, p_label2; string p_label = p.getLabel();
	p_label1 = points[pid1].getLabel(); p_label2 = points[pid2].getLabel();
	int var_id, var_id1, var_id2; var_id = vars[p_label]; var_id1 = vars[p_label1]; var_id2 = vars[p_label2];
	string equation = "(x" + to_string(var_id2) + " - x" + to_string(var_id1) +
		")*(y" + to_string(var_id) + " - y" + to_string(var_id1) + ") - (y" +
		to_string(var_id2) + " - y" + to_string(var_id1) + ")*(x" + to_string(var_id)
		+ " - x" + to_string(var_id1) + " = 0";
	return equation;
}

string parallelL2S(line_class l1, line_class l2, map<string, int> vars, vector<point_class> &points)
{
	int var_id1 = vars[points[l1.getPt1Id()].getLabel()];
	int var_id2 = vars[points[l1.getPt2Id()].getLabel()];
	int var_id3 = vars[points[l2.getPt1Id()].getLabel()];
	int var_id4 = vars[points[l2.getPt2Id()].getLabel()];
	string equation = "(y" + to_string(var_id4) + " - y" + to_string(var_id3) + ")*(x"
		+ to_string(var_id2) + " - x" + to_string(var_id1) + ") - (y" + to_string(var_id2) + " - y"
		+ to_string(var_id1) + ") * (x" + to_string(var_id4) + " - x" + to_string(var_id3) + ") = 0";
	return equation;

}

string perperdicularL2S(line_class l1, line_class l2, map<string, int> vars, vector<point_class> &points)
{
	int var_id1 = vars[points[l1.getPt1Id()].getLabel()];
	int var_id2 = vars[points[l1.getPt2Id()].getLabel()];
	int var_id3 = vars[points[l2.getPt1Id()].getLabel()];
	int var_id4 = vars[points[l2.getPt2Id()].getLabel()];
	string equation = "(y" + to_string(var_id4) + " - y" + to_string(var_id3) + ")*(y"
		+ to_string(var_id2) + " - y" + to_string(var_id1) + ") + (x" + to_string(var_id2) + " - x"
		+ to_string(var_id1) + ") * (x" + to_string(var_id4) + " - x" + to_string(var_id3) + ") = 0";
	return equation;
}

string sameLenL2S(line_class l1, line_class l2, map<string, int> vars, vector<point_class> &points)
{
	int var_id1 = vars[points[l1.getPt1Id()].getLabel()];
	int var_id2 = vars[points[l1.getPt2Id()].getLabel()];
	int var_id3 = vars[points[l2.getPt1Id()].getLabel()];
	int var_id4 = vars[points[l2.getPt2Id()].getLabel()];
	string equation = "(y" + to_string(var_id4) + " - y" + to_string(var_id3) + ")^2 + (x"
		+ to_string(var_id4) + " - x" + to_string(var_id3) + ")^2 - (y" + to_string(var_id2) + " - y"
		+ to_string(var_id1) + ")^2 - (x" + to_string(var_id2) + " - x" + to_string(var_id1) + ")^2 = 0";
	return equation;
}
string sameLenL2S(int pid, line_class l, map<string, int> vars, vector<point_class> &points)
{
	int var_id = vars[points[pid].getLabel()];
	int var_id1 = vars[points[l.getPt1Id()].getLabel()];
	int var_id2 = vars[points[l.getPt2Id()].getLabel()];
	string equation = "(y" + to_string(var_id2) + " - y" + to_string(var_id) + ")^2 + (x"
		+ to_string(var_id2) + " - x" + to_string(var_id) + ")^2 - (y" + to_string(var_id1) + " - y"
		+ to_string(var_id) + ")^2 - (x" + to_string(var_id1) + " - x" + to_string(var_id) + ")^2 = 0";
	return equation;
}

string lSlope2str(line_class l, map<string, int> vars, vector<point_class> &points)
{
	int var_id1 = vars[points[l.getPt1Id()].getLabel()];
	int var_id2 = vars[points[l.getPt2Id()].getLabel()];
	string equation = "atan2( (y" + to_string(var_id2) + " - y" + to_string(var_id1) + ")/(x"
		+ to_string(var_id2) + " -x" + to_string(var_id1) + ") ) - " + to_string(l.getSlope_r()) + " = 0";
	return equation;
}

string lLen2str(line_class l, map<string, int> vars, vector<point_class> &points)
{
	int var_id1 = vars[points[l.getPt1Id()].getLabel()];
	int var_id2 = vars[points[l.getPt2Id()].getLabel()];
	string equation = "(y" + to_string(var_id2) + " - y" + to_string(var_id1) + ")^2" +
		"(x" + to_string(var_id2) + " -x" + to_string(var_id1) + ")^2 - (" + to_string(l.getLen(points)) + ")^2 = 0";
	return equation;
}