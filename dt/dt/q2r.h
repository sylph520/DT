#pragma once
#include "stdafx.h"
#include "diagram.h"

bool same_len_line(line_class a, line_class b, vector<point_class> &pcs);

bool _same_len_line(Vec2i pt, Vec2i pt1, Vec2i pt2, vector<point_class> &pcs);

bool same_len_line(point_class a, point_class b, point_class c, vector<point_class> &pcs);

bool same_ang_lines(line_class a, line_class b, line_class c, line_class d, vector<point_class> &pcs);

bool pt_line_mid(line_class l, point_class p, vector<point_class> &pcs);


//equation part
string pOnCircle2string(point_class p, circle_class c, map<string, int> vars);

string pInLine2string(point_class p, line_class l, map<string, int> vars, vector<point_class> &points);

string parallelL2S(line_class l1, line_class l2, map<string, int> vars , vector<point_class> &points);

string perperdicularL2S(line_class l1, line_class l2, map<string, int> vars, vector<point_class> &points);

string sameLenL2S(line_class l1, line_class l2, map<string, int> vars, vector<point_class> &points);

string sameLenL2S(int pid, line_class l, map<string, int> vars, vector<point_class> &points);

string lSlope2str(line_class l, map<string, int> vars, vector<point_class> &points);

string lLen2str(line_class l, map<string, int> vars, vector<point_class> &points);