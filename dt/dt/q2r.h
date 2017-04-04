#pragma once
#include "stdafx.h"
#include "diagram.h"

bool same_len_line(line_class a, line_class b, vector<point_class> &pcs);

bool _same_len_line(Vec2i pt, Vec2i pt1, Vec2i pt2, vector<point_class> &pcs);

bool same_len_line(point_class a, point_class b, point_class c, vector<point_class> &pcs);

bool same_ang_lines(line_class a, line_class b, line_class c, line_class d, vector<point_class> &pcs);

bool pt_line_mid(line_class l, point_class p, vector<point_class> &pcs);

bool isPerpendicular();