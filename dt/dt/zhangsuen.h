#ifndef __ZHANGSUEN_H__
#define __ZHANGSUEN_H__

#include "stdafx.h"


typedef std::pair<int, int> zPoint;
typedef unsigned char uchar_t;

int num_one_pixel_neighbours(const cv::Mat& image, const zPoint& zPoint);

int num_zero_pixel_neighbours(const cv::Mat& image, const zPoint& zPoint);

int connectivity(const cv::Mat& image, const zPoint& zPoint);

int yokoi_connectivity(const cv::Mat& image, const zPoint& zPoint);

void delete_pixels(const cv::Mat& image, const std::set<zPoint>& zPoints);

void remove_staircases(cv::Mat& image);

void zhangsuen_thin(cv::Mat& img);

void thin(cv::Mat& img, bool need_boundary_smoothing, 
          bool need_acute_angle_emphasis, bool destair);

void boundary_smooth(cv::Mat& image);

void acute_angle_emphasis(cv::Mat& image);

bool match(const cv::Mat& image, const std::vector<zPoint>& zPoints,
           const std::vector<uchar_t>& values);

bool match_templates(const cv::Mat& image, const zPoint& zPoint, int k);

#endif
