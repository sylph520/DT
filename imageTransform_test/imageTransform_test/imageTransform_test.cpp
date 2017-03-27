// imageTransform_test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <opencv.hpp>

using  namespace  std;
using namespace  cv;
int _tmain(int argc, _TCHAR* argv[])
{
	int count = 0;
	for (auto i = 1; i <= 12; i++)
	{
		for (auto j = 1; j <= 50;j++)
		{
			cout <<i<<"   "<< j << endl;

			char in_filePath[256],out_filePath[256],in_fileName[256],out_fileName[256],in_fileUNC[256],out_fileUNC[256];
			char outName[256];
			sprintf_s(in_filePath,"cated/%d", i);
			sprintf_s(in_fileName, "%s/%d", in_filePath,j);
			sprintf_s(in_fileUNC, "%s.png",in_fileName);
			sprintf_s(out_filePath, "cated_tiff/%d", i);
			sprintf_s(out_fileName, "%s/%d", out_filePath, j);
			sprintf_s(out_fileUNC, "%s.png", out_fileName);

			Mat tmp = imread(in_fileUNC);
			imwrite(out_fileUNC, tmp);
			sprintf_s(outName, "%d.tiff", count++); 
			imwrite(outName, tmp);

		}
	}
	return 0;
}

