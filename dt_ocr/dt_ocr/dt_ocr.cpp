// dt_ocr.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "ocr_method.h"
#include <core/mat.hpp>
#include <imgcodecs.hpp>
#include <fstream>
#include <iostream>
using namespace std;
int _tmain(int argc, _TCHAR* argv[])
{
	// the image path
	char* base_path = "D:\\data\\graph-DB\\nt7\\charImgs";
	char gtFilePath[100]; sprintf_s(gtFilePath, "%s\\mess_labels.txt", base_path);
	std::fstream gtFile(gtFilePath); fstream compareFile("D:\\data\\graph-DB\\nt7\\charImgs\\compare.txt");
	string gtStr; gtFile >> gtStr; 
//	compareFile << "test begin" << endl;
	char ocrCmd3[100];
	int all_num, checkout_num, right_num;
	all_num = checkout_num = right_num = 0;
	for (auto i = 0; i < 527; ++i)
	{
		all_num++;
		char gtChar = gtStr[i];
		if (gtChar != '`')
		{
			checkout_num++;
			sprintf_s(ocrCmd3, "gocr049.exe -i %s\\charImg-%d.pnm -o tmpResult.txt", base_path, i);
			system(ocrCmd3);
			cout << i<<"  ocring" << endl;
			std::fstream singleCharFile;
			singleCharFile.open("tmpResult.txt", std::ios::in);
			char singleResult;
			singleCharFile >> singleResult;
			singleCharFile.close();
			compareFile << singleResult << " vs " << gtChar << endl;
			if (gtChar == singleResult || abs(gtChar-singleResult)==32)
			{
				right_num++;
			}
		}
	}
	compareFile.close();
	cout << right_num << " : " << checkout_num << " : " << all_num << endl;
	return 0;
}

