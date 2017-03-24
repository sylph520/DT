// createDirtest.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Windows.h"
#include <iostream>
using namespace std;

int _tmain(int argc, _TCHAR* argv[])
{
	char test[100] = "D:\\data\\graph-DB\\nt6\\charImgs\\2";
//	char cmd[100]; sprintf_s(cmd, "mkdir %s", test);
//	WinExec(cmd, WM_SHOWWINDOW);
//	system(cmd);
	cout<< CreateDirectory((LPCTSTR)(test), NULL);
	WaitMessage();
	return 0;
}

