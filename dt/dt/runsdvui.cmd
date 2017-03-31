cd /d "D:\Document\Visual Studio 2013\Projects\DT\dt\dt" &msbuild "dt.vcxproj" /t:sdvViewer /p:configuration="Debug" /p:platform=Win32
exit %errorlevel% 