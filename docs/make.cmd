setlocal
set PYTHONPATH=..
rd /S /Q html
::cd ..
::X:\Python310\Scripts\sphinx-build -b spelling . ./build
X:\Python310\Scripts\sphinx-build -b html . ./html
endlocal
start html/index.html
pause