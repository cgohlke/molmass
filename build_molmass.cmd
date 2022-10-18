setlocal

rd /S /Q molmass.egg-info

rd /Q /S build
call py.exe setup.py sdist
start X:\Python310\Scripts\restview.exe README.rst
pause

rd /Q /S build
call py.exe setup.py bdist_wheel --python-tag=py3

echo Running Tests
call py38x32 -m molmass --test
call py38x32 -m molmass --doctest

call py311 -X dev -m molmass --test
call py311 -X dev -m molmass --doctest
::call py -m pytest molmass --doctest-modules

call py311 -X dev molmass\elements.py > _tmp.txt

call py -X dev -m molmass --web --url=http://127.0.0.1:5005

endlocal
pause
