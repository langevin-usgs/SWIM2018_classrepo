rem install Miniconda
rmdir %cd%\..\Miniconda3 /s /q
..\software\Miniconda3-latest-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%cd%\..\Miniconda3
..\Miniconda3\python -V

rem install python packages
..\Miniconda3\python install_packages.py

rem pause to see if any errors
pause
