@echo off
call :normalizepath %cd%\..\Miniconda3
set mc3=%retval%
echo Miniconda will be installed in %mc3%


rem remove Miniconda3 directory if it already exists
if exist %mc3% (
  echo Existing %mc3% directory will be deleted.
  rmdir %mc3% /s /q
)


rem install Miniconda
echo Installing Miniconda to %mc3%
cd ..\software
Miniconda3-latest-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%mc3%
cd ..\installation
rem start /wait "" /D ..\software Miniconda3-latest-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /AddToPath=0 /S /D=%cd%\..\Miniconda3
if not exist %mc3% (
  echo Miniconda could not be installed.
  exit /B
)

echo on
rem verify the miniconda install
%mc3%\python -V


rem pause so screen will not go away
pause


:: ========== FUNCTIONS ==========
exit /B

:normalizepath
  set retval=%~dpfn1
  exit /B
