@echo off
REM Run Maximin test with GWM-VI

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial

echo Running GWM-VI.exe with maximin_pll.gwm
cd ..
call StartRunners 2
cd maximin
..\..\Bin\GWM-VI.exe ..\data\maximin_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with maximin_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\maximin_serial.gwm
goto end


:end
pause
