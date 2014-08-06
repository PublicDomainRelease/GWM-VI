@echo off
REM Run Supply2 test with GWM-VI 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial


echo Running GWM-VI.exe with supply2_pll.gwm
cd ..
call StartRunners 2
cd supply2
..\..\Bin\GWM-VI.exe ..\data\supply2_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with supply2_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\supply2_serial.gwm
goto end



:end
pause
