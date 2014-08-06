@echo off
REM  Run Seawater test with GWM-VI (parallel or 
REM  serial) 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial


echo Running GWM-VI.exe with seawater_pll.gwm
cd ..
call StartRunners 2
cd seawater
..\..\Bin\GWM-VI.exe ..\data\seawater_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with seawater_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\seawater_serial.gwm
goto end


:end
pause
