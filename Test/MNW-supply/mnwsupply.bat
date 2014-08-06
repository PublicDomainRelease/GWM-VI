@echo off
REM  Run mnwsupply test with GWM-VI (parallel or serial) 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial


echo Running GWM-VI.exe with mnwsupply_pll.gwm
cd ..
call StartRunners 2
cd mnw-supply
..\..\Bin\GWM-VI.exe ..\data\mnwsupply_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with mnwsupply_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\mnwsupply_serial.gwm
goto end



:end
pause
