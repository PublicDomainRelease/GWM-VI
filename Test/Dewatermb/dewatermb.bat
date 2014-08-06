@echo off
REM  Run Dewatermb test with GWM-VI 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial

echo Running GWM-VI.exe with dewatermb_pll.gwm
cd ..
call StartRunners 2
cd dewatermb
..\..\Bin\GWM-VI.exe ..\data\dewatermb_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with dewatermb_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\dewatermb_serial.gwm
goto end



:end
pause
