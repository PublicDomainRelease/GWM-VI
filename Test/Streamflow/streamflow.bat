@echo off
REM Run Streamflow test with GWM-VI 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial


echo Running GWM-VI.exe with seawaterdryQ_pll.gwm
cd ..
call StartRunners 2
cd streamflow
..\..\Bin\GWM-VI.exe ..\data\streamflow_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with streamflow_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\streamflow_serial.gwm
goto end


:end
pause
