@echo off
REM Run drain_volume test with GWM-VI 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial

echo Running GWM-VI.exe with drain_volume_pll.gwm
cd ..
call StartRunners 2
cd drain
..\..\Bin\GWM-VI.exe ..\data\drain_volume_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with drain_volume_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\drain_volume_serial.gwm
goto end



:end
pause
