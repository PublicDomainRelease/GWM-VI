@echo off
REM Run drain_rate test with GWM-VI 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial

echo Running GWM-VI.exe with drain_rate_pll.gwm
cd ..
call StartRunners 2
cd drain
..\..\Bin\GWM-VI.exe ..\data\drain_rate_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with drain_rate_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\drain_rate_serial.gwm
goto end


:end
pause
