@echo off
REM  Run Dewater test with GWM-VI (parallel or 
REM  serial)

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial

echo Running GWM-VI.exe with dewater_pll.gwm
cd ..
call StartRunners 2
cd dewater
..\..\Bin\GWM-VI.exe ..\data\dewater_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with dewater_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\dewater_serial.gwm
goto end



:end
pause
