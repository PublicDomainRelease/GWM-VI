@echo off
REM Run Storage test with GWM-VI 

if "%1"=="serial" goto serial
if "%1"=="SERIAL" goto serial


echo Running GWM-VI.exe with storage_pll.gwm
cd ..
call StartRunners 2
cd storage
..\..\Bin\GWM-VI.exe ..\data\storage_pll.gwm
goto end

:serial
echo Running GWM-VI.exe with storage_serial.gwm
..\..\Bin\GWM-VI.exe ..\data\storage_serial.gwm
goto end


:end
pause
