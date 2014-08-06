@echo off

REM This is StartRunners.bat
REM StartRunners takes one argument, an integer, 
REM which is the number of runners to be started.

if "%1"=="" goto NoRunners
set N=%1

if %N%==0 goto NoRunners

cd Runner1
Echo Starting Runner 1
Start "Runner 1" /min ..\..\Bin\JRunnerM

if %N%==1 goto GoHome

cd ..\Runner2
Echo Starting Runner 2
Start "Runner 2" /min ..\..\Bin\JRunnerM

if %N%==2 goto GoHome

:GoHome
cd ..\
goto End

:NoRunners
Echo No Runners started

:End
