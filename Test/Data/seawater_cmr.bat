@echo off
REM  SEAWATER Composite Model Run
REM  Invoke MMProc to preprocess, run seawater simulation,
REM  and postprocess 
echo seawater_cmr.bat is deleting old Modflow output files...
if exist seawater.lst del seawater.lst
if exist seawaterhd.bin del seawaterhd.bin 
echo seawater_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo seawater_cmr.bat has finished.
