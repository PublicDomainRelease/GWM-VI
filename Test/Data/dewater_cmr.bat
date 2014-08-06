@echo off
REM  DEWATER Composite Model Run
REM  Invoke MMProc to preprocess, run dewater simulation,
REM  and postprocess 
echo dewater_cmr.bat is deleting old Modflow output files...
if exist dewater.lst del dewater.lst
if exist dewaterhd.bin del dewaterhd.bin 
if exist dewatercbc.bin del dewatercbc.bin
echo dewater_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo dewater_cmr.bat has finished.
