@echo off
REM  DEWATER Composite Model Run
REM  Invoke MMProc to preprocess, run dewatermb simulation,
REM  and postprocess 
echo dewatermb_cmr.bat is deleting old Modflow output files...
if exist dewatermb.lst del dewatermb.lst
if exist dewatermbhd.bin del dewatermbhd.bin 
echo dewatermb_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo dewatermb_cmr.bat has finished.
