@echo off
REM  drain_rate Composite Model Run
REM  Invoke MMProc to preprocess, run drain_rate simulation,
REM  and postprocess 
echo drain_rate_cmr.bat is deleting old Modflow output files...
if exist drain_rate.lst del drain_rate.lst
if exist drainstate.bin del drainstate.bin
if exist drainheads.bin del drainheads.bin
echo drain_rate_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo drain_rate_cmr.bat has finished.
