@echo off
REM  drain_volume Composite Model Run
REM  Invoke MMProc to preprocess, run drain_volume simulation,
REM  and postprocess 
echo drain_volume_cmr.bat is deleting old Modflow output files...
if exist drain_volume.lst del drain_volume.lst
if exist drainstate.bin del drainstate.bin
if exist drainheads.bin del drainheads.bin
echo drain_volume_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo drain_volume_cmr.bat has finished.
