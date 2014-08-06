@echo off
REM  STREAMFLOW Composite Model Run
REM  Invoke MMProc to preprocess, run streamflow simulation,
REM  and postprocess 
echo streamflow_cmr.bat is deleting old Modflow output files...
if exist streamflow.lst del streamflow.lst
if exist streamflowstate.bin del streamflowstate.bin 
if exist streamflowheads.bin del streamflowheads.bin
echo streamflow_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo streamflow_cmr.bat has finished.
