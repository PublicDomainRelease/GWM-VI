@echo off
REM  mnwsupply Composite Model Run
REM  Invoke MMProc to preprocess, run mnwsupply simulation,
REM  and postprocess 
echo mnwsupply_cmr.bat is deleting old Modflow output files...
if exist mnwsupply.lst del mnwsupply.lst
if exist mnwsupplyheads.bin del mnwsupplyheads.bin
if exist mnwsupply.gwmwell del mnwsupply.gwmwell
echo mnwsupply_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo mnwsupply_cmr.bat has finished.
