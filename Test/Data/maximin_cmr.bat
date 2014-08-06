@echo off
REM  MAXIMIN Composite Model Run
REM  Invoke MMProc to preprocess, run maximin simulation,
REM  and postprocess 
echo maximin_cmr.bat is deleting old Modflow output files...
if exist maximin.lst del maximin.lst
if exist maximinheads.bin del maximinheads.bin
if exist maximin.gwmwell del maximin.gwmwell
echo maximin_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo maximin_cmr.bat has finished.
