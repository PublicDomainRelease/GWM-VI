@echo off
REM  SUPPLY2 Composite Model Run
REM  Invoke MMProc to preprocess, run supply2 simulation,
REM  and postprocess 
echo supply2_cmr.bat is deleting old Modflow output files...
if exist supply2.lst del supply2.lst
if exist supply2sf.bin del supply2sf.bin 
if exist supply2hd.bin del supply2hd.bin
echo supply2_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo supply2_cmr.bat has finished.
