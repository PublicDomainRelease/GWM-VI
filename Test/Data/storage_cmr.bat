@echo off
REM  STORAGE Composite Model Run
REM  Invoke MMProc to preprocess, run storage simulation,
REM  and postprocess 
echo storage_cmr.bat is deleting old Modflow output files...
REM if exist storage.lst del storage.lst
REM if exist storageheads.bin del storageheads.bin 
REM if exist storagestate.bin del storagestate.bin
echo storage_cmr.bat is invoking MMProc.exe...
..\..\Bin\MMProc.exe
echo storage_cmr.bat has finished.
