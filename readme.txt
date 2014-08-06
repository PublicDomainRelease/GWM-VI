readme.txt

                 GWM-VI Suite of Programs - Version: 1.0.1
              Groundwater-Management with Parallel Processing 
                       for Multiple MODFLOW Versions


NOTE: Any use of trade, product or firm names is for descriptive purposes 
      only and does not imply endorsement by the U.S. Government.

This version of the GWM-VI suite of programs is packaged for personal computers 
using the Microsoft Windows Vista, XP, 7, or 8 operating systems.  Executable 
files for personal computers are provided as well as the source code.  The 
source code can be compiled to run on other computers.

IMPORTANT: Users should also review the file release.txt, which describes 
changes that have been introduced into GWM-VI with each official release; these 
changes may substantially affect users.

Instructions for installation, execution, and testing of GWM-VI are provided 
below.


                            TABLE OF CONTENTS

                         A. DISTRIBUTION FILE
                         B. INSTALLING
                         C. EXECUTING GWM-VI
                         D. TESTING
                         E. COMPILING
                         F. USING EXISTING MODFLOW AND GWM-2005 FILES

A. DISTRIBUTION FILE

The following compressed file contains the GWM-VI distribution:

         GWMVI_1_0_1.zip

The distribution file contains:

         Compiled runfiles and source code for GWM-VI
         Test data sets
         Documentation files

The distribution file can be uncompressed by any utility or program capable of 
uncompressing a zip file and extracting its contents.  The Microsoft Windows 
operating system includes such a utility.  An alternative is to use 7Zip 
(http://www.7-zip.org/).

Uncompressing the distribution file creates a directory containing a number of 
subdirectories and individual files.  The installation instructions assume that 
the files are restored into folder C:\WRDAPP.  If so, the following folder 
structure will be created in C:\WRDAPP:

| GWMVI_1_0_1
   | Bin              ; Executable files for personal computers
   | Doc              ; Documentation files
   | Src              ; Source code for use on any computer
      | GWM-VI        ; Source for GWM-VI program
      | JRunnerM      ; Source for JRunnerM program
      | JUPITER_lib   ; Source for JUPITER API library
      | MMProc        ; Source for MMProc program
      | PrecUtls_lib  ; Source for PrecUtls library
   | Test             ; Examples/verification tests
      | Data          ; Input files for all tests
      | Dewater       ; Files for the Dewater test
      | Dewatermb     ; Files for the Dewatermb test
      | Drain         ; Files for the Drain test
      | Maximin       ; Files for the Maximin test
      | MNW-supply    ; Files for the MNW-supply test
      | Runner1       ; Runner folder for parallel processing
      | Runner2       ; Runner folder for parallel processing
      | Seawater      ; Files for the Seawater test
      | Storage       ; Files for the Storage test
      | Streamflow    ; Files for the Streamflow test
      | Supply2       ; Files for the Supply2 test
      | Test-out      ; Output files for verifying all tests
      | Test-run      ; Folder for running all tests at once


It is recommended that no user files are kept in the GWMVI_1_0_1 folder 
structure.  

Included in the Doc folder are various documentation files.  Some of them are 
Portable Document Format (PDF) files. The PDF files are readable and printable 
on various computer platforms using Acrobat Reader from Adobe. The Acrobat 
Reader is freely available from the following World Wide Web site:

      http://www.adobe.com/


B. INSTALLING

To make the executable files of the GWM-VI suite accessible from any folder, 
the folder containing the executables (Bin) should be included in the PATH 
environment variable.  Also, if a prior release of GWM-VI is installed on your 
system, the folder containing the executables for the prior release should be 
removed from the PATH environment variable.

As an alternative, the executable files in the Bin folder can be copied into 
a folder already included in the PATH environment variable. 

       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
                 WINDOWS 7 SYSTEMS
             
From the Start menu, select Control Panel.  Select System and Security, and 
within that screen choose the System option. Then select the Advanced System 
Settings option.  Select the Environment Variables button.  In the System 
Variables pane, select the PATH variable followed by Edit.  In the Edit window, 
add ";C:\WRDAPP\GWMVI_1_0_1\Bin" to the end of the Variable Value (ensure that 
the current contents of the User Value are not deleted) and click OK. Click OK 
in the Environment Variables window and then exit from the control panel 
windows.  Initiate and use a new Windows Command window.

       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
                  WINDOWS XP SYSTEMS
             
From the Start menu, select Settings and then Control Panel.  Double click 
System and select the Advanced tab.  Click on Environment Variables.  If a PATH 
user variable already is defined, click on it in the User Variables pane, then 
click Edit.  In the Edit User Variable window, add ";C:\WRDAPP\GWMVI_1_0_1\Bin" 
to the end of the Variable Value (ensure that the current contents of the User 
Value are not deleted) and click OK.  If a PATH user variable is not already 
defined in the User variables pane of the Environment Variables window, click 
New.  In the New User Variable window, define a new variable PATH as shown 
above.  Click OK.  Click OK in the Environment Variables window and again in 
the System Properties window.  Initiate and use a new Windows Command window.

       HOW TO ADD TO THE PATH ENVIRONMENT VARIABLE
                 WINDOWS VISTA SYSTEMS
             
From the Start menu, select Settings and then Control Panel.  Select System & 
Maintenance followed by System.  Choose the Advanced System option.  Select the 
Settings Task, and then select the Environmental Variables button.  In the 
System Variables pane, select the PATH variable followed by Edit.  In the Edit 
window, add ";C:\WRDAPP\GWMVI_1_0_1\Bin" to the end of the Variable Value 
(ensure that the current contents of the User Value are not deleted) and click 
OK. Click OK in the Environment Variables window and then exit from the control 
panel windows.  Initiate and use a new Windows Command window.


C. EXECUTING GWM-VI

After the executable files in the Bin folder are installed in a folder that is 
included in your PATH, GWM-VI is initiated in a Windows Command-Prompt window 
using the following command: 

          GWM-VI [Fname]

The optional Fname argument is the input file.  If no argument is used, the 
user is prompted to enter the input-file name.  The input file is the GWM file 
for the problem.  In addition to the entries that define the management 
problem, the input file contains entries for the MODFLOW name file and the 
GWM-VI control file.  All other input files (such as SimulatedValues.jif, 
Modflow_status.jif, MMProc.in.jtf, mnw2_input.jtf, and MMProc.in) are 
automatically created by GWM-VI.


D. TESTING

Test data sets are provided to verify that GWM-VI is correctly installed and 
running on the system.  The tests may also be looked at as examples of how to 
use the program.  The Test folder contains subfolders for running each test.  
The Test\Test-out folder contains the output files from running each test.  
Each of the following folders contains a batch file that can be used to run one 
test in either serial or parallel mode:

  Dewater
  Dewatermb
  Drain
  Maximin
  MNW-supply
  Seawater
  Storage
  Streamflow
  Supply2

Additional information on each of the tests is provided in the readme_tests.txt
file in the Test folder.

For convenience, batch files for all tests also are provided in the Test-run 
folder.  Each test can be run from this folder, and multiple tests can be run 
without changing the current folder.

Each test can be run in a command-prompt window by making one of the folders 
listed above the current folder and entering the name of the corresponding 
batch file.  For example, to run the Dewater test, make Test\Dewater the 
current folder and type the name of the batch file in that folder (dewater.bat) 
at the command prompt.  To make a serial-mode run, add "serial" as a command-
line option:

dewater serial

If your computer has multiple processor cores, a parallel-mode run of GWM-VI 
can be started by omitting the "serial" command-line option.

MODFLOW input that is not manipulated by GWM-VI is read from the Test\Data 
folder.  Prior to each MODFLOW run, GWM-VI writes additional MODFLOW input to 
one or more files in the current folder.  Output generated by GWM-VI also is 
written to the current folder and can then be compared to output for the test 
of interest in the Test\Test-out folder.  When parallel mode is used, the 
Runner1 and Runner2 folders under the Test folder are used for making some 
MODFLOW runs.


E. COMPILING

The executable files provided in the Bin folder were created using the Intel
Visual Fortran 13.0 compiler.  Although executable versions of the programs 
are provided, the source code is provided in the Src folder so that the 
programs can be recompiled if necessary.  However, the USGS cannot provide 
assistance to those compiling GWM-VI.  In general, the requirements are a 
Fortran compiler and the knowledge of using the compilers.  Additional guidance 
follows.

Fortran-90 source code is provided to build five projects: two are intended to 
be built as static libraries (JUPITER_lib and PrecUtls_lib), and three are to 
be built as console-mode executable programs (GWM-VI, JRunnerM, and MMProc).  
Each project has, under the Src folder, a corresponding subfolder containing 
source-code files.  Source-code files ending in the extension ".f" or ".for" 
use the Fortran-90 fixed source form; whereas files ending in the extension 
".f90" use the free source form.  GWM is designed to use double-precision 
variables for floating-point numbers; however, some of the source code declares 
some variables as "REAL," which by default would be single-precision.  To 
successfully compile the GWM-VI suite of programs, a compiler option needs to 
be used when compiling the GWM-VI and MMProc projects to promote all "REAL" 
variables to "DOUBLE PRECISION."  The JUPITER_lib, PrecUtls_lib, and JRunnerM 
projects are designed to be compiled without the compiler option set to promote 
variables declared as "REAL."  Because of dependencies among projects, certain 
projects need to be built before others.  The following build order, with 
dependencies as noted, can be used:

1. JUPITER_lib
2. PrecUtls_lib (dependent on JUPITER_lib)
3. JRunnerM (dependent on JUPITER_lib)
4. GWM-VI (dependent on JUPITER_lib)
5. MMProc (dependent on JUPITER_lib and PrecUtls_lib)

Each of the project folders contains source-code files which, together with 
any library(ies) on which the project depends as indicated above, enable the 
corresponding project to be built.


F. USING EXISTING MODFLOW AND GWM-2005 FILES

GWM-VI imposes several requirements on MODFLOW and GWM input files that may 
require modification of existing versions of those files used either with 
MODFLOW or with GWM-2005.

GWM-VI is generally backwards-compatible with prior versions of GWM and 
includes all GWM capabilities except those few specifically excluded in the 
documentation report, specifically, use of input variable CRITMFC, references 
to heads at MNW wells, use of multi-grid flow processes and constraints on 
streams modeled with the STR package.  The user needs to ensure that the names 
of decision variables and head, streamflow or drain constraints follow the 
JUPITER naming convention as described in the documentation report.

MODFLOW input files must adhere to several requirements:  

1. The Fortran file unit numbers specified in the MODFLOW NAME file are 
limited.  Unit numbers 96-99 are reserved by MODFLOW.  In addition, GWM-VI 
reserves unit numbers 980-1000. 

2. In the BAS input file, the Options in data set 1 must include the FREE 
option.  This ensures that stresses can be read from and written to files with 
adequate precision.  Including the FREE option may require reformatting other 
MODFLOW input files.

3. GWM-VI reads binary output files from MODFLOW to determine heads and flows 
needed to evaluate constraints and state variables.  It is preferred that these 
binary files be double precision, however, GWM-VI can operate with single-
precision binary files.  GWM-VI will determine if a given binary file is single 
or double precision.
