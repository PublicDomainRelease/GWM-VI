readme_tests.txt

                      GWM-VI Sample Problems

Data files for 9 test problems are provided to confirm that GWM-VI is 
correctly installed and running on the system.  The tests also may be used as 
examples of how to formulate optimization problems and how to organize data 
for serial and parallel runs.
  
The input data are organized by sample problem:

Test name        Description of test
------------     ------------------------------------------------------------
Dewater          Dewater problem, linear formulation
Dewatermb        Dewater problem, mixed-binary linear formulation
Drain            Drain problem with state variables (rate and volume version)
Maximin          Maximin problem with state variables
MNW-supply       MNW-supply problem, demonstrating use of MNW2 Package
Seawater         Seawater problem, nonlinear formulation
Storage          Storage problem with state variables
Streamflow       Streamflow problem with state variables
Supply2          Supply problem, nonlinear formulation


The Dewater, Dewatermb, and Seawater problems are described in Ahlfeld and 
others (2005). The Supply2 problem was initially the Supply problem as 
described in Ahlfeld and others (2005). However, experience with the sample 
problem indicated that there could be numerical-stability issues associated 
with its solution. Therefore, the problem was revised and renamed Supply2, as 
described in Ahlfeld and others (2009).

The Maximin, Storage, and Streamflow problems are described in Ahlfeld and 
others (2011).

The Drain sample problem is described in the document 'Drains.pdf,' which is 
provided in the doc subdirectory.  It has two variations in two batch files: 
(1) drain_rate.bat, in which the instantaneous drain flow rate at a specified 
time is constrained, and (2) drain_volume.bat, in which the cumulative drain 
volume over a specified time is constrained.

The MNW-supply sample problem is described in Ahlfeld and Barlow (2013).

The test problems can be run by executing one of the batch files in directory 
test-run.  Executation is transferred to the appropriate test problem sub-
directory in the Data directory.  Existing output for each test problem can 
be found in directory test-out and can be compared with the corresponding 
output in the test problem sub-directory.

Each test problem can be run in either parallel or serial mode.  For parallel
mode, just execute the batch file.  For serial mode, enter the name of the 
test problem as a command in a DOS command-prompt window with the word 
'serial' following the command.  For example to run the dewater test problem
in serial mode, from a DOS window, manuever to the test-run directory and 
type 'dewater serial'.


References:

Ahlfeld, D.P., Baker, K.M., and Barlow, P.M., 2009, GWM-2005--A Groundwater-
Management Process for MODFLOW-2005 with Local Grid Refinement (LGR)
capability: U.S. Geological Survey Techniques and Methods, book 6, chap. A33,
65 p.

Ahlfeld, D.P., and Barlow, P.M., 2013, Use of multi-node wells in the 
Groundwater-Management Process of MODFLOW–2005 (GWM–2005): U.S. Geological 
Survey Techniques and Methods, book 6, chap. A47, 26 p., 
http://pubs.usgs.gov/tm/06/a47/.

Ahlfeld, D.P., Barlow, P.M., and Baker, K.M., 2011, Documentation for the
State Variables Package for the Groundwater-Management Process of
MODFLOW-2005 (GWM-2005): U.S. Geological Survey Techniques and Methods, 
book 6, chap. A36, 45 p.

Ahlfeld, D.P., Barlow, P.M., and Mulligan, A.E., 2005, GWM--A ground-water
management process for the U.S. Geological Survey modular ground-water
model (MODFLOW-2000): U.S. Geological Survey Open-File Report 2005-1072,
124 p.
