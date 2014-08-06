      MODULE GWM1RMS3
C     VERSION: 16JULY2009
      IMPLICIT NONE
      PRIVATE
      PUBLIC::IRM,AMAT,CST,BNDS,RHS,NRMC,NV,NDV,NCON,NCONF,NVF,NSIGDIG,
     1        NPGNMX,IPGNA,NPGNA,NINFMX,DINIT,DMIN,DSC,DELTA,AFACT,
     2        PGFACT,LPITMAX,VARBASE,CRITMFC,MFCNVRG,IBASE,IREF,DEWATER,
     3        DEWATERQ,DELINC,NONLIN,SLPITCNT,SLPITMAX,LASTLP,NRESET,
     4        SLPITPRT,SLPINFCNT,SLPVCRIT,SLPZCRIT,OBJOLD,FVOLD,BBITMAX,
     5        BBITPRT,RHSREL,RHSRLL,RHSREU,RHSRLU,CSTREL,CSTRLL,CSTREU,
     6        CSTRLU,RHSRLB,RHSRUB,RHSROR,CSTROR,CSTRLB,CSTRUB,CONTYP,
     7        RANGENAME,RANGENAMEF,RHSIN,RHSINF,RANGEFLG,OBJ,HCLOSEG
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
C
C-----GENERAL RMS VARIABLES
      LOGICAL(LGT),SAVE::NONLIN
      INTEGER(I4B),SAVE::IBASE,IREF
      REAL(DP),SAVE::OBJ
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MFCNVRG,DEWATER,DEWATERQ
      INTEGER(I4B),SAVE::NRMC,NV,NDV,NCON,NCONF,NVF
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::RANGENAME(:),RANGENAMEF(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::CONTYP(:)
      REAL(DP),SAVE,ALLOCATABLE::RHSIN(:),RHSINF(:)
C
C       NRMC     -number of constraints in the response matrix
C       NCON     -total number of constraints equations
C       NV       -total number of variables, including slacks/surpluses
C       NDV      -total number of decision variables + 1
C       NONLIN   -indicator that problem will be solved as a nonlinear problem 
C       IBASE    -indicator that base values have been read (1-yes, 0-no)
C       IREF     -indicator that reference simulation is needed
C       OBJ      -current value of objective function
C       DEWATER  -status of presence of dewatering at constraints in most recent simulation
C       DEWATERQ -status of presence of dewatering at wells in most recent simulation
C       MFCNVRG  -status of flow-process convergence for most recent simulation
C       RANGENAME-mapping from a variable number to a name; used for output 
C       RHSIN    -original input value of the rhs for each constraint;used for output
C       RANGENAMEF,RHSINF-final values of RANGENAME and RHSINF;used for output
C       NVF,NCONF-final values of NVF and NCONF; used for output
C       CONTYP   -indicates if constraint equality, inequality or transformed inequality
C
C-----VARIABLES FOR RESPONSE MATRIX AND PERTURBATION CALCULATIONS
      INTEGER(I4B),SAVE::IRM
      REAL(DP),SAVE::DINIT,DMIN,DSC,DELTA,AFACT,PGFACT,CRITMFC
      REAL(DP),SAVE,ALLOCATABLE::DELINC(:)
      INTEGER(I4B),SAVE::NSIGDIG,NPGNMX,NRESET
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPGNA,NPGNA
      REAL(DP),SAVE::HCLOSEG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VARBASE
C
C       IRM      -indicator for input/output of response matrix
C       DINIT,DMIN,DSC -input parameters that control perturbation parameter
C       DELTA    -computed value of perturbation parameter
C       DELINC   -array for storing individual perturbation values
C       AFACT    -relaxation parameter on the base solution
C       NRESET   -counter for number of base solution resets
C       PGFACT   -perturbation step length adjustment factor
C       CRITMFC  -criteria used to determine acceptability of GWF solution
C       NSIGDIG  -required number of significant digits in response coefficient
C       HCLOSEG  -flow process convergence criteria used for perturbation testing
C       NPGNMX   -maximum number of perturbation generation attempts before failure
C       IPGNA    -index of perturbation generation attempts
C       NPGNA    -number of perturbation generation attempts for this variable
C       VARBASE  -storage for base value of perturbed flow-rate variable
C   
C-----VARIABLES RELATED TO SOLVING THE LINEAR PROGRAM
      REAL(DP),SAVE,ALLOCATABLE::AMAT(:,:),CST(:),BNDS(:),RHS(:)
      INTEGER(I4B),SAVE::LPITMAX
C
C       LPITMAX  -maximum number of iterations allowed by the LP solver
C       AMAT     -constraint matrix with dense storage
C       CST      -storage for cost coefficients and LP solution
C       BNDS     -storage for upper bounds
C       RHS      -storage for right hand side and shadow prices
C
C-----VARIABLES RELATED TO SOLVING SEQUENTIAL LINEAR PROGRAM
      INTEGER(I4B),SAVE::SLPITCNT,SLPITMAX,SLPITPRT,SLPINFCNT,NINFMX
      REAL(DP),SAVE::SLPVCRIT,SLPZCRIT,OBJOLD
      REAL(DP),SAVE,ALLOCATABLE::FVOLD(:)
C
C       SLPITCNT -SLP iteration counter 
C       SLPITMAX -maximum number of SLP iterations
C       SLPINFCNT-counter for number of iterations that are infeasible
C       NINFMX   -maximum number of SLP iterations that are infeasible before failure
C       SLPVCRIT -SLP convergence criteria on difference in variables
C       SLPZCRIT -SLP convergence criteria on difference in objective
C       OBJOLD   -objective value from prior iteration
C       FVOLD    -flow variable values from prior iteration
C       SLPITPRT -SLP print flag for detailed output (1-yes,0-no)
C
C-----VARIABLES RELATED TO BRANCH AND BOUND ALGORITHM
      INTEGER(I4B),SAVE::BBITMAX,BBITPRT
      LOGICAL(LGT),SAVE::LASTLP
C
C       BBITMAX  -maximum number of branch and bound subproblems
C       BBITPRT  -print flag for detailed output (1-yes,0-no)
C       LASTLP   -indicator that this LP call is the last one
C
C-----ARRAYS THAT STORE LINEAR PROGRAM RANGE ANALYSIS RESULTS
      INTEGER(I4B),SAVE::RANGEFLG
      INTEGER(I4B),SAVE,ALLOCATABLE::RHSREL(:),RHSRLL(:),RHSREU(:),
     &                               RHSRLU(:),CSTREL(:),CSTRLL(:),
     &                               CSTREU(:),CSTRLU(:)
      REAL(DP),SAVE,ALLOCATABLE::RHSRLB(:),RHSRUB(:),RHSROR(:),
     &                           CSTROR(:),CSTRLB(:),CSTRUB(:)
C
C       RANGEFLG -indicator for performing range analysis (1-yes, 0-no)
C       RHSRLB   -lower bound on RHS
C       RHSRUB   -upper bound on RHS
C       RHSROR   -original RHS
C       RHSREL   -entering variable if RHS lower bound is exceeded
C       RHSRLL   -leaving variable if RHS lower bound is exceeded
C       RHSREU   -entering variable if RHS upper bound is exceeded
C       RHSRLU   -leaving variable if RHS upper bound is exceeded
C       CSTRLB   -lower bound on cost coefficient
C       CSTRUB   -upper bound on cost coefficient
C       CSTROR   -original cost coefficient
C       CSTREL   -entering variable is cost lower bound is exceeded
C       CSTRLL   -leaving variable if cost lower bound is exceeded
C       CSTREU   -entering variable if cost upper bound is exceeded
C       CSTRLU   -leaving variable if cost upper bound is exceeded
C
C
      END MODULE GWM1RMS3
