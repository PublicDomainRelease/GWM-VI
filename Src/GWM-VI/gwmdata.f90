MODULE GWMDATA
  USE DATATYPES
  USE GLOBAL_DATA, ONLY: LENDNAM, MAX_STRING_LEN
  IMPLICIT NONE
  SAVE
  PRIVATE
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  !
  !   Public data
  PUBLIC INCTRL, INDIS, IOUTMESS, &
         AVET, BINC, DNPP, &
         IPERTURBMETHOD, IPTR, IUDIS,   &
         MAXRECL, MCVCAT, MCVUSE, MODELVAL, NCATS, NDINC, NNEGT,   &
         NOBS, NOMIT, NOPNT, NPOST, NRUNS, NW, OBSNAM, OUTNAME, PADJ,   &
         PARNAM, PRECISPRO, PVAL, PVALTMP, RSQ, RSQP,   &
         TEMPVAL, HCLOSET, HDRY, HNOFLO, WELLFILE
  PUBLIC HOBDRY, HOBDRYSET, NDEP, NPTC, NDDEP, NEDEP, NTDEP, NDEPSTAT
  PUBLIC :: OPTHEAD, MCLHEAD, MIFHEAD, MOFHEAD, PLLHEAD, RUNHEAD, &
            SIMHEAD, TAIL, PRECISION, DECPOINT, DELB, PARGP, &
            PASSIGNED, PERTURBMIN, PERTURBRAT, PVALFWD, PVALMAX,  &
            PVALMIN, WTVEC
  !
  !
  TYPE (LLIST), POINTER :: OPTHEAD ! Pointer to head of Options input-block list
  TYPE (LLIST), POINTER :: SIMHEAD ! Pointer to head of Simulation input-block list
  TYPE (LLIST), POINTER :: MCLHEAD ! Pointer to head of Model_Command_Lines input-block list
  TYPE (LLIST), POINTER :: MIFHEAD ! Pointer to head of Model_Input_Files input-block list
  TYPE (LLIST), POINTER :: MOFHEAD ! Pointer to head of Model_Output_Files input-block list
  TYPE (LLIST), POINTER :: PLLHEAD ! Pointer to head of Parallel_Control input-block list
  TYPE (LLIST), POINTER :: RUNHEAD ! Pointer to head of Parallel_Runners input-block list
  TYPE (LLIST), POINTER :: TAIL
  !
  !   Model-calculated values categories
  INTEGER                            :: NCATS
  PARAMETER (NCATS=1)
  CHARACTER(LEN=6), DIMENSION(NCATS) :: MCVCAT
  LOGICAL, DIMENSION(NCATS)          :: MCVUSE
  DATA MCVCAT/'OBS   '/
  DATA MCVUSE/.FALSE./
  !
  !
  INTEGER                            :: IUDIS = 24   ! IUNIT element where INDIS is stored
  INTEGER                            :: MAXRECL ! Max. record length for GWMOUT
  INTEGER                            :: NOPNT   ! Decimal point protocol
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NW      ! Minimum word length of a parameter
  INTEGER                            :: PRECISPRO  ! Precision protocol
  !
  !   Variables related to a full parameter set (parameters to be estimated or
  !   analyzed, plus parameters that are fixed or not to be analyzed)
  CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PARNAM   ! Parameter names
  CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: PARGP    ! Parameter groups
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: DELB     ! Delta-B amount, as written to model input file
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: PVAL     ! Current parameter values
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: PVALFWD  ! Forward-perturbed parameter value
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: PVALTMP  ! Parameter temporary values
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: PVALMIN  ! Minimum parameter value
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: PVALMAX  ! Maximum parameter value
  LOGICAL, ALLOCATABLE, DIMENSION(:)           :: PADJ       ! Is parameter adjustable?
  LOGICAL, ALLOCATABLE, DIMENSION(:)           :: PASSIGNED  ! Has PVAL been assigned?
  INTEGER                                      :: NPARGPCOLSAR, NPARGPS
  CHARACTER (LEN=8)                            :: PRECISION='single' ! Precision protocol
  CHARACTER (LEN=8)                            :: DECPOINT='point'   ! Decimal-point protocol
  !
  !   Variables related to parameter subset (only parameters to be estimated or
  !   analyzed)
  INTEGER, ALLOCATABLE, DIMENSION(:)      :: IPTR ! Position of parameter in full set
  !
  !   For PARAMETER_GROUPS input, do not define a default order
  INTEGER, PARAMETER                      :: NPARGPCOLS = 0
  CHARACTER(LEN=40), DIMENSION(1), TARGET :: PARGPCOL = (/' '/)
  !
  !   For PARAMETER_DATA input, do not define a default column order
  INTEGER, PARAMETER                      :: NPARCOLS = 0
  CHARACTER(LEN=40), DIMENSION(1), TARGET :: PARCOL = (/' '/)
  !
  !   For PARAMETERVALUES input, define default column order
  INTEGER, PARAMETER                            :: NPVCOLS = 2
  CHARACTER(LEN=40), DIMENSION(NPVCOLS), TARGET :: PVCOL =    &
      (/ 'PARAMNAME   ','STARTVALUE  '/)
  !
  !   Variables related to dependents
  INTEGER               :: NOBS    ! Number of observations
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: MODELVAL      ! Simulated equivalents
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: TEMPVAL       ! Temporary storage for simulated equivalent
  CHARACTER(LEN=LENDNAM), ALLOCATABLE, DIMENSION(:) :: OBSNAM
  DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:) :: WTVEC         ! (1/variance)
  !
  !   Variables related to calculation of sensitivities
  DOUBLE PRECISION                            :: HDRY, HNOFLO
  REAL                                        :: HCLOSET
  INTEGER, ALLOCATABLE, DIMENSION(:)          :: IPERTURBMETHOD
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PERTURBRAT ! Perturbation ratio for each parameter
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PERTURBMIN ! Minimum (absolute) perturbation for each parameter
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BINC       ! Perturbation increment for each parameter
  !
  !   Other application variables
  CHARACTER(LEN=MAX_STRING_LEN) :: OUTNAME, WELLFILE
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DNPP
  ! NPTC is total number of parameters, including HCLOSE
  INTEGER :: NDINC, NNEGT, NOMIT, NPOST, NRUNS, NPTC, NDEP, NDDEP, NEDEP, &
             NTDEP, NDEPSTAT
  INTEGER :: INDIS=0, INCTRL=0, IOUTMESS=0
  DOUBLE PRECISION :: AVET, RSQ, RSQP, WTRL
  DOUBLE PRECISION :: HOBDRY
  LOGICAL :: HOBDRYSET
END MODULE GWMDATA