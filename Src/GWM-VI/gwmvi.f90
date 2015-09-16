PROGRAM GWMVI
!     ******************************************************************
!     MAIN CODE FOR GWM-VI
!         GROUNDWATER MANAGEMENT - VERSION INDEPENDENT
!     ******************************************************************
!        SPECIFICATIONS:
!     ------------------------------------------------------------------
USE GWM_STOP, ONLY: GSTOP
USE GLOBAL, ONLY: NCOL, NROW, NLAY, NPER, PERLEN, IUNIT, NBOTM, NCNFBD, &
                  ITMUNI, LENUNI, ITRSS, IXSEC, INBAS, IFREFM, NODES, &
                  IOUT, NIUNIT, NSTP
USE GWM_SUBS, ONLY: MNWCOPYOK, MNW2FILE, MNW2FILEORIG
USE GWMDATA, ONLY:  &
  !   DATA
  INCTRL, IOUTMESS, AVET, BINC, DNPP, HCLOSET, IPERTURBMETHOD, IPTR, &
  MAXRECL, MCVCAT, MCVUSE, MODELVAL, NCATS, NDINC, NNEGT, &
  NOMIT, NOPNT, NPOST, NRUNS, NW, OBSNAM, OUTNAME, PADJ, PARNAM, &   
  PRECISPRO, PVAL, PVALTMP, RSQ, RSQP, & 
  NPTC, NDEP, NDEPSTAT
USE GWMMOD, ONLY:  &
  !   SUBROUTINES
  MAIN_INI, MAIN_INI_DEP, MAIN_INI_EQN, MAIN_INI_PARS, &   
  MAIN_DEF, MAIN_GEN, MAIN_UEV, MAIN_CLN, DECVARS_TO_PVAL, &
  MAIN_WRITE_MIF_BLOCK, MAIN_WRITE_MOF_BLOCK, MAIN_WRITE_MF_STATUS_JIF
USE GWM_UTLS, ONLY: APPEND_MESSAGE, COPY_FILE, TEXT_FILE_COPY
USE GWM_SUBS, ONLY: CLEAN_UP
USE PLLM, ONLY: PLLM_STOP_RUNNERS
USE GWM1BAS3SUBS, ONLY: GWM1BAS3AR
USE GWM1RMS3SUBS, ONLY: GWM1RMS3PL,GWM1RMS3PP,GWM1RMS3FP, &
                       GWM1RMS3FM,GWM1RMS3AP,GWM1RMS3OT
USE GWM1DCV3, ONLY: NFVAR, FVBASE, FVNCELL, FVMIN, FVMAX, FVNAME
USE GWM1DCV3, ONLY: MNW2HERE
USE GWM1BAS3, ONLY: GWMOUT, GWMWFILE
USE GWM1RMS3, ONLY: DELINC, CST, MFCNVRG, DEWATER, DEWATERQ, IPGNA, HCLOSEG
USE GWM1HDC3, ONLY: HDCNAME,NHB,NDD,NDF,NGD
USE GWM1STC3, ONLY: STCNAME,NSF,NSD,STCNUM
USE GWM1STA3, ONLY: SVNAME,STANUM,NHVAR,NRVAR,NSVAR,NDVAR,NSTADEP,SVILOC
USE MKMMPROC, ONLY: WRMMPROCIN,GWMWFILENAME
USE MKMMPROC_RDNAM, ONLY: READ_NAME_FILE_DIS
USE GWFBASMODULE,   ONLY: HDRY,HNOFLO
USE GWM1OBJ3, ONLY : SOLNTYP
!
!   Use modules of the JUPITER API
USE DATATYPES
USE GLOBAL_DATA, ONLY: AMESSAGE, HYPHENS, IVERB, LENDNAM, MAX_STRING_LEN
USE UTILITIES
USE BASIC, ONLY:   & 
  COVMATARR, DERIV_INTERFACE, BAS_INI_COVMAT, BAS_INI_GETOPTIONS, &   
  BAS_INI_MODELEXEC, BAS_GEN, BAS_EXE_SELECT, BAS_EXE, BAS_CLN
USE MODEL_IO, ONLY:   &   
  MIO_INI_ALLOC, MIO_INI_ARRAYS, MIO_INI_DIMENSION, MIO_INI_INPUTFILES, &   
  MIO_INI_OUTPUTFILES, MIO_INI_INSTRUCT1, MIO_INI_INSTRUCTALLOC, &   
  MIO_INI_INSTRUCT2, MIO_INI_TEMPLATE, MIO_ADA_WRITEFILES, MIO_EXT
USE DEPENDENTS, ONLY:   &   
  OBSVAL, WTCORRELATED, DEP_DX_WRITE_GM, DEP_INI_ALLOC, DEP_INI_READ, &   
  DEP_INI_STORE, DEP_UEV_DX_WRITE_OS, DEP_UEV_DX_WRITE_P,    &   
  DEP_UEV_DX_WRITE_R, DEP_UEV_RESIDUALS, DEP_UEV_WRITEOBSTABLE, DEP_CLN
USE PLLM, ONLY:   &   
  PLLM_INI_DISPATCHER, PLLM_CLN, PLLM_RMS2FP
USE PLL_GWM, ONLY: PLLGWM_MAKE_RUNS
!
IMPLICIT NONE
!
! Fortran parameters
!
! ASSIGN VERSION NUMBER AND DATE
CHARACTER*40 VERSION
CHARACTER(LEN=10)::MFVNAM      
PARAMETER (VERSION=' 1.0.2')! Insert version number here
PARAMETER (MFVNAM='VI')     ! Do not change - this is a flag for GWM1BAS3AR     
INTEGER JOBDIM, JOBLEN
PARAMETER (JOBDIM=5, JOBLEN=20)
!
!   Declare variables
INTEGER               :: INUNIT = 9  ! Unit number for input file
INTEGER, DIMENSION(8) :: IBDT        ! Begin time
CHARACTER(LEN=10)     :: CHDATE, CHTIME, CHZONE
CHARACTER(LEN=MAX_STRING_LEN) :: FILEIN=' '  ! Input data file name
CHARACTER(LEN=MAX_STRING_LEN) :: FILEOUT     ! Output data file name
INTEGER                       :: ICOMMAND    ! Number of model command line to execute
INTEGER                       :: IL, IFAIL
CHARACTER(LEN=180)            :: INSTRUCTION
INTEGER                       :: ISTAT       ! Status: 0 for success
INTEGER                       :: NPT         ! Total number of parameters
INTEGER                       :: NPE         ! Number of parameters with adjustable=yes
INTEGER                       :: NPD         ! Number of derived parameters
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XOBS   ! Sensitivity of simulated equivalents to parameters
INTEGER                       :: ICTRL       ! Control-loop job pointer
LOGICAL                       :: CTRLDONE
LOGICAL                       :: SENSDONE
CHARACTER(LEN=JOBLEN), DIMENSION(JOBDIM) :: CTRLJOB  ! Sequence of jobs required of Control loop
LOGICAL :: NEEDSENS = .TRUE.  ! True unless a forward-only run is specifically needed (other than for perturbation)
!
!   Parallel processing
INTEGER :: KLOOP    ! Counter for iterations potentially parallel outer loop
INTEGER :: KLPTR
INTEGER :: KPI
INTEGER :: KPPL     ! Counter for iterations of potentially parallel inner loops
INTEGER :: LCIS           ! size of CINSTSET array
INTEGER :: MAXRUNSPLL = 1
INTEGER :: NINSTRUCT      ! size of LOCINS array
INTEGER :: NMIFILE        ! Number of model-input files
INTEGER :: NMOFILE        ! Number of model-output files
INTEGER :: NRUNSPLL = 1
INTEGER :: NUMLADV
INTEGER :: NUMLOOP  ! Number of iterations for potentiall parallel outer loop
INTEGER :: NUMPPL   ! Number of iterations for potentially parallel inner loops
LOGICAL :: DO_PARALLEL
LOGICAL :: PARALLEL_ACTIVE
CHARACTER (LEN=MAX_STRING_LEN) :: COMMANDPLL
CHARACTER(LEN=6),   ALLOCATABLE, DIMENSION(:) :: CATGOR    ! Model-calculated value category for each model-output file
CHARACTER(LEN=1),   ALLOCATABLE, DIMENSION(:) :: CINSTSET  ! holds compressed instruction set
CHARACTER(LEN=MAX_STRING_LEN), ALLOCATABLE, DIMENSION(:) :: INSTRUCTFILE
INTEGER,            ALLOCATABLE, DIMENSION(:) :: KPE
INTEGER,            ALLOCATABLE, DIMENSION(:) :: LOCINS    ! pointer to instructions
CHARACTER(LEN=MAX_STRING_LEN), ALLOCATABLE, DIMENSION(:) :: MIFILE    ! Names of model-input files
CHARACTER(LEN=MAX_STRING_LEN), ALLOCATABLE, DIMENSION(:) :: MOFILE    ! Names of model-output files
CHARACTER(LEN=1),   ALLOCATABLE, DIMENSION(:) :: MRKDL
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PARVALSETS
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PARVALSETSE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DEPVALSETS
CHARACTER(LEN=MAX_STRING_LEN), ALLOCATABLE, DIMENSION(:) :: TFILE    ! Names of template files
!
INTEGER, ALLOCATABLE, DIMENSION(:) :: OMIT
DOUBLE PRECISION :: DTLA
LOGICAL :: DONOTEVALRUNS = .TRUE.
!
DATA CTRLJOB /'FORWARD','SENSITIVITY','SOLVE_GWM','REPEAT','STOP'/
!
CHARACTER*80 HEADNG(2)
CHARACTER(LEN=200) :: FNAME,LINE,GWMSTRG
! 
!  GWM loop controls and variables
LOGICAL :: FIRSTSIM,LASTSIM,FINISH,GWMCNVRG,PERTDONE
INTEGER :: IPERT,NPERT,INGWM
!
!  Other variables
LOGICAL :: LEX, OK
INTEGER :: I,NCOVMAT,OPTIONFLAG,ISCRATCH, NUNOPN
CHARACTER*2000 SIMVALJIF_FILE, MMPROCJTF_FILE, MNW2JTF_FILE, NAMFILE
REAL(KIND(1.0)) :: HCLOSETSP
TYPE (CDMATRIX) :: WTMAT    ! Weight matrix, not used
TYPE (CDMATRIX) :: WTMATSQR ! Square-root of weight matrix, not used

! for debugging
!integer :: ii,jj

INCLUDE 'openspec_f90.inc'
!
!   Format statements
    1 FORMAT (/,34X,'GWM-VI',/,  &
      4X,'U.S. GEOLOGICAL SURVEY GROUNDWATER MANAGEMENT',  &
      ' VERSION INDEPENDENT PROGRAM',/,17X,'Version ',A/)
    2 FORMAT(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ', &
      I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2,/)
    5 FORMAT(1X,'Enter name of input file:')
   10 FORMAT(A)
   50 FORMAT(20X,'Output from program ',A)
  210 FORMAT(/,1X,'Error in instruction: "',A,'"')
  250 FORMAT(1X,'Control-loop job = "',A,'" and KLPTR = ',I3)
  260 FORMAT(1X,'GWM solver iteration: ',I4)
  300 FORMAT(/,1X,'For perturbed parameter "',A, &
          '", model output follows:')
  400 FORMAT(/,'Warning: No file associated with unit ',i4,' in name file.',/, &
             'Flow-rate decision variables not written to Well-Type file.',/)
  900 FORMAT(/,1X,'Normal termination of ',A)
!
! ******************************************************************************
! ****************************** INITIALIZE ************************************
! ******************************************************************************
!
IFAIL = 0
PARALLEL_ACTIVE = .FALSE.
NCOVMAT = 0
IPERT = 0
OK = .TRUE.
MNWCOPYOK = .FALSE.
! Add "STATUSFILE modflow.status 1" line to jdispatch.rdy 
! to activate model status option in JRunnerM
OPTIONFLAG = 1
WRITE (*,1) TRIM(VERSION)
CALL DATE_AND_TIME(CHDATE,CHTIME,CHZONE,IBDT)
WRITE(*,2) (IBDT(I),I=1,3),(IBDT(I),I=5,7)
!
! Open main GWM input file
FNAME = ' '
FNAME = UTL_GETARG(1)      ! Obtain filename from screen prompt
! Not all compilers support command-line arguments, so provide 
! another way to get name of input file
IF (FNAME == ' ') THEN
  WRITE(*,5)
  READ(*,'(A)') FNAME
ENDIF
INQUIRE(FILE=FNAME,EXIST=LEX)
ISTAT = 0
IF (LEX) THEN
  INGWM = UTL_GETUNIT(980,1000)   ! Unit number for GWM input file
  OPEN(INGWM,FILE=FNAME,IOSTAT=ISTAT)
ELSE
  WRITE(*,22)TRIM(FNAME)
  22 FORMAT(1X,'ERROR: FILE "',A,'" DOES NOT EXIST')
  CALL UTL_STOP('')
ENDIF
IF (ISTAT .NE. 0) THEN
  WRITE(*,44)TRIM(FNAME)
  44 FORMAT(1X,'ERROR OPENING FILE: ',A)
  CALL UTL_STOP('')
ENDIF
! 
NUNOPN = UTL_GETUNIT(980,1000)    ! NUNOPN is temporary file to hold outpuf from MAIN_INI
OPEN(UNIT=NUNOPN,STATUS='SCRATCH')       
CALL SGWF2BAS7PNT(1)       ! GWM-VI allows only single grid problems
ALLOCATE(NCOL,NROW,NLAY,NPER,NBOTM,NCNFBD,ITMUNI,LENUNI,ITRSS)    
ALLOCATE(IXSEC,INBAS,IFREFM,NODES,IOUT,HDRY,HNOFLO)

ALLOCATE(IUNIT(NIUNIT))
IUNIT = 0
IUNIT(56) = INGWM          ! GWM1BAS3AR accesses IUNIT(56) to obtain unit number for GWM file
!
CALL MAIN_INI(INGWM,NUNOPN,NAMFILE) ! Read CTRL file and NAMFILE name; rewind INGWM
CALL SGWF2BAS7PSV(1)       ! Save pointers to data for grid 1
HCLOSETSP = HCLOSET        ! GWM1BAS3AR expects a single precision variable for HCLOSE
!
CALL READ_NAME_FILE_DIS(NAMFILE) ! Open NAMFILE to get DIS file and read GLOBAL data
!
CALL GWM1BAS3AR(IOUTMESS,NPER,PERLEN,HCLOSETSP,VERSION,MFVNAM,1) ! Read all GWM input; open GWMOUT; skip CTRL file
! 
REWIND(NUNOPN)             ! Append the output written by MAIN_INI to output file GWMOUT
99 CONTINUE
READ(NUNOPN,'(1X,A)',END=100)LINE
WRITE(GWMOUT,'(1X,A)')TRIM(LINE)
GO TO 99                                 
100 CONTINUE 
CLOSE(NUNOPN)                 
!
!   Initialize model-execution data (read MODEL_COMMAND_LINES input block)
CALL BAS_INI_MODELEXEC(INCTRL,IOUTMESS)
SIMVALJIF_FILE = 'SimulatedValues.jif'
MMPROCJTF_FILE = 'MMProc.in.jtf'
MNW2JTF_FILE = 'mnw2_input.jtf'
!
!   Write MMProc.in.jtf file, SimulatedValues.jif and mnw2_input.jtf files
CALL WRMMPROCIN(NAMFILE,SIMVALJIF_FILE,MMPROCJTF_FILE,MNW2JTF_FILE,GWMOUT)      
!
! Write modflow_status.jif
CALL MAIN_WRITE_MF_STATUS_JIF()
!
!   Initialize GWM loop controls
FIRSTSIM = .TRUE.
GWMCNVRG = .FALSE.
FINISH   = .FALSE.
LASTSIM  = .FALSE.
DO I=1,NFVAR
  MFCNVRG(I) = .TRUE.
ENDDO
!
!   Initialize flow-rate and derived parameter data 
CALL MAIN_INI_PARS(INGWM,NFVAR,FVMIN,FVMAX,FVBASE,FVNCELL,FVNAME, &
                   NPT,NPE,NPD,OK)
IF (.NOT. OK) GOTO 999
!
!   Allocate constraint-related arrays.  MAIN_INI_DEP invokes DEP_INI_READ
CALL MAIN_INI_DEP(NCOVMAT,NHB,NDD,NDF,NGD,HDCNAME,STCNUM,STCNAME,&
                  STANUM,SVNAME,NHVAR,NRVAR,NSVAR,NDVAR,NSTADEP,SVILOC)
!
!   Read, store, and echo constraint-related data (part 2).
CALL DEP_INI_ALLOC()
CALL DEP_INI_STORE(2,IOUTMESS,NCOVMAT,NDEPSTAT,COVMATARR,0, &
    OBSNAM,DTLA,WTMAT,WTMATSQR)
!
!   Initialize Model_IO module.
CALL MIO_INI_ALLOC(IFAIL,NPTC)
!
! Write a Model_Input_Files input block on scratch file and read it
ISCRATCH = UTL_GETUNIT(7,1000)
OPEN(ISCRATCH,STATUS='SCRATCH')
CALL MAIN_WRITE_MIF_BLOCK(ISCRATCH,MNW2FILE)
REWIND(ISCRATCH)
!
! Read, store, and echo file names for model input and output
CALL MIO_INI_INPUTFILES(ISCRATCH,IOUTMESS)  ! (read MODEL_INPUT_FILES block)
CLOSE(ISCRATCH)
!
! Write a Model_Output_Files input block on scratch file and read it
OPEN(ISCRATCH,STATUS='SCRATCH')
CALL MAIN_WRITE_MOF_BLOCK(ISCRATCH)
REWIND(ISCRATCH)
CALL MIO_INI_OUTPUTFILES(ISCRATCH,IOUTMESS) ! (read MODEL_OUTPUT_FILES block)
CLOSE(ISCRATCH)
!
!   Do initialization related to model-input files
CALL MIO_INI_TEMPLATE(IFAIL,NPTC,PARNAM,NW)
IF (IFAIL .NE. 0) THEN
  OK = .FALSE.
  GOTO 999
ENDIF
!
!   Initialize and populate instruction arrays
CALL MIO_INI_INSTRUCT1(IFAIL)
IF (IFAIL .NE. 0) THEN
  OK = .FALSE.
  GOTO 999
ENDIF
!
!   Allocate memory for storage of instructions
CALL MIO_INI_INSTRUCTALLOC(NDEPSTAT,IFAIL)
IF (IFAIL .NE. 0) THEN
  OK = .FALSE.
  GOTO 999
ENDIF
CALL MIO_INI_INSTRUCT2(IFAIL,NCATS,MCVCAT)
IF (IFAIL .NE. 0) THEN
  OK = .FALSE.
  GOTO 999
ENDIF
!
!   Initialize parallel-processing data
CALL MIO_INI_DIMENSION(LCIS,NINSTRUCT,NMIFILE,NMOFILE,NUMLADV)
ALLOCATE(CATGOR(NMOFILE),CINSTSET(LCIS),INSTRUCTFILE(NMOFILE), &
    LOCINS(NINSTRUCT),MIFILE(NMIFILE),MOFILE(NMOFILE), &
    MRKDL(NMOFILE),TFILE(NMIFILE))
CALL MIO_INI_ARRAYS(IFAIL,LCIS,NINSTRUCT,NMIFILE,NMOFILE,CATGOR, &
    CINSTSET,INSTRUCTFILE,LOCINS,MIFILE,MOFILE,MRKDL,TFILE)
IF (IFAIL .NE. 0) THEN
  CALL APPEND_MESSAGE('Programming error: MIO_INI_ARRAYS reports failure')
  OK = .FALSE.
  GOTO 999
ENDIF
!
!   (read PARALLEL_CONTROL and PARALLEL_RUNNERS blocks)
MAXRUNSPLL = NPE
CALL PLLM_INI_DISPATCHER(INCTRL,IOUTMESS,LCIS,NCATS,NINSTRUCT, &    
    NMIFILE,NMOFILE,NOPNT,NPT,NUMLADV,NDEPSTAT,PRECISPRO,CATGOR, &
    CINSTSET,OBSNAM,INSTRUCTFILE,LOCINS,MCVUSE,MIFILE,MOFILE, &    
    MRKDL,NW,PARNAM,TFILE,PARALLEL_ACTIVE,MAXRUNSPLL, &
    MMPROCJTF_FILE,MNW2JTF_FILE,MNW2HERE,OK,OPTIONFLAG)
IF (.NOT. OK) GOTO 999
!
!   Allocate arrays declared in main program unit
ALLOCATE(PARVALSETS(NPTC,MAXRUNSPLL),PARVALSETSE(NPE,MAXRUNSPLL), &
    DEPVALSETS(NDEPSTAT,MAXRUNSPLL),KPE(MAXRUNSPLL))
ALLOCATE(XOBS(NPE,NDEPSTAT),OMIT(NDEPSTAT))
PARVALSETS = 0.0D0
PARVALSETSE = 0.0D0
DEPVALSETS = 0.0D0
KPE = 0
XOBS = 0.0D0
OMIT = 0
!
IF (MNW2HERE) THEN
  ! Save copy of original MNW2 input file, before it gets modified
  MNWCOPYOK = TEXT_FILE_COPY(MNW2FILE,MNW2FILEORIG)
ENDIF
!
CTRLDONE = .FALSE.
SENSDONE = .TRUE.  ! Necessary initial condition
!
! ******************************************************************************
! ************************** TOP OF CONTROL LOOP *******************************
! ******************************************************************************
!
ICTRL = 0
CONTROL: DO WHILE (.NOT. CTRLDONE)
  !
  ! ****************************************************************************
  ! ************** DEFINE JOB OF CURRENT ITERATION OF CONTROL LOOP *************
  ! ****************************************************************************
  !
  CALL MAIN_DEF(GWMOUT,JOBDIM,JOBLEN,NPE,CTRLJOB,NEEDSENS,FIRSTSIM, &
                LASTSIM, FINISH,CTRLDONE,ICTRL,SENSDONE,KPPL,NUMPPL)
  IF (CTRLDONE) EXIT CONTROL
  IF (CTRLJOB(ICTRL) .EQ. 'FORWARD') THEN
    CALL GWM1RMS3PL(IPERT,NPERT,FIRSTSIM,LASTSIM)
  ENDIF
  IF (NUMPPL>0) THEN
    !
    ! **************************************************************************
    ! ************* GENERATE PARAMETER VALUES--Populate PARVALSETS *************
    ! **************************************************************************
    !
    DO KPPL=1,NUMPPL
      !   Populate one column of PARVALSETS with one set of values of
      !   adjustable parameters
      IF (CTRLJOB(ICTRL)=='FORWARD') THEN
        CALL BAS_GEN(NOPNT,NPT,NW,PRECISPRO,PVAL)
      ELSEIF (CTRLJOB(ICTRL)=='SENSITIVITY') THEN
        !   Populate a full set of parameter values, including perturbed parameter
        CALL MAIN_GEN(CTRLJOB(ICTRL),KPPL,NPE,NPT,NFVAR,IOUTMESS,DELINC,IPTR,  &
                      NOPNT,NW,PARNAM,FVBASE,KPE(KPPL),PARVALSETSE(:,KPPL:KPPL), &
                      OK)
        IF (.NOT. OK) GOTO 999
        !   Assign flows to cells where flow-rate decision variable KPPL controls
        !   flow at multiple cells
        CALL DECVARS_TO_PVAL(NFVAR,NPT,PARVALSETSE(:,KPPL:KPPL),  &
                             FVNCELL,PVALTMP)
        DO IL=1,NPTC
          PARVALSETS(IL,KPPL) = PVALTMP(IL)
        ENDDO
      ENDIF
    ENDDO
    !
    !   Establish iterative structure within potentially parallel loop
    IF (PARALLEL_ACTIVE .AND. NUMPPL>1) THEN
      DO_PARALLEL = .TRUE.
      !   Execute outer (PPL) loop once and inner loops multiple times
      NUMLOOP = 1
      NRUNSPLL = NUMPPL
    ELSE
      DO_PARALLEL = .FALSE.
      !   Execute outer loop multiple times and inner loops once
      NUMLOOP = NUMPPL
      NRUNSPLL = 1
    ENDIF
    !
    ! **************************************************************************
    ! ********************* TOP OF POTENTIALLY PARALLEL LOOP *******************
    ! **************************************************************************
    !
    PPL: DO KLOOP=1,NUMLOOP
      PPL1: DO KPPL=1,NRUNSPLL
        IF (DO_PARALLEL) THEN
          KLPTR = KPPL
        ELSE
          KLPTR = KLOOP
        ENDIF
        IF (CTRLJOB(ICTRL)=='FORWARD') THEN
          PVALTMP = PVAL
        ELSEIF (CTRLJOB(ICTRL)=='SENSITIVITY' ) THEN
          !   Populate PVALTMP with current parameter values, with substitution
          !   of adjustable parameters as needed for current iteration.
          IF (.NOT. DO_PARALLEL) THEN
            PVALTMP = RESHAPE(PARVALSETS(:,KLPTR:KLPTR),SHAPE(PVALTMP))
          ENDIF
        ENDIF
      ENDDO PPL1
      !
      IF (.NOT. DO_PARALLEL) THEN
        !-------BEGIN PERTURBATION LOOP FOR GWM PROCESS
        PERTDONE = .FALSE.
        DOPERT: DO WHILE (.NOT. PERTDONE)
          !
          ! **********************************************************************
          ! ************************ ADAPT PARAMETER VALUES **********************
          ! **********************************************************************
          !
          !---------PREPARE CALCULATION OF SIMULATION RESPONSE
          ! GWM1RMS3PP may adjust FVBASE and/or DELINC
          CALL GWM1RMS3PP(IOUTMESS,IPERT,FIRSTSIM,LASTSIM,NPERT,gwmstrg)
          CALL DECVARS_TO_PVAL(NFVAR,NPT,FVBASE,FVNCELL,PVALTMP)
          !
          CALL MIO_ADA_WRITEFILES(IFAIL,NPTC,PARNAM,NOPNT,NW,PRECISPRO,PVALTMP)
          IF (IFAIL .NE. 0) THEN
            OK = .FALSE.
            GOTO 999
          ENDIF
          !
          ! GWM: PREPARE FLOW PROCESS FILES FOR NEXT SIMULATION
          !
          ! **********************************************************************
          ! ********************** EXECUTE APPLICATION MODEL *********************
          ! **********************************************************************
          !
          IF (CTRLJOB(ICTRL) .EQ. 'FORWARD') THEN
            CALL BAS_EXE_SELECT(IOUTMESS,CTRLJOB(ICTRL),ICOMMAND)
          ENDIF
          IF (ICOMMAND>0) THEN
          !
            CALL BAS_EXE(ICOMMAND,-1,KLOOP,NUMPPL)
          ! 
          ENDIF
          !
          ! **********************************************************************
          ! ******************* EXTRACT VALUES FROM MODEL OUTPUT *****************
          ! **********************************************************************
          !
          IF (KPE(KLPTR)>0) THEN
            IF (IVERB .GT. 3) THEN
              ! Identify perturbed decision variable for extracted values
              WRITE(IOUTMESS,300) TRIM(PARNAM(IPTR(KPE(KLPTR))))
            ENDIF
          ENDIF
          !   Extract simulated equivalents to constraints
          CALL MIO_EXT(IFAIL,IOUTMESS,NCATS,NDEPSTAT,OBSNAM,MCVUSE, &
              DEPVALSETS(:,KLOOP:KLOOP),INSTRUCTION)
          IF (IFAIL .NE. 0) THEN
            WRITE(IOUTMESS,210) TRIM(AMESSAGE)
            IF (INSTRUCTION .NE. ' ') WRITE(IOUTMESS,210) TRIM(INSTRUCTION)
            OK = .FALSE.
            GOTO 999
          ENDIF
          ! Use extracted values to populate HDCSTATE, STCSTATE, and STASTATE arrays
          ! and run-status variables
          CALL MAIN_UEV(NDEPSTAT,DEPVALSETS(:,KLOOP:KLOOP),IPERT,NHB,NDD,NDF, &
                        NGD,NSF,NSD,NFVAR,DEWATERQ,MFCNVRG,HCLOSEG,HDRY,GWMOUT, &
                        OK)
          IF (.NOT. OK) GOTO 999
          IF (MFCNVRG(IPERT) .AND. .NOT. DEWATER(IPERT)   &
              .AND. .NOT. DEWATERQ(IPERT)) THEN
            PERTDONE = .TRUE.
          ENDIF
          IF (SOLNTYP.EQ.'FR') THEN
            LASTSIM = .TRUE.
          ENDIF
          CALL GWM1RMS3FP(IPERT,NPERT,FIRSTSIM,LASTSIM)
          IF (IPERT.GE.1 .AND. IPERT.LE.NPERT) THEN
            IF (IPGNA(IPERT)==2) PERTDONE = .FALSE.
          ENDIF
          IF (LASTSIM) PERTDONE = .TRUE.
        ENDDO DOPERT
        FIRSTSIM = .FALSE.       ! SIMULATION COMPLETED, NO LONGER FIRSTSIM
      ELSE
        !
        ! **********************************************************************
        ! *************** MAKE PERTURBATION RUNS IN PARALLEL *******************
        ! **********************************************************************
        !
        !   Assign command to be run in parallel
        CALL BAS_EXE_SELECT(IOUTMESS,'FORWARD',ICOMMAND,COMMANDPLL)
        !
        !   Adapt, execute, and extract in parallel
        CALL PLLGWM_MAKE_RUNS(IOUTMESS,LCIS,NCATS,NINSTRUCT, &
                NMIFILE,NMOFILE,NOPNT,NPT,   &
                NUMLADV,NRUNSPLL,NDEPSTAT,PRECISPRO,CATGOR,CINSTSET,OBSNAM,  &
                INSTRUCTFILE,LOCINS,MCVUSE,MIFILE,MOFILE,MRKDL,NW,PARNAM,   &
                PARVALSETS,TFILE,DEPVALSETS,HDRY,OK,COMMANDPLL,OPTIONFLAG)
      ENDIF
      !
      PPL2: DO KPPL=1,NRUNSPLL
        IF (DO_PARALLEL) THEN
          KLPTR = KPPL
        ELSE
          KLPTR = KLOOP
        ENDIF
        !
        IF (CTRLJOB(ICTRL) .EQ. 'FORWARD') THEN
          MODELVAL = RESHAPE(DEPVALSETS(:,KLPTR:KLPTR),(/NDEPSTAT/))
        ELSEIF (CTRLJOB(ICTRL) .EQ. 'SENSITIVITY') THEN
          !   Save simulated equivalents for model run with perturbed parameters
!          IF (KPE(KLPTR)>0) THEN
!            IF (IVERB .GT. 3) THEN
!              WRITE(IOUTMESS,300) TRIM(PARNAM(IPTR(KPE(KLPTR))))
!            ENDIF
!          ENDIF
        ENDIF
        !
        ! GWM: END OF GWM PERTURBATION LOOP
        !
        ! **********************************************************************
        ! ******************** END OF POTENTIALLY PARALLEL LOOP ****************
        ! **********************************************************************
        !
      ENDDO PPL2
    ENDDO PPL
  ENDIF
  ! for debugging
!  write(gwmout,*)CTRLJOB(ICTRL)//' contents of depvalsets:'
!  do jj=1,maxrunspll
!    write(gwmout,'(a,1x,i10)')'jj = ',jj
!    write(gwmout,'(8g15.8)')(depvalsets(ii,jj),ii=1,ndepstat)
!  enddo      
  !
  ! ****************************************************************************
  ! ************************* USE EXTRACTED VALUES *****************************
  ! ****************************************************************************
  !
  IF (CTRLJOB(ICTRL) .EQ. 'SENSITIVITY') THEN
      SENSDONE=.TRUE.
  ENDIF
  IF (CTRLJOB(ICTRL) .EQ. 'SOLVE_GWM') THEN
      CALL GWM1RMS3FM               ! FORMULATE THE GWM PROBLEM
      CALL GWM1RMS3AP(GWMCNVRG)     ! SOLVE THE GWM PROBLEM
      CALL DECVARS_TO_PVAL(NFVAR,NPT,CST,FVNCELL,PVAL) !   Assign updated flow values to PVAL
      IF(GWMCNVRG)THEN
        CALL GWM1RMS3OT(8,1)        ! WRITE GWM SOLUTION OUTPUT 
        LASTSIM = .TRUE.            ! EXECUTE ONE LAST SIMULATION
      ELSEIF(.NOT.GWMCNVRG)THEN
        CALL GWM1RMS3OT(-8,1)       ! WRITE GWM SOLUTION OUTPUT 
      ENDIF
  ELSEIF (CTRLJOB(ICTRL) .EQ. 'FORWARD') THEN
    IF (LASTSIM)THEN
      FINISH = .TRUE.
      IF (GWMWFILE(1).GT.0) THEN    ! Write the GWM managed files into a WEL file format
        IF (GWMWFILENAME .NE. '') THEN
          OPEN(GWMWFILE(1),FILE=GWMWFILENAME,STATUS='REPLACE')
          CALL GWM1RMS3OT(7,1)         
          CLOSE(GWMWFILE(1))      
        ELSE
          WRITE(*,400)GWMWFILE(1)
        ENDIF
      ENDIF
    ENDIF
  ENDIF
!
! GWM: END OF GWM ITERATION LOOP
!
! ******************************************************************************
! ************************** END OF CONTROL LOOP *******************************
! ******************************************************************************
!
ENDDO CONTROL
!
! ******************************************************************************
! ******************************* CLEAN UP *************************************
! ******************************************************************************
!
999 CONTINUE
CALL CLEAN_UP()

IF (GWMCNVRG .OR. SOLNTYP.EQ.'FR') THEN
  WRITE(*,*) ' Normal termination of GWM-VI'
ELSE
  WRITE(*,*) ' GWM-VI failed to converge'
ENDIF
!
IF (PARALLEL_ACTIVE) CALL PLLM_STOP_RUNNERS()
CALL BAS_CLN()
CALL DEP_CLN()
CALL PLLM_CLN()
CALL MAIN_CLN()
!
CALL UTL_ENDTIME(IBDT,GWMOUT)
CALL UTL_STOP(' ',GWMOUT)
END PROGRAM GWMVI
