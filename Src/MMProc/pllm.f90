MODULE PLLM
  USE DATATYPES
  USE GLOBAL_DATA, ONLY: MAX_STRING_LEN
  USE DATATYPES
  PRIVATE
  !
  !   Selected data are PUBLIC:
  PUBLIC :: AUTOSTOPRUNNERS, COMMAND, FNDISFIN, FNDISPAR, FNDISRDY, FNRUNDEP,  &
            FNRUNFAIL, FNRUNFIN, FNRUNRDY, LLPTRPLLCTRL, LLPTRPLLRUN,   &
            SNAME, WAIT, WAITRUNNERS, PARALLEL_ACTIVE_GWM, &
            OS_PLL, OS_TEMP, DOPLL, RUNRECFIL, NRUNNERCOLS, RUNNERCOLS, BTIME, &
            RUNNERSTAT, KOD, RUNTIME, RUNNERDIR, RUNNERNAME, IRUNNERSTAT, &
            IVERBRUNNER, TIMEOUTFAC, OS_DIS, RENAME_DIS, SRENAME, ERRRUNNER, &
            NUMRUNNERS, FLNM, FMTARG, ACCARG, FILACT, USESLEEP
  !
  !   PUBLIC data
  LOGICAL                       :: AUTOSTOPRUNNERS=.TRUE.
  LOGICAL                       :: PARALLEL_ACTIVE_GWM=.FALSE.
  CHARACTER(LEN=MAX_STRING_LEN) :: COMMAND = ' '
  !   Names of signal files
  CHARACTER(LEN=13), PARAMETER  :: FNDISFIN = 'jdispatch.fin'
  CHARACTER(LEN=11), PARAMETER  :: FNDISPAR = 'jdispar.rdy'
  CHARACTER(LEN=13), PARAMETER  :: FNDISRDY = 'jdispatch.rdy'
  CHARACTER(LEN=11), PARAMETER  :: FNRUNDEP = 'jrundep.rdy'
  CHARACTER(LEN=12), PARAMETER  :: FNRUNFAIL = 'jrunfail.fin'
  CHARACTER(LEN=11), PARAMETER  :: FNRUNFIN = 'jrunner.fin'
  CHARACTER(LEN=11), PARAMETER  :: FNRUNRDY = 'jrunner.rdy'
  !
  TYPE (LLIST),    POINTER      :: LLPTRPLLCTRL !  Parallel control data
  TYPE (LLIST),    POINTER      :: LLPTRPLLRUN  !  Parallel runners data
  CHARACTER(LEN=20)             :: SNAME = ' '
  DOUBLE PRECISION              :: WAIT = 0.001D0, WAITRUNNERS = 0.001D0
  !
  !   PRIVATE data
  INTEGER,                     ALLOCATABLE, DIMENSION(:,:) :: BTIME       ! Run begin times
  LOGICAL                                                  :: DOPLL
  CHARACTER(LEN=MAX_STRING_LEN)                            :: ERRRUNNER = ' '
  INTEGER,                       ALLOCATABLE, DIMENSION(:) :: IRUNNERSTAT, KOD
  INTEGER                                                  :: IVERBRUNNER = 3
  INTEGER, PARAMETER                                       :: NRUNNERCOLS = 3
  INTEGER                                                  :: NUMRUNNERS
  CHARACTER(LEN=20)                                        :: OS_DIS    ! Operating system on which dispatcher runs
  CHARACTER(LEN=20)                                        :: OS_PLL    ! Operating system of both dispatcher and runner
  CHARACTER(LEN=20)                                        :: OS_TEMP
  CHARACTER(LEN=6)                                         :: RENAME_DIS  ! "rename" command on dispatcher
  CHARACTER(LEN=40),                DIMENSION(NRUNNERCOLS) :: RUNNERCOLS
  CHARACTER(LEN=MAX_STRING_LEN), ALLOCATABLE, DIMENSION(:) :: RUNNERDIR
  CHARACTER(LEN=20),             ALLOCATABLE, DIMENSION(:) :: RUNNERNAME
  CHARACTER(LEN=MAX_STRING_LEN)                            :: RUNRECFIL
  DOUBLE PRECISION,              ALLOCATABLE, DIMENSION(:) :: RUNTIME
  CHARACTER(LEN=6)                                         :: SRENAME = 'ren' ! "rename" command on runner
  DOUBLE PRECISION                                         :: TIMEOUTFAC = 3.0D0
  CHARACTER(LEN=16), DIMENSION(-1:11) :: RUNNERSTAT
  !
  !   Added for GWM
  CHARACTER(LEN=200)::FLNM
  CHARACTER(LEN=20)::FILACT,FMTARG,ACCARG,FILSTAT
  ! Set USESLEEP = .FALSE. if compiler does not support SLEEPQQ subroutine
  LOGICAL :: USESLEEP = .TRUE.
  !
  DATA RUNNERSTAT/'Nonresponsive   ','Not active      ','Active          ',   &
        'Active          ','Running model   ','Active          ',   &
        'unused          ','unused          ','unused          ',   &
        'Notified to stop','Terminated      ','Stopped -- Error',   &
        'Overdue         '/
  !
  DATA RUNNERCOLS /'RUNNERNAME','RUNNERDIR ','RUNTIME  '/
  !
  !   IRUNNERSTAT flags:
  !     <0 -- Runner nonresponsive; absolute value indicates last
  !              status before nonresponse
  !      0 -- Initial status
  !      1 or greater -- Runner recognized as activated
  !      2 -- jdispatch.rdy file (FNDISRDY) copied to runner directory
  !      3 -- jdispar.rdy file (FNDISPAR) copied to runner directory;
  !              presumably, model is running
  !      4 -- jrundep.rdy file (FNRUNDEP) present in runner directory;
  !              model run done
  !      8 -- jdispatch.fin file (FNDISFIN) copied to runner directory;
  !              signal to stop j_runner executable
  !      9 -- jrunner.fin file (FNRUNFIN) found, indicating j_runner
  !              executable has stopped execution
  !     10 -- jrunfail.fin file (FNRUNFAIL) found, indicating j_runner
  !              executable has stopped execution due to error
  !     11 -- runner is overdue
  !
  !   Selected subprograms are PUBLIC:
  PUBLIC :: PLLM_INI_DISPATCHER, PLLM_INI_RUNNER_DIM, PLLM_INI_RUNNER_POP,  &
            PLLM_ADA, PLLM_EXE, PLLM_EXT, PLLM_CLN,   &
            PLLM_READ_DISPAR, PLLM_RUNNER_STOP,    &
            PLLM_WRITE_RUNDEP, PLLM_RMS2FP, PLLM_USTOP, PLLM_STOP_RUNNERS, &
            PLLM_WAIT, PLLM_DONE, PLLM_READ_RUNDEP, PLLM_WRITE_DISRDY, &
            PLLM_WRITE_DISPAR
  !
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_INI_DISPATCHER(INUNIT,IOUT,LCIS,NCATS,NINSTRUCT,NMIFILE, &
                            NMOFILE,NOPNT,NPT,NUMLADV,NVEXT,PRECISPRO, &
                            CATGOR,CINSTSET,EXTNAM,INSTRUCTFILE,LOCINS, &
                            MCVUSE,MIFILE,MOFILE,MRKDL,NW,PARNAM,TFILE, &
                            PARACTIVE,MAXRUNSPLL,MMPROCJTF_FILE, &
                            MNW2JTF_FILE,MNW2HERE,OK,OPTIONFLAG)
    !   Initialize PARALLEL module: Read parallel control and runner data and
    !   allocate and populate arrays.
    USE DATATYPES
    USE GLOBAL_DATA, ONLY: AMESSAGE, ERRSUB, IVERB, LENDNAM, MAX_STRING_LEN
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                                             INTENT(IN)  :: INUNIT
    INTEGER,                                             INTENT(IN)  :: IOUT
    INTEGER,                                             INTENT(IN)  :: LCIS
    INTEGER,                                             INTENT(IN)  :: NCATS
    INTEGER,                                             INTENT(IN)  :: NINSTRUCT
    INTEGER,                                             INTENT(IN)  :: NMIFILE
    INTEGER,                                             INTENT(IN)  :: NMOFILE
    INTEGER,                                             INTENT(IN)  :: NOPNT
    INTEGER,                                             INTENT(IN)  :: NPT
    INTEGER,                                             INTENT(IN)  :: NUMLADV
    INTEGER,                                             INTENT(IN)  :: NVEXT
    INTEGER,                                             INTENT(IN)  :: PRECISPRO
    CHARACTER(LEN=6),              DIMENSION(NMOFILE),   INTENT(IN)  :: CATGOR ! Model-calc'd val. cat. for each model-output file
    CHARACTER(LEN=1),              DIMENSION(LCIS),      INTENT(IN)  :: CINSTSET  ! holds compressed instruction set
    CHARACTER(LEN=LENDNAM),        DIMENSION(NVEXT),     INTENT(IN)  :: EXTNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE),   INTENT(IN)  :: INSTRUCTFILE ! Instruction file names
    INTEGER,                       DIMENSION(NINSTRUCT), INTENT(IN)  :: LOCINS    ! pointer to instructions
    LOGICAL,                       DIMENSION(NCATS),     INTENT(IN)  :: MCVUSE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE),   INTENT(IN)  :: MIFILE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE),   INTENT(IN)  :: MOFILE ! Model output file names
    CHARACTER(LEN=1),              DIMENSION(NMOFILE),   INTENT(IN)  :: MRKDL  ! Marker delimiters
    INTEGER,                       DIMENSION(NPT),       INTENT(IN)  :: NW
    CHARACTER(LEN=12),             DIMENSION(NPT),       INTENT(IN)  :: PARNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE),   INTENT(IN)  :: TFILE
    LOGICAL,                                             INTENT(OUT) :: PARACTIVE
    INTEGER,                                             INTENT(IN)  :: MAXRUNSPLL    !
    CHARACTER(LEN=*),                                    INTENT(IN)  :: MMPROCJTF_FILE
    CHARACTER(LEN=*),                                    INTENT(IN)  :: MNW2JTF_FILE
    LOGICAL,                                             INTENT(IN)  :: MNW2HERE
    LOGICAL,                                           INTENT(INOUT) :: OK
    INTEGER,                                   OPTIONAL, INTENT(IN)  :: OPTIONFLAG
    !   Local variables
    INTEGER :: I, IERR, ISTAT, IU, KACTIVE, KCTRL, KLOOPS, KTRY, LD, MORE
    CHARACTER(LEN=40), DIMENSION(0) :: COLNAMES
    CHARACTER(LEN=MAX_STRING_LEN)   :: FNAME, FRUNNERNAME
    TYPE (LLIST),      POINTER      :: TAIL
    LOGICAL :: LEX, LOP, OKLOCAL
    INTEGER :: OPTIONFLAGLOCAL
    !
    !   Format statements
    30 FORMAT(1X,'Runner "',A,'" is active')
    40 FORMAT(1X,'File ',A,' exists on runner ',A,   &
              ' but access to the file is prevented. ',/,   &
              ' Runner ',A,' is assumed to be inactive.')
    50 FORMAT(/,1X,I4,' runner(s) are active')
    60 FORMAT(/,1X,'*** ERROR: No runners are active')
    100 FORMAT(//,1X,'PARALLEL PROCESSING option is enabled.')
    150 FORMAT(1X,'Operating system under which dispatcher program runs: ',A)
    200 FORMAT(1X,'Number of runners = ',I4)
    250 FORMAT(1X,'WARNING: Unrecognized operating system: "',A,'"',/,   &
               1X,'OPERATINGSYSTEM defaults to: "',A,'"')
    300 FORMAT(/,1X,'Runner number ',I4,':  Name = ',A,/,   &
               1X,'Directory = ',A,/,   &
               1X,'Estimated runtime = ',G11.4,' seconds')
    350 FORMAT(/,1X,'ERROR: The following directory name does not end',   &
                  ' in a back slash or forward slash:',/,1X,A)
    400 FORMAT(A)
    !
    NULLIFY(LLPTRPLLCTRL,LLPTRPLLRUN,TAIL)
    !   Assign defaults
    OS_PLL = ' '
    OS_TEMP = 'WINDOWS'
    PARACTIVE = .FALSE.
    DOPLL = .FALSE.
    KCTRL = 0
    NUMRUNNERS = 0
    RUNRECFIL = ' '
    WAIT = 0.001D0
    WAITRUNNERS = 0.001D0
    OKLOCAL = .TRUE.
    !
    IF (PRESENT(OPTIONFLAG)) THEN
      OPTIONFLAGLOCAL = OPTIONFLAG
    ELSE
      OPTIONFLAGLOCAL = 0
    ENDIF
    !
    CALL UTL_READBLOCK(0,'PARALLEL_CONTROL',COLNAMES,INUNIT,IOUT,'*',   &
                       .FALSE.,LLPTRPLLCTRL,TAIL,KCTRL)
    IF (KCTRL>0) THEN
      CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'PARALLEL',DOPLL)
    ENDIF
    CALL UTL_READBLOCK(NRUNNERCOLS,'PARALLEL_RUNNERS',RUNNERCOLS,INUNIT,   &
                       IOUT,'RUNNERNAME',DOPLL,LLPTRPLLRUN,TAIL,NUMRUNNERS)
    IF (KCTRL>0) THEN
      PARACTIVE = DOPLL
      IF (DOPLL) THEN
        IF (NUMRUNNERS==0) THEN
          AMESSAGE = 'PARALLEL is TRUE, but zero runners have been defined'
          CALL UTL_SUBERROR(ERRSUB)
        ENDIF
        ALLOCATE(BTIME(8,NUMRUNNERS),IRUNNERSTAT(NUMRUNNERS),   &
                 KOD(NUMRUNNERS), RUNTIME(NUMRUNNERS),   &
                 RUNNERDIR(NUMRUNNERS),RUNNERNAME(NUMRUNNERS))
        !
        !   Populate arrays with defaults as appropriate
        IRUNNERSTAT = 0
        KOD = 0
        RUNTIME = 10.0D0
        !
        !   Populate arrays and assign variables
        CALL UTL_FILTERLIST(LLPTRPLLRUN,IOUT,'RUNNERNAME',NUMRUNNERS,IERR,   &
                            RUNNERNAME,MORE)
        CALL UTL_FILTERLIST(LLPTRPLLRUN,IOUT,'RUNNERDIR',NUMRUNNERS,IERR,   &
                            RUNNERDIR,MORE)
        CALL UTL_FILTERLIST(LLPTRPLLRUN,IOUT,'RUNTIME',NUMRUNNERS,IERR,   &
                            RUNTIME,MORE)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'WAIT',WAIT)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'WAITRUNNERS',WAITRUNNERS)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'RUNRECORD',RUNRECFIL)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'VERBOSERUNNER',IVERBRUNNER)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'AUTOSTOPRUNNERS',   &
                        AUTOSTOPRUNNERS)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'OSDISPATCHER',OS_TEMP)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'OPERATINGSYSTEM',OS_PLL)
        CALL UTL_FILTER(IERR,LLPTRPLLCTRL,IOUT,'TIMEOUTFACTOR',TIMEOUTFAC)
        !
        IF (OS_PLL .NE. ' ') OS_TEMP = OS_PLL
        CALL UTL_CASETRANS(OS_TEMP,'hi')
        OS_DIS = OS_TEMP
        !   Echo input
        IF (IVERB>2) THEN
          WRITE(IOUT,100)
          WRITE(IOUT,150)TRIM(OS_DIS)
          WRITE(IOUT,200)NUMRUNNERS
        ENDIF
        !
        !   Assign variables dependent on operating system of dispatcher program
        IF (OS_TEMP == 'WINDOWS') THEN
          RENAME_DIS = 'ren'
        ELSEIF (OS_TEMP == 'DOS') THEN
          OS_DIS = 'WINDOWS'
          RENAME_DIS = 'ren'
        ELSEIF (OS_TEMP == 'UNIX') THEN
          RENAME_DIS = 'mv'
        ELSEIF (OS_TEMP == 'LINUX') THEN
          OS_DIS = 'UNIX'
          RENAME_DIS = 'mv'
        ELSE
          OS_DIS = 'WINDOWS'
          IF (IVERB>0) WRITE(IOUT,250)TRIM(OS_TEMP),TRIM(OS_DIS)
          RENAME_DIS = 'ren'
        ENDIF

        DO I=1,NUMRUNNERS
          LD = LEN_TRIM(RUNNERDIR(I))
          IF (IVERB>2) WRITE(IOUT,300)I,TRIM(RUNNERNAME(I)),   &
              TRIM(RUNNERDIR(I)),RUNTIME(I)
          IF (RUNNERDIR(I)(LD:LD) .NE. '/' .AND.   &
              RUNNERDIR(I)(LD:LD) .NE. '\') THEN
            WRITE(IOUT,350)TRIM(RUNNERDIR(I))
            CALL PLLM_STOP_RUNNERS()
          ENDIF
          !   Delete jrunner.rdy file
          FRUNNERNAME = TRIM(RUNNERDIR(I))//FNRUNRDY
          INQUIRE(FILE=FRUNNERNAME,EXIST=LEX)
          IF (LEX) THEN
            IU = UTL_GETUNIT(7,1000)
            OPEN(IU,FILE=FRUNNERNAME,ERR=310,IOSTAT=ISTAT)
            310 CONTINUE
            IF (ISTAT==0) THEN
              ISTAT = PLL_CLOSE_AND_DELETE(IU,FRUNNERNAME)
            ENDIF
          ENDIF
        ENDDO
        IF (IVERB>2) WRITE(IOUT,400)' '
        !
        KACTIVE = 0
        KLOOPS = 0
        340 CONTINUE
        CALL PLLM_WAIT(1.0D0)   ! Allow time for runners to recreate FNRUNRDY file
        !   Determine status of runners
        DO I=1,NUMRUNNERS
          IF (IRUNNERSTAT(I)==0) THEN
            FRUNNERNAME = TRIM(RUNNERDIR(I))//FNRUNRDY
            INQUIRE(FILE=FRUNNERNAME,EXIST=LEX)
            IF (LEX) THEN
              IU = UTL_GETUNIT(7,1000)
              KTRY = 0
              405 CONTINUE
              ISTAT = 0
              INQUIRE(UNIT=IU,OPENED=LOP)
              IF (.NOT. LOP) OPEN(IU,FILE=FRUNNERNAME,ERR=410,IOSTAT=ISTAT)
              410 CONTINUE
              IF (ISTAT==0) THEN
                ISTAT = PLL_CLOSE_AND_DELETE(IU,FRUNNERNAME)
              ENDIF
              IF (ISTAT==0) THEN
                IRUNNERSTAT(I) = 1
                KOD(I) = 0
                IF (IVERB>2) WRITE(*,30)TRIM(RUNNERNAME(I))
                KACTIVE = KACTIVE+1
              ELSE
                IF (KTRY > 10) THEN
                  WRITE(*,40)TRIM(FNRUNRDY),TRIM(RUNNERNAME(I)),TRIM(RUNNERNAME(I))
                  IRUNNERSTAT(I) = -1
                ELSE
                  CALL PLLM_WAIT()
                  KTRY = KTRY + 1
                  GOTO 405
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF (KACTIVE==0 .AND. KLOOPS<20) THEN
          KLOOPS = KLOOPS+1
          GOTO 340
        ENDIF
        !
        IF (IVERB>2) WRITE(*,50)KACTIVE
        IF (KACTIVE==0) THEN
          WRITE(*,60)
          WRITE(IOUT,60)
          IF (IVERB>3) THEN
            WRITE(*,*)' PLLM_INI_DISPATCHER calling PLLM_STOP_RUNNERS'
            WRITE(IOUT,*)' PLLM_INI_DISPATCHER calling PLLM_STOP_RUNNERS'
          ENDIF
          CALL PLLM_STOP_RUNNERS()
          OKLOCAL = .FALSE.
          GOTO 999
        ENDIF
        !
        ! Copy template file mmproc.in.jtf to all runner directories
        CALL COPY_TO_RUNNERS(MMPROCJTF_FILE,OKLOCAL)
        IF (.NOT. OKLOCAL) GOTO 999
        IF (MNW2HERE) THEN
          ! Copy template file mnw2_input.jtf to all runner directories
          CALL COPY_TO_RUNNERS(MNW2JTF_FILE,OKLOCAL)
          IF (.NOT. OKLOCAL) GOTO 999
        ENDIF
        !
        !   Create jdispatch.rdy files to pass model-interaction data to runners
        DO I=1,NUMRUNNERS
          IF (IRUNNERSTAT(I)==1) THEN
            !   Create "jdispatch.rdy" file for this runner
            CALL PLLM_WRITE_DISRDY(I,LCIS,NCATS,NINSTRUCT,NMIFILE,NMOFILE,   &
                                  NOPNT,NPT,NUMLADV,NVEXT,PRECISPRO,   &
                                  CATGOR,CINSTSET,EXTNAM,INSTRUCTFILE,   &
                                  LOCINS,MCVUSE,MIFILE,MOFILE,MRKDL,NW,   &
                                  PARNAM,TFILE,OPTIONFLAGLOCAL)
          ENDIF
        ENDDO
        !
        FNAME = 'jparallel.interrupt'
        INQUIRE(FILE=FNAME,EXIST=LEX)
        IF (LEX) THEN
          IU = UTL_GETUNIT(7,1000)
          OPEN(IU,FILE=FNAME)
          CLOSE(IU,STATUS='DELETE',IOSTAT=ISTAT)
        ENDIF
      ENDIF
    ENDIF
    PARALLEL_ACTIVE_GWM = PARACTIVE
    !
    999 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE PLLM_INI_DISPATCHER
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_INI_RUNNER_DIM(IOUNIT,RESET,LCIS,NCATS,NINSTRUCT,NMIFILE,   &
                               NMOFILE,NOPNT,NPT,NUMLADV,NVEXT,   &
                               PRECISPRO,INIMESS)
    !   Look for jdispatch.rdy file.  When found, read scalar values from it.
    USE GLOBAL_DATA, ONLY: AMESSAGE, IVERB
    USE UTILITIES
    USE MODEL_IO, ONLY: MIO_INI_ALLOC
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                    INTENT(IN)    :: IOUNIT
    LOGICAL,                    INTENT(INOUT) :: RESET
    INTEGER,                    INTENT(OUT)   :: LCIS
    INTEGER,                    INTENT(OUT)   :: NCATS
    INTEGER,                    INTENT(OUT)   :: NINSTRUCT
    INTEGER,                    INTENT(OUT)   :: NMIFILE ! Number of model-input files
    INTEGER,                    INTENT(OUT)   :: NMOFILE ! Number of model-output files
    INTEGER,                    INTENT(OUT)   :: NOPNT
    INTEGER,                    INTENT(OUT)   :: NPT
    INTEGER,                    INTENT(OUT)   :: NUMLADV
    INTEGER,                    INTENT(OUT)   :: NVEXT   ! Number of dependent values to extract
    INTEGER,                    INTENT(OUT)   :: PRECISPRO
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)    :: INIMESS
    !
    !   Local variables
    INTEGER :: IFAIL, ISTAT, ISTAT1, ISTAT2, ISTAT3, K, KFND, KTRY, KTRY1, IU
    LOGICAL :: LEX
    CHARACTER(LEN=5) :: FINMESS
    CHARACTER(LEN=MAX_STRING_LEN) :: FILERR
    CHARACTER(LEN=20) :: VARERR = ' '
    !
    !   Format statements
    100 FORMAT(1X,A,/)
    150 FORMAT(A)
    170 FORMAT(1X,'File "',A,'" has been created.')
    200 FORMAT(F20.0)
    250 FORMAT(I10)
    300 FORMAT(1X,'File "',A,'" found...')
    320 FORMAT(1X,'VERBOSERUNNER = ',I2)
    340 FORMAT(1X,'This is runner "',A,'".  The rename command is "',A,'"')
    700 FORMAT(1X,'Error reading from file "',A,'"')
    720 FORMAT(1X,'Error reading variable "',A,'" from file "',A,'"')
    740 FORMAT(1X,'Closing file and making attempt: ',I2)
    800 FORMAT(/,1X,'Runner "',A,'" received signal to stop execution')
    810 FORMAT(/,1X,'Runner received signal to stop execution')
    820 FORMAT(/,1X,'Runner "',A,   &
        '" resetting to continue execution (PLLM_INI_RUNNER_DIM)',/)
    830 FORMAT(/,1X,   &
        'Runner resetting to continue execution (PLLM_INI_RUNNER_DIM)',/)
    840 FORMAT(/,1X,'Runner cannot open file ',A)
    850 FORMAT(/,1X,'Runner cannot close file ',A)
    !
    IF (IVERB>2) THEN
      IF (PRESENT(INIMESS)) WRITE(*,100)TRIM(INIMESS)
    ENDIF
    !
    !   Write FNRUNRDY file to signal that runner is active
    KTRY = 0
    3 CONTINUE
    OPEN(IOUNIT,FILE=FNRUNRDY,STATUS='REPLACE',IOSTAT=ISTAT)
    IF (ISTAT .NE. 0) THEN
      KTRY = KTRY+1
      CALL PLLM_WAIT(0.05D0)
      IF (KTRY<10) GOTO 3
      WRITE(*,840)TRIM(FNRUNRDY)
      RESET = .TRUE.
      RETURN
    ENDIF
    IF (IVERB>4) WRITE(*,170)TRIM(FNRUNRDY)
    !
    KTRY = 0
    6 CONTINUE
    CLOSE(IOUNIT,STATUS='KEEP',IOSTAT=ISTAT)
    IF (ISTAT .NE. 0) THEN
      KTRY = KTRY+1
      CALL PLLM_WAIT(0.05D0)
      IF (KTRY<10) GOTO 6
      WRITE(*,850)TRIM(FNRUNRDY)
      RESET = .TRUE.
      RETURN
    ENDIF
    !
    !   Look for FNDISRDY here, and initialize arrays
    !
    K = 0
    KFND = 0
    10 CONTINUE
    INILOOP: DO WHILE (.TRUE.)
      !   Ensure that FNRUNRDY file exists, to show runner is alive
      INQUIRE(FILE=FNRUNRDY,EXIST=LEX)
      IF (.NOT. LEX) THEN
        IU = UTL_GETUNIT(7,1000)
        KTRY1 = 0
        20 CONTINUE
        OPEN(IU,FILE=FNRUNRDY,IOSTAT=ISTAT)
        IF (ISTAT>0 .AND. KTRY1<11) THEN
          KTRY1 = KTRY1+1
          CALL PLLM_WAIT()
          GOTO 20
        ENDIF
        KTRY1 = 0
        30 CONTINUE
        CLOSE(IU,STATUS='KEEP',IOSTAT=ISTAT)
        IF (ISTAT>0 .AND. KTRY1<11) THEN
          KTRY1 = KTRY1+1
          CALL PLLM_WAIT()
          GOTO 30
        ENDIF
        IF (IVERB>4) WRITE(*,170)TRIM(FNRUNRDY)
      ENDIF
      !   Look for FNDISRDY file, indicating dispatcher is in communication
      !   with runner
      INQUIRE(FILE=FNDISRDY,EXIST=LEX)
      IF (LEX) THEN
        IF (IVERB>2 .AND. KFND==0) WRITE(*,300)FNDISRDY
        !   Read initialization data provided by dispatcher in jdispatch.rdy
        OPEN(IOUNIT,FILE=FNDISRDY,STATUS='OLD',IOSTAT=ISTAT)
        IF (ISTAT .NE. 0) THEN
          KFND = KFND+1
          GOTO 50
        ENDIF
        FILERR=FNDISRDY
        VARERR='SNAME'
        READ(IOUNIT,150,ERR=990,END=990)SNAME
        VARERR='SRENAME'
        READ(IOUNIT,150,ERR=990,END=990)SRENAME
        VARERR='WAIT'
        READ(IOUNIT,200,ERR=990,END=990)WAIT
        VARERR='IVERB'
        READ(IOUNIT,250,ERR=990,END=990)IVERB
        VARERR='LCIS'
        READ(IOUNIT,250,ERR=990,END=990)LCIS
        VARERR='NINSTRUCT'
        READ(IOUNIT,250,ERR=990,END=990)NINSTRUCT
        VARERR='NMIFILE'
        READ(IOUNIT,250,ERR=990,END=990)NMIFILE
        VARERR='NMOFILE'
        READ(IOUNIT,250,ERR=990,END=990)NMOFILE
        VARERR='NPT'
        READ(IOUNIT,250,ERR=990,END=990)NPT
        VARERR='NOPNT'
        READ(IOUNIT,250,ERR=990,END=990)NOPNT
        VARERR='PRECISPRO'
        READ(IOUNIT,250,ERR=990,END=990)PRECISPRO
        VARERR='NCATS'
        READ(IOUNIT,250,ERR=990,END=990)NCATS
        VARERR='NVEXT'
        READ(IOUNIT,250,ERR=990,END=990)NVEXT
        VARERR='NUMLADV'
        READ(IOUNIT,250,ERR=990,END=990)NUMLADV
        !
        IF (IVERB>2) WRITE(*,320) IVERB
        IF (IVERB>3) WRITE(*,340) TRIM(SNAME),TRIM(SRENAME)
        !
        !   Initialize the MIO module for model-input files
        CALL MIO_INI_ALLOC(IFAIL,NPT)
        IF (IFAIL .NE. 0) THEN
          ERRRUNNER = 'Error: runner cannot initialize the Model_IO module.'
          WRITE(*,*) TRIM(ERRRUNNER)
          CALL PLLM_RUNNER_STOP(1)
        ENDIF
        !
        RETURN
      ENDIF
      50 CONTINUE
      !
      !   Look for FNDISFIN file, indicating runner
      !   should stop execution or reset
      INQUIRE(FILE=FNDISFIN,EXIST=LEX)
      IF (LEX) THEN
        IF (IVERB>2)WRITE(*,300)TRIM(FNDISFIN)
        KTRY = 0
        560 CONTINUE
        OPEN(IOUNIT,FILE=FNDISFIN,IOSTAT=ISTAT1)
        IF (ISTAT1 .NE. 0) THEN
          KTRY = KTRY+1
          IF (KTRY<10) THEN
            CALL PLLM_WAIT()
            GOTO 560
          ELSE
            FINMESS=' '
            WRITE(*,'(A)')' Warning: WAIT time may be too small'
          ENDIF
        ELSE
          READ(IOUNIT,150,END=590,ERR=590,IOSTAT=ISTAT2)FINMESS
          590 CONTINUE
          KTRY = 0
          595 CONTINUE
          CLOSE(IOUNIT,STATUS='KEEP',IOSTAT=ISTAT3)
          IF (ISTAT3 .NE. 0) THEN
            KTRY = KTRY+1
            IF (KTRY<10) THEN
              CALL PLLM_WAIT()
              GOTO 595
            ELSE
              WRITE(*,'(A)')' Warning: WAIT time may be too small'
            ENDIF
          ENDIF
          !
          IF (FINMESS=='RESET') THEN
            RESET=.TRUE.
            IF (SNAME==' ') THEN
              WRITE(*,830)
            ELSE
              WRITE(*,820)TRIM(SNAME)
            ENDIF
            !
            RETURN
            !
          ELSE
            IF (SNAME==' ') THEN
              WRITE(*,810)
            ELSE
              WRITE(*,800)TRIM(SNAME)
            ENDIF
          ENDIF
          CALL PLLM_RUNNER_STOP()
        ENDIF
      ENDIF
    ENDDO INILOOP
    !
    990 CONTINUE
    WRITE(*,700)TRIM(FILERR)
    IF (VARERR .NE. ' ') THEN
      WRITE(*,720)TRIM(VARERR),TRIM(FNDISRDY)
      ERRRUNNER = ' '
      WRITE(ERRRUNNER,720)TRIM(VARERR),TRIM(FNDISRDY)
      K = K+1
      IF (K<11) THEN
        CLOSE(IOUNIT,IOSTAT=ISTAT)
        CALL PLLM_WAIT(WAIT)
        WRITE(*,740)K
        GOTO 10
      ELSE
        WRITE(*,'(A)')' Warning: WAIT time may be too small'
      ENDIF
    ENDIF
    AMESSAGE = 'Error in subroutine PLLM_INI_RUNNER_DIM'
    CALL PLLM_RUNNER_STOP(1)
    !
    RETURN
  END SUBROUTINE PLLM_INI_RUNNER_DIM
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_INI_RUNNER_POP(IOUNIT,LCIS,NCATS,NINSTRUCT,NMIFILE,NMOFILE,   &
                               NPT,NUMLADV,NVEXT,CATGOR,CINSTSET,EXTNAM,   &
                               INSTRUCTFILE,LOCINS,MCVUSE,MIFILE,   &
                               MOFILE,MRKDL,NW,PARNAM,TFILE)
    !   Read array values from jdispatch.rdy file and populate arrays,
    !   then close and delete the file
    USE GLOBAL_DATA, ONLY: AMESSAGE, IVERB, LENDNAM, MAX_STRING_LEN
    USE MODEL_IO, ONLY: MIO_INI_INPUTFILES_RUNNER, MIO_INI_INSTRUCT_RUNNER,   &
                        MIO_INI_OUTPUTFILES_RUNNER, MIO_INI_TEMPLATE
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                                            INTENT(IN)  :: IOUNIT
    INTEGER,                                            INTENT(IN)  :: LCIS
    INTEGER,                                            INTENT(IN)  :: NCATS
    INTEGER,                                            INTENT(IN)  :: NINSTRUCT
    INTEGER,                                            INTENT(IN)  :: NMIFILE
    INTEGER,                                            INTENT(IN)  :: NMOFILE
    INTEGER,                                            INTENT(IN)  :: NPT
    INTEGER,                                            INTENT(IN)  :: NUMLADV
    INTEGER,                                            INTENT(IN)  :: NVEXT
    CHARACTER(LEN=6),              DIMENSION(NMOFILE),  INTENT(OUT) :: CATGOR ! Model-calc'd vaL. cat. for each model-output file
    CHARACTER(LEN=1),              DIMENSION(LCIS),     INTENT(OUT) :: CINSTSET ! holds compressed instruction set
    CHARACTER(LEN=LENDNAM),        DIMENSION(NVEXT),    INTENT(OUT) :: EXTNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE),  INTENT(OUT) :: INSTRUCTFILE ! Instruction file names
    INTEGER,                       DIMENSION(NINSTRUCT), INTENT(OUT) :: LOCINS ! pointer to instructions
    LOGICAL,                       DIMENSION(NCATS),    INTENT(OUT) :: MCVUSE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE),  INTENT(OUT) :: MIFILE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE),  INTENT(OUT) :: MOFILE ! Model output file names
    CHARACTER(LEN=1),              DIMENSION(NMOFILE),  INTENT(OUT) :: MRKDL  ! Marker delimiters
    INTEGER,                       DIMENSION(NPT),      INTENT(OUT) :: NW
    CHARACTER(LEN=12),             DIMENSION(NPT),      INTENT(OUT) :: PARNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE),  INTENT(OUT) :: TFILE
    !
    !   Local variables
    INTEGER :: I, IFAIL, ISTAT, K
    !
    !   Format statements
    100 FORMAT(A)
    200 FORMAT(6(1X,A12))
    300 FORMAT(6(1X,I4,8X))
    360 FORMAT(3(1X,A))
    380 FORMAT(10(1X,L1))
    420 FORMAT(10(1X,A6))
    440 FORMAT(40(1X,A1))
    450 FORMAT(1X,'Deleting file "',A,'"')
    460 FORMAT(80A1)
    480 FORMAT(8I10)
    500 FORMAT(1X,'Error reading from file "',A,'"')
    !
    READ(IOUNIT,100,ERR=990,END=990)(MIFILE(I),I=1,NMIFILE)
    READ(IOUNIT,100,ERR=990,END=990)(TFILE(I),I=1,NMIFILE)
    READ(IOUNIT,200,ERR=990,END=990)(PARNAM(I),I=1,NPT)
    READ(IOUNIT,300,ERR=990,END=990)(NW(I),I=1,NPT)
    READ(IOUNIT,420,ERR=990,END=990)(CATGOR(I),I=1,NMOFILE)
    READ(IOUNIT,360,ERR=990,END=990)(EXTNAM(I),I=1,NVEXT)
    READ(IOUNIT,100,ERR=990,END=990)(INSTRUCTFILE(I),I=1,NMOFILE)
    READ(IOUNIT,100,ERR=990,END=990)(MOFILE(I),I=1,NMOFILE)
    READ(IOUNIT,440,ERR=990,END=990)(MRKDL(I),I=1,NMOFILE)
    READ(IOUNIT,380,ERR=990,END=990)(MCVUSE(I),I=1,NCATS)
    READ(IOUNIT,460,ERR=990,END=990)(CINSTSET(I),I=1,LCIS)
    READ(IOUNIT,480,ERR=990,END=990)(LOCINS(I),I=1,NINSTRUCT)
    IF(IVERB>3)WRITE(*,450)FNDISRDY
    K = 0
    50 CONTINUE
    CLOSE(IOUNIT,STATUS='DELETE',IOSTAT=ISTAT)
    IF (ISTAT .NE. 0) THEN
      K = K+1
      IF (K<11) THEN
        CALL PLLM_WAIT()
        GOTO 50
      ENDIF
    ENDIF
    !
    !   Initialize model-input and template file arrays
    CALL MIO_INI_INPUTFILES_RUNNER(IFAIL,NMIFILE,MIFILE,TFILE)
    IF (IFAIL .NE. 0) THEN
      ERRRUNNER = 'Error: Runner cannot initialize model-input file data.'
      CALL PLLM_RUNNER_STOP(1)
    ENDIF
    CALL MIO_INI_TEMPLATE(IFAIL,NPT,PARNAM,NW)
    IF (IFAIL .NE. 0) THEN
      ERRRUNNER = 'Error: Runner cannot initialize template file data.'
      CALL PLLM_RUNNER_STOP(1)
    ENDIF
    !
    !   Initialize model-output and instruction file arrays
    CALL MIO_INI_OUTPUTFILES_RUNNER(IFAIL,NMOFILE,CATGOR,INSTRUCTFILE,MOFILE,MRKDL)
    IF (IFAIL .NE. 0) THEN
      ERRRUNNER = 'Error: Runner cannot initialize model-output file data.'
      CALL PLLM_RUNNER_STOP(1)
    ENDIF
    CALL MIO_INI_INSTRUCT_RUNNER(IFAIL,LCIS,NINSTRUCT,NUMLADV,NVEXT,CINSTSET,   &
                                LOCINS)
    IF (IFAIL .NE. 0) THEN
      ERRRUNNER = 'Error: Runner cannot initialize instruction data.'
      CALL PLLM_RUNNER_STOP(1)
    ENDIF
    !
    RETURN
    !
    990 CONTINUE
    WRITE(*,500)FNDISRDY
    ERRRUNNER = ' '
    WRITE(ERRRUNNER,500)FNDISRDY
    AMESSAGE = ' '
    WRITE(AMESSAGE,500)FNDISRDY
    CALL PLLM_RUNNER_STOP(1)
  END SUBROUTINE PLLM_INI_RUNNER_POP
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_ADA(NPT,NOPNT,PARNAM,PRECISPRO,NW,PVAL)
    !   Write model-input file(s).  This subroutine is intended to be
    !   called by a runner.
    USE GLOBAL_DATA, ONLY: AMESSAGE, IVERB
    USE MODEL_IO, ONLY: MIO_ADA_WRITEFILES
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                                 INTENT(IN)    :: NPT     ! Number of parameters
    INTEGER,                                 INTENT(IN)    :: NOPNT   ! Decimal point protocol
    CHARACTER(LEN=12),   DIMENSION(NPT),     INTENT(IN)    :: PARNAM    ! Parameter names
    INTEGER,                                 INTENT(IN)    :: PRECISPRO ! Precision protocol
    INTEGER,             DIMENSION(NPT),     INTENT(INOUT) :: NW    ! Minimum word length of a parameter
    DOUBLE PRECISION,    DIMENSION(NPT),     INTENT(INOUT) :: PVAL  ! parameter values
    !
    !   Local variables
    INTEGER :: IFAIL
    !
    !   Write the model-input files
    CALL MIO_ADA_WRITEFILES(IFAIL,NPT,PARNAM,NOPNT,NW,PRECISPRO,PVAL)
    IF (IFAIL .NE. 0) THEN
      ERRRUNNER = 'Error: Runner cannot write model-input file(s).'
      CALL PLLM_RUNNER_STOP(1)
    ENDIF
    RETURN
  END SUBROUTINE PLLM_ADA
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_EXE(IRUNRUNNER,NRUNRUNNER)
    USE GLOBAL_DATA, ONLY: IVERB
    USE UTILITIES
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER, INTENT(IN) :: IRUNRUNNER
    INTEGER, INTENT(IN) :: NRUNRUNNER
    !
    !   Local variables
    !
    !   Format statements
    100 FORMAT(/,1X,'Running system command: ',1X,A,/)
    110 FORMAT(1X,'(Run ',I5,' of ',I5,')')
    120 FORMAT(1X,'Finished system command:',1X,A)
    !
    IF (IVERB>3) WRITE(*,100) TRIM(COMMAND)
    IF (IVERB>2) WRITE(*,110) IRUNRUNNER,NRUNRUNNER
    !
    CALL UTL_SYSTEM(TRIM(COMMAND))
    !
    IF (IVERB>3) WRITE(*,120) TRIM(COMMAND)
    !
    RETURN
  END SUBROUTINE PLLM_EXE
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_EXT(NCATS,NVEXT,EXTNAM,MCVUSE,EXTVAL)
    !   Extract model-calculated dependents
    USE GLOBAL_DATA, ONLY: MAX_STRING_LEN
    USE MODEL_IO, ONLY: MIO_EXT
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                             INTENT(IN)  :: NCATS  ! Number of categories of model-calculated values
    INTEGER,                             INTENT(IN)  :: NVEXT  ! number of values to extract
    CHARACTER(LEN=*), DIMENSION(NVEXT),  INTENT(IN)  :: EXTNAM ! extracted-value names
    LOGICAL,          DIMENSION(NCATS),  INTENT(IN)  :: MCVUSE ! Do extractions for this category?
    DOUBLE PRECISION, DIMENSION(NVEXT),  INTENT(OUT) :: EXTVAL ! extracted values
    !
    !   Local variables
    INTEGER :: IFAIL
    CHARACTER(LEN=MAX_STRING_LEN) :: INSTRUCTION ! instruction of error
    !
    !   Format statements
    100 FORMAT(/,1X,'Error encountered in extraction instruction:',/,1X,A)
    !
    CALL MIO_EXT(IFAIL,-1,NCATS,NVEXT,EXTNAM,MCVUSE,EXTVAL,INSTRUCTION)
    IF (IFAIL .NE. 0) THEN
      WRITE(*,100) TRIM(INSTRUCTION)
      ERRRUNNER = 'Error encountered in extraction instruction: '//  &
          TRIM(INSTRUCTION)
      CALL PLLM_RUNNER_STOP(1)
    ENDIF
    !
    RETURN
  END SUBROUTINE PLLM_EXT
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_CLN()
    !   Deallocate all arrays in the Parallel module
    USE GLOBAL_DATA, ONLY: IVERB
    IMPLICIT NONE
    !
    IF (ALLOCATED(BTIME)) DEALLOCATE(BTIME)
    IF (ALLOCATED(IRUNNERSTAT)) DEALLOCATE(IRUNNERSTAT)
    IF (ALLOCATED(KOD)) DEALLOCATE(KOD)
    IF (ALLOCATED(RUNTIME)) DEALLOCATE(RUNTIME)
    IF (ALLOCATED(RUNNERDIR)) DEALLOCATE(RUNNERDIR)
    IF (ALLOCATED(RUNNERNAME)) DEALLOCATE(RUNNERNAME)
    CALL TYP_DEALLOC(LLPTRPLLCTRL)
    CALL TYP_DEALLOC(LLPTRPLLRUN)
    RETURN
  END SUBROUTINE PLLM_CLN
  !-----------------------------------------------------------------------------
  !----------------------- PARALLEL-PROCESSING UTILITIES -----------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_DONE(NRUNSPLL,IRUNSTAT,DONE,NUMDONE)
    !   Return TRUE only if all elements of IRUNSTAT >= 3.
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                      INTENT(IN)  :: NRUNSPLL
    INTEGER, DIMENSION(NRUNSPLL), INTENT(IN)  :: IRUNSTAT
    LOGICAL,                      INTENT(OUT) :: DONE
    INTEGER,                      INTENT(OUT) :: NUMDONE
    !
    !   Local variables
    INTEGER :: I
    !
    NUMDONE = 0
    DONE = .FALSE.
    DO I=1,NRUNSPLL
      IF (IRUNSTAT(I)==3) NUMDONE = NUMDONE+1
    ENDDO
    IF (NUMDONE==NRUNSPLL) DONE = .TRUE.
    RETURN
  END SUBROUTINE PLLM_DONE
  !-----------------------------------------------------------------------------
  INTEGER FUNCTION PLLM_FIND_RUNNER() RESULT (IDRUNNER)
    !   Return the number of a runner that is available to make a model run.
    !   Return -1 if no runners are available.
    IMPLICIT NONE
    !
    !   Local variables
    INTEGER :: I
    CHARACTER(LEN=MAX_STRING_LEN) :: FNAME
    LOGICAL :: LEX
    !
    IDRUNNER = -1
    DO I=1,NUMRUNNERS
      FNAME=TRIM(RUNNERDIR(I))//FNRUNRDY
      INQUIRE(FILE=FNAME,EXIST=LEX)
      IF (LEX) THEN
        !   Double check that file still exists
        CALL PLLM_WAIT()
        INQUIRE(FILE=FNAME,EXIST=LEX)
        IF (LEX) THEN
          IDRUNNER = I
          RETURN
        ENDIF
      ENDIF
    ENDDO
    RETURN
  END FUNCTION PLLM_FIND_RUNNER
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_READ_DISPAR(IFAIL,NPT,IRUNRUNNER,NRUNRUNNER,PVAL)
    !   Read the contents of a "jdispar.rdy" file and delete it
    USE GLOBAL_DATA, ONLY: IVERB
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                          INTENT(OUT) :: IFAIL
    INTEGER,                          INTENT(IN)  :: NPT
    INTEGER,                          INTENT(OUT) :: IRUNRUNNER
    INTEGER,                          INTENT(OUT) :: NRUNRUNNER
    DOUBLE PRECISION, DIMENSION(NPT), INTENT(OUT) :: PVAL
    !
    !   Local variables
    INTEGER :: I, IOUNIT, ISTAT, K, KD, KTRY
    CHARACTER(LEN=20) :: VARERR = ' '
    LOGICAL :: LEX, LOP
    !
    !   Format statements
    100 FORMAT(A)
    120 FORMAT(I10)
    200 FORMAT(3(1X,G24.16))
    270 FORMAT(1X,'Deleting file "',A,'"')
    280 FORMAT(1X,'WARNING: Unable to delete file ',A)
    300 FORMAT(1X,'Error opening file "',A,'"')
    720 FORMAT(1X,'Error reading variable "',A,'" from file "',A,'"')
    740 FORMAT(1X,'Closing file and making attempt: ',I2)
    !
    IFAIL = 0
    IOUNIT = UTL_GETUNIT(7,1000)
    K = 0
    KD = 0
    KTRY = 0
    10 CONTINUE
    INQUIRE(UNIT=IOUNIT,OPENED=LOP)
    IF (LOP) CLOSE(IOUNIT,IOSTAT=ISTAT)
    IF (ISTAT>0 .AND. KTRY<11) THEN
      KTRY = KTRY+1
      CALL PLLM_WAIT()
      GOTO 10
    ENDIF
    OPEN(IOUNIT,FILE=FNDISPAR,STATUS='OLD',ERR=500)
    VARERR = 'COMMAND'
    READ(IOUNIT,100,ERR=400,END=400) COMMAND
    VARERR = 'IRUNRUNNER'
    READ(IOUNIT,120,ERR=400,END=400) IRUNRUNNER
    VARERR = 'NRUNRUNNER'
    READ(IOUNIT,120,ERR=400,END=400) NRUNRUNNER
    VARERR = 'PVE'
    READ(IOUNIT,200,ERR=400,END=400)(PVAL(I),I=1,NPT)
    IF (IVERB>3) WRITE(*,270) FNDISPAR
    KTRY = 0
    20 CONTINUE
    ISTAT = PLL_CLOSE_AND_DELETE(IOUNIT,FNDISPAR,1000)
    IF (ISTAT>0 .AND. KTRY<11) THEN
      KTRY = KTRY+1
      CALL PLLM_WAIT()
      GOTO 20
    ENDIF
    IF (ISTAT>0) WRITE(*,280)FNDISPAR
    RETURN
    !
    400 CONTINUE
    IF (VARERR .NE. ' ') THEN
      WRITE(*,720)TRIM(VARERR),TRIM(FNDISPAR)
      ERRRUNNER = ' '
      WRITE(ERRRUNNER,720)TRIM(VARERR),TRIM(FNDISPAR)
      K = K+1
      CLOSE(IOUNIT)
      IF (K<11) THEN
        CALL PLLM_WAIT(WAIT)
        INQUIRE(FILE=FNDISPAR,EXIST=LEX)
        IF (LEX) THEN
          WRITE(*,740)K
          GOTO 10
        ELSE
          IFAIL = 1
          RETURN
        ENDIF
      ELSE
        WRITE(*,'(A)')' Warning: WAIT time may be too small'
        CALL PLLM_RUNNER_STOP(1)
      ENDIF
    ENDIF
    !
    500 CONTINUE
    KD = KD+1
    IF (KD<11) THEN
      CALL PLLM_WAIT(WAIT)
      GOTO 10
    ENDIF
    WRITE(*,300) FNDISPAR
    ERRRUNNER = ' '
    WRITE(ERRRUNNER,300) FNDISPAR
    CALL PLLM_RUNNER_STOP(1)
  END SUBROUTINE PLLM_READ_DISPAR
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_READ_RUNDEP(IERR,IRUN,IRUNNER,NVEXT,EXTVAL)
    !   Read the contents of a "jrundep.rdy" file, then delete it
    USE GLOBAL_DATA, ONLY: IVERB, MAX_STRING_LEN
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER, INTENT(OUT) :: IERR
    INTEGER, INTENT(IN)  :: IRUN     ! Run number
    INTEGER, INTENT(IN)  :: IRUNNER  ! Runner number
    INTEGER, INTENT(IN)  :: NVEXT
    DOUBLE PRECISION, DIMENSION(NVEXT), INTENT(OUT) :: EXTVAL
    !
    !   Local variables
    INTEGER :: I, ISTAT, IRUNRUNNER, IUNIT, KTRY
    CHARACTER(LEN=MAX_STRING_LEN) :: FNAME
    LOGICAL :: LEX
    !
    !   Format statements
    100 FORMAT(3(1X,G24.16))
    150 FORMAT(I10)
    180 FORMAT(1X,'Run number (',I5,') read from ',A,' does not match run number ',   &
               I5,/,' -- dependent values not read')
    200 FORMAT(1X,'ERROR opening file "',A,'"',/,1X,'for run number ',I5,   &
               ' on runner number ',I4)
    250 FORMAT(1X,'ERROR reading from file "',A,'"',/,1X,'for run number ',I5)
    300 FORMAT(1X,'ERROR reading from file "',A,'"',/,1X,'for run number ',I5, &
               ' on runner number ',I4)
    350 FORMAT(1X,'ERROR: Cant''t delete file "',A,'"',/,1X,'for run number ',I5, &
               ' on runner number ',I4)
    360 FORMAT(1X,'Successfully deleted file "',A,'"',/,1X,'for run number ',I5, &
               ' on runner number ',I4)
    !
    IERR = 0
    FNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNDEP
    IUNIT = UTL_GETUNIT(7,1000)
    KTRY = 0
    10 CONTINUE
    OPEN(IUNIT,FILE=FNAME,STATUS='OLD',ERR=20,IOSTAT=ISTAT)
    20 CONTINUE
    IF (ISTAT .NE. 0) THEN
      IF (KTRY<11) THEN
        KTRY = KTRY+1
        CALL PLLM_WAIT()
        GOTO 10
      ELSE
        IF (IVERB>2) WRITE(*,200)TRIM(FNAME),IRUN,IRUNNER
        IERR = 1
        RETURN
      ENDIF
    ENDIF
    !
    !   Read run number and dependent values from jrundep.rdy file
    KTRY = 0
    30 CONTINUE
    READ(IUNIT,150,IOSTAT=ISTAT)IRUNRUNNER
    IF (ISTAT .NE. 0) THEN
      IF (KTRY>=11) THEN
        IF (IVERB>2) WRITE(*,250)TRIM(FNAME),IRUN
        IERR = 1
        WRITE(*,'(A)')' Warning: WAIT time may be too small'
        35 CONTINUE
        CLOSE(IUNIT,STATUS='KEEP',IOSTAT=ISTAT)
        IF (ISTAT>0 .AND. KTRY<21) THEN
          KTRY = KTRY+1
          CALL PLLM_WAIT()
          GOTO 35
        ENDIF
        RETURN
      ELSE
        KTRY = KTRY+1
        CALL PLLM_WAIT()
        GOTO 30
      ENDIF
    ENDIF
    IF (IRUNRUNNER==IRUN) THEN
      READ(IUNIT,100,ERR=40,END=40,IOSTAT=ISTAT)(EXTVAL(I),I=1,NVEXT)
      40 CONTINUE
      IF (ISTAT .NE. 0) THEN
        CLOSE(IUNIT,STATUS='KEEP')
        IF (KTRY==10) THEN
          IF (IVERB>2) WRITE(*,300)TRIM(FNAME),IRUN,IRUNNER
          IERR = 1
          WRITE(*,'(A)')' Warning: WAIT time may be too small'
          RETURN
        ELSE
          KTRY = KTRY+1
          CALL PLLM_WAIT()
          GOTO 10
        ENDIF
      ENDIF
    ELSE
      IERR = 2 ! Signifies that IRUNRUNNER does not match IRUN
      IF (IVERB>3) WRITE(*,180)IRUNRUNNER,FNRUNDEP,IRUN
    ENDIF
    !
    IF (IERR==0 .OR. IERR==2) THEN
      KTRY = 0
      60 CONTINUE
      LEX = .FALSE.
      INQUIRE(FILE=FNAME,EXIST=LEX,IOSTAT=ISTAT)
      IF (LEX) THEN
        ISTAT = PLL_CLOSE_AND_DELETE(IUNIT,FNAME,1000)
        IF (ISTAT .NE. 0) THEN
          IF (KTRY < 20) THEN
            KTRY = KTRY+1
            CALL PLLM_WAIT()
            GOTO 60
          ELSE
            IF (IVERB>2) WRITE(*,350)TRIM(FNAME),IRUN,IRUNNER
            IERR = 3 ! Signifies that run completed and model-simulated values 
                     ! were extracted, but jrundep.rdy could not be deleted
          ENDIF
        ENDIF
      ENDIF
      IF (ISTAT==0 .AND. IVERB>3) THEN
        WRITE(*,360)TRIM(FNAME),IRUN,IRUNNER
      ENDIF
    ENDIF
    RETURN
  END SUBROUTINE PLLM_READ_RUNDEP
  ! ----------------------------------------------------------------------------
  LOGICAL FUNCTION PLLM_RUNNER_READY(IRUNNER) RESULT (READY)
    IMPLICIT NONE
    !
    !   Argument-list variable
    INTEGER, INTENT(IN) :: IRUNNER
    !
    !   Local variables
    CHARACTER(LEN=MAX_STRING_LEN) :: FNAME
    LOGICAL :: LEX
    !
    READY = .FALSE.
    FNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNRDY
    INQUIRE(FILE=FNAME,EXIST=LEX)
    READY = LEX
    RETURN
  END FUNCTION PLLM_RUNNER_READY
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_RUNNER_STOP(IFLAG)
    !   To be executed by runner when it needs to stop.  If IFLAG=1, create a
    !   jrunfail.fin file, indicating that the runner is stopping due to an
    !   error.
    USE GLOBAL_DATA, ONLY: AMESSAGE, IVERB
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER, OPTIONAL, INTENT(IN) :: IFLAG
    !
    !   Local variables
    INTEGER :: IFL, ISTAT, IU
    LOGICAL :: LEX
    !
    !   Format statements
    100 FORMAT(1X,'Runner "',A,'" stopping')
    120 FORMAT(1X,A)
    150 FORMAT(1X,'Runner stopping')
    200 FORMAT(1X,'Runner "',A,'" stopping due to error')
    250 FORMAT(1X,'Runner stopping due to error')
    270 FORMAT(1X,'Deleting file "',A,'"')
    !
    IU = UTL_GETUNIT(7,1000)
    IF (PRESENT(IFLAG)) THEN
      IFL = IFLAG
    ELSE
      IFL = 0
    ENDIF
    !
    INQUIRE(FILE=FNRUNRDY,EXIST=LEX)
    IF (LEX) THEN
      OPEN(IU,FILE=FNRUNRDY,IOSTAT=ISTAT)
      IF (ISTAT==0) THEN
        IF (IVERB>3) WRITE(*,270)TRIM(FNRUNRDY)
        ISTAT = PLL_CLOSE_AND_DELETE(IU,FNRUNRDY)
      ENDIF
    ENDIF
    !
    INQUIRE(FILE=FNRUNDEP,EXIST=LEX)
    IF (LEX) THEN
      OPEN(IU,FILE=FNRUNDEP,IOSTAT=ISTAT)
      IF (ISTAT==0) THEN
        IF (IVERB>3) WRITE(*,270)TRIM(FNRUNDEP)
        ISTAT = PLL_CLOSE_AND_DELETE(IU,FNRUNDEP)
      ENDIF
    ENDIF
    !
    IF (IFL==1) THEN
      OPEN(IU,FILE=FNRUNFAIL,STATUS='REPLACE')
      WRITE(IU,120)TRIM(AMESSAGE)//'  '//TRIM(ERRRUNNER)
      CLOSE(IU,STATUS='KEEP',IOSTAT=ISTAT)
      IF (SNAME .NE. ' ') THEN
        WRITE(*,200)TRIM(SNAME)
      ELSE
        WRITE(*,250)
      ENDIF
    ENDIF
    !
    INQUIRE(FILE=FNDISFIN,EXIST=LEX)
    IF (LEX) THEN
      OPEN(IU,FILE=FNDISFIN,IOSTAT=ISTAT)
      IF (ISTAT==0) THEN
        IF (IVERB>3) WRITE(*,270)TRIM(FNDISFIN)
        ISTAT = PLL_CLOSE_AND_DELETE(IU,FNDISFIN)
        500 CONTINUE
      ENDIF
    ENDIF
    !
    OPEN(IU,FILE=FNRUNFIN,STATUS='REPLACE')
    CLOSE(IU,STATUS='KEEP',IOSTAT=ISTAT)
    IF (SNAME .NE. ' ') THEN
      WRITE(*,100)TRIM(SNAME)
    ELSE
      WRITE(*,150)
    ENDIF
    CALL UTL_STOP(' ')
    !
    RETURN
  END SUBROUTINE PLLM_RUNNER_STOP
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_WRITE_DISPAR(IFAIL,IRUN,IRUNNER,NPT,NRUNSPLL,PVAL)
    !   Write a "jdispar.rdy" file
    USE GLOBAL_DATA, ONLY: IVERB, MAX_STRING_LEN
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument list variables
    INTEGER, INTENT(OUT) :: IFAIL
    INTEGER, INTENT(IN)  :: IRUN        ! Run number
    INTEGER, INTENT(IN)  :: IRUNNER     ! Runner number
    INTEGER, INTENT(IN)  :: NPT         ! Number of  parameters
    INTEGER, INTENT(IN)  :: NRUNSPLL
    DOUBLE PRECISION, DIMENSION(NPT), INTENT(IN) :: PVAL ! Parameter values
    !
    !   Local variables
    CHARACTER(LEN=MAX_STRING_LEN) :: RENAM, TEMPNAME
    INTEGER :: I, IRDY, ISTAT
    CHARACTER(LEN=10) :: CHDATE, CHTIME, CHZONE
    INTEGER, DIMENSION(8) :: IEDT
    !
    !   Format statements
    100 FORMAT(A)
    120 FORMAT(I10)
    200 FORMAT(3(1X,G24.16))
    270 FORMAT(1X,'Deleting file "',A,'"')
    300 FORMAT(1X,'Executing command: "',A,'" for run ',I5)
    320 FORMAT(1X,'Signal sent to start run ',I3,' on runner ',A,   &
               ' at ',I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)
    !
    TEMPNAME = TRIM(RUNNERDIR(IRUNNER))//'j_temp.rdy'
    IFAIL = 0
    IRDY = UTL_GETUNIT(7,1000)
    OPEN(IRDY,FILE=TEMPNAME,STATUS='REPLACE',IOSTAT=ISTAT)
    IF (ISTAT .NE. 0) THEN
      IFAIL = 1
      RETURN
    ENDIF
    !   Write data to local, temporary file
    WRITE(IRDY,100)COMMAND
    WRITE(IRDY,120)IRUN
    WRITE(IRDY,120)NRUNSPLL
    WRITE(IRDY,200)(PVAL(I),I=1,NPT)
    CLOSE(IRDY,STATUS='KEEP')
    !
    IF (OS_DIS == 'WINDOWS') THEN
      RENAM=TRIM(RENAME_DIS)//' '//TRIM(TEMPNAME)//' '//FNDISPAR
    ELSE
      RENAM=TRIM(RENAME_DIS)//' '//TRIM(TEMPNAME)//' '//   &
            TRIM(RUNNERDIR(IRUNNER))//FNDISPAR
    ENDIF
    IF (IVERB>3) WRITE(*,300)TRIM(RENAM),IRUN
    CALL UTL_SYSTEM(RENAM)
    OPEN(IRDY,FILE=TEMPNAME)
    IF (IVERB>3) WRITE(*,270)TRIM(TEMPNAME)
    ISTAT = PLL_CLOSE_AND_DELETE(IRDY,TEMPNAME)
    CALL DATE_AND_TIME(CHDATE,CHTIME,CHZONE,IEDT)
    IF (IVERB>2) WRITE(*,320)IRUN,TRIM(RUNNERNAME(IRUNNER)),   &
                             (IEDT(I),I=1,3),(IEDT(I),I=5,7)
    !
    RETURN
  END SUBROUTINE PLLM_WRITE_DISPAR
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_WRITE_DISRDY(IRUNNER,LCIS,NCATS,NINSTRUCT,NMIFILE,NMOFILE,   &
                              NOPNT,NPT,NUMLADV,NVEXT,PRECISPRO,CATGOR,   &
                              CINSTSET,EXTNAM,INSTRUCTFILE,LOCINS,   &
                              MCVUSE,MIFILE,MOFILE,MRKDL,NW,PARNAM,   &
                              TFILE,OPTIONFLAG)
    !   Indicate dispatcher is ready, by creating file "jdispatch.rdy" in a runner
    !   directory.  jdispatch.rdy contains command to start model, and other
    !   data needed by the runner.
    !
    !   Contents of jdispatch.rdy:
    !      RUNNERNAME
    !      COMMAND to start model
    !      WAIT (seconds)
    !      IVERBRUNNER
    !      Number of parameters(?)
    !      Number of model-input files
    !      Number of dependents to be extracted
    !      Number of model-output files
    !      Length of instructions array
    !      Extraction instructions
    !      ... (more data from MODEL_IO module?)
    !
    USE UTILITIES
    USE GLOBAL_DATA, ONLY: IVERB, LENDNAM, MAX_STRING_LEN
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                                    INTENT(IN) :: IRUNNER
    INTEGER,                                    INTENT(IN) :: LCIS
    INTEGER,                                    INTENT(IN) :: NCATS
    INTEGER,                                    INTENT(IN) :: NINSTRUCT
    INTEGER,                                    INTENT(IN) :: NMIFILE
    INTEGER,                                    INTENT(IN) :: NMOFILE
    INTEGER,                                    INTENT(IN) :: NOPNT
    INTEGER,                                    INTENT(IN) :: NPT
    INTEGER,                                    INTENT(IN) :: NUMLADV
    INTEGER,                                    INTENT(IN) :: NVEXT
    INTEGER,                                    INTENT(IN) :: PRECISPRO
    CHARACTER(LEN=6),              DIMENSION(NMOFILE), INTENT(IN) :: CATGOR ! Model-calc'd val. cat. for each model-output file
    CHARACTER(LEN=1),              DIMENSION(LCIS),    INTENT(IN) :: CINSTSET ! holds compressed instruction set
    CHARACTER(LEN=LENDNAM),        DIMENSION(NVEXT),   INTENT(IN) :: EXTNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE), INTENT(IN) :: INSTRUCTFILE ! Instruction file names
    INTEGER,                       DIMENSION(NINSTRUCT), INTENT(IN) :: LOCINS ! pointer to instructions
    LOGICAL,                       DIMENSION(NCATS),   INTENT(IN) :: MCVUSE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE), INTENT(IN) :: MIFILE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE), INTENT(IN) :: MOFILE ! Model output file names
    CHARACTER(LEN=1),              DIMENSION(NMOFILE), INTENT(IN) :: MRKDL  ! Marker delimiters
    INTEGER,                       DIMENSION(NPT),     INTENT(IN) :: NW
    CHARACTER(LEN=12),             DIMENSION(NPT),     INTENT(IN) :: PARNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE), INTENT(IN) :: TFILE
    INTEGER, OPTIONAL, INTENT(IN) :: OPTIONFLAG
    !
    !   Local variables
    INTEGER :: IRDY, ISTAT, ISTAT2, J, KTRY
    LOGICAL :: LEX, LOP
    CHARACTER(LEN=MAX_STRING_LEN) :: FNAME, TEMPNAME, RENAM
    INTEGER :: OPTIONFLAGLOCAL
    !
    !   Format statements
    100 FORMAT(A)
    150 FORMAT(I10)
    200 FORMAT(G12.5)
    300 FORMAT(1X,'Executing command: "',A,'"')
    320 FORMAT(6(1X,A12))
    330 FORMAT(6(1X,I4,8X))
    360 FORMAT(3(1X,A))
    380 FORMAT(10(1X,L1))
    420 FORMAT(10(1X,A6))
    440 FORMAT(40(1X,A1))
    460 FORMAT(80A1)
    480 FORMAT(8I10)
    500 FORMAT('STATUSFILE  modflow.status  1')
    !
    IF (.NOT. DOPLL) RETURN
    !
    IF (PRESENT(OPTIONFLAG)) THEN
      OPTIONFLAGLOCAL = OPTIONFLAG
    ELSE
      OPTIONFLAGLOCAL = 0
    ENDIF
    !
    TEMPNAME = TRIM(RUNNERDIR(IRUNNER))//'j_temp.rdy'
    IRDY = UTL_GETUNIT(7,1000)
    KTRY = 0
    10 CONTINUE
    OPEN(IRDY,FILE=TEMPNAME,STATUS='REPLACE',IOSTAT=ISTAT)
    INQUIRE(UNIT=IRDY,OPENED=LOP,IOSTAT=ISTAT2)
    IF (.NOT. LOP .AND. (ISTAT.NE.0 .OR. ISTAT2.NE.0)) THEN
      KTRY = KTRY + 1
      CALL PLLM_WAIT()
      IF (KTRY<20) GOTO 10
      IRUNNERSTAT(IRUNNER) = -1
      RETURN
    ENDIF
    !
    WRITE(IRDY,100)TRIM(RUNNERNAME(IRUNNER))
    WRITE(IRDY,100)RENAME_DIS
    WRITE(IRDY,200)WAITRUNNERS
    WRITE(IRDY,150)IVERBRUNNER
    WRITE(IRDY,150)LCIS
    WRITE(IRDY,150)NINSTRUCT
    WRITE(IRDY,150)NMIFILE
    WRITE(IRDY,150)NMOFILE
    WRITE(IRDY,150)NPT
    WRITE(IRDY,150)NOPNT
    WRITE(IRDY,150)PRECISPRO
    WRITE(IRDY,150)NCATS
    WRITE(IRDY,150)NVEXT
    WRITE(IRDY,150)NUMLADV
    WRITE(IRDY,100)(MIFILE(J),J=1,NMIFILE)
    WRITE(IRDY,100)(TFILE(J),J=1,NMIFILE)
    WRITE(IRDY,320)(PARNAM(J),J=1,NPT)
    WRITE(IRDY,330)(NW(J),J=1,NPT)
    WRITE(IRDY,420)(CATGOR(J),J=1,NMOFILE)
    WRITE(IRDY,360)(EXTNAM(J),J=1,NVEXT)
    WRITE(IRDY,100)(INSTRUCTFILE(J),J=1,NMOFILE)
    WRITE(IRDY,100)(MOFILE(J),J=1,NMOFILE)
    WRITE(IRDY,440)(MRKDL(J),J=1,NMOFILE)
    WRITE(IRDY,380)(MCVUSE(J),J=1,NCATS)
    WRITE(IRDY,460)(CINSTSET(J),J=1,LCIS)
    WRITE(IRDY,480)(LOCINS(J),J=1,NINSTRUCT)
    IF (OPTIONFLAGLOCAL==1) THEN
      WRITE(IRDY,500)
    ENDIF
    CLOSE(IRDY)
    !
    FNAME = TRIM(RUNNERDIR(IRUNNER))//FNDISRDY
    INQUIRE(FILE=FNAME,EXIST=LEX)
    IF (LEX) THEN
      OPEN(IRDY,FILE=FNAME)
      ISTAT = PLL_CLOSE_AND_DELETE(IRDY,FNAME)
    ENDIF
    !
    IF (OS_DIS == 'WINDOWS') THEN
      RENAM=TRIM(RENAME_DIS)//' '//TRIM(TEMPNAME)//' '//FNDISRDY
    ELSE
      RENAM=TRIM(RENAME_DIS)//' '//TRIM(TEMPNAME)//' '//   &
            TRIM(RUNNERDIR(IRUNNER))//FNDISRDY
    ENDIF
    IF (IVERB>3) WRITE(*,300)TRIM(RENAM)
    CALL UTL_SYSTEM(RENAM)
    IRUNNERSTAT(IRUNNER) = 2
    !
    RETURN
  END SUBROUTINE PLLM_WRITE_DISRDY
  ! ----------------------------------------------------------------------------
  SUBROUTINE PLLM_WRITE_RUNDEP(IRUN,NVEXT,EXTVAL)
    !   Write dependent values to "jrundep.rdy" file
    USE GLOBAL_DATA, ONLY: IVERB, MAX_STRING_LEN
    USE UTILITIES
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                            INTENT(IN) :: IRUN
    INTEGER,                            INTENT(IN) :: NVEXT
    DOUBLE PRECISION, DIMENSION(NVEXT), INTENT(IN) :: EXTVAL
    !
    !   Local variables
    INTEGER :: I, IUNIT
    CHARACTER(LEN=MAX_STRING_LEN) :: SYSCOM
    CHARACTER(LEN=11) :: FNAME
    !
    !   Format statements
    100 FORMAT(3(1X,G24.16))
    150 FORMAT(I10)
    200 FORMAT(1X,'File "',A,'" has been written')
    !
    IUNIT = UTL_GETUNIT(7,1000)
    FNAME = 'temp_rundep'
    OPEN(IUNIT,FILE=FNAME,STATUS='REPLACE')
    WRITE(IUNIT,150)IRUN
    WRITE(IUNIT,100)(EXTVAL(I),I=1,NVEXT)
    CLOSE(IUNIT,STATUS='KEEP')
    SYSCOM = TRIM(SRENAME)//' '//FNAME//' '//FNRUNDEP
    CALL UTL_SYSTEM(SYSCOM)
    IF (IVERB>3) WRITE(*,200)FNRUNDEP
    RETURN
  END SUBROUTINE PLLM_WRITE_RUNDEP
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_STOP_RUNNERS()
    !   Delete "jdispatch.rdy" file in directories of all runners
    !   Write "jdispatch.fin" file in directories of all runners
    USE GLOBAL_DATA, ONLY: AMESSAGE, IVERB, MAX_STRING_LEN
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    !
    !   Argument-list variables
    !
    !   Local variables
    LOGICAL :: LEX
    INTEGER :: I, IRSTAT, ISTAT, IU, K, KACTIVE, KTRY
    CHARACTER(LEN=MAX_STRING_LEN) :: FNAME, SYSCOM
    !
    !   Format statements
    120 FORMAT(/,1X,'Status of runners:')
    140 FORMAT(1X,I5,2X,A20,2X,A)
    270 FORMAT(1X,'Deleting file "',A,'"')
    300 FORMAT(1X,'Executing PLLM_STOP_RUNNERS')
    320 FORMAT(1X,'Executing system command: ',A)
    650 FORMAT(1X,'ERROR opening file: ',A)
    !
    IF (IVERB>3) WRITE(*,300)
    IU = UTL_GETUNIT(7,1000)
    DO I=1,NUMRUNNERS
      !   First, delete no-longer relevant signal file(s)
      !
      FNAME = TRIM(RUNNERDIR(I))//FNDISPAR
      INQUIRE(FILE=FNAME,EXIST=LEX)
      IF (LEX) THEN
        OPEN(IU,FILE=FNAME,IOSTAT=ISTAT)
        IF (IVERB>3) WRITE(*,270)TRIM(FNAME)
        ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME)
      ENDIF
      !
      FNAME = TRIM(RUNNERDIR(I))//FNDISRDY
      5 CONTINUE
      INQUIRE(FILE=FNAME,EXIST=LEX)
      IF (LEX) THEN
        OPEN(IU,FILE=FNAME,IOSTAT=ISTAT)
        IF (ISTAT .NE. 0) THEN
          CALL PLLM_WAIT()
          GOTO 5
        ENDIF
        IF (IVERB>3) WRITE(*,270)TRIM(FNAME)
        KTRY = 0
        7 CONTINUE
        ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME)
        IF (ISTAT>0 .AND. KTRY<11) THEN
          KTRY = KTRY+1
          CALL PLLM_WAIT()
          GOTO 7
        ENDIF
      ENDIF
      !
      !   Write "jdispatch.fin" file
      !
      FNAME = TRIM(RUNNERDIR(I))//FNDISFIN
      INQUIRE(FILE=FNAME,IOSTAT=ISTAT,EXIST=LEX)
      IF (LEX .AND. ISTAT==0) THEN
        OPEN(IU,FILE=FNAME,IOSTAT=ISTAT)
        ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME)
      ENDIF
      !
      FNAME = TRIM(RUNNERDIR(I))//'j_temp.fin'
      OPEN(IU,FILE=FNAME,STATUS='REPLACE',ERR=550,IOSTAT=ISTAT)
      550 CONTINUE
      IF (ISTAT .NE. 0) THEN
        WRITE(*,650)TRIM(FNAME)
      ELSE
        IF (AUTOSTOPRUNNERS) THEN
          WRITE(IU,'(A)')'STOP'
        ELSE
          WRITE(IU,'(A)')'RESET'
        ENDIF
        KTRY = 0
        20 CONTINUE
        CLOSE(IU,STATUS='KEEP',IOSTAT=ISTAT)
        IF (ISTAT>0 .AND. KTRY<11) THEN
          KTRY = KTRY+1
          CALL PLLM_WAIT()
          GOTO 20
        ENDIF
        IF (OS_DIS == 'WINDOWS') THEN
          SYSCOM = TRIM(RENAME_DIS)//' '//TRIM(FNAME)//' '//FNDISFIN
        ELSE
          SYSCOM = TRIM(RENAME_DIS)//' '//TRIM(FNAME)//' '//   &
                   TRIM(RUNNERDIR(I))//FNDISFIN
        ENDIF
        IF (IVERB>3) WRITE(*,320)TRIM(SYSCOM)
        CALL UTL_SYSTEM(SYSCOM)
      ENDIF
      !
    ENDDO
    !
    ! Find, close, and delete jrunner.fin files
    !
    K = 0
    LOOP1: DO WHILE (K<10)
      K = K+1
      DO I=1,NUMRUNNERS
        FNAME = TRIM(RUNNERDIR(I))//FNRUNFIN
        INQUIRE(FILE=FNAME,EXIST=LEX)
        IF (LEX) THEN
          OPEN(IU,FILE=FNAME)
          IF (IVERB>3) WRITE(*,270)TRIM(FNAME)
          KTRY = 0
          10 CONTINUE
          ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME,5)
          IF (ISTAT .NE. 0) THEN
            IF (KTRY==10) THEN
              AMESSAGE = 'File access conflict for file "'//TRIM(FNAME)//'"'
              IF(K<3)WRITE(*,'(A)')AMESSAGE
              CYCLE LOOP1
            ELSE
              KTRY = KTRY+1
              CALL PLLM_WAIT()
              GOTO 10
            ENDIF
          ENDIF
          IF (IRUNNERSTAT(I) .NE. 10) IRUNNERSTAT(I) = 9
          FNAME = TRIM(RUNNERDIR(I))//FNDISFIN
          INQUIRE(FILE=FNAME,EXIST=LEX)
          IF (LEX) THEN
            OPEN(IU,FILE=FNAME,IOSTAT=ISTAT)
            IF (IVERB>3) WRITE(*,270)TRIM(FNAME)
            ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME)
          ENDIF
        ENDIF
      ENDDO
      !
      KACTIVE = 0
      DO I=1,NUMRUNNERS
        IF (IRUNNERSTAT(I)>-1 .AND. IRUNNERSTAT(I)<9) KACTIVE = KACTIVE+1
      ENDDO
      !
      IF (KACTIVE==0) GOTO 800
      CALL PLLM_WAIT()
    ENDDO LOOP1
    !
    800 CONTINUE
    !   Write status of runners to screen
      WRITE(*,120)
      DO I=1,NUMRUNNERS
        IF (IRUNNERSTAT(I)>0) THEN
          IRSTAT = IRUNNERSTAT(I)
        ELSE
          IRSTAT = -1
        ENDIF
        WRITE(*,140)I,RUNNERNAME(I),RUNNERSTAT(IRSTAT)
      ENDDO
    !
    RETURN
  END SUBROUTINE PLLM_STOP_RUNNERS
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_WAIT(TSECS)
    !   Wait a while
    ! Comment-out following line if compiler does not support SLEEPQQ
    USE IFPORT, ONLY: SLEEPQQ
    IMPLICIT NONE
    !
    !   Argument-list variable
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: TSECS
    !
    !   Local variables
    CHARACTER(LEN=8)  :: DATE
    CHARACTER(LEN=10) :: TIME
    CHARACTER(LEN=5)  :: ZONE
    INTEGER, DIMENSION(8) :: VALUES0, VALUES1
    DOUBLE PRECISION :: DIFF, SECS, YRSEC0, YRSEC1
    DOUBLE PRECISION, PARAMETER :: SPD = 86400.D0
    DOUBLE PRECISION, PARAMETER :: SPH = 3600.D0
    DOUBLE PRECISION, PARAMETER :: SPM = 60.D0
    INTEGER(4) :: MSEC    
    !
    IF (PRESENT(TSECS)) THEN
      SECS = TSECS
    ELSE
      SECS = WAIT
    ENDIF
    !
    IF (USESLEEP) THEN
      MSEC = 1000*SECS
      ! If compiler does not support SLEEPQQ or equivalent functionality, 
      ! comment-out following line and set USESLEEP = .FALSE. in data
      ! section of this module.
      CALL SLEEPQQ(MSEC)
    ELSE
      10 CONTINUE
      CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES0)
      !   Calculate seconds since beginning of month
      YRSEC0 = (VALUES0(3)-1)*SPD + VALUES0(5)*SPH + VALUES0(6)*SPM   &
               + VALUES0(7)*1.D0 + VALUES0(8)*1.0D-3
      DO WHILE (.TRUE.)
        CALL DATE_AND_TIME(DATE,TIME,ZONE,VALUES1)
        IF (VALUES0(2) .NE. VALUES1(2)) GOTO 10
        !   Again calculate seconds since beginning of month
        YRSEC1 = (VALUES1(3)-1)*SPD + VALUES1(5)*SPH + VALUES1(6)*SPM   &
                 + VALUES1(7)*1.D0 + VALUES1(8)*1.0D-3
        DIFF = YRSEC1-YRSEC0
        IF (DIFF .GE. SECS) RETURN
      ENDDO
    ENDIF
    !
    RETURN
  END SUBROUTINE PLLM_WAIT
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLM_USTOP(STOPMESS,IOUT)
    USE UTILITIES, ONLY: UTL_STOP
    IMPLICIT NONE
    !
    ! Argument-list variables
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: STOPMESS
    INTEGER, OPTIONAL, INTENT(IN) :: IOUT
    !
    ! Local variables
    CHARACTER(LEN=200) :: MESS
    !
    IF (PARALLEL_ACTIVE_GWM) CALL PLLM_STOP_RUNNERS()
    !
    MESS = ' '
    IF (PRESENT(STOPMESS)) MESS = STOPMESS
    IF (PRESENT(IOUT)) THEN
      CALL UTL_STOP(MESS,IOUT)
    ELSE
      CALL UTL_STOP(MESS)
    ENDIF
    STOP
  END SUBROUTINE PLLM_USTOP
  !-----------------------------------------------------------------------------
  SUBROUTINE COPY_TO_RUNNERS(FNAME,OK)
    ! Copy a file to all runner directories
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, COPY_FILE
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_SYSTEM_FUNC
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, K, KTRY, L1, L2, NUNIT
    LOGICAL :: LEX, OK1, OKLOCAL
    CHARACTER(LEN=1) :: SEPARATOR
    CHARACTER(LEN=100) :: ERRORMSG
    CHARACTER(LEN = 2000) :: COMMAND, PATH
    !
    OKLOCAL = .TRUE.
    ERRORMSG = ''
    IF (OS_DIS == 'WINDOWS') THEN
      SEPARATOR = '\'
    ELSE
      SEPARATOR = '/'
    ENDIF
    OKLOCAL = .TRUE.
    NUNIT = UTL_GETUNIT(7,10000)
    ! Delete file from all runner directories
    K = 0
    LOOP0: DO I=1,NUMRUNNERS
      IF (IRUNNERSTAT(I)>0) THEN
        K = K + 1
        PATH  = TRIM(RUNNERDIR(I)) // FNAME
        L1 = PLL_CLOSE_AND_DELETE(NUNIT,PATH,50,0.1D0)
        IF (L1.NE.0) THEN
          ! Runner is nonresponsive
          IRUNNERSTAT(I) = -IRUNNERSTAT(I)
          K = K - 1
        ENDIF
      ENDIF
    ENDDO LOOP0
    !
    IF (K>0) THEN
      K = 0
      LOOP1: DO I=1,NUMRUNNERS
        IF (IRUNNERSTAT(I)>0) THEN
          K = K + 1
          PATH  = TRIM(RUNNERDIR(I)) // FNAME
          KTRY = 0
          10 CONTINUE
          CALL COPY_FILE_TO_DIR(FNAME,RUNNERDIR(I),OK1)
          IF (.NOT. OK1) THEN
            INQUIRE(FILE=PATH,EXIST=LEX)
            IF (LEX) THEN 
              CYCLE LOOP1
            ELSE
              IF (KTRY<20) THEN
                KTRY = KTRY + 1
                CALL PLLM_WAIT(0.1D0)
                GOTO 10
              ELSE
                ! Runner is nonresponsive
                IRUNNERSTAT(I) = -IRUNNERSTAT(I)
                K = K - 1
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO LOOP1
      !
    ENDIF
    !
    IF (K==0) THEN
      OKLOCAL = .FALSE.
      ERRORMSG = 'Error: All runner directories are nonresponsive'
    ENDIF
    !
    IF (.NOT. OKLOCAL) THEN
      CALL APPEND_MESSAGE(ERRORMSG)
      IF (OK) OK = OKLOCAL
    ENDIF
    RETURN
  END SUBROUTINE COPY_TO_RUNNERS
  !-----------------------------------------------------------------------------
  SUBROUTINE COPY_FILE_TO_DIR(FNAME,DIR,OK)
    ! Copy file FNAME to directory DIR
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, COPY_FILE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    CHARACTER(LEN=*), INTENT(IN) :: DIR
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    CHARACTER(LEN=2000) :: FNAME2
    CHARACTER(LEN=200) :: MESS
    CHARACTER(LEN=1) :: SEPARATOR
    LOGICAL :: OKLOCAL
    !
    IF (OS_DIS == 'WINDOWS') THEN
      SEPARATOR = '\'
    ELSE
      SEPARATOR = '/'
    ENDIF
    OKLOCAL = .TRUE.
    FNAME2 = TRIM(DIR)//SEPARATOR//TRIM(FNAME)
    OK = COPY_FILE(FNAME,FNAME2)
    !
    ! Normal return
    IF (OK) RETURN
    !
    ! Error handling
    MESS = 'Error: Unable to copy file "' // TRIM(FNAME) // '"' &
           // ' to "' // TRIM(FNAME2) // '"'
    CALL APPEND_MESSAGE(MESS)
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE COPY_FILE_TO_DIR
  !-----------------------------------------------------------------------------
END MODULE PLLM
