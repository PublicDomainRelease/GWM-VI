PROGRAM JRUNNERM
  USE GLOBAL_DATA, ONLY: LENDNAM, MAX_STRING_LEN, IVERB
  USE UTILITIES
  USE PARALLEL_PROCESSING, ONLY:   &
      !   Data
      COMMAND, FNDISFIN, FNDISPAR, FNRUNFAIL, FNRUNRDY, SNAME, WAIT, &
      USESLEEP, &
      !   Subprograms
      PLL_INI_RUNNER_DIM, PLL_INI_RUNNER_POP, PLL_ADA, PLL_EXE, PLL_EXT,   &
      PLL_CLN,  PLL_READ_DISPAR, PLL_RUNNER_STOP, PLL_WAIT, &
      PLL_WRITE_RUNDEP, PLL_CLOSE_AND_DELETE
  !
  USE PARALLEL_PROCESSING, ONLY: STATUSFILE, ISTATFILETYPE
  !
  USE MODEL_IO, ONLY: MIO_CLN_DEALLOC
  USE JRUNNERMMOD, ONLY:   &
      !   Data (scalars)
      FINISH, FINMESS, IFAIL, IOUNIT, IRUNRUNNER, KWAIT, LCIS, LEX,   &
      NCATS, NINSTRUCT, NMIFILE, NMOFILE, NOPNT, NPT, NRUNRUNNER,   &
      NUMLADV, NVEXT, PRECISPRO, RESET,   &
      !   Data (arrays)
      CATGOR, CINSTSET, EXTNAM, EXTVAL, INSTRUCTFILE, LOCINS, MCVUSE,   &
      MIFILE, MOFILE, MRKDL, NW, PARNAM, PVAL, TFILE,   &
      !   Subprograms
      JRN_INI, JRN_INI_ALLOC, JRN_UEV_WRITE, JRN_CLN, JRN_POPEXTVAL
  IMPLICIT NONE
  !
  CHARACTER(LEN=40) :: PROGNAM
  PARAMETER (PROGNAM='JRUNNERM version 1.0.1 2/20/2014')
  INTEGER :: IEXTOUT, ISTAT, IU, KFIN, KTRY
  CHARACTER(LEN=8) :: FNEXT
  CHARACTER(LEN=5) :: CRUN
  !   Added for JRunnerM
  INTEGER :: DWSTAT, IMFSTAT, MFSTAT
  DOUBLE PRECISION :: HCLOSE, HDRY, HNOFLO
  CHARACTER(LEN=20) :: CDUM
  !
  !   Format statements
  50 FORMAT(1X,A)
  100 FORMAT(A)
  120 FORMAT(I5.5)
  150 FORMAT(1X,'File "',A,'" has been created.')
  160 FORMAT(1X,'File "',A,'" not found after model run on runner ',A)
  270 FORMAT(1X,'Deleting file "',A,'"')
  300 FORMAT(1X,'File "',A,'" found...')
  600 FORMAT(/,1X,'Runner "',A,'" received signal to stop execution')
  610 FORMAT(/,1X,'Runner received signal to stop execution')
  620 FORMAT(/,1X,'Runner "',A,'" resetting to continue execution (JRUNNER)',/)
  630 FORMAT(/,1X,'Runner resetting to continue execution (JRUNNER)',/)
  !
  WRITE(*,50)TRIM(PROGNAM)
  10 CONTINUE
  RESET = .FALSE.
  IVERB = 3
  IOUNIT = 7
  IU = 8
  IMFSTAT = 0
  IEXTOUT = 0
  KFIN = 0
  FNEXT = 'run     '
  SNAME = ' '
  COMMAND = ' '
  WAIT = 0.001D0
  !
  !   Delete spurious signal files
  CALL JRN_INI()
  !   Create jrunner.rdy file, indicating runner is active.  Look for
  !   jdispatch.rdy file; when found, read scalars for initialization.
  CALL PLL_INI_RUNNER_DIM(IOUNIT,RESET,LCIS,NCATS,NINSTRUCT,NMIFILE,NMOFILE,   &
                         NOPNT,NPT,NUMLADV,NVEXT,PRECISPRO,   &
                         'Program jrunnerm is running...')
  IF (RESET) THEN
    CALL PLL_CLN()
    CALL MIO_CLN_DEALLOC()
    GOTO 10
  ENDIF
  !   Allocate arrays
  CALL JRN_INI_ALLOC()
  !   Populate arrays that are constant for all runs
  CALL PLL_INI_RUNNER_POP(IOUNIT,LCIS,NCATS,NINSTRUCT,NMIFILE,NMOFILE,   &
                         NPT,NUMLADV,NVEXT,CATGOR,CINSTSET,EXTNAM,   &
                         INSTRUCTFILE,LOCINS,MCVUSE,MIFILE,MOFILE,   &
                         MRKDL,NW,PARNAM,TFILE)
  FINISH = .FALSE.
  LOOP1: DO WHILE (.NOT. FINISH)
    !   Ensure that FNRUNRDY file exists, to show runner is alive
    INQUIRE(FILE=FNRUNRDY,EXIST=LEX)
    IF (.NOT. LEX) THEN
      KTRY = 0
      15 CONTINUE
      OPEN(IU,FILE=FNRUNRDY,IOSTAT=ISTAT)
      IF (ISTAT>0 .AND. KTRY<11) THEN
        KTRY = KTRY+1
        CALL PLL_WAIT()
        GOTO 15
      ENDIF
      KTRY = 0
      20 CONTINUE
      CLOSE(IU,STATUS='KEEP',IOSTAT=ISTAT)
      IF (ISTAT>0 .AND. KTRY<11) THEN
        KTRY = KTRY+1
        CALL PLL_WAIT()
        GOTO 20
      ENDIF
      IF (IVERB>4) WRITE(*,150)TRIM(FNRUNRDY)
    ENDIF
    !
    !   Look for FNDISPAR file, indicating runner should make a model run
    INQUIRE(FILE=FNDISPAR,EXIST=LEX)
    IF (LEX) THEN
      IF (IVERB>3) WRITE(*,300)FNDISPAR
      INQUIRE(FILE=FNRUNRDY,EXIST=LEX)
      IF (LEX) THEN
        OPEN(IOUNIT,FILE=FNRUNRDY)
        IF(IVERB>3)WRITE(*,270)FNRUNRDY
        ISTAT = PLL_CLOSE_AND_DELETE(IOUNIT,FNRUNRDY)
      ENDIF
      !
      !   Read parameter values from "jdispar.rdy" file into PVAL array
      CALL PLL_READ_DISPAR(IFAIL,NPT,IRUNRUNNER,NRUNRUNNER,PVAL)
      IF (IFAIL .NE. 0) GOTO 550
      !
      !   Adapt parameter values -- Write required model-input files
      CALL PLL_ADA(NPT,NOPNT,PARNAM,PRECISPRO,NW,PVAL)
      !
      IF (IVERB>3) THEN
        !   Open 'run#####._ext' to hold extracted values, where #####
        !   is the run number
        IF (IRUNRUNNER<100000) THEN
          WRITE(CRUN,120) IRUNRUNNER
          FNEXT(4:8) = CRUN
          IEXTOUT = UTL_DX_OPEN(FNEXT,'_ext','REPLACE')
        ENDIF
      ENDIF
      !
      !   Execute model
      CALL PLL_EXE(IRUNRUNNER,NRUNRUNNER)
      !
      IF (ISTATFILETYPE==0) THEN
        MFSTAT = 0
      ELSEIF (ISTATFILETYPE==1) THEN
        !   Open model status file
        IMFSTAT = UTL_GETUNIT(7,2000)
        OPEN(IMFSTAT,FILE=STATUSFILE,STATUS='OLD',IOSTAT=ISTAT)
        IF (ISTAT==0) THEN
          !   Read variables from model status file
          READ(IMFSTAT,*)CDUM,CDUM,MFSTAT
          IF (MFSTAT==0 .OR. MFSTAT==1) THEN
            READ(IMFSTAT,*)CDUM,CDUM,DWSTAT
          ELSE
            DWSTAT = -1
          ENDIF
          CLOSE(IMFSTAT)
        ELSE
          OPEN(IMFSTAT,FILE=FNRUNFAIL)
          WRITE(IMFSTAT,160)TRIM(STATUSFILE),TRIM(SNAME)
          CLOSE(IMFSTAT)
          CALL UTL_STOP()
        ENDIF
      ENDIF
      !
      IF (MFSTAT==0) THEN
        !   Extract dependents
        CALL PLL_EXT(NCATS,NVEXT,EXTNAM,MCVUSE,EXTVAL)
      ELSE
        !   Populate EXTVAL with bogus dependents and values read from modflow.status
        CALL JRN_POPEXTVAL(MFSTAT,DWSTAT,HCLOSE,HDRY,HNOFLO)
      ENDIF
      !
      !   Write dependent values to "jrundep.rdy" file
      CALL PLL_WRITE_RUNDEP(IRUNRUNNER,NVEXT,EXTVAL)
      !
      IF (IVERB>3 .AND. IEXTOUT>0) THEN
        !   Write dependent values to _ext file
        CALL JRN_UEV_WRITE(IEXTOUT,IFAIL)
        IEXTOUT = UTL_DX_CLOSE('_ext','KEEP')
      ENDIF
      !
    ENDIF
    550 CONTINUE
    KWAIT = 0
    570 CONTINUE
    !   Look for FNDISFIN file, indicating runner should stop execution
    INQUIRE(FILE=FNDISFIN,EXIST=LEX)
    IF (LEX) THEN
      CALL JRN_CLN()
      CALL PLL_CLN()
      CALL MIO_CLN_DEALLOC()
      IF (IVERB>2 .AND. KFIN==0)WRITE(*,300)TRIM(FNDISFIN)
      KFIN = KFIN+1
      OPEN(IOUNIT,FILE=FNDISFIN,ERR=580,IOSTAT=ISTAT)
      580 CONTINUE
      IF (ISTAT==0) THEN
        READ(IOUNIT,100,ERR=700,IOSTAT=ISTAT)FINMESS
        700 CONTINUE
        IF (ISTAT==0) THEN
          ISTAT = PLL_CLOSE_AND_DELETE(IOUNIT,FNDISFIN)
        ELSE
          CLOSE(IOUNIT,STATUS='KEEP',IOSTAT=ISTAT)
          CALL PLL_WAIT()
          KWAIT = KWAIT+1
          IF (KWAIT<5) GOTO 570
        ENDIF
        !
        IF (FINMESS=='STOP') THEN
          IF (SNAME==' ') THEN
            WRITE(*,610)
          ELSE
            WRITE(*,600)TRIM(SNAME)
          ENDIF
          FINISH = .TRUE.
        ELSEIF (FINMESS=='RESET') THEN
          IF (SNAME==' ') THEN
            WRITE(*,630)
          ELSE
            WRITE(*,620)TRIM(SNAME)
          ENDIF
          GOTO 10
        ENDIF
      ENDIF
    ENDIF
    IF (USESLEEP) CALL PLL_WAIT()
  ENDDO LOOP1
  CALL PLL_RUNNER_STOP()
END PROGRAM JRUNNERM
