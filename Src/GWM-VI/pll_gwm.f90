MODULE PLL_GWM
  USE PLLM
  USE DATATYPES
  USE GLOBAL_DATA, ONLY: MAX_STRING_LEN
  !
  PUBLIC :: PLLGWM_MAKE_RUNS
  !
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE PLLGWM_MAKE_RUNS(IOUT,LCIS,NCATS,NINSTRUCT, &
                            NMIFILE,NMOFILE,NOPNT,   &
                            NPT,NUMLADV,NRUNSPLL,NVEXT,PRECISPRO,   &
                            CATGOR,CINSTSET,EXTNAM,INSTRUCTFILE,   &
                            LOCINS,MCVUSE,MIFILE,MOFILE,MRKDL,NW,PARNAM,   &
                            PARVALSETS,TFILE,DEPVALSETS,HDRY,OK,COMMANDPLL, &
                            OPTIONFLAG)
    !   Use runners to perform ADApt, EXEcute, and EXTract tasks in parallel
    USE GLOBAL_DATA, ONLY: AMESSAGE, IVERB, LENDNAM, MAX_STRING_LEN
    USE UTILITIES
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    USE GWMMOD, ONLY: DECVARS_TO_PVAL, MAIN_UEV
    USE GWM_UTLS, ONLY: APPEND_MESSAGE
    USE GWM1DCV3, ONLY : NFVAR, FVBASE, FVNCELL, FVCURRENT
    USE GWMDATA, ONLY: PVALTMP, NDEPSTAT
    USE GWM1RMS3, ONLY: MFCNVRG, DEWATER, DEWATERQ, IPGNA, HCLOSEG
    USE GWM1RMS3SUBS, ONLY: GWM1RMS3PP, GWM1RMS3FP
    USE GWM1HDC3, ONLY: NHB,NDD,NDF,NGD
    USE GWM1STC3, ONLY: NSF,NSD
    USE GWM1BAS3, ONLY: GWMOUT
    IMPLICIT NONE
    !
    !   Argument-list variables
    INTEGER,                                     INTENT(IN)  :: IOUT
    INTEGER,                                     INTENT(IN)  :: LCIS
    INTEGER,                                     INTENT(IN)  :: NCATS
    INTEGER,                                     INTENT(IN)  :: NINSTRUCT
    INTEGER,                                     INTENT(IN)  :: NMIFILE
    INTEGER,                                     INTENT(IN)  :: NMOFILE
    INTEGER,                                     INTENT(IN)  :: NOPNT
    INTEGER,                                     INTENT(IN)  :: NPT
    INTEGER,                                     INTENT(IN)  :: NUMLADV
    INTEGER,                                     INTENT(IN)  :: NRUNSPLL
    INTEGER,                                     INTENT(IN)  :: NVEXT
    INTEGER,                                     INTENT(IN)  :: PRECISPRO
    CHARACTER(LEN=6),              DIMENSION(NMOFILE),  INTENT(IN)  :: CATGOR ! Model-calc'd val. cat. for each model-output file
    CHARACTER(LEN=1),              DIMENSION(LCIS),     INTENT(IN)  :: CINSTSET ! holds compressed instruction set
    CHARACTER(LEN=LENDNAM),        DIMENSION(NVEXT),    INTENT(IN)  :: EXTNAM
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE),  INTENT(IN)  :: INSTRUCTFILE ! Instruction file names
    INTEGER,                       DIMENSION(NINSTRUCT), INTENT(IN) :: LOCINS ! pointer to instructions
    LOGICAL,                       DIMENSION(NCATS),    INTENT(IN)  :: MCVUSE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE),  INTENT(IN)  :: MIFILE
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMOFILE),  INTENT(IN)  :: MOFILE ! Model output file names
    CHARACTER(LEN=1),              DIMENSION(NMOFILE),  INTENT(IN)  :: MRKDL  ! Marker delimiters
    INTEGER,                       DIMENSION(NPT),      INTENT(IN)  :: NW
    CHARACTER(LEN=12),             DIMENSION(NPT),      INTENT(IN)  :: PARNAM
    DOUBLE PRECISION,        DIMENSION(NPT,NRUNSPLL), INTENT(INOUT) :: PARVALSETS
    CHARACTER(LEN=MAX_STRING_LEN), DIMENSION(NMIFILE),  INTENT(IN)  :: TFILE
    DOUBLE PRECISION,        DIMENSION(NVEXT,NRUNSPLL), INTENT(OUT) :: DEPVALSETS
    DOUBLE PRECISION, INTENT(IN) :: HDRY
    LOGICAL, INTENT(INOUT) :: OK
    CHARACTER(LEN=MAX_STRING_LEN),   OPTIONAL,   INTENT(IN)  :: COMMANDPLL
    INTEGER, OPTIONAL, INTENT(IN) :: OPTIONFLAG
    !
    !   Local variables
    CHARACTER(LEN=200)::GWMSTRG
    LOGICAL :: DONE, LEX, OKLOCAL
    INTEGER :: IERR, IPERT, IRUN, IRUNNER, ISTAT, IU, KDONE, KND, KNS, KSL1,   &
               KSTART, KTRY, KWAIT, NUMDONE, NUMDONELAST
    INTEGER                       :: NR2ST    ! Number of run to be started
    INTEGER, DIMENSION(NRUNSPLL)  :: IRUNSTAT ! Status of each run that's needed
    INTEGER, DIMENSION(NUMRUNNERS) :: IRUNNUM ! Number of run in progress on runner
    INTEGER, DIMENSION(8)         :: VALUES
    CHARACTER(LEN=MAX_STRING_LEN) :: FNAME, FRUNNERNAME
    CHARACTER(LEN=10)             :: CHDATE, CHTIME, CHZONE
    DOUBLE PRECISION              :: ELMIN, ELSECS, RTIME, TLIM
    CHARACTER(LEN=1) :: OK2INTERRUPT
    INTEGER :: IFIL
    !   Added for GWM
    INTEGER :: DWSTAT, MFSTAT, I, IPERTTEMP, OPTIONFLAGLOCAL
    LOGICAL :: PERTDONEPLL
    INTEGER :: II
    !
    !   Format statements
    13 FORMAT(' in PLLGWM_MAKE_RUNS, PLLM_WRITE_DISPAR reports',   &
              ' failure on runner ',A)
    20 FORMAT(1X,'Making ',I5,' runs of command "',A,'" on runners...')
    30 FORMAT(1X,'Runner "',A,'" is active')
    40 FORMAT(1X,'File inquiry on runner ',A,   &
              ' failed--runner is assumed inactive.')
    50 FORMAT(1X,'Runner "',A,'" reports failure and needs to be restarted.')
    60 FORMAT(A)
    70 FORMAT(1X,'WARNING: Cannot read error message from file "',A,'"')
    80 FORMAT(/,1X,'Runner "',A,'" reports the following error:')
    100 FORMAT(1X,I5,' of ',I5,' model runs completed.')
    120 FORMAT(1X,'No runners active -- Stopping')
    130 FORMAT(1X,'Runner ',A,' has successfully completed run number ',I5,   &
        ' in ',F11.5,' min.')
    131 FORMAT(1X,'Runner ',A,' has successfully completed run number ',I5,   &
        ' in ',F11.5,' min. but ',A,' file could not be deleted')
    135 FORMAT(1X,'Runtime for runner "',A,'" changed to ',G12.5,/,   &
        ' because elapsed time < expected runtime.')
    136 FORMAT(1X,'Runtime for runner "',A,'" changed to ',G12.5,/,   &
        ' because elapsed time > expected runtime.')
    150 FORMAT(/,1X,'Runner "',A,'" is overdue for run number ',I5,/   &
        1X,'Expected run time = ',G9.3,' seconds; elapsed time = ',G9.3,   &
        ' seconds.',/)
    160 FORMAT(/,1X,'Runner "',A,'" is overdue for run number ',I5,/)
    170 FORMAT(/,1X,'Runner "',A,'" is overdue for run number ',I5,/   &
        1X,'and needs to be restarted.',/)
    180 FORMAT(1X,'Increasing expected runtime for runner "',A,'" to ',G10.3,   &
        1X,'seconds.')
    200 FORMAT(1X,'Failed to read dependent values from "jrundep.rdy" file',   &
               /,' for run number ',I5,' and runner number ',I5)
    260 FORMAT(/,1X,'Runner "',A,'" has been terminated at run number ',I5,/   &
        1X,'and needs to be restarted.',/)
    300 FORMAT(' File ',A,' found -- OK to interrupt main program (Y or N)?')
    425 FORMAT(' Warning: WAIT time may be too small')
    !
    !   IRUNSTAT flags for model runs:
    !      0 -- Run not started
    !      1 -- Run in progress
    !      2 -- Run done
    !      3 -- Extraction done--job complete
    !     -1 -- Run failed or was terminated before completion
    !     -2 -- Run failed to converge
    !
    IF (.NOT. DOPLL) RETURN
    OKLOCAL = .TRUE.
    IF (PRESENT(COMMANDPLL)) COMMAND = COMMANDPLL
    IF (COMMAND == ' ') THEN
      CALL PLLM_STOP_RUNNERS()
      CALL APPEND_MESSAGE('Error in PLLGWM_MAKE_RUNS -- No command to execute')
      OKLOCAL = .FALSE.
      GOTO 999
    ENDIF
    !
    IF (PRESENT(OPTIONFLAG)) THEN
      OPTIONFLAGLOCAL = OPTIONFLAG
    ELSE
      OPTIONFLAGLOCAL = 0
    ENDIF
    !
    IF (IVERB>0) WRITE(*,20) NRUNSPLL,TRIM(COMMAND)
    !
    DONE = .FALSE.
    NUMDONE = 0
    NUMDONELAST = 0
    IRUNSTAT = 0
    IRUNNUM = 0
    LOOP1: DO WHILE (.NOT. DONE)
      !
      !#################################################################
      !   Determine number of runs not yet started
      !#################################################################
      KNS = 0
      LOOP2: DO IRUN=1,NRUNSPLL
        IF (IRUNSTAT(IRUN)==0 .OR. IRUNSTAT(IRUN)==-1 .OR. IRUNSTAT(IRUN)==-2) THEN
          KNS = KNS+1
        ENDIF
      ENDDO LOOP2
      !
      !#########################################################################
      !   If any model runs need to be started, start as many runs as possible
      !#########################################################################
      IF (KNS>0) THEN
        KSTART = 0
        KSL1 = 0
        KWAIT = 0
        320 CONTINUE
        LOOP3: DO IRUNNER=1,NUMRUNNERS
          !   Look for jrunfail.fin file, indicating runner has failed
          FRUNNERNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNFAIL
          INQUIRE(FILE=FRUNNERNAME,EXIST=LEX,IOSTAT=ISTAT)
          IF (ISTAT>0) THEN
            !   File inquiry has failed.  Assume runner is inactive.
            LEX = .FALSE.
            IF (IRUNNERSTAT(IRUNNER)>0) THEN
              IRUNNERSTAT(IRUNNER) = -IRUNNERSTAT(IRUNNER)
              IF (IVERB>2) WRITE(*,40)TRIM(RUNNERNAME(IRUNNER))
            ELSE
              IRUNNERSTAT(IRUNNER) = -1
            ENDIF
            IF (IRUNNUM(IRUNNER)>0) THEN
              IRUNSTAT(IRUNNUM(IRUNNER)) = -1
            ENDIF
            IRUNNUM(IRUNNER) = 0
          ENDIF
          IF (LEX) THEN
            !   jrunfail.fin file has been found.
            IRUNNERSTAT(IRUNNER) = 10
            IF (IRUNNUM(IRUNNER)>0) THEN
              IRUNSTAT(IRUNNUM(IRUNNER)) = -1
            ENDIF
            IRUNNUM(IRUNNER) = 0
            IU = UTL_GETUNIT(7,1000)
            OPEN(IU,FILE=FRUNNERNAME)
            IF (IVERB>0) THEN
              WRITE(*,50)TRIM(RUNNERNAME(IRUNNER))
              ERRRUNNER = ' '
              READ(IU,60,IOSTAT=ISTAT)ERRRUNNER
              IF (ISTAT==0) THEN
                IF (ERRRUNNER .NE. ' ') THEN
                  WRITE(*,80)TRIM(RUNNERNAME(IRUNNER))
                  AMESSAGE = ERRRUNNER
                  CALL UTL_WRITE_MESSAGE()
                  WRITE(*,60)' '
                  AMESSAGE = ' '
                ENDIF
              ELSE
                WRITE(*,70)TRIM(FRUNNERNAME)
              ENDIF
            ENDIF
            KTRY = 0
            500 CONTINUE
            ISTAT = PLL_CLOSE_AND_DELETE(IU,FRUNNERNAME)
            IF (ISTAT>0 .AND. KTRY<11) THEN
              KTRY = KTRY+1
              CALL PLLM_WAIT()
              GOTO 500
            ENDIF
            !
            !   Look for jrunner.rdy, indicating runner is 
            !   ready to make a model run
            FRUNNERNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNRDY
            INQUIRE(FILE=FRUNNERNAME,EXIST=LEX)
            IF (LEX) THEN
              KTRY = 0
              510 CONTINUE
              OPEN(IU,FILE=FRUNNERNAME,IOSTAT=ISTAT)
              IF (ISTAT>0 .AND. KTRY<11) THEN
                KTRY = KTRY+1
                CALL PLLM_WAIT()
                GOTO 510
              ENDIF
              IF (ISTAT==0) THEN
                KTRY = 0
                520 CONTINUE
                ISTAT = PLL_CLOSE_AND_DELETE(IU,FRUNNERNAME)
                IF (ISTAT>0 .AND. KTRY<11) THEN
                  KTRY = KTRY+1
                  CALL PLLM_WAIT()
                  GOTO 520
                ENDIF
              ENDIF
            ENDIF
          ENDIF
          !
          !   Check for jrunner.fin, indicating runner has been terminated
          FNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNFIN
          INQUIRE(FILE=FNAME,EXIST=LEX,IOSTAT=ISTAT)
          IF (ISTAT>0) THEN
            !   File inquiry has failed.  Assume runner is inactive.
            LEX = .FALSE.
            IF (IRUNNERSTAT(IRUNNER)>0) THEN
              IRUNNERSTAT(IRUNNER) = -IRUNNERSTAT(IRUNNER)
              IF (IVERB>2) WRITE(*,40)TRIM(RUNNERNAME(IRUNNER))
            ELSE
              IRUNNERSTAT(IRUNNER) = -1
            ENDIF
            IF (IRUNNUM(IRUNNER)>0) THEN
              IRUNSTAT(IRUNNUM(IRUNNER)) = -1
            ENDIF
            IRUNNUM(IRUNNER) = 0
          ENDIF
          IF (LEX) THEN
            !   jrunner.fin has been found
            IF (IRUNNUM(IRUNNER)>0) THEN
              IF (IVERB>1) WRITE(*,260)TRIM(RUNNERNAME(IRUNNER)),   &
                         IRUNNUM(IRUNNER)
              IRUNSTAT(IRUNNUM(IRUNNER)) = -1
            ENDIF
            IRUNNERSTAT(IRUNNER) = 9
            IRUNNUM(IRUNNER) = 0
            IU = UTL_GETUNIT(7,1000)
            KTRY = 0
            530 CONTINUE
            OPEN(IU,FILE=FNAME,IOSTAT=ISTAT)
            IF (ISTAT>0 .AND. KTRY<11) THEN
              KTRY = KTRY+1
              CALL PLLM_WAIT()
              GOTO 530
            ENDIF
            IF (ISTAT==0) THEN
              KTRY = 0
              540 CONTINUE
              ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME)
              IF (ISTAT>0 .AND. KTRY<11) THEN
                KTRY = KTRY+1
                CALL PLLM_WAIT()
                GOTO 540
              ENDIF
            ENDIF
          ENDIF
          !
          IF (IRUNNERSTAT(IRUNNER)==2) THEN
            KSL1 = KSL1 + 1
            !   Select a run that needs to be started
            LOOP4: DO IRUN=1,NRUNSPLL
              IF (IRUNSTAT(IRUN)==0 .OR. IRUNSTAT(IRUN)==-1 .OR.  &
                  IRUNSTAT(IRUN)==-2) THEN
                NR2ST = IRUN
                IPERT = IRUN
                EXIT LOOP4
              ENDIF
            ENDDO LOOP4
            !
            !   Adjust FVBASE and/or DELINC
            DO II=1,NFVAR
              FVBASE(II) = FVCURRENT(II)
            ENDDO
            CALL GWM1RMS3PP(IOUT,IPERT,.FALSE.,.FALSE.,NRUNSPLL,gwmstrg)
            IF (.NOT. OKLOCAL) GOTO 999
            !
            CALL DECVARS_TO_PVAL(NFVAR,NPT,FVBASE,  &
                                 FVNCELL,PARVALSETS(:,NR2ST))
            !   Create "jdispar.rdy" file, which signals to runner to ADApt
            !   parameters, EXEcute model, and EXTract dependent values
            CALL PLLM_WRITE_DISPAR(IERR,NR2ST,IRUNNER,NPT,NRUNSPLL,   &
                                  PARVALSETS(:,NR2ST))
            IF (IERR==0) THEN
              CALL DATE_AND_TIME(CHDATE,CHTIME,CHZONE,VALUES)
              BTIME(1:8,IRUNNER) = VALUES(1:8)
              IRUNSTAT(NR2ST) = 1
              IRUNNERSTAT(IRUNNER) = 3
              IRUNNUM(IRUNNER) = NR2ST
              KSTART = KSTART+1
            ELSE
              WRITE(*,13)TRIM(RUNNERNAME(IRUNNER))
              IRUNNERSTAT(IRUNNER) = -2
            ENDIF
            IF (KSTART==KNS) EXIT LOOP3
          ELSEIF (IRUNNERSTAT(IRUNNER)==3 .OR. IRUNNERSTAT(IRUNNER)==4) THEN
            KSL1 = KSL1+1
          ELSEIF (IRUNNERSTAT(IRUNNER)<1 .OR. IRUNNERSTAT(IRUNNER)==9   &
                  .OR. IRUNNERSTAT(IRUNNER)==10   &
                  .OR. IRUNNERSTAT(IRUNNER)==11) THEN
            !   See if non-active runner is now available
            FRUNNERNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNRDY
            INQUIRE(FILE=FRUNNERNAME,EXIST=LEX)
            IF (LEX) THEN
                IRUNNERSTAT(IRUNNER) = 1
                KOD(IRUNNER) = 0
                KSL1 = KSL1+1
                IF (IVERB>2) WRITE(*,30)TRIM(RUNNERNAME(IRUNNER))
                !   Write "jdispatch.rdy" file and set IRUNNERSTAT=2
                CALL PLLM_WRITE_DISRDY(IRUNNER,LCIS,NCATS,NINSTRUCT,NMIFILE, &
                          NMOFILE,NOPNT,NPT,NUMLADV,NVEXT,PRECISPRO,CATGOR, &
                          CINSTSET,EXTNAM,INSTRUCTFILE,LOCINS,MCVUSE,MIFILE, &
                          MOFILE,MRKDL,NW,PARNAM,TFILE,OPTIONFLAGLOCAL)
            ENDIF
          ENDIF
        ENDDO LOOP3
        IF (KSL1==0) THEN
          IF (KWAIT<21) THEN
            KWAIT = KWAIT+1
            CALL PLLM_WAIT()
            GOTO 320
          ELSE
            WRITE(*,425)
          ENDIF
          WRITE(*,120)
          CALL PLLM_STOP_RUNNERS()
          OKLOCAL = .FALSE.
          GOTO 999
        ENDIF
      ENDIF
      !
      !#################################################################
      !   Determine number of model runs for which dependents are needed
      !#################################################################
      KND = 0
      LOOP5: DO IRUN=1,NRUNSPLL
        IF (IRUNSTAT(IRUN)==1) THEN
          KND = KND+1
        ENDIF
      ENDDO LOOP5
      !
      !#########################################################################
      !   If any model runs are in progress, look to see if jrundep.rdy exists
      !   or if any model runs are long overdue
      !#########################################################################
      IF (KND>0) THEN
        KDONE = 0
        LOOP6: DO IRUNNER=1,NUMRUNNERS
          IF (IRUNNERSTAT(IRUNNER)==3) THEN
            !   When available, copy dependent values from runners and store them
            FNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNDEP
            INQUIRE(FILE=FNAME,EXIST=LEX,IOSTAT=ISTAT)
            IF (ISTAT>0) THEN
              !   File inquiry has failed.  Assume runner is inactive.
              LEX = .FALSE.
              IF (IRUNNERSTAT(IRUNNER)>0) THEN
                IRUNNERSTAT(IRUNNER) = -IRUNNERSTAT(IRUNNER)
                IF (IVERB>2) WRITE(*,40)TRIM(RUNNERNAME(IRUNNER))
              ELSE
                IRUNNERSTAT(IRUNNER) = -1
              ENDIF
              IF (IRUNNUM(IRUNNER)>0) THEN
                IRUNSTAT(IRUNNUM(IRUNNER)) = -1
              ENDIF
              IRUNNUM(IRUNNER) = 0
            ENDIF
            IF (LEX) THEN
              IRUNNERSTAT(IRUNNER) = 4
              IRUN = IRUNNUM(IRUNNER)
              IPERT = IRUN
              IRUNSTAT(IRUN) = 2
              !   jrundep.rdy is available
              CALL PLLM_READ_RUNDEP(IERR,IRUN,IRUNNER,NVEXT,   &
                                   DEPVALSETS(:,IRUN))
              IF (IERR==0) THEN
                IF (IVERB>2) THEN
                  ELMIN = UTL_ELAPSED_TIME(BTIME(:,IRUNNER))/60.0D0
                  WRITE(*,130)TRIM(RUNNERNAME(IRUNNER)),IRUN,ELMIN
                ENDIF
                !
                !   Use extracted values to assign MFCNVRG, DEWATER, and DEWATERQ and
                !   populate HDCSTATE and STCSTATE arrays
                CALL MAIN_UEV(NDEPSTAT,DEPVALSETS(:,IRUN:IRUN),IRUN,NHB,NDD,NDF, &
                              NGD,NSF,NSD,NFVAR,DEWATERQ,MFCNVRG,HCLOSEG,HDRY, &
                              GWMOUT,OKLOCAL)
                IF (.NOT. OKLOCAL) GOTO 999
                !
                IF (MFCNVRG(IRUN)   &
                    .AND. .NOT. DEWATER(IRUN)  &
                    .AND. .NOT. DEWATERQ(IRUN)) THEN
                  PERTDONEPLL = .TRUE.
                ELSE
                  PERTDONEPLL = .FALSE.
                ENDIF
                IRUNNERSTAT(IRUNNER) = 2
                IRUNNUM(IRUNNER) = 0
                ELSECS = UTL_ELAPSED_TIME(BTIME(:,IRUNNER))
                RTIME = RUNTIME(IRUNNER)
                IF (ELSECS<RTIME) THEN
                  RUNTIME(IRUNNER) = (2.0D0*RTIME + ELSECS)/3.0D0
                  IF (IVERB>4) WRITE(*,135)TRIM(RUNNERNAME(IRUNNER)),RUNTIME(IRUNNER)
                ELSEIF (ELSECS>RTIME) THEN
                  RUNTIME(IRUNNER) = ELSECS
                  IF (IVERB>4) WRITE(*,136)TRIM(RUNNERNAME(IRUNNER)),RUNTIME(IRUNNER)
                ENDIF
                !   Use simulation results to compute response matrix and 
                !   augmented right hand side
                ! GWM1RMS3FP increments IPERT argument, so use another variable
                IPERTTEMP = IPERT
                CALL GWM1RMS3FP(IPERTTEMP,NRUNSPLL,.FALSE.,.FALSE.)
                IF (.NOT. OKLOCAL) GOTO 999
!
                IF (IPGNA(IPERT)==2)PERTDONEPLL = .FALSE. ! 2 means Response Precision Inadequate
                IF (PERTDONEPLL) THEN
                  IRUNSTAT(IRUN) = 3
                  KDONE = KDONE+1
                ELSE
                  IRUNSTAT(IRUN) = -2
                ENDIF
              ELSEIF (IERR==2) THEN
                !   Run number read from jrundep.rdy file does not match
                !   expected run number.  Assume runner is ready to accept
                !   another job and the run needs to be rerun.
                IRUNNERSTAT(IRUNNER) = 2
                IRUNNUM(IRUNNER) = 0
                IRUNSTAT(IRUN) = 0
              ELSEIF (IERR==3) THEN
                ! File jrundep.rdy could not be deleted, but extracted values were successfully read
                IF (IVERB>2) THEN
                  ELMIN = UTL_ELAPSED_TIME(BTIME(:,IRUNNER))/60.0D0
                  WRITE(*,131)TRIM(RUNNERNAME(IRUNNER)),IRUN,ELMIN,TRIM(FNRUNDEP)
                ENDIF
                !
                !   Use extracted values to assign MFCNVRG, DEWATER, and DEWATERQ and
                !   populate HDCSTATE and STCSTATE arrays
                CALL MAIN_UEV(NDEPSTAT,DEPVALSETS(:,IRUN:IRUN),IRUN,NHB,NDD,NDF, &
                              NGD,NSF,NSD,NFVAR,DEWATERQ,MFCNVRG,HCLOSEG,HDRY, &
                              GWMOUT,OKLOCAL)
                IF (.NOT. OKLOCAL) GOTO 999
                !
                IF (MFCNVRG(IRUN)   &
                    .AND. .NOT. DEWATER(IRUN)  &
                    .AND. .NOT. DEWATERQ(IRUN)) THEN
                  PERTDONEPLL = .TRUE.
                ELSE
                  PERTDONEPLL = .FALSE.
                ENDIF
                IRUNNERSTAT(IRUNNER) = -4  ! Runner is nonresponsive
                IRUNNUM(IRUNNER) = 0
                ELSECS = UTL_ELAPSED_TIME(BTIME(:,IRUNNER))
                RTIME = RUNTIME(IRUNNER)
                IF (ELSECS<RTIME) THEN
                  RUNTIME(IRUNNER) = (2.0D0*RTIME + ELSECS)/3.0D0
                  IF (IVERB>4) WRITE(*,135)TRIM(RUNNERNAME(IRUNNER)),RUNTIME(IRUNNER)
                ELSEIF (ELSECS>RTIME) THEN
                  RUNTIME(IRUNNER) = ELSECS
                  IF (IVERB>4) WRITE(*,136)TRIM(RUNNERNAME(IRUNNER)),RUNTIME(IRUNNER)
                ENDIF
                !   Use simulation results to compute response matrix and 
                !   augmented right hand side
                ! GWM1RMS3FP increments IPERT argument, so use another variable
                IPERTTEMP = IPERT
                CALL GWM1RMS3FP(IPERTTEMP,NRUNSPLL,.FALSE.,.FALSE.)
                IF (.NOT. OKLOCAL) GOTO 999
                !
                IF (IPGNA(IPERT)==2)PERTDONEPLL = .FALSE. ! 2 means Response Precision Inadequate
                IF (PERTDONEPLL) THEN
                  IRUNSTAT(IRUN) = 3
                  KDONE = KDONE+1
                ELSE
                  IRUNSTAT(IRUN) = -2
                ENDIF
              ELSE
                WRITE(*,200)IRUNNUM(IRUNNER),IRUNNER
                CALL PLLM_STOP_RUNNERS()
              ENDIF
            ELSE
              !   Check for long overdue model run
              ELSECS = UTL_ELAPSED_TIME(BTIME(:,IRUNNER))
              TLIM = RUNTIME(IRUNNER)*TIMEOUTFAC
              IF (ELSECS>TLIM) THEN
                !   Model run is overdue
                IF (IVERB>2) THEN
                  WRITE(*,150)TRIM(RUNNERNAME(IRUNNER)),IRUNNUM(IRUNNER),   &
                      RUNTIME(IRUNNER),ELSECS
                ELSE
                  WRITE(*,160)TRIM(RUNNERNAME(IRUNNER)),IRUNNUM(IRUNNER)
                ENDIF
                IF (IRUNNUM(IRUNNER)>0) THEN
                  IRUNSTAT(IRUNNUM(IRUNNER)) = -1
                ENDIF
                IF (KOD(IRUNNER)>2) THEN
                  !  Too many overdue runs; assume runner is no longer running
                  WRITE(*,170)TRIM(RUNNERNAME(IRUNNER))
                  IRUNNERSTAT(IRUNNER) = -1
                  IRUNNUM(IRUNNER) = 0
                ELSE
                  !   Increase expected run time for next run
                  RUNTIME(IRUNNER) = RUNTIME(IRUNNER)*1.5D0
                  IF (IVERB>3) WRITE(*,180)TRIM(RUNNERNAME(IRUNNER)),RUNTIME(IRUNNER)
                  KOD(IRUNNER) = KOD(IRUNNER) + 1
                  IRUNNERSTAT(IRUNNER) = 11
                ENDIF
              ENDIF
              !   Check for jrunner.fin, indicating runner has been terminated
              FNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNFIN
              INQUIRE(FILE=FNAME,EXIST=LEX,IOSTAT=ISTAT)
              IF (ISTAT>0) THEN
                !   File inquiry has failed.  Assume runner is inactive.
                LEX = .FALSE.
                IF (IRUNNERSTAT(IRUNNER)>0) THEN
                  IRUNNERSTAT(IRUNNER) = -IRUNNERSTAT(IRUNNER)
                  IF (IVERB>2) WRITE(*,40)TRIM(RUNNERNAME(IRUNNER))
                ELSE
                  IRUNNERSTAT(IRUNNER) = -1
                ENDIF
                IF (IRUNNUM(IRUNNER)>0) THEN
                  IRUNSTAT(IRUNNUM(IRUNNER)) = -1
                ENDIF
                IRUNNUM(IRUNNER) = 0
              ENDIF
              IF (LEX) THEN
                !   jrunner.fin has been found
                IF (IRUNNUM(IRUNNER)>0) THEN
                  IF (IVERB>1) WRITE(*,260)TRIM(RUNNERNAME(IRUNNER)),   &
                             IRUNNUM(IRUNNER)
                  IRUNSTAT(IRUNNUM(IRUNNER)) = -1
                ENDIF
                IRUNNERSTAT(IRUNNER) = -1
                IRUNNUM(IRUNNER) = 0
              ENDIF
              !   Check for jrunfail.fin, indicating runner reports error
              FNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNFAIL
              INQUIRE(FILE=FNAME,EXIST=LEX,IOSTAT=ISTAT)
              IF (ISTAT>0) THEN
                !   File inquiry has failed.  Assume runner is inactive.
                LEX = .FALSE.
                IF (IRUNNERSTAT(IRUNNER)>0) THEN
                  IRUNNERSTAT(IRUNNER) = -IRUNNERSTAT(IRUNNER)
                  IF (IVERB>2) WRITE(*,40)TRIM(RUNNERNAME(IRUNNER))
                ELSE
                  IRUNNERSTAT(IRUNNER) = -1
                ENDIF
                IF (IRUNNUM(IRUNNER)>0) THEN
                  IRUNSTAT(IRUNNUM(IRUNNER)) = -1
                ENDIF
                IRUNNUM(IRUNNER) = 0
              ENDIF
              IF (LEX) THEN
                !   jrunfail.fin has been found
                IRUNNERSTAT(IRUNNER) = 10
                IF (IRUNNUM(IRUNNER)>0) THEN
                  IRUNSTAT(IRUNNUM(IRUNNER)) = -1
                ENDIF
                IRUNNUM(IRUNNER) = 0
                IU = UTL_GETUNIT(7,1000)
                KTRY = 0
                545 CONTINUE
                OPEN(IU,FILE=FNAME,IOSTAT=ISTAT)
                IF (ISTAT>0 .AND. KTRY<11) THEN
                  KTRY = KTRY+1
                  CALL PLLM_WAIT()
                  GOTO 545
                ENDIF
                IF (IVERB>0) THEN
                  WRITE(*,50)TRIM(RUNNERNAME(IRUNNER))
                  ERRRUNNER = ' '
                  READ(IU,60)ERRRUNNER
                  IF (ERRRUNNER .NE. ' ') THEN
                    WRITE(*,80)TRIM(RUNNERNAME(IRUNNER))
                    AMESSAGE = ERRRUNNER
                    CALL UTL_WRITE_MESSAGE()
                    WRITE(*,60)' '
                    AMESSAGE = ' '
                  ENDIF
                ENDIF
                IF (ISTAT==0) THEN
                  KTRY = 0
                  560 CONTINUE
                  ISTAT = PLL_CLOSE_AND_DELETE(IU,FNAME)
                  IF (ISTAT>0 .AND. KTRY<11) THEN
                    KTRY = KTRY+1
                    CALL PLLM_WAIT()
                    GOTO 560
                  ENDIF
                ENDIF
                !
                !   Delete jrunner.rdy file for failed runner
                FRUNNERNAME = TRIM(RUNNERDIR(IRUNNER))//FNRUNRDY
                INQUIRE(FILE=FRUNNERNAME,EXIST=LEX)
                IF (LEX) THEN
                  KTRY = 0
                  580 CONTINUE
                  OPEN(IU,FILE=FRUNNERNAME,IOSTAT=ISTAT)
                  IF (ISTAT>0 .AND. KTRY<11) THEN
                    KTRY = KTRY+1
                    CALL PLLM_WAIT()
                    GOTO 580
                  ENDIF
                  IF (ISTAT==0) THEN
                    KTRY = 0
                    620 CONTINUE
                    ISTAT = PLL_CLOSE_AND_DELETE(IU,FRUNNERNAME)
                    IF (ISTAT>0 .AND. KTRY<11) THEN
                      KTRY = KTRY+1
                      CALL PLLM_WAIT()
                      GOTO 620
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO LOOP6
        !
        !   Provide capability to interrupt main program while parallel runs
        !   are being made.
        FNAME = 'jparallel.interrupt'
        INQUIRE(FILE=FNAME,EXIST=LEX)
        IF (LEX) THEN
          WRITE(*,300)TRIM(FNAME)
          READ(*,*)OK2INTERRUPT
          IFIL = UTL_GETUNIT(7,1000)
          KTRY = 0
          640 CONTINUE
          OPEN(IFIL,FILE=FNAME,IOSTAT=ISTAT)
          IF (ISTAT>0 .AND. KTRY<11) THEN
            KTRY = KTRY+1
            CALL PLLM_WAIT()
            GOTO 640
          ENDIF
          IF (ISTAT==0) THEN
            KTRY = 0
            660 CONTINUE
            CLOSE(IFIL,STATUS='DELETE',IOSTAT=ISTAT)
            IF (ISTAT>0 .AND. KTRY<11) THEN
              KTRY = KTRY+1
              CALL PLLM_WAIT()
              GOTO 660
            ENDIF
          ENDIF
          IF (UTL_SAMENAME(OK2INTERRUPT,'Y')) THEN
            CALL PLLM_STOP_RUNNERS()
            OKLOCAL = .FALSE.
            GOTO 999
          ENDIF
        ENDIF
        !
      ENDIF
      !
      CALL PLLM_DONE(NRUNSPLL,IRUNSTAT,DONE,NUMDONE)
      IF (.NOT. DONE) THEN
        IF (NUMDONE .NE. NUMDONELAST) THEN
          NUMDONELAST = NUMDONE
          IF (IVERB>2) WRITE(*,100)NUMDONE,NRUNSPLL
        ENDIF
      ELSE
        IF (IVERB>2) WRITE(*,100)NRUNSPLL,NRUNSPLL
      ENDIF
      !
      ! When SLEEPQQ is supported, waiting may result in less CPU usage
      IF (USESLEEP) CALL PLLM_WAIT()   
    ENDDO LOOP1
    !
    999 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE PLLGWM_MAKE_RUNS
  ! ----------------------------------------------------------------------------
END MODULE PLL_GWM
