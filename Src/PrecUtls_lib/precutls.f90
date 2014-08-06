!
!  This Module must be compiled with default real KIND = 4 bytes 
!                    and default double precision KIND = 8 bytes
!

MODULE PRECUTLSMOD
  PRIVATE
  PUBLIC :: GET_BINARY_HEAD_DATASET, GET_BUDGET_3D_ARRAY, &
            EXTRACT_BUD_DRN_DEPS, EXTRACT_BUD_SFR_DEPS, &
            EXTRACT_QNET_MNW2, &
            GET_BINARY_FILE_LENGTH, BUDGETPRECISION, &
            HEADPRECISION
  PUBLIC :: AUXDIM
  INTEGER, PARAMETER :: AUXDIM = 20
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IBUF
  REAL, ALLOCATABLE, DIMENSION(:) :: BUF
  REAL, ALLOCATABLE, DIMENSION(:,:) :: BUFF2D
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BUFF3D
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DBUF
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DBUF3D
  REAL, DIMENSION(6) :: QV
  DOUBLE PRECISION, DIMENSION(6) :: QVD
CONTAINS
  !-------------------------------------------------------------------
  INTEGER FUNCTION BUDGETPRECISION(IU)
    ! Determine single or double precision file type for a MODFLOW
    ! budget file:  0=unrecognized, 1=single, 2=double.
    IMPLICIT NONE
    ! Argument-list variables
    INTEGER, INTENT(IN) :: IU
    ! Local variables
    INTEGER :: ICELL, ICODE, IERR, IPREC, KPER, KSTP, &
               N, NC, NCOL, NL, NLAY, NODES, NLST, NR, NROW
    REAL :: DELT, PERTIM, TOTIM, VAL
    DOUBLE PRECISION :: DELTD, PERTIMD, TOTIMD, VALD
    CHARACTER(LEN=16) :: TEXT1, TEXT2
    LOGICAL :: MOREDATA
    !
    !  Default is unrecognized file
    IPREC=0
    !
    !  SINGLE check
    READ(IU,ERR=100,END=100) KSTP,KPER,TEXT1,NCOL,NROW,NLAY
    ICODE=0
    IF (NLAY.LT.0) THEN
      NLAY=-NLAY
      READ(IU,ERR=50,END=50) ICODE,DELT,PERTIM,TOTIM
    ENDIF
    IF (NCOL.LT.1 .OR. NROW.LT.1 .OR. NLAY.LT.1) GO TO 100
    IF (NCOL.GT.100000000 .OR. NROW.GT.100000000 .OR. &
                        NLAY.GT.100000000) GO TO 100
    IF (NCOL*NROW.GT.100000000 .OR. NCOL*NLAY.GT.100000000 .OR. &
                    NROW*NLAY.GT.100000000) GO TO 100
    ALLOCATE(BUFF3D(NCOL,NROW,NLAY))
    ALLOCATE(DBUF3D(NCOL,NROW,NLAY))
    NODES=NCOL*NROW*NLAY
    !
    !  Read data depending on ICODE.  
    IF (ICODE.EQ.0 .OR. ICODE.EQ.1) THEN
      READ(IU,ERR=50,END=50) BUFF3D
    ELSEIF (ICODE.EQ.2) THEN
      READ(IU,ERR=50,END=50) NLST
      IF (NLST.LT.0) GO TO 50
      IF (NLST.GT.0) THEN
        DO N=1,NLST
          READ(IU,END=50,ERR=50) ICELL,VAL
          IF (ICELL.LE.0 .OR. ICELL.GT.NODES) GO TO 50
        ENDDO
      ENDIF
    ELSEIF (ICODE==3) THEN
      CALL READ_2D_ARRAY_LAYER(IU,1,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 50
    ELSEIF (ICODE==4) THEN
      CALL READ_2D_ARRAY_L1(IU,1,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 50
    ELSEIF (ICODE==5) THEN
      CALL READ_LIST_AUX(IU,1,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 50
    ELSE
      GO TO 100
    ENDIF
    !
    !  Read 2nd header and check for valid type.
    MOREDATA = .FALSE.
    READ(IU,ERR=50,END=40) KSTP,KPER,TEXT2
    MOREDATA = .TRUE.
    IF (VALID_BUDGET_TEXT(TEXT1) .AND. VALID_BUDGET_TEXT(TEXT2)) THEN
      IPREC = 1
    ENDIF
    40 CONTINUE
    IF (.NOT. MOREDATA) IPREC = 1 ! EOF indicates file contains a single, valid data set
    IF (IPREC==1) GOTO 100
    !
    !  DOUBLE check
    50 REWIND(IU)
    READ(IU,ERR=100,END=100) KSTP,KPER,TEXT1,NC,NR,NL
    ICODE=0
    IF (NL.LT.0) THEN
      NL=-NL
      READ(IU,ERR=100,END=100) ICODE,DELTD,PERTIMD,TOTIMD
    ENDIF
    !
    !  Read data depending on ICODE.
    IF (ICODE.EQ.0 .OR. ICODE.EQ.1) THEN
      READ(IU,ERR=100,END=100) DBUF3D
    ELSEIF (ICODE.EQ.2) THEN
      READ(IU,ERR=100,END=100) NLST
      IF (NLST.LT.0) GO TO 100
      IF (NLST.GT.0) THEN
        DO N=1,NLST
          READ(IU,END=100,ERR=100) ICELL,VALD
          IF (ICELL.LE.0 .OR. ICELL.GT.NODES) GO TO 100
        ENDDO
      ENDIF
    ELSEIF (ICODE==3) THEN
      CALL READ_2D_ARRAY_LAYER(IU,2,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 100
    ELSEIF (ICODE==4) THEN
      CALL READ_2D_ARRAY_L1(IU,2,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 100
    ELSEIF (ICODE==5) THEN
      CALL READ_LIST_AUX(IU,2,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 100
    ELSE
      GO TO 100
    ENDIF
    !
    !  Read 2nd header and check for valid type.
    MOREDATA = .FALSE.
    READ(IU,ERR=100,END=90) KSTP,KPER,TEXT2
    MOREDATA = .TRUE.
    IF (VALID_BUDGET_TEXT(TEXT1) .AND. VALID_BUDGET_TEXT(TEXT2)) THEN
      IPREC = 2
    ENDIF
    90 CONTINUE
    IF (.NOT. MOREDATA) IPREC = 2
    !
    100 REWIND(IU)
    BUDGETPRECISION = IPREC
    IF (ALLOCATED(BUFF3D)) DEALLOCATE(BUFF3D)
    IF (ALLOCATED(DBUF3D)) DEALLOCATE(DBUF3D)
    RETURN
  END FUNCTION BUDGETPRECISION
  !-------------------------------------------------------------------
  INTEGER FUNCTION HEADPRECISION(IU)
    ! Return 1 if head or drawdown file opened on unit IU is single
    ! precision.  Return 2 if it's double precision.  Return value
    ! of -1 means precision cannot be determined.
    IMPLICIT NONE
    ! Argument
    INTEGER, INTENT(IN) :: IU
    ! Local variables
    INTEGER :: I, ILAY, ISTAT, J, KP, KS, NCOL, NCOL2, NROW, NROW2
    CHARACTER(LEN=16) :: TEXT
    REAL :: PERTIMS, TOTIMS
    REAL, ALLOCATABLE, DIMENSION(:,:) :: VALS
    DOUBLE PRECISION :: PERTIMD, TOTIMD
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VALD
    !
    HEADPRECISION = -1 ! Unknown
    REWIND(IU)
    !
    ! First assume single precision
    READ(IU,ERR=20,END=20)KS,KP,PERTIMS,TOTIMS,TEXT,NCOL,NROW,ILAY
    IF (NCOL>0 .AND. NROW>0) THEN
      IF (NCOL*NROW<100000000) THEN
        ALLOCATE(VALS(NCOL,NROW),STAT=ISTAT)
        IF (ISTAT==0) THEN
          DO I=1,NROW
            READ(IU,ERR=20,END=20)(VALS(J,I),J=1,NCOL)
          ENDDO
          READ(IU,ERR=20,END=10)KS,KP,PERTIMS,TOTIMS,TEXT,NCOL2,NROW2,ILAY
          IF (NCOL2==NCOL .AND. NROW2==NROW) THEN
            HEADPRECISION = 1
          ELSE 
            GOTO 20
          ENDIF
        ELSE
          GOTO 20
        ENDIF
      ELSE
        GOTO 20
      ENDIF
    ELSE 
      GOTO 20
    ENDIF
    !
    10 CONTINUE
    ! If binary head file contains only 1 2-D head array, READ of 2nd header 
    ! record will fail with EOF error.  Since other READs must have been 
    ! successful, file is single precision
    HEADPRECISION = 1
    !
    ! Error reading as single precision
    20 CONTINUE
    !
    IF (HEADPRECISION<0) THEN
      ! Retry, assuming double precision
      REWIND(IU)
      READ(IU,ERR=40,END=40)KS,KP,PERTIMD,TOTIMD,TEXT,NCOL,NROW,ILAY
      IF (NCOL>0 .AND. NROW>0) THEN
        IF (NCOL*NROW<100000000) THEN
          ALLOCATE(VALD(NCOL,NROW),STAT=ISTAT)
          IF (ISTAT==0) THEN
            DO I=1,NROW
              READ(IU,ERR=40,END=40)(VALD(J,I),J=1,NCOL)
            ENDDO
            READ(IU,ERR=40,END=30)KS,KP,PERTIMD,TOTIMD,TEXT,NCOL2,NROW2,ILAY
            IF (NCOL2==NCOL .AND. NROW2==NROW) THEN
              HEADPRECISION = 2
            ELSE
              GOTO 40
            ENDIF
          ELSE
            GOTO 40
          ENDIF
        ELSE
          GOTO 40
        ENDIF
      ELSE
        GOTO 40
      ENDIF
      30 CONTINUE
      ! If binary head file contains only 1 2-D head array, READ of 2nd header 
      ! record will fail.  Since other READs must have been successful, file
      ! is double precision
      HEADPRECISION = 2
    ENDIF      
    !
    ! Deallocate arrays and exit
    40 CONTINUE
    IF (ALLOCATED(VALS)) DEALLOCATE(VALS)
    IF (ALLOCATED(VALD)) DEALLOCATE(VALD)
    REWIND(IU)
    RETURN
  END FUNCTION HEADPRECISION
  !-------------------------------------------------------------------
  SUBROUTINE GET_BINARY_HEAD_DATASET(IU,PRECISION,KSTP,KPER,TEXT, &
                                      NCOL,NROW,ILAY,DHEADS,OK)
    ! Read binary head data written by ULASAV
    USE GWM_UTLS, ONLY: APPEND_MESSAGE
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, NCOL, NROW, PRECISION
    INTEGER, INTENT(OUT) :: KSTP, KPER, ILAY
    CHARACTER(LEN=16), INTENT(OUT) :: TEXT
    DOUBLE PRECISION, DIMENSION(NCOL,NROW), INTENT(OUT) :: DHEADS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, J, NC, NR
    REAL :: PERTIMS, TOTIMS
    DOUBLE PRECISION :: PERTIMD, TOTIMD
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    DHEADS = 0.0D0
    SELECT CASE (PRECISION)
    CASE (1)
      IF (.NOT. ALLOCATED(BUFF2D)) ALLOCATE(BUFF2D(NCOL,NROW))
      READ(IU,ERR=900,END=950)KSTP,KPER,PERTIMS,TOTIMS,TEXT,NC,NR,ILAY
      DO I=1,NR
        READ(IU,ERR=900,END=950)(BUFF2D(J,I),J=1,NC)
        DO J=1,NC
          DHEADS(J,I) = BUFF2D(J,I)
        ENDDO
      ENDDO
    CASE (2)
      READ(IU,ERR=900,END=950)KSTP,KPER,PERTIMD,TOTIMD,TEXT,NC,NR,ILAY
      DO I=1,NR
        READ(IU,ERR=900,END=950)(DHEADS(J,I),J=1,NC)
      ENDDO
    END SELECT
    TEXT = ADJUSTL(TEXT)
    !
    ! Normal return
    RETURN
    !
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error extracting simulated head value(s)')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    950 CONTINUE
    CALL APPEND_MESSAGE('Binary head file missing simulated value(s)')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE GET_BINARY_HEAD_DATASET
  !-------------------------------------------------------------------
  SUBROUTINE GET_BUDGET_3D_ARRAY(IU,IPREC,NCOL,NROW,NLAY, &
                                 QUANTITY,NCB,KPER,KSTP,DBUF3D,OK)
    ! Read a data set from binary cell-by-cell budget file written by 
    ! (1) UBUDSV (3D array), or (2) UBDSV1 (3D array with time data), 
    ! or (3) UBDSV2 and UBDSVA (list), or (4) UBDSV3 (2D array with 
    ! optional 2D array of layer numbers), or (5) UBDSV4 and UBDSVB 
    ! (list).  Return data in 3D array DBUF3D.  If IPREC==1, the 
    ! budget file is assumed to be single precision; if IPREC==2, it 
    ! is assumed to be double precision.
    USE GWM_UTLS, ONLY: APPEND_MESSAGE
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NCOL, NLAY, NROW
    CHARACTER(LEN=*), INTENT(OUT) :: QUANTITY ! TEXT associated with quantity of interest
    INTEGER, INTENT(OUT) :: KPER, KSTP
    DOUBLE PRECISION, DIMENSION(NCOL,NROW,NLAY), INTENT(OUT) :: DBUF3D
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ICRL, IERR, IST, J, K, L, N, NC, NCB, NL, NLIST, &
               NR, NT, NVAL
    REAL :: DELTLOCAL, PERTIM, Q, TOTIM
    DOUBLE PRECISION :: DELTDLOCAL, PERTIMD, QD, TOTIMD
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=16) :: TEXT
    !
    OKLOCAL = .TRUE.
    DBUF3D = 0.0D0
    !
    ! Read a header record from unformatted budget file
    READ(IU,ERR=900,END=950,IOSTAT=IST) KSTP,KPER,TEXT,NC,NR,NL
    QUANTITY = ADJUSTL(TEXT)
    IF (IST.NE.0) THEN
      OKLOCAL = .FALSE.
      GOTO 900
    ENDIF
    ! Negative NLAY indicates COMPACT BUDGET
    IF (NL.LT.0) THEN
      NL=-NL
      ! Read the extra record that defines the compact type
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) NCB,DELTLOCAL,PERTIM,TOTIM
      CASE (2)
        READ(IU,END=950,ERR=900) NCB,DELTDLOCAL,PERTIMD,TOTIMD
      END SELECT
    ELSE
      NCB = 1
    ENDIF
    !
    SELECT CASE (NCB)
    CASE (1)
      ! Full 3-D array (many packages when budget file is not COMPACT)
      CALL READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (2)
      !  List without auxiliary values (SFR)
      CALL READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (3)
      !  2-D array plus a layer indicator array (ETS, EVT, RCH)
      CALL READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (4)
      !  2-D array for layer 1 (ETS, EVT, RCH)
      CALL READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (5)
      !  List with auxiliary values (WEL, DRN, RIV, GHB)
      CALL READ_LIST_AUX(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    END SELECT
    IF (IERR==1) THEN
      GOTO 900
    ELSEIF (IERR==2) THEN
      GOTO 950
    ENDIF
    !
    ! Normal return
    RETURN
    !
    ! Error handling
    900 CONTINUE
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    CALL APPEND_MESSAGE('Error extracting data from cell-by-cell budget file')
    RETURN
    950 CONTINUE
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    CALL APPEND_MESSAGE('Binary cell-by-cell budget file missing simulated value(s)')
    RETURN
  END SUBROUTINE GET_BUDGET_3D_ARRAY
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_BUD_DRN_DEPS(NCOL,NROW,NLAY,NPER,NSV,MXACTD, &
                                  MAXTS,NSTP,FNAME,AUXVAR,DEPS,DELT,OK)
    ! Extract model-simulated DRN-FLOW and/or DRN-VOLUME values from
    ! binary budget file and store in DEPS()%SIMVAL
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, FINDCELL, OPEN_FILE
    USE MMPTYPES, ONLY: DEPENDENT
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NCOL, NROW, NLAY, NPER
    CHARACTER(LEN=*), INTENT(IN) :: FNAME, AUXVAR
    INTEGER, INTENT(IN) :: NSV, MXACTD, MAXTS
    TYPE (DEPENDENT), DIMENSION(NSV),     INTENT(INOUT) :: DEPS
    DOUBLE PRECISION, DIMENSION(MAXTS,NPER), INTENT(IN) :: DELT
    INTEGER, DIMENSION(NPER), INTENT(IN) :: NSTP
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables - scalars
    INTEGER :: I, ID, II, IERR, IPREC, ISPEC, IST, IU, &
               J, JD, JJ, K, KD, KPER, KSTP, &
               MAXAUXVALUE, &
               NAUX, NC, NCB, NL, NLIST, NR
    INTEGER :: AUXINDX, IAUX
    REAL :: DELTLOCAL, PERTIM, TOTIM
    DOUBLE PRECISION :: DELTDLOCAL, PERTIMD, TOTIMD
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=16) :: DRNTEXT, TEXT
    CHARACTER(LEN=15) :: CTYPE
    ! Local variables - arrays
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ICRLS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: AUXVALS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DBUFF3D
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QVALS
    CHARACTER(LEN=16), DIMENSION(AUXDIM) :: AUXTXT
    !
    OKLOCAL = .TRUE.
    DRNTEXT = '          DRAINS'
    ! AUXINDX is index 1 in AUXVALS for aux-variable specified in mmproc.in item 20
    AUXINDX = 0 
    IPREC = 0
    MAXAUXVALUE = -2
    DO I=1,NSV
      IF (DEPS(I)%PKGID=='DRNF' .OR. DEPS(I)%PKGID=='DRNV') THEN
        IF (DEPS(I)%INDX(4)>MAXAUXVALUE) MAXAUXVALUE = DEPS(I)%INDX(4)
      ENDIF
    ENDDO
    ALLOCATE(ICRLS(MXACTD))
    ALLOCATE(AUXVALS(AUXDIM,MXACTD))
    ALLOCATE(DBUFF3D(NCOL,NROW,NLAY))
    ALLOCATE(QVALS(MXACTD))
    IF (FNAME=='') OKLOCAL = .FALSE.
    IF (OKLOCAL) THEN
      IU = UTL_GETUNIT(7,9000)
      OKLOCAL = OPEN_FILE(IU,FNAME,'OLD',.TRUE.)
      IF (OKLOCAL) THEN
        IPREC = BUDGETPRECISION(IU)
        ! Determine precision of binary budget file
        IF (IPREC<1) THEN
          CALL APPEND_MESSAGE('Error: Unable to determine precision ' // &
                              'of budget file "' // TRIM(FNAME) // '"')
          OKLOCAL = .FALSE.
        ENDIF
      ENDIF
    ENDIF
    !
    ! Iterate through data sets in budget file associated with IDRNCB
    DATASET: DO
      !
      ! Read a header record from unformatted budget file
      READ(IU,ERR=900,END=1000,IOSTAT=IST) KSTP,KPER,TEXT,NC,NR,NL
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NC*NR*NL))
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NC*NR*NL))
      IF (IST.NE.0) THEN
        OKLOCAL = .FALSE.
        GOTO 900
      ENDIF
      ! Negative NLAY indicates COMPACT BUDGET
      IF (NL.LT.0) THEN
        NL=-NL
        ! Read the extra record that defines the compact type
        SELECT CASE (IPREC)
        CASE (1)
          READ(IU,END=950,ERR=900) NCB,DELTLOCAL,PERTIM,TOTIM
        CASE (2)
          READ(IU,END=950,ERR=900) NCB,DELTDLOCAL,PERTIMD,TOTIMD
        END SELECT
      ELSE
        NCB = 1
      ENDIF
      !
      ! Read data set from budget file.
      ! When COMPACT BUDGET is specified, DRN writes 
      ! list data using UBDSV4 and UBDSVB (NCB=5)
      SELECT CASE (NCB)
      CASE (1) ! Full 3-D array (many packages when budget file is not COMPACT)
        CALL READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (2)
        ! Read list without auxiliary values and ignore.
        CALL READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (3) ! 2-D array plus a layer indicator array (ETS, EVT, RCH)
        IF (.NOT. ALLOCATED(IBUF)) ALLOCATE(IBUF(NC*NR*NL))
        CALL READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (4) ! 2-D array for layer 1 (ETS, EVT, RCH)
        CALL READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (5) ! List with auxiliary values (WEL, DRN, RIV, GHB)
               ! Data set may contain DRN flows
        CALL READ_NCB_5(IU,IPREC,MXACTD,AUXDIM,AUXTXT, &
                        AUXVALS,MAXAUXVALUE,AUXVAR,FNAME, &
                        AUXINDX,NAUX,NLIST,ICRLS,QVALS,OKLOCAL)
        IF (.NOT. OKLOCAL) GOTO 900
        IF (TEXT==DRNTEXT) THEN
          ! Data set contains DRN data
          DO II=1,NSV
            ID = DEPS(II)%INDX(2) ! ROW
            JD = DEPS(II)%INDX(1) ! COLUMN
            KD = DEPS(II)%INDX(3) ! LAYER
            CTYPE = DEPS(II)%CTYPE
            IF (CTYPE=='DRN-FLOW') THEN
              IF (KPER==DEPS(II)%INDX(5) .AND. KSTP==NSTP(KPER)) THEN
                ! Dataset is for last time step of stress period of interest
                DO JJ=1,NLIST
                  CALL FINDCELL(ICRLS(JJ),NROW,NCOL,NLAY,I,J,K)
                  IF (I==ID .AND. J==JD .AND. K==KD) THEN
                    ! Cell indices match
                    IF (AUXINDX>0) THEN
                      ! Cell may contain multiple drains. Check 
                      ! to see if this is the specified drain.
                      IAUX = NINT(AUXVALS(AUXINDX,JJ))
                      ISPEC = DEPS(II)%INDX(4)
                      IF (IAUX==ISPEC) THEN
                        DEPS(II)%SIMVAL = QVALS(JJ)
                      ENDIF
                    ELSE
                      ! Simulated equivalent combines all drains in cell
                      DEPS(II)%SIMVAL = DEPS(II)%SIMVAL + QVALS(JJ)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ELSEIF (CTYPE=='DRN-VOLUME') THEN
              IF (KPER>=DEPS(II)%INDX(5) .AND. KPER<=DEPS(II)%INDX(6)) THEN
                ! Stress period is in range (BEGSP:ENDSP) of interest
                DO JJ=1,NLIST
                  CALL FINDCELL(ICRLS(JJ),NROW,NCOL,NLAY,I,J,K)
                  IF (I==ID .AND. J==JD .AND. K==KD) THEN
                    ! Cell indices match
                    IF (AUXINDX>0) THEN
                      ! Cell may contain multiple drains. Check 
                      ! to see if this is the specified drain.
                      IAUX = NINT(AUXVALS(AUXINDX,JJ))
                      ISPEC = DEPS(II)%INDX(4)
                      IF (IAUX==ISPEC) THEN
                        DEPS(II)%SIMVAL = DEPS(II)%SIMVAL + QVALS(JJ)*DELT(KSTP,KPER)
                      ENDIF
                    ELSE
                      ! Simulated equivalent combines all drains in cell
                      DEPS(II)%SIMVAL = DEPS(II)%SIMVAL + QVALS(JJ)*DELT(KSTP,KPER)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      END SELECT
      IF (IERR==1) THEN
        GOTO 900
      ELSEIF (IERR==2) THEN
        GOTO 950
      ENDIF
    ENDDO DATASET
    !
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error reading DRN budget data')
    OKLOCAL = .FALSE.
    GOTO 1000
    !
    950 CONTINUE
    CALL APPEND_MESSAGE('Error: End of file encountered while ' // &
                        'reading DRN budget data')
    OKLOCAL = .FALSE.
    GOTO 1000
    !
    1000 CONTINUE
    !
    DEALLOCATE(ICRLS,AUXVALS,DBUFF3D,QVALS)
    IF (ALLOCATED(IBUF)) DEALLOCATE(IBUF)
    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
    IF (ALLOCATED(DBUF)) DEALLOCATE(DBUF)
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_BUD_DRN_DEPS
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_BUD_SFR_DEPS(NCOL,NROW,NLAY,NPER,NSTRM,NSV, &
                                  FNAME,ITYP,DEPS,NSTP,DBUFF3D,OK)
    ! Extract model-simulated SFR-FLOW (instream flow) or SFR-LEAK 
    ! (SW-GW exchange) values from binary budget file and store 
    ! in DEP()%SIMVAL
    !
    ! Argument ITYP meaning:
    !   1: Extract only SFR-FLOW values
    !   2: Extract only SFR-LEAK values
    !
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, OPEN_FILE
    USE MMPTYPES, ONLY: DEPENDENT
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NCOL, NROW, NLAY, NPER
    INTEGER, INTENT(IN) :: NSTRM
    INTEGER, INTENT(IN) :: NSV
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    INTEGER, INTENT(IN) :: ITYP
    TYPE (DEPENDENT), DIMENSION(NSV), INTENT(INOUT) :: DEPS
    INTEGER, DIMENSION(NPER), INTENT(IN) :: NSTP
    DOUBLE PRECISION, DIMENSION(NCOL,NROW,NLAY), INTENT(INOUT) :: DBUFF3D
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, IERR, IPREC, IST, IU, IUSFR, KPER, KSTP, NC, NCB, &
               NL, NR, NREACH
    INTEGER :: IGRID = 1
    INTEGER :: KPERSFR
    REAL :: DELTLOCAL, PERTIM, TOTIM
    DOUBLE PRECISION :: DELTDLOCAL, PERTIMD, TOTIMD
    CHARACTER(LEN=2000) :: SFRFFILE
    CHARACTER(LEN=16) :: TEXT, QUANTITY
    LOGICAL :: OKLOCAL
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QD
    CHARACTER(LEN=16), DIMENSION(2) :: BUDTEXT
    CHARACTER(LEN=15), DIMENSION(2) :: DEPTEXT
    DATA BUDTEXT/'STREAMFLOW OUT  ','  STREAM LEAKAGE'/
    DATA DEPTEXT/'SFR-FLOW       ','SFR-LEAK       '/
    !
    ALLOCATE(QD(NSTRM))
    OKLOCAL = .TRUE.
    IPREC = 0
    NREACH = 0
    IF (FNAME=='') OKLOCAL = .FALSE.
    IF (OKLOCAL) THEN
      IU = UTL_GETUNIT(7,9000)
      OKLOCAL = OPEN_FILE(IU,FNAME,'OLD',.TRUE.)
      IF (OKLOCAL) THEN
        IPREC = BUDGETPRECISION(IU)
        ! Determine precision of binary budget file
        IF (IPREC<1) THEN
          CALL APPEND_MESSAGE('Error: Unable to determine precision ' // &
                              'of budget file "' // TRIM(FNAME) // '"')
          OKLOCAL = .FALSE.
        ENDIF
      ENDIF
    ENDIF
    !
    ! Iterate through data sets in budget file associated with ISTCB2
    DATASET: DO
      !
      ! Read a header record from unformatted budget file
      READ(IU,ERR=900,END=1000,IOSTAT=IST) KSTP,KPER,TEXT,NC,NR,NL
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NC*NR*NL))
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NC*NR*NL))
      QUANTITY = ADJUSTL(TEXT)
      IF (IST.NE.0) THEN
        OKLOCAL = .FALSE.
        GOTO 900
      ENDIF
      ! Negative NLAY indicates COMPACT BUDGET
      IF (NL.LT.0) THEN
        NL=-NL
        ! Read the extra record that defines the compact type
        SELECT CASE (IPREC)
        CASE (1)
          READ(IU,END=950,ERR=900) NCB,DELTLOCAL,PERTIM,TOTIM
        CASE (2)
          READ(IU,END=950,ERR=900) NCB,DELTDLOCAL,PERTIMD,TOTIMD
        END SELECT
      ELSE
        NCB = 1
      ENDIF
      !
      ! Read data set from budget file
      ! When COMPACT BUDGET is specified, SFR writes 
      ! list data using UBDSV2 and UBDSVA (NCB=2)
      SELECT CASE (NCB)
      CASE (1) ! Full 3-D array (many packages when budget file is not COMPACT)
        CALL READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (2) ! Data set may contain SFR streamflows
        IF (TEXT==BUDTEXT(ITYP)) THEN
          ! Data set contains SFR data of indicated type.
          ! Get list of SFR data from binary file
          CALL READ_LIST_SFR(IU,IPREC,NSTRM,QD,IERR)
          ! Loop through all dependents to find SFR dependents of
          ! indicated type, and see if stress period matches.
          DO I=1,NSV
            IF (DEPS(I)%CTYPE==DEPTEXT(ITYP)) THEN
              IF (DEPS(I)%INDX(3)==KPER .AND. KSTP==NSTP(KPER)) THEN
                DEPS(I)%SIMVAL = QD(DEPS(I)%INDX(4))
              ENDIF
            ENDIF
          ENDDO
        ELSE
          ! Data set does not contain SFR data of indicated type.  
          ! Read list without auxiliary values and ignore.
          CALL READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
        ENDIF
      CASE (3) ! 2-D array plus a layer indicator array (ETS, EVT, RCH)
        IF (.NOT. ALLOCATED(IBUF)) ALLOCATE(IBUF(NC*NR*NL))
        CALL READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (4) ! 2-D array for layer 1 (ETS, EVT, RCH)
        CALL READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (5) ! List with auxiliary values (WEL, DRN, RIV, GHB)
        CALL READ_LIST_AUX(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      END SELECT
      IF (IERR==1) THEN
        GOTO 900
      ELSEIF (IERR==2) THEN
        GOTO 950
      ENDIF
    ENDDO DATASET
    !
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error reading SFR budget data')
    OKLOCAL = .FALSE.
    GOTO 1000
    !
    950 CONTINUE
    CALL APPEND_MESSAGE('Error: End of file encountered while ' // &
                        'reading SFR budget data')
    OKLOCAL = .FALSE.
    GOTO 1000
    !
    1000 CONTINUE
    DEALLOCATE(QD)
    IF (ALLOCATED(IBUF)) DEALLOCATE(IBUF)
    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
    IF (ALLOCATED(DBUF)) DEALLOCATE(DBUF)
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_BUD_SFR_DEPS
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_QNET_MNW2(NCOL,NROW,NLAY,NPER,NMQMNW,MAXDIM, &
                               MAXTS,NSTP,FNAME,AUXVAR,MNW2WELLS, &
                               DELT,OK)
                               
    ! Extract model-simulated node-by-node MNW2 Qnet (Qact) values from
    ! binary budget file, then sum by MNW well and compare to MNW2WELLS()%QDES
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, FINDCELL, OPEN_FILE
    USE MMPMNWMOD, ONLY: MNW2WELL
    USE GWM_SUBS, ONLY: FUZZY_EQUALS, EPSQNET
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NCOL, NROW, NLAY, NPER, MAXDIM
    CHARACTER(LEN=*), INTENT(IN) :: FNAME, AUXVAR
    INTEGER, INTENT(IN) :: NMQMNW, MAXTS !, MMNWMAX
    TYPE (MNW2WELL), DIMENSION(NMQMNW),     INTENT(INOUT) :: MNW2WELLS
    DOUBLE PRECISION, DIMENSION(MAXTS,NPER), INTENT(IN) :: DELT
    INTEGER, DIMENSION(NPER), INTENT(IN) :: NSTP
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables - scalars
    INTEGER :: I, ID, II, IERR, IPREC, ISPEC, IST, IU, &
               J, JD, JJ, K, KD, KPER, KSTP, &
               MAXAUXVALUE, &
               NAUX, NC, NCB, NL, NLIST, NR
    INTEGER :: AUXINDX, AUXVAL, IAUX, ICELL, NCELLS
    REAL :: DELTLOCAL, PERTIM, TOTIM
    DOUBLE PRECISION :: DELTDLOCAL, EPS, PERTIMD, QNET, TOTIMD
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=16) :: MNWTEXT, TEXT
    CHARACTER(LEN=15) :: CTYPE
    ! Local variables - arrays
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ICRLS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: AUXVALS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DBUFF3D
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QVALS
    CHARACTER(LEN=16), DIMENSION(AUXDIM) :: AUXTXT
    !
    IERR = 0
    OKLOCAL = .TRUE.
    EPS = 1.0D-5
    MNWTEXT = '            MNW2'
    ! AUXINDX is index 1 in AUXVALS for aux-variable specified in mmproc.in item 20
    AUXINDX = 0 
    IPREC = 0
    MAXAUXVALUE = -2
    ALLOCATE(ICRLS(MAXDIM))
    ALLOCATE(AUXVALS(AUXDIM,MAXDIM))
    ALLOCATE(DBUFF3D(NCOL,NROW,NLAY))
    ALLOCATE(QVALS(MAXDIM))
    IF (FNAME=='') OKLOCAL = .FALSE.
    IF (OKLOCAL) THEN
      IU = UTL_GETUNIT(7,9000)
      OKLOCAL = OPEN_FILE(IU,FNAME,'OLD',.TRUE.)
      IF (OKLOCAL) THEN
        IPREC = BUDGETPRECISION(IU)
        ! Determine precision of binary budget file
        IF (IPREC<1) THEN
          CALL APPEND_MESSAGE('Error: Unable to determine precision ' // &
                              'of budget file "' // TRIM(FNAME) // '"')
          OKLOCAL = .FALSE.
        ENDIF
      ENDIF
    ENDIF
    !
    ! Iterate through data sets in budget file associated with IWL2CB
    DATASET: DO
      !
      ! Read a header record from unformatted budget file
      READ(IU,ERR=900,END=1000,IOSTAT=IST) KSTP,KPER,TEXT,NC,NR,NL
      IF (IST.NE.0) THEN
        OKLOCAL = .FALSE.
        GOTO 900
      ENDIF
      ! Negative NLAY indicates COMPACT BUDGET
      IF (NL.LT.0) THEN
        NL=-NL
        ! Read the extra record that defines the compact type
        SELECT CASE (IPREC)
        CASE (1)
          READ(IU,END=950,ERR=900) NCB,DELTLOCAL,PERTIM,TOTIM
        CASE (2)
          READ(IU,END=950,ERR=900) NCB,DELTDLOCAL,PERTIMD,TOTIMD
        END SELECT
      ELSE
        NCB = 1
      ENDIF
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NC*NR*NL))
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NC*NR*NL))
      !
      ! Read data set from budget file.
      ! When COMPACT BUDGET is specified, DRN writes 
      ! list data using UBDSV4 and UBDSVB (NCB=5)
      SELECT CASE (NCB)
      CASE (1) ! Full 3-D array (many packages when budget file is not COMPACT)
        CALL READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (2)
        ! Read list without auxiliary values and ignore.
        CALL READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
        ! Ned TODO: for MNW2, can't ignore list without aux variables.  
        ! 2nd thought: maybe aux var. is needed, to identify well.
      CASE (3) ! 2-D array plus a layer indicator array (ETS, EVT, RCH)
        IF (.NOT. ALLOCATED(IBUF)) ALLOCATE(IBUF(NC*NR*NL))
        CALL READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (4) ! 2-D array for layer 1 (ETS, EVT, RCH)
        CALL READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUFF3D)
      CASE (5) ! List with auxiliary values (WEL, DRN, RIV, GHB)
               ! Data set may contain DRN flows
        CALL READ_NCB_5(IU,IPREC,MAXDIM,AUXDIM,AUXTXT, &
                        AUXVALS,MAXAUXVALUE,AUXVAR,FNAME, &
                        AUXINDX,NAUX,NLIST,ICRLS,QVALS,OKLOCAL)
        IF (.NOT. OKLOCAL) GOTO 900
        IF (TEXT==MNWTEXT) THEN
          ! Data set contains MNW2 data
          DO II=1,NMQMNW
            IF (MNW2WELLS(II)%IACTIVE(KPER)==1) THEN
              ! Dataset is for a time step in a stress period in which the MNW well is active,
              ! so need to process all cells in this MNW well to compute Qnet.
              !
              ! AUXVAL is an integer that uniquely identifies this MNW well
              AUXVAL = MNW2WELLS(II)%AUXVALUE
              ! Initialize Qnet, which will be summed over the cells included in this MNW well
              QNET = 0.0D0
              NCELLS = MNW2WELLS(II)%NC
              DO ICELL=1,NCELLS
                ID = MNW2WELLS(II)%INDX(ICELL,2)  ! DEPS(II)%INDX(2) ! ROW
                JD = MNW2WELLS(II)%INDX(ICELL,3)  ! DEPS(II)%INDX(1) ! COLUMN
                KD = MNW2WELLS(II)%INDX(ICELL,1)  ! DEPS(II)%INDX(3) ! LAYER
                ! Iterate through data written by NLIST calls to UBDSVB
                DO JJ=1,NLIST
                  CALL FINDCELL(ICRLS(JJ),NROW,NCOL,NLAY,I,J,K)
                  IF (I==ID .AND. J==JD .AND. K==KD) THEN
                    ! Cell indices match
!                    IF (AUXINDX>0) THEN
                      ! Cell may contain multiple MNW nodes. Check 
                      ! to see if this node belongs to the specified MNW well.
                      IAUX = NINT(AUXVALS(1,JJ))
                      IF (IAUX==AUXVAL) THEN
                        QNET = QNET + QVALS(JJ)
                      ENDIF
!                    ELSE
!                      QNET = QNET + QVALS(JJ)
!                    ENDIF
                  ENDIF
                ENDDO ! End of list of nodes included in dataset
              ENDDO ! End of current cell
              IF (.NOT. FUZZY_EQUALS(QNET,MNW2WELLS(II)%QDES,EPSQNET)) THEN
                MNW2WELLS(II)%DRY = .TRUE.
              ENDIF
            ENDIF 
          ENDDO ! End of current MNW well
        ENDIF
      END SELECT
      IF (IERR==1) THEN
        GOTO 900
      ELSEIF (IERR==2) THEN
        GOTO 950
      ENDIF
    ENDDO DATASET
    !
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error reading MNW2 budget data')
    OKLOCAL = .FALSE.
    GOTO 1000
    !
    950 CONTINUE
    CALL APPEND_MESSAGE('Error: End of file encountered while ' // &
                        'reading MNW2 budget data')
    OKLOCAL = .FALSE.
    GOTO 1000
    !
    1000 CONTINUE
    !
    DEALLOCATE(ICRLS,AUXVALS,DBUFF3D,QVALS)
    IF (ALLOCATED(IBUF)) DEALLOCATE(IBUF)
    IF (ALLOCATED(BUF)) DEALLOCATE(BUF)
    IF (ALLOCATED(DBUF)) DEALLOCATE(DBUF)
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_QNET_MNW2
  !-------------------------------------------------------------------
  SUBROUTINE READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, J, K, N, NRCL
    !
    IERR = 0
    NRCL = NC*NR*NL
    SELECT CASE (IPREC)
    CASE (1)
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NC*NR*NL))
      READ(IU,END=950,ERR=900) (BUF(I),I=1,NRCL)
    CASE (2)
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NC*NR*NL))
      READ(IU,END=950,ERR=900) (DBUF(I),I=1,NRCL)
    CASE DEFAULT
      GOTO 900
    END SELECT
    IF (.TRUE.) THEN
      ! Copy values from BUF to 3-D array
      N = 0
      DO K=1,NL
        DO I=1,NR
          DO J=1,NC
            N = N+1
            SELECT CASE (IPREC)
            CASE (1)
              DBUF3D(J,I,K) = REAL(BUF(N),8)
            CASE (2)
              DBUF3D(J,I,K) = DBUF(N)
            END SELECT
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_3D_ARRAY
  !-------------------------------------------------------------------
  SUBROUTINE READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    USE GWM_UTLS, ONLY: FINDCELL
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, ICRL, J, K, L, NLIST
    REAL :: Q
    DOUBLE PRECISION :: QD
    !
    IERR = 0
    READ(IU,END=950,ERR=900) NLIST
    DO L=1,NLIST
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) ICRL,Q
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + REAL(Q,8)
      CASE (2)
        READ(IU,END=950,ERR=900) ICRL,QD
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + QD
      END SELECT
    ENDDO
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_LIST
  !-------------------------------------------------------------------
  SUBROUTINE READ_LIST_SFR(IU,IPREC,NSTRM,QD,IERR)
    ! Read NLIST (written by UBDSV2) and reach-by-reach data
    ! written by SFR using UBDSVA
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NSTRM
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NSTRM), INTENT(OUT) :: QD
    ! Local variables
    INTEGER :: I, ICRL, J, K, L, NLIST
    REAL :: Q
    !
    QD = 0.0D0
    IERR = 0
    READ(IU,END=950,ERR=900) NLIST
    DO L=1,NLIST
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) ICRL,Q
        QD(L) = REAL(Q,8)
      CASE (2)
        READ(IU,END=950,ERR=900) ICRL,QD(L)
      END SELECT
    ENDDO
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_LIST_SFR
  !-------------------------------------------------------------------
  SUBROUTINE READ_LIST_AUX(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    USE GWM_UTLS, ONLY: FINDCELL
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, ICRL, J, K, L, N, NLIST, NVAL
    REAL :: Q
    DOUBLE PRECISION :: QD
    CHARACTER(LEN=16) :: TEXT
    !
    IERR = 0
    READ(IU,END=950,ERR=900) NVAL
    IF (NVAL.GT.1) THEN
      READ(IU,END=950,ERR=900) (TEXT,N=2,NVAL)
    ENDIF
    READ(IU,END=950,ERR=900) NLIST
    DO L=1,NLIST
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) ICRL,(QV(N),N=1,NVAL)
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + REAL(QV(1),8)
      CASE (2)
        READ(IU,END=950,ERR=900) ICRL,(QVD(N),N=1,NVAL)
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + QVD(1)
      END SELECT
    ENDDO
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_LIST_AUX
  !-------------------------------------------------------------------
  SUBROUTINE READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    USE GWM_UTLS, ONLY: FINDCELL
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, J, K, N, NRC
    !
    IERR = 0
    NRC = NC*NR
    IF (.NOT. ALLOCATED(IBUF)) ALLOCATE(IBUF(NRC))
    READ(IU,END=950,ERR=900) (IBUF(N),N=1,NRC)
    SELECT CASE (IPREC)
    CASE (1)
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NRC))
      READ(IU,END=950,ERR=900) (BUF(N),N=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,IBUF(N)) = REAL(BUF(N),8)
      ENDDO
    CASE (2)
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NRC))
      READ(IU,END=950,ERR=900) (DBUF(I),I=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,IBUF(N)) = DBUF(N)
        ENDDO
    END SELECT
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_2D_ARRAY_LAYER
  !-------------------------------------------------------------------
  SUBROUTINE READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    USE GWM_UTLS, ONLY: FINDCELL
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, J, K, N, NRC
    !
    IERR = 0
    NRC = NC*NR
    SELECT CASE (IPREC)
    CASE (1)
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NRC))
      READ(IU,END=950,ERR=900) (BUF(I),I=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,1) = REAL(BUF(N),8)
      ENDDO
    CASE (2)
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NRC))
      READ(IU,END=950,ERR=900) (DBUF(I),I=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,1) = DBUF(N)
      ENDDO
    END SELECT
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_2D_ARRAY_L1
  !-------------------------------------------------------------------
  SUBROUTINE READ_NCB_5(IU,IPREC,MAXDIM,AUXDIM,AUXTXT,AUXVALS, &
                        MAXAUXVALUE,AUXVAR,FNAME,AUXINDX,NAUX, &
                        NLIST,ICRLS,QVALS,OK)
    ! Invoke when NCB = 5 (data written by UBDSV4 and UBDSVB).
    ! Read list of values with or without auxiliary values, starting
    ! with NAUX+1 (written by UBDSV4), followed (optionally) by
    ! AUXTXT (names of AUX variables), followed by NLIST and NLIST 
    ! repetitions of ICRL, Q, and (optionally) AUX values (written
    ! by UBDSVB).  Names of AUX variables are returned in AUXTXT, 
    ! ICRL values are returned in ICRLS, and Q values are returned 
    ! in QVALS.  NAUX read from binary file is returned.
    USE GWM_UTLS, ONLY: APPEND_MESSAGE
    USE UTILITIES, ONLY: UTL_SAMENAME
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, MAXDIM, AUXDIM
    CHARACTER(LEN=16) :: AUXTXT(AUXDIM)
    DOUBLE PRECISION, DIMENSION(AUXDIM,MAXDIM) :: AUXVALS
    ! MAXAUXVALUE is the maximum aux-value read from mmproc.in.
    INTEGER, INTENT(IN) :: MAXAUXVALUE
    ! AUXVAR IS name of AUX variable that identifies boundary 
    ! (e.g. drain) of interest.
    CHARACTER(LEN=*), INTENT(IN) :: AUXVAR, FNAME
    ! AUXINDX is the index in AUXVALS that
    ! corresponds to AUX variable AUXVAR.
    INTEGER, INTENT(OUT) :: AUXINDX
    INTEGER, INTENT(OUT) :: NAUX, NLIST
    INTEGER, DIMENSION(MAXDIM), INTENT(INOUT) :: ICRLS
    DOUBLE PRECISION, DIMENSION(MAXDIM),   INTENT(INOUT) :: QVALS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, II, J, K, N, NAUXP1
    REAL :: AUXVAL, Q
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=400) :: ERRMSG
    !
    OKLOCAL = .TRUE.
    AUXINDX = 0
    !
    IF (OKLOCAL) THEN
      ! Read remaining data written by UBDSV4
      READ(IU) NAUXP1
      NAUX = NAUXP1 - 1
!      READ(IU) NAUX
      IF (NAUX>0) READ(IU) (AUXTXT(N),N=1,NAUX)
      IF (MAXAUXVALUE>-1) THEN
        ! Assign AUXINDX
        IF (NAUX>0) THEN
          DO N=1,NAUX
            IF (UTL_SAMENAME(AUXTXT(N),AUXVAR)) THEN
              AUXINDX = N
              EXIT
            ENDIF
          ENDDO
        ENDIF
        IF (AUXINDX==0) THEN
          ERRMSG = 'Required AUX variable ("' // TRIM(AUXVAR) // &
              '") not found in binary file "' // TRIM(FNAME) // &
              '". Are you missing "COMPACT BUDGET AUXILIARY"' // &
              ' in the Output Control input file?'
          CALL APPEND_MESSAGE(ERRMSG)
          OKLOCAL = .FALSE.
          GOTO 900
        ENDIF
      ENDIF
      READ(IU) NLIST
      !
      IF (NLIST>MAXDIM) THEN
        ERRMSG = 'MAXDIM dimension too small in READ_NCB_5'
        CALL APPEND_MESSAGE(ERRMSG)
        OKLOCAL = .FALSE.
        GOTO 900
      ENDIF
      !
      ! Read data written by UBDSVB
      IF (IPREC==1) THEN
        DO II=1,NLIST
          READ(IU) ICRLS(II)
          READ(IU) Q
          QVALS(II) = REAL(Q,8)
          IF(NAUX.GT.0) THEN
             DO N=1,NAUX
               READ(IU)AUXVAL
               AUXVALS(N,II) = REAL(AUXVAL,8)
             ENDDO
          ELSE
          ENDIF
        ENDDO
      ELSEIF (IPREC==2) THEN
        DO II=1,NLIST
          READ(IU) ICRLS(II),QVALS(II)
          IF(NAUX.GT.0) THEN
            READ(IU) (AUXVALS(N,II),N=1,NAUX)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
    900 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_NCB_5
  !-------------------------------------------------------------------
  INTEGER FUNCTION GET_BINARY_FILE_LENGTH(IU)
    ! Return length of a binary file, in bytes.
    ! Assume unit IU has been opened (return -1 if not).
    ! Rewind IU before return.
    IMPLICIT NONE
    ! Argument
    INTEGER, INTENT(IN) :: IU
    ! Local variables
    INTEGER :: K
    CHARACTER(LEN=1) :: BYT
    LOGICAL :: LOP, OK
    !
    INQUIRE(UNIT=IU,OPENED=LOP)
    IF (LOP) THEN
      REWIND(IU)
      K = 0
      OK = .TRUE.
      DO WHILE (OK)
        READ(IU,END=10)BYT
        K = K + 1
      ENDDO
      10 CONTINUE
      REWIND(IU)
      GET_BINARY_FILE_LENGTH = K
    ELSE
      GET_BINARY_FILE_LENGTH = -1
    ENDIF
    !
    RETURN
  END FUNCTION GET_BINARY_FILE_LENGTH
  !-------------------------------------------------------------------
  LOGICAL FUNCTION VALID_BUDGET_TEXT(TEXT)
    ! Determine if TEXT is a valid text entry in a cell-by-cell budget
    ! file.  Try to make it support all MODFLOW packages.
    USE UTILITIES, ONLY: UTL_SAMENAME
    IMPLICIT NONE
    ! Argument
    CHARACTER(LEN=*), INTENT(IN) :: TEXT
    ! Local variables
    CHARACTER(LEN=16) :: LOCALTEXT, LEFTTEXT
    LOGICAL :: VALID
    INTEGER :: I, IFAIL, TEXTSDIM, TEXTLEN
    PARAMETER (TEXTSDIM=30)
    CHARACTER(LEN=16), DIMENSION(TEXTSDIM) :: TEXTS
    DATA TEXTS/'CONSTANT HEAD   ', &
               'FLOW RIGHT FACE ', &
               'FLOW FRONT FACE ', &
               'FLOW LOWER FACE ', &
               'STORAGE         ', &
               'DRAINS          ', &
               'DRAINS (DRT)    ', &
               'ET SEGMENTS     ', &
               'ET              ', &
               'SPECIFIED FLOWS ', &
               'HEAD DEP BOUNDS ', &
               'INTERBED STORAGE', &
               'LAKE SEEPAGE    ', &
               'MNW             ', &
               'MNW2            ', &
               'RECHARGE        ', &
               'RESERV. LEAKAGE ', &
               'RIVER LEAKAGE   ', &
               'STREAM LEAKAGE  ', &
               'STREAMFLOW OUT  ', &
               'STREAM LEAKAGE  ', &
               'STREAM FLOW OUT ', &
               'INST. IB STORAGE', &
               'DELAY IB STORAGE', &
               'SWR LEAKAGE     ', &
               'SWR GWET        ', &
               'UZF RECHARGE    ', &
               'GW ET           ', &
               'SURFACE LEAKAGE ', &
               'WELLS           '/               
    !
    ! If TEXT is blank, it's invalid
    IF (TEXT==' ') THEN
      VALID_BUDGET_TEXT = .FALSE.
      RETURN
    ENDIF
    !
    VALID = .FALSE.
    TEXTLEN = LEN(TEXT)
    !
    LOCALTEXT = TEXT
    ! Left-justify text
    LEFTTEXT = ADJUSTL(LOCALTEXT)
    !
    ! Iteratively compare left-justified text with valid values
    COMPARE: DO I=1,TEXTSDIM
      IF (UTL_SAMENAME(LEFTTEXT,TEXTS(I))) THEN
        VALID = .TRUE.
        EXIT COMPARE
      ENDIF
    ENDDO COMPARE
    !
    VALID_BUDGET_TEXT = VALID
    RETURN
  END FUNCTION VALID_BUDGET_TEXT
  !-------------------------------------------------------------------
END MODULE PRECUTLSMOD