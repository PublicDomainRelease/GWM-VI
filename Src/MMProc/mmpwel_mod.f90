MODULE MMPWELMOD
  PRIVATE
  ! Public variables
  PUBLIC :: ACTIVEPARS, IWELCB, MXPW, MXPWNAUX, NAUX, NPAR, NPWEL, &
            NWELVL, WELLP, WELOPTIONLIST
  ! Public type and program units
  PUBLIC :: WELTYPE, WELTYPE_COMBINE, WELTYPE_DEALLOC, &
            WELTYPE_POPULATE, READ_WEL_FILE
  !
  ! Data related to parameter wells
  INTEGER :: NPWEL  ! Number of WEL parameters
  INTEGER :: MXPW   ! Max number of parameter wells
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NPAR   ! Number of parameters used each stress period
  REAL, ALLOCATABLE, DIMENSION(:,:) :: WELLP   ! Use only for parameter data
  CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:,:) :: ACTIVEPARS
  !
  ! Data related to non-parameter wells
  INTEGER :: IWELCB, IPRWEL
  TYPE WELTYPE
    ! Holds data related to non-parameter wells
    INTEGER  :: MXNPW  ! Max number of non-parameter wells
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ITMP 
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NNPWEL ! Number of non-parameter wells used each stress period
    ! WELL array extended to store data for all stress periods.
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WELLNP
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: DRY
  END TYPE WELTYPE
  !
  ! Data that are independent of parameter/non-parameter status
  INTEGER :: NAUX, NWELVL
  ! WELAUX and WELL arrays extended to 
  ! store data for all stress periods.
  CHARACTER(LEN=16), DIMENSION(20)  :: WELAUX
  CHARACTER(LEN=200) :: WELOPTIONLIST
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE WELTYPE_CONSTRUCT(WT,NWELVALUES,NSPER,MXNPW)
    IMPLICIT NONE
    ! Arguments
    TYPE (WELTYPE), INTENT(INOUT) :: WT
    INTEGER, INTENT(IN) :: NWELVALUES, NSPER, MXNPW
    ! Local variables
    INTEGER :: I
    !
    WT%MXNPW = MXNPW
    ! Allocate arrays to hold non-parameter data
    IF (ALLOCATED(WT%ITMP)) DEALLOCATE(WT%ITMP)
    IF (ALLOCATED(WT%NNPWEL)) DEALLOCATE(WT%NNPWEL)
    IF (ALLOCATED(WT%WELLNP)) DEALLOCATE(WT%WELLNP)
    IF (ALLOCATED(WT%DRY)) DEALLOCATE(WT%DRY)
    ALLOCATE(WT%ITMP(NSPER))
    ALLOCATE(WT%NNPWEL(NSPER))
    ALLOCATE(WT%WELLNP(NWELVALUES,WT%MXNPW,NSPER))
    ALLOCATE(WT%DRY(WT%MXNPW))
    WT%ITMP = 0
    WT%NNPWEL = 0
    WT%WELLNP = 0.0
    WT%DRY = .FALSE.
    !
    RETURN
  END SUBROUTINE WELTYPE_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE WELTYPE_POPULATE(NMQWEL,NPER,IACTIVEWEL,INDXWEL, &
                              QWEL,WT)
    ! Populate WELTYPE structure WT with managed WEL data
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NMQWEL, NPER
    INTEGER, INTENT(IN), DIMENSION(NPER,NMQWEL) :: IACTIVEWEL
    INTEGER, INTENT(IN),DIMENSION(3,NMQWEL) :: INDXWEL
    DOUBLE PRECISION, INTENT(IN), DIMENSION(NMQWEL) :: QWEL
    TYPE (WELTYPE), INTENT(INOUT) :: WT
    ! Local variables
    INTEGER :: I, K, N
    !
    CALL WELTYPE_CONSTRUCT(WT,NWELVL,NPER,NMQWEL)
    DO K=1,NPER
      N = 0
      DO I=1,NMQWEL
        IF (IACTIVEWEL(K,I)>0) THEN
          WT%WELLNP(1,I,K) = INDXWEL(3,I) ! Layer index
          WT%WELLNP(2,I,K) = INDXWEL(2,I) ! Row index
          WT%WELLNP(3,I,K) = INDXWEL(1,I) ! Column index
          WT%WELLNP(4,I,K) = QWEL(I)      ! Well pump rate (<0 for discharge)
          N = N  + 1
        ENDIF
      ENDDO
      WT%NNPWEL(K) = N
    ENDDO
    !
    RETURN
  END SUBROUTINE WELTYPE_POPULATE
  !-------------------------------------------------------------------
  SUBROUTINE WELTYPE_COMBINE(WT1,WT2,WT3)
    ! Combine WEL data stored in two WELTYPE structures.  WT3 is 
    ! generated as the combination of WT1 and WT2.
    IMPLICIT NONE
    ! Arguments
    TYPE (WELTYPE), INTENT(IN) :: WT1, WT2
    TYPE (WELTYPE), INTENT(INOUT) :: WT3
    ! Local variables
    INTEGER :: I, II, J, K, MXNPW3, NSP1, NSP2, NWVL1, NWVL2
    !
    MXNPW3 = WT1%MXNPW + WT2%MXNPW
    NSP1 = SIZE(WT1%NNPWEL)
    NSP2 = SIZE(WT2%NNPWEL)
    IF (NSP1 .NE. NSP2) THEN
      STOP 'ERROR IN WELTYPE_COMBINE: NNPWEL DIMENSION MISMATCH'
    ENDIF
    NWVL1 = SIZE(WT1%WELLNP,1)
    NWVL2 = SIZE(WT2%WELLNP,1)
    IF (NWVL1 .NE. NWVL2) THEN
      STOP 'ERROR IN WELTYPE_COMBINE: WELLNP DIMENSION MISMATCH'
    ENDIF
    CALL WELTYPE_CONSTRUCT(WT3,NWVL1,NSP1,MXNPW3)
    DO K=1,NSP1
      WT3%NNPWEL(K) = WT1%NNPWEL(K) + WT2%NNPWEL(K)
      DO J=1,NWVL1
        II = 0
        DO I=1,WT1%MXNPW
          II = II + 1
          WT3%WELLNP(J,II,K) = WT1%WELLNP(J,I,K)
        ENDDO
        DO I=1,WT2%MXNPW
          II = II + 1
          WT3%WELLNP(J,II,K) = WT2%WELLNP(J,I,K)
        ENDDO
      ENDDO
      II = 0
      DO I=1,WT1%MXNPW
        II = II + 1
        WT3%DRY(II) = WT1%DRY(I)
      ENDDO
      DO I=1,WT2%MXNPW
        II = II + 1
        WT3%DRY(II) = WT2%DRY(I)
      ENDDO
    ENDDO
    !
    RETURN
  END SUBROUTINE WELTYPE_COMBINE
  !-------------------------------------------------------------------
  SUBROUTINE READ_WEL_FILE(WELFILE,IOUT,NPER,IFREFM,NCOL,NROW,NLAY, &
                           ORIG_WELLS,OK)
    ! Open file and read WEL file.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: WELFILE
    INTEGER, INTENT(IN) :: IOUT, NPER, IFREFM, NCOL, NROW, NLAY
    TYPE (WELTYPE), INTENT(INOUT) :: ORIG_WELLS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: INWEL, KPER, MXACTW, MXNPW
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    INWEL = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(INWEL,WELFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      ! Read and store SIMULATION data
      CALL READ_WEL_SIM_DATA(INWEL,IOUT,IFREFM,NCOL,NROW,NLAY,MXACTW,OKLOCAL)
    ENDIF
    ALLOCATE(NPAR(NPER), ACTIVEPARS(NPWEL,NPER))
    !
    ! Find max number of non-parameter wells
    MXNPW = MXACTW - MXPW
    CALL WELTYPE_CONSTRUCT(ORIG_WELLS,NWELVL,NPER,MXNPW)
    !
    IF (OKLOCAL) THEN
      ! Read and store STRESS PERIOD data 
      DO KPER=1,NPER
        IF (OKLOCAL) THEN
          CALL READ_WEL_PER_DATA(INWEL,IOUT,IFREFM,KPER,NCOL,NROW,NLAY, &
                               ORIG_WELLS,OKLOCAL)
        ENDIF
      ENDDO
    ENDIF
    !
    CLOSE(INWEL)
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_WEL_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_WEL_SIM_DATA(INWEL,IOUT,IFREFM,NCOL,NROW,NLAY, &
                               MXACTW,OK)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(INOUT) :: INWEL, MXACTW
    INTEGER, INTENT(IN) :: IOUT, IFREFM, NCOL,NROW,NLAY
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, IP, ISTART, ISTOP, K, LLOC, LSTBEG, LSTSUM, N, &
               NLST, NUMINST, NINLST
    INTEGER :: IWELPB, MXWELL
    REAL :: R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .FALSE.
    WELOPTIONLIST = ''
    CALL URDCOM(INWEL,IOUT,LINE)
    CALL UPARLSTAL(INWEL,IOUT,LINE,NPWEL,MXPW)
    IF(IFREFM.EQ.0) THEN
      READ(LINE,'(2I10)') MXACTW,IWELCB
      LLOC=21
    ELSE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTW,R,IOUT,INWEL)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IWELCB,R,IOUT,INWEL)
    END IF
    !C
    !C3------READ AUXILIARY VARIABLES AND PRINT FLAG.
    NAUX=0
    IPRWEL=1
    10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INWEL)
    WELOPTIONLIST = LINE(ISTART:LEN_TRIM(LINE))
    IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR. &
             LINE(ISTART:ISTOP).EQ.'AUX') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INWEL)
      IF(NAUX.LT.20) THEN
        NAUX=NAUX+1
        WELAUX(NAUX)=LINE(ISTART:ISTOP)
        WRITE(IOUT,12) WELAUX(NAUX)
        12 FORMAT(1X,'AUXILIARY WELL VARIABLE: ',A)
      END IF
      GO TO 10
    ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
      WRITE(IOUT,13)
      13 FORMAT(1X,'LISTS OF WELL CELLS WILL NOT BE PRINTED')
      IPRWEL = 0
      GO TO 10
    END IF
    !C3A-----THERE ARE FOUR INPUT VALUES PLUS ONE LOCATION FOR
    !C3A-----CELL-BY-CELL FLOW.
    NWELVL=5+NAUX
    !C
    !C4------ALLOCATE SPACE FOR THE WELL DATA.
    ! In MMProc, only need to store parameter wells in WELLP
    IWELPB=1
    MXWELL=MXPW
    IF(MXACTW.LT.1) THEN
      WRITE(IOUT,17)
      17 FORMAT(1X, &
          'In Well Package MXACTW=0')
    END IF
    ALLOCATE (WELLP(NWELVL,MXPW)) ! Only needed for parameter wells
    !C
    !C5------READ NAMED PARAMETERS.
    WRITE(IOUT,18) NPWEL
    18 FORMAT(1X,//1X,I5,' Well parameters')
    IF(NPWEL.GT.0) THEN
      LSTSUM=IWELPB
      DO 120 K=1,NPWEL
        LSTBEG=LSTSUM
        CALL UPARLSTRP(LSTSUM,MXWELL,INWEL,IOUT,IP,'WEL','Q',1, &
                        NUMINST)
        NLST=LSTSUM-LSTBEG
        IF(NUMINST.EQ.0) THEN
          !C5A-----READ PARAMETER WITHOUT INSTANCES.
          CALL ULSTRD(NLST,WELLP,LSTBEG,NWELVL,MXWELL,1,INWEL, &
             IOUT,'WELL NO.  LAYER   ROW   COL   STRESS FACTOR', &
             WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
        ELSE
          !C5B-----READ INSTANCES.
          NINLST=NLST/NUMINST
          DO 110 I=1,NUMINST
            CALL UINSRP(I,INWEL,IOUT,IP,IPRWEL)
            CALL ULSTRD(NINLST,WELLP,LSTBEG,NWELVL,MXWELL,1,INWEL, &
                    IOUT, &
                    'WELL NO.  LAYER   ROW   COL   STRESS FACTOR', &
                    WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
            LSTBEG=LSTBEG+NINLST
          110 CONTINUE
        END IF
      120 CONTINUE
    END IF
    OKLOCAL = .TRUE.
    IF (OK) OK = OKLOCAL
    !
    RETURN
  END SUBROUTINE READ_WEL_SIM_DATA
  !-------------------------------------------------------------------
  SUBROUTINE READ_WEL_PER_DATA(INWEL,IOUT,IFREFM,KPER,NCOL,NROW, &
                               NLAY,WT,OK)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: INWEL, IOUT, IFREFM, KPER, NCOL, NROW, &
                           NLAY
    TYPE (WELTYPE), INTENT(INOUT) :: WT
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, IOUTU, J, N, NREAD
    CHARACTER(LEN=6) :: CWELL
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    !C1----READ NUMBER OF WELLS (OR FLAG SAYING REUSE WELL DATA).
    !C1----AND NUMBER OF PARAMETERS
    IF(NPWEL.GT.0) THEN
      IF(IFREFM.EQ.0) THEN
        READ(INWEL,'(2I10)') WT%ITMP(KPER),NPAR(KPER)
      ELSE
        READ(INWEL,*) WT%ITMP(KPER),NPAR(KPER)
      END IF
    ELSE
      NPAR(KPER)=0
      IF(IFREFM.EQ.0) THEN
        READ(INWEL,'(I10)') WT%ITMP(KPER)
      ELSE
        READ(INWEL,*) WT%ITMP(KPER)
      END IF
    END IF
    !C
    !C------Calculate some constants.
    NAUX=NWELVL-5
    IOUTU = IOUT
    IF (IPRWEL.EQ.0) IOUTU=-IOUTU
    !C
    !C1A-----IF ITMP LESS THAN ZERO REUSE NON-PARAMETER DATA. PRINT MESSAGE.
    !C1A-----IF ITMP=>0, SET NUMBER OF NON-PARAMETER WELLS EQUAL TO ITMP.
    IF(WT%ITMP(KPER).LT.0) THEN
      WRITE(IOUT,6)
      6 FORMAT(1X,/ &
             1X,'REUSING NON-PARAMETER WELLS FROM LAST STRESS PERIOD')
    ELSE
      WT%NNPWEL(KPER) = WT%ITMP(KPER)
    END IF
    !C
    !C1B-----IF THERE ARE NEW NON-PARAMETER WELLS, READ THEM.
    IF (WT%ITMP(KPER).GT.0) THEN
      CALL ULSTRD(WT%NNPWEL(KPER),WT%WELLNP(:,:,KPER),1,NWELVL,WT%MXNPW,1,INWEL,IOUT, &
                  'WELL NO.  LAYER   ROW   COL   STRESS RATE', &
                  WELAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,4,IPRWEL)
    ELSEIF (WT%ITMP(KPER).LT.0) THEN
      ! Copy data from previous stress period
      WT%NNPWEL(KPER) = WT%NNPWEL(KPER-1)
      DO J=1,WT%MXNPW
        DO I=1,NWELVL
          WT%WELLNP(I,J,KPER) = WT%WELLNP(I,J,KPER-1)
        ENDDO
      ENDDO
    END IF
    !
    ! If parameters are active this stress period, save the item-7 lines
    DO N=1,NPAR(KPER)
      READ(INWEL,'(A80)')ACTIVEPARS(N,KPER)
    ENDDO
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_WEL_PER_DATA
  !-------------------------------------------------------------------
  SUBROUTINE WELTYPE_DEALLOC(WT)
    IMPLICIT NONE
    ! Argument
    TYPE (WELTYPE), INTENT(INOUT) :: WT
    IF (ALLOCATED(WT%ITMP)) DEALLOCATE(WT%ITMP)
    IF (ALLOCATED(WT%NNPWEL)) DEALLOCATE(WT%NNPWEL)
    IF (ALLOCATED(WT%WELLNP)) DEALLOCATE(WT%WELLNP)
    IF (ALLOCATED(WT%DRY)) DEALLOCATE(WT%DRY)
    RETURN
  END SUBROUTINE WELTYPE_DEALLOC
  !-------------------------------------------------------------------
END MODULE MMPWELMOD
