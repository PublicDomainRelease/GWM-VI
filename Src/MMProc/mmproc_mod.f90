MODULE MMPROCMOD
  USE MMPTYPES
  USE MMPWELMOD, ONLY: WELTYPE
  USE GWM_UTLS, ONLY: APPEND_MESSAGE
  PRIVATE
  ! PUBLIC VARIABLES
  PUBLIC :: COMMAND, DEPS, IACTIVEMNW, IACTIVEWEL, INDXWEL, IOUT, &
            NAMEFILE, NIUNIT, NMQWEL, MAXDEPPKGS, NSV, &
            QMNW, QWEL, WELLIDMNW, WELLIDWEL, NPER, NSTP, IFREFM, &
            NCOL, NROW, NLAY, WELFILE, NHEDDEP, NSFRFDEP, &
            NSFRLDEP, NSTORDEP, NSTRFDEP, NSTRLDEP, NDRNFDEP, &
            NDRNVDEP, MXACTD
  PUBLIC :: AUXVAR, DRNFILE, MNW2FILE
  ! PUBLIC PROCEDURES
  PUBLIC :: EXTRACT_HEADS, EXTRACT_BUD_GRID_DEPS, FIND_BINARY_FILES, &
            MMP_ALLOC, MMP_INI, READ_DEPENDENTS, READ_INPUT_1, &
            MMP_READ_DRN_FILE, READ_WEL_DATA, RESTORE_WEL_FILE, &
            SAVE_WEL_FILE, WRITE_STATUS, WRITE_WEL_FILE, &
            WRITE_DEPENDENT_DATA, INITIALIZE_MMPMOD, READ_FILENAMES, &
            SFR_IDENTIFY, EXTRACT_SFR_DEPS, MMP_READ_NAME_FILE, &
            EXTRACT_DRN_DEPS, MMP_READ_MNW2_FILE, CHECK_MNW_QNET
  INTEGER, PARAMETER :: IDDIM = 20
  INTEGER, PARAMETER :: MAXDEPPKGS = 20
  INTEGER, PARAMETER :: MAXUNITS = 100
  ! Scalars
  INTEGER :: IOUT, ICHFLG, IFREFM, IPRTIM, ITRSS, IXSEC, KCP, NBOTM, &
             NCOL, NLAY, NMQWEL, NPER, NROW, NSV
  INTEGER :: IBCFCB, IBOUUN, IHEDUN, IHUFCB, ILPFCB, IUPWCB, ISTCB1, &
             ISTCB2, IFLOWCB, IDRNCB, IMNWCB, MNWMAX=0, NODTOT, MXACTD, MAXTS
  INTEGER :: NHEDDEP, NSFRFDEP, NSFRLDEP, NSTORDEP, NSTRFDEP, &
             NSTRLDEP
  REAL :: HCLOSE, HDRY, HNOFLO
  CHARACTER(LEN=2000) :: SIMVALOUT = 'SimulatedValues.out'
  CHARACTER(LEN=2000) :: COMMAND, NAMEFILE
  CHARACTER(LEN=2000) :: BASFILE, BCFFILE, DISFILE, DE4FILE, DRNFILE, &
                         LMGFILE, GMGFILE, HUFFILE, LISTFILE, &
                         LPFFILE, MNW2FILE, NWTFILE, OCFILE, PCGFILE, &
                         PCGNFILE, SFRFILE, SIPFILE, SORFILE, STRFILE, &
                         UPWFILE, WELFILE, WELFILECOPY
  CHARACTER(LEN=2000) :: HEADFILE, FLOWBUDFILE, SFRLEAKBUDFILE, &
                         SFRFLOWBUDFILE, DRNBUDFILE
  CHARACTER(LEN=200) :: BASOPTIONLIST
  CHARACTER(LEN=16) :: AUXVAR
  ! Arrays
  INTEGER, DIMENSION(MAXUNITS)           :: IUNIT
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: IBUF
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: NSTP
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: LBOTM, LAYCBD, ISSFLG
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: IACTIVEMNW
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: IACTIVEWEL
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: INDXWEL  ! (C,R,L)
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: IOFLG
  INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: ISTRMSFR
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IBOUND
  REAL, ALLOCATABLE, DIMENSION(:) :: BUF
  REAL, ALLOCATABLE, DIMENSION(:) :: DELR, DELC 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PERLEN, TSMULT
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BOTM
  REAL, ALLOCATABLE, DIMENSION(:,:) :: BUFF2D
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BUFF3D
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DBUF
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DBUFF2D
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DELT  ! Dimension: (MAXTS,NPER)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DBUFF3D
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     :: QMNW, QWEL
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)   :: STRMSFR
  CHARACTER(LEN=200), DIMENSION(MAXUNITS)          :: FNAME, FTYPE
  CHARACTER (LEN=IDDIM), ALLOCATABLE, DIMENSION(:) :: WELLIDMNW
  CHARACTER (LEN=IDDIM), ALLOCATABLE, DIMENSION(:) :: WELLIDWEL
  TYPE (DEPENDENT), ALLOCATABLE, DIMENSION(:) :: DEPS
  TYPE (DEPENDENT_PKG), DIMENSION(MAXDEPPKGS) :: DEPPKGS
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE MMP_INI()
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    IMPLICIT NONE
    ! Local variables
    INTEGER :: ISTAT, IU, MAXTRIES
    DOUBLE PRECISION :: WAITSECS
    !    
    MAXTRIES = 1000
    WAITSECS = 0.1
    IU = UTL_GETUNIT(7,9000)
    ISTAT = PLL_CLOSE_AND_DELETE(IU,SIMVALOUT,MAXTRIES,WAITSECS)
    RETURN
  END SUBROUTINE MMP_INI
  !-------------------------------------------------------------------
  SUBROUTINE MMP_READ_DRN_FILE(OK)
    ! Open DRN file and read IDRNCB
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_RWORD
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(OUT) :: OK
    ! Local variables
    INTEGER :: ICOL, ISTART, ISTOP, ISTAT, IU, N
    DOUBLE PRECISION :: R
    LOGICAL :: LEX
    CHARACTER(LEN=2000) :: ERRMSG, LINE, WORD
    !
    OK = .TRUE.
    ERRMSG = ''
    INQUIRE(FILE=DRNFILE,EXIST=LEX)
    IF (LEX) THEN
      IU = UTL_GETUNIT(7,9000)
      OPEN(IU,FILE=DRNFILE,STATUS='OLD',IOSTAT=ISTAT)
      IF (ISTAT==0) THEN
        10 CONTINUE
        READ(IU,'(A2000)')LINE
        IF (LINE(1:1)=='#') GOTO 10
        ICOL = 1
        CALL UTL_RWORD(0,1,.FALSE.,ICOL,LINE,ISTART,ISTOP,N,R)
        IF (LINE(ISTART:ISTOP)=='PARAMETER') GOTO 10
        ICOL = 1
        CALL UTL_RWORD(0,2,.FALSE.,ICOL,LINE,ISTART,ISTOP,MXACTD,R)
        CALL UTL_RWORD(0,2,.FALSE.,ICOL,LINE,ISTART,ISTOP,IDRNCB,R)
      ELSE
        ERRMSG = 'Error opening file: ' // TRIM(DRNFILE)
        OK = .FALSE.
      ENDIF
    ELSE
      ERRMSG = 'File does not exist: ' // TRIM(DRNFILE)
      OK = .FALSE.
    ENDIF
    IF (.NOT. OK) CALL APPEND_MESSAGE(ERRMSG)
    RETURN
  END SUBROUTINE MMP_READ_DRN_FILE
  !-------------------------------------------------------------------
  SUBROUTINE MMP_READ_MNW2_FILE(OK)
    ! Open MNW2 file and read MNWMAX,[NODTOT] and IWL2CB (IMNWCB)
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_RWORD
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(OUT) :: OK
    ! Local variables
    INTEGER :: ICOL, ISTART, ISTOP, ISTAT, IU, N
    DOUBLE PRECISION :: R
    LOGICAL :: LEX
    CHARACTER(LEN=2000) :: ERRMSG, LINE, WORD
    !
    OK = .TRUE.
    ERRMSG = ''
    INQUIRE(FILE=MNW2FILE,EXIST=LEX)
    IF (LEX) THEN
      IU = UTL_GETUNIT(7,9000)
      OPEN(IU,FILE=MNW2FILE,STATUS='OLD',IOSTAT=ISTAT)
      IF (ISTAT==0) THEN
        10 CONTINUE
        READ(IU,'(A2000)')LINE
        IF (LINE(1:1)=='#') GOTO 10
        ICOL = 1
!        CALL UTL_RWORD(0,1,.FALSE.,ICOL,LINE,ISTART,ISTOP,N,R)
!        IF (LINE(ISTART:ISTOP)=='PARAMETER') GOTO 10
!        ICOL = 1
        CALL UTL_RWORD(0,2,.FALSE.,ICOL,LINE,ISTART,ISTOP,MNWMAX,R)
        ! If MNWMAX<0 NODTOT is present 
        IF(MNWMAX.LT.0)CALL UTL_RWORD(0,2,.FALSE.,ICOL,LINE,ISTART,ISTOP,NODTOT,R) 
        CALL UTL_RWORD(0,2,.FALSE.,ICOL,LINE,ISTART,ISTOP,IMNWCB,R)
      ELSE
        ERRMSG = 'Error opening file: ' // TRIM(MNW2FILE)
        OK = .FALSE.
      ENDIF
    ELSE
      ERRMSG = 'File does not exist: ' // TRIM(MNW2FILE)
      OK = .FALSE.
    ENDIF
    IF (.NOT. OK) CALL APPEND_MESSAGE(ERRMSG)
    RETURN
  END SUBROUTINE MMP_READ_MNW2_FILE
  !-------------------------------------------------------------------
  SUBROUTINE MMP_READ_SFR_FILE(OK)
    ! Open file and read ISTCB1, ISTCB2, and list of reaches
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    USE GWFSFRMODULE, ONLY: NSTRM
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, LLOC, ISTART, ISTOP, NSS, NSFRPAR, &
               NPARSEG, IN, IU
    INTEGER :: ISFROPT, IUZT, IRTFLG, NUMTIM, NSTRAIL, ISUZN, &
               NSFRSETS, NSSAR, NSTRMAR, NSEGDIM, NSFRAUX, MAXPTS, &
               IUNITGWT, NUZST, NSTOTRL, NUMAVE, MSTRMAR, NSEG, NREACH, &
               II, KRCH, IRCH, JRCH, JSEG, IREACH
    DOUBLE PRECISION :: STRIN, STROUT, FXLKOT
    REAL :: CONST, DLEAK, FLWTOL, R, WEIGHT 
    DOUBLE PRECISION :: HSTRM, QSTRM, HWDTH, HWTPRM, &
                        DLKOTFLW, SLKOTFLW, DLKSTAGE, SFRQ, CONCQ, &
                        CONCRUN, CONCPPT, SGOTFLW, DVRSFLW
    CHARACTER(LEN=200) :: LINE
    CHARACTER(LEN=16), DIMENSION(20) :: SFRAUX
    INTEGER, ALLOCATABLE, DIMENSION(:) :: IOTSG, NSEGCK
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISEG, IDIVAR
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: THTS, THTR, THTI, EPS, &
                                                    UHC, SFRUZBD
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: XSEC, QSTAGE, SEG
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    IN = IU
    OKLOCAL = OPEN_FILE(IU,SFRFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU, IOUT, LINE)
      LLOC = 1
      ISFROPT = 0
      IUZT = 0
      IRTFLG = 0
      NUMTIM = 1
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NSTRM, R, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NSS, R, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NSFRPAR, R, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NPARSEG, R, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, CONST, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, DLEAK, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, ISTCB1, R, IOUT, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, ISTCB2, R, IOUT, IU)
      !
      !C3------READ ISFROPT FLAGS WHEN NSTRM IS LESS THAN ZERO.
      IF ( NSTRM.LT.0 ) THEN
        NSTRM = ABS(NSTRM)
        CALL URWORD(line, lloc, istart, istop, 2, ISFROPT, r, IOUT, In)
        !C
        !C4------READ UNSATURATED FLOW VARIABLES WHEN ISFROPT GREATER THAN 1.
        IF ( ISFROPT.GE.2 ) THEN
          IUZT = 1
          CALL URWORD(line, lloc, istart, istop, 2, NSTRAIL, r, IOUT,  &
                       In)
          CALL URWORD(line, lloc, istart, istop, 2, ISUZN, r, IOUT, In)
          CALL URWORD(line, lloc, istart, istop, 2, NSFRSETS, r, IOUT, &
                       In)
        END IF
        !4b-----Data read for transient routing.
        CALL URWORD(line, lloc, istart, istop, 2, IRTFLG, r, IOUT, In)
        IF ( IRTFLG .GT. 0 ) THEN
          NUMTIM = 1
          WEIGHT = 1.0
          FLWTOL = 1.0D-6
          CALL URWORD(line, lloc, istart, istop, 2, NUMTIM, r, IOUT, In)
          CALL URWORD(line, lloc, istart, istop, 3, i, WEIGHT, IOUT, In)
          CALL URWORD(line, lloc, istart, istop, 3, i, FLWTOL, IOUT, In)
          IF ( NUMTIM.LT.1 ) NUMTIM = 1
          IF ( WEIGHT.LT.0.0 .OR. WEIGHT.GT.1.0 ) WEIGHT=1.0
          IF ( FLWTOL.LT.0.0 ) FLWTOL=0.0
        ELSE
          NUMTIM = 1
          WEIGHT = 1.0
          FLWTOL = 0.0
        END IF
      END IF
      IF ( NSS.LT.0 ) NSS = 0
      IF ( NSFRPAR.LE.0 ) THEN
        NSFRPAR = 0
        nparseg = 0
      END IF
      IF ( nparseg.LT.0 ) nparseg = 0
      nssar = 1
      IF (NSS.GT.0) nssar = NSS
      nstrmar = 1
      IF (NSTRM.GT.0) nstrmar = NSTRM
      nsegdim = NSS + nparseg
      IF (nsegdim.LT.1) nsegdim = 1
      !C
      !C4C-----READ AUXILLARY VARIABLES
      NSFRAUX=0
      10 CONTINUE
      CALL URWORD(line,lloc,istart,istop,1,i,r,IOUT,In)
      IF(line(istart:istop).EQ.'AUXILIARY' .OR. &
               line(istart:istop).EQ.'AUX') THEN
        CALL URWORD(line,lloc,istart,istop,1,i,r,IOUT,In)
        IF(NSFRAUX.LT.20) THEN
          NSFRAUX=NSFRAUX+1
          SFRAUX(NSFRAUX)=line(istart:istop)
        END IF
        GO TO 10
      ENDIF
      !C
      !C5------CALCULATE SPACE NEEDED FOR TABULATED DISCHARGE VERSUS FLOW
      !C         AND WIDTH RELATIONS.
      MAXPTS = 3*50
      !C
      !C     ALLOCATE ARRAY STORAGE FOR STREAMS
      STRIN = 0.0
      STROUT = 0.0
      FXLKOT = 0.0
      ALLOCATE (ISTRMSFR(5,nstrmar))
      ALLOCATE (STRMSFR(30,nstrmar))
      STRMSFR = 0.0
      HSTRM = 0.0
      QSTRM = 0.0
      HWDTH = 0.0
      HWTPRM = 0.0
      ISTRMSFR = 0
      ALLOCATE (SEG(27+NSFRAUX,nsegdim), ISEG(4,nsegdim),  &
               IDIVAR(2,nsegdim))
      SEG = 0.0
      ISEG = 0
      IDIVAR = 0
      DLKOTFLW = 0.0D0
      DLKSTAGE = 0.0D0
      SLKOTFLW = 0.0D0
      ALLOCATE (IOTSG(nsegdim))
      IOTSG = 0
      SFRQ = 0.0
      !C
      IF ( IUZT.EQ.1 ) THEN
        IF ( NSTRAIL.LT.0 ) THEN
          NSTRAIL = ABS(NSTRAIL)
        END IF
        IF ( NSTRAIL.EQ.0 ) THEN
          IUZT = 0
        END IF
      END IF
      IF ( DLEAK.LE.0.0 ) THEN
        DLEAK = 0.00001
      END IF
      !C
      IF ( IUZT.EQ.1 ) THEN
        !C
        !C8------ALLOCATE SPACE FOR UNSATURATED FLOW.
        NUZST = NSTRM
        NSTOTRL = ISUZN*NSTRAIL*NSFRSETS
        NUMAVE = 21
        mstrmar = nstrmar
      ELSE
      !C
      !C9------ALLOCATE ONLY ONE ARRAY ELEMENT IF UNSATURATED FLOW IS INACTIVE.
        NUZST = 1
        NSTOTRL = 1
        NUMAVE = 1
        ISUZN = 1
        NSTRAIL = 1
        NSFRSETS = 1
        mstrmar = 1
      END IF
      !C
      !C ALLOCATE AND INITIALIZE ARRAYS
      !C
      ALLOCATE (THTS(NUZST), THTR(NUZST), THTI(NUZST), EPS(NUZST))
      THTS = 0.0D0
      THTR = 0.0D0
      THTI = 0.0D0
      EPS = 0.0D0
      ALLOCATE (UHC(NUZST))
      UHC = 0.0
      ALLOCATE (XSEC(16, nsegdim), QSTAGE(MAXPTS,nsegdim)) 
      XSEC = 0.0
      QSTAGE = 0.0
      ALLOCATE (NSEGCK(nssar))
      NSEGCK = 0
      SGOTFLW = 0.0
      DVRSFLW = 0.0
      ALLOCATE (SFRUZBD(10))
      SFRUZBD = 0.0
      !C
      !C11-----READ AND WRITE DATA FOR EACH REACH ON BASIS OF ISFROPT.
      nseg = 0
      nreach = 0
      ! READ REACH INFO (ITEM 2), LINES 353-478 IN GWF2SFR7AR
      DO ii = 1, NSTRM
        IF ( ISFROPT.EQ.0 ) THEN
          READ (In, *) krch, irch, jrch, jseg, ireach, STRMSFR(1, ii)
        ELSE IF ( ISFROPT.EQ.1 ) THEN
          READ (In, *) krch, irch, jrch, jseg, ireach, STRMSFR(1, ii), &
                        STRMSFR(3, ii), STRMSFR(2, ii), STRMSFR(8, ii), &
                        STRMSFR(6, ii)
          STRMSFR(4, ii) = STRMSFR(3, ii) - STRMSFR(8, ii)
        ELSE IF ( ISFROPT.EQ.2 ) THEN
          READ (In, *) krch, irch, jrch, jseg, ireach, STRMSFR(1, ii), &
                        STRMSFR(3, ii), STRMSFR(2, ii), STRMSFR(8, ii), &
                        STRMSFR(6, ii), THTS(ii), THTI(ii), EPS(ii)
          STRMSFR(4, ii) = STRMSFR(3, ii) - STRMSFR(8, ii)
        ELSE IF ( ISFROPT.EQ.3 ) THEN
          READ (In, *) krch, irch, jrch, jseg, ireach, STRMSFR(1, ii), &
                        STRMSFR(3, ii), STRMSFR(2, ii), STRMSFR(8, ii), &
                        STRMSFR(6, ii), THTS(ii), THTI(ii), EPS(ii), UHC(ii)
          STRMSFR(4, ii) = STRMSFR(3, ii) - STRMSFR(8, ii)
        ELSE IF ( ISFROPT.EQ.4 .OR. ISFROPT.EQ.5 ) THEN
          READ (In, *) krch, irch, jrch, jseg, ireach, STRMSFR(1, ii)
        END IF
        !C
        !C13-----CHECK RANGE AND ORDER FOR SEGMENTS AND REACHES.
        IF ( jseg.NE.nseg ) THEN
          nseg = nseg + 1
          nreach = 0
        END IF
        nreach = nreach + 1
        ISTRMSFR(1, ii) = krch
        ISTRMSFR(2, ii) = irch
        ISTRMSFR(3, ii) = jrch
        ISTRMSFR(4, ii) = jseg
        ISTRMSFR(5, ii) = ireach
        SEG(1, ISTRMSFR(4, ii)) = SEG(1, ISTRMSFR(4, ii)) + STRMSFR(1, ii)
        !C       Number of reaches in segment added to ISEG
        ISEG(4, jseg) = ireach
      END DO
      !
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE MMP_READ_SFR_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_DEPENDENTS(IU,OK)
    ! Read dependents: mmproc.in items 17 and 18.
    ! Allocate DEPS array and read all information about simulated
    ! values needed for dependents.  Also find unit number and file
    ! name for binary file from which simulated values are to be 
    ! extracted for each unique DEPENDENT type.
    USE MMPTYPES, ONLY: DEPENDENT, DEPENDENT_READ
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I
    CHARACTER(LEN=40) :: MESS
    LOGICAL :: OKLOCAL
    ! Format statements
    100 FORMAT('NSV: ',I6)
    !
    OKLOCAL = .TRUE.
    !
    MESS = 'reading NSV'
    READ(IU,*,ERR=950)NSV
    WRITE(IOUT,100)NSV
    !
    ALLOCATE(DEPS(NSV))
    KCP = 0     ! Count unique DEPENDENT packages
    DO I=1,NSV
      CALL DEPENDENT_READ(IU,NCOL,NROW,NLAY,DEPS(I),OKLOCAL)
      IF (.NOT. OKLOCAL) GOTO 900
      SELECT CASE (DEPS(I)%CTYPE)
      CASE ('HEAD')
        NHEDDEP = NHEDDEP + 1
      CASE ('DRN-FLOW')
        NDRNFDEP = NDRNFDEP + 1
      CASE ('DRN-VOLUME')
        NDRNVDEP = NDRNVDEP + 1
      CASE ('SFR-FLOW')
        NSFRFDEP = NSFRFDEP + 1
      CASE ('SFR-LEAK')
        NSFRLDEP = NSFRLDEP + 1
      CASE ('STORAGE-CHANGE')
        NSTORDEP = NSTORDEP + 1
      CASE DEFAULT
      END SELECT
      IF (IS_NEW_DEPPKG(DEPS(I))) THEN
        KCP = KCP + 1
        DEPPKGS(KCP)%PKGID = DEPS(I)%PKGID
        DEPPKGS(KCP)%UNIT = GET_BINARY_UNIT(DEPPKGS(KCP)%PKGID)
        CALL GET_FILENAME(DEPPKGS(KCP)%UNIT,DEPPKGS(KCP)%FILENAME)
      ENDIF
    ENDDO
    900 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
    !
    ! Error handling
    950 CONTINUE
    CALL APPEND_MESSAGE('Error ' // TRIM(MESS))
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    
  END SUBROUTINE READ_DEPENDENTS
  !-------------------------------------------------------------------
  SUBROUTINE READ_INPUT_1(IU,ORIG_WELLS,MANAGED_WELLS,NMQMNW,OK)
    ! Read, store, and echo items 1-14 from main input file
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_RWORD
    USE MMPWELMOD, ONLY: READ_WEL_FILE
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU
    TYPE (WELTYPE), INTENT(INOUT) :: ORIG_WELLS, MANAGED_WELLS
    INTEGER, INTENT(OUT) :: NMQMNW
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ICOL, ISTART, ISTOP, N, NABF
    REAL :: R
    CHARACTER(LEN=40) :: FTYPE, MESS
    CHARACTER(LEN=2000) :: FNAME, LINE
    LOGICAL :: OKLOCAL
    ! Format statements
    20 FORMAT('COMMAND: ',A)
    40 FORMAT('NAME FILE: ',A)
    45 FORMAT('LIST FILE: ',A)
    50 FORMAT('NLAY: ',I3,/, &
              'NROW: ',I6,/, &
              'NCOL: ',I6,/, &
              'NPER: ',I6)
    55 FORMAT('HNOFLO: ',G14.7)
    60 FORMAT('NMQWEL: ',I5)
    80 FORMAT('NMQMNW: ',I5)
    90 FORMAT('MNW2FILE: 'A)
    100 FORMAT('NSV: ',I6)
    110 FORMAT('IFREFM: ',I1)
    120 FORMAT('HDRY: ',G14.7)
    130 FORMAT('HCLOSE: ',G14.7)
    140 FORMAT('HEADFILE: ',A)
    150 FORMAT('FLOWBUDFILE: ',A)
    160 FORMAT('NABF: ',I1)
    170 FORMAT('WELFILE: ',A)
    180 FORMAT('Binary budget file containing data of type ',A,':  ',A)
    !
    OKLOCAL = .TRUE.
    NMQMNW = 0
    !
    ! Item 1
    MESS = 'reading COMMAND'
    READ(IU,*,ERR=900)COMMAND
    WRITE(IOUT,20)TRIM(COMMAND)
    !
    ! Item 2
    MESS = 'reading NAMEFILE'
    READ(IU,*,ERR=900)NAMEFILE
    WRITE(IOUT,40)TRIM(NAMEFILE)
    !
    ! Item 2a
    MESS = 'reading LISTFILE'
    READ(IU,*,ERR=900)LISTFILE
    WRITE(IOUT,45)TRIM(LISTFILE)
    !
    ! Item 2b
    MESS = 'reading NLAY NROW NCOL NPER'
    READ(IU,*,ERR=900)NLAY,NROW,NCOL,NPER
    WRITE(IOUT,50)NLAY,NROW,NCOL,NPER
    !
    ! Item 2c
    ALLOCATE(NSTP(NPER))
    MESS = 'reading NSTP'
    READ(IU,*,ERR=900)(NSTP(I),I=1,NPER)
    WRITE(IOUT,'(A)')'NSTP values have been read'
    !
    ! Item 3
    MESS = 'reading HNOFLO'
    READ(IU,*,ERR=900)HNOFLO
    WRITE(IOUT,55)HNOFLO
    !
    ! Item 4
    MESS = 'reading IFREFM'
    READ(IU,*,ERR=900)IFREFM
    WRITE(IOUT,110)IFREFM
    !
    ! Item 5
    MESS = 'reading HDRY'
    READ(IU,*,ERR=900)HDRY
    WRITE(IOUT,120)HDRY
    !
    ! Item 6
    MESS = 'reading HCLOSE'
    READ(IU,*,ERR=900)HCLOSE
    WRITE(IOUT,130)HCLOSE
    !
    ! Item 7
    MESS = 'reading HEADFILE'
    READ(IU,*,ERR=900)HEADFILE
    WRITE(IOUT,140)TRIM(HEADFILE)
    !
    ! Item 8
    MESS = 'reading FLOWBUDFILE'
    READ(IU,*,ERR=900)FLOWBUDFILE
    WRITE(IOUT,150)TRIM(FLOWBUDFILE)
    !
    ! Item 9
    MESS = 'reading NABF'
    READ(IU,*,ERR=900)NABF
    WRITE(IOUT,160)NABF
    !
    ! Item(s) 10
    IF (NABF>0) THEN
      ! Read names of additional binary budget files
      DO I=1,NABF
        READ(IU,'(A2000)',ERR=900)LINE
        ICOL = 1
        ! Read budget file type and name
        CALL UTL_RWORD(IOUT,1,.TRUE.,ICOL,LINE,ISTART,ISTOP,N,R)
        FTYPE = LINE(ISTART:ISTOP)
        CALL UTL_RWORD(IOUT,0,.TRUE.,ICOL,LINE,ISTART,ISTOP,N,R)
        FNAME = LINE(ISTART:ISTOP)
        SELECT CASE (FTYPE)
          CASE ('SFR-LEAK')
            SFRLEAKBUDFILE = FNAME
          CASE ('SFR-FLOW')
            SFRFLOWBUDFILE = FNAME
          CASE ('DRN')
            DRNBUDFILE = FNAME
          CASE DEFAULT
            MESS = 'Unknown binary budget file type: ' // TRIM(FTYPE)
            GOTO 900
        END SELECT
        WRITE(IOUT,180)TRIM(FTYPE),TRIM(FNAME)
      ENDDO
    ENDIF
    !
    ! Item 11
    MESS = 'reading NMQWEL'
    READ(IU,*,ERR=900)NMQWEL
    WRITE(IOUT,60)NMQWEL
    !
    IF (NMQWEL>0) THEN
      ! Item 12
      MESS = 'reading WELFILE'
      READ(IU,*,ERR=900)WELFILE
      WRITE(IOUT,170)TRIM(WELFILE)
      ! Read WELL input file
      CALL READ_WEL_FILE(WELFILE,IOUT,NPER,IFREFM,NCOL,NROW,NLAY,ORIG_WELLS,OKLOCAL)
      IF (.NOT. OKLOCAL) GOTO 900
      ! Items 13a and 13b
      CALL READ_WEL_DATA(IU,MANAGED_WELLS,OKLOCAL)
      IF (.NOT. OKLOCAL) GOTO 900
    ENDIF
    !
    ! Item 14
    MESS = 'reading NMQMNW'
    READ(IU,*,ERR=900)NMQMNW
    WRITE(IOUT,80)NMQMNW
    !
!    IF (NMQMNW>0) THEN
!      ! Item 15
!      MESS = 'reading MNW2 input file name'
!      READ(IU,*,ERR=900)MNW2FILE
!      WRITE(IOUT,90)TRIM(MNW2FILE)
!    ENDIF
    !
    WRITE(IOUT,'(1X)')
    IF (OK) OK = OKLOCAL
    RETURN
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error ' // TRIM(MESS))
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_INPUT_1
  !-------------------------------------------------------------------
  SUBROUTINE READ_WEL_DATA(IU,MANAGED_WELLS,OK)
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE MMPWELMOD, ONLY: WELTYPE, WELTYPE_POPULATE
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU
    TYPE (WELTYPE), INTENT(INOUT) :: MANAGED_WELLS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, J
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    ALLOCATE(WELLIDWEL(NMQWEL),IACTIVEWEL(NPER,NMQWEL), &
             INDXWEL(3,NMQWEL),QWEL(NMQWEL))
    DO I=1,NMQWEL
      READ(IU,*,ERR=900)WELLIDWEL(I),(INDXWEL(J,I),J=3,1,-1),QWEL(I)
      READ(IU,*,ERR=920)(IACTIVEWEL(J,I),J=1,NPER)
    ENDDO
    !
    ! Store managed-well data in WELTYPE structure 
    CALL WELTYPE_POPULATE(NMQWEL,NPER,IACTIVEWEL,INDXWEL,QWEL,MANAGED_WELLS)
    !
    IF (OK) OK = OKLOCAL
    RETURN
    !
    ! Error handling
    900 CONTINUE
    WRITE(AMESSAGE,910)I
    910 FORMAT('Error reading WEL data item ',I5)
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
    920 CONTINUE
    WRITE(AMESSAGE,930)TRIM(WELLIDWEL(I))
    930 FORMAT('Error reading WEL data for well with ID: ',A)
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
  END SUBROUTINE READ_WEL_DATA
  !-------------------------------------------------------------------
  SUBROUTINE INITIALIZE_MMPMOD()
    USE PARAMMODULE, ONLY: MXPAR, MXCLST, MXINST, B, IPLOC, IPCLST, &
                           PARNAM, PARTYP, INAME, IACTIVE, IPSUM, &
                           ICLSUM, INAMLOC, NPVAL
    USE GWFSFRMODULE, ONLY: NSTRM
    IMPLICIT NONE
    ! Local variable
    INTEGER :: I
    !
    ! Integer
    IBCFCB = 0
    IHEDUN = 0
    IHUFCB = 0
    ILPFCB = 0
    IUPWCB = 0
    IFLOWCB = 0
    ISTCB1 = 0
    ISTCB2 = 0
    IDRNCB = 0
    IUNIT = 0
    NDRNFDEP = 0 
    NDRNVDEP = 0
    NHEDDEP = 0
    NSFRFDEP = 0
    NSFRLDEP = 0
    NSTORDEP = 0
    NSTRFDEP = 0
    NSTRLDEP = 0
    ! Character
    FTYPE = ''
    FNAME = ''
    BASFILE = ''
    BCFFILE = ''
    DISFILE = ''
    DE4FILE = ''
    DRNFILE = ''
    LMGFILE = ''
    GMGFILE = ''
    HUFFILE = ''
    LISTFILE = ''
    LPFFILE = ''
    MNW2FILE = ''
    NWTFILE = ''
    OCFILE = ''
    PCGFILE = ''
    PCGNFILE = ''
    SFRFILE = ''
    SIPFILE = ''
    SORFILE = ''
    STRFILE = ''
    UPWFILE = ''
    WELFILE = ''
    HEADFILE = ''
    FLOWBUDFILE = ''
    SFRLEAKBUDFILE = ''
    SFRFLOWBUDFILE = ''
    DRNBUDFILE = ''
    ! Arrays related to parameters
    ALLOCATE(B(MXPAR))
    ALLOCATE(IACTIVE(MXPAR))
    ALLOCATE(IPLOC(4,MXPAR))
    ALLOCATE(IPCLST(14,MXCLST))
    ALLOCATE(PARNAM(MXPAR))
    ALLOCATE(PARTYP(MXPAR))
    ALLOCATE(INAME(MXINST))
    B = 0.0D0
    IACTIVE = 0
    IPLOC = 0
    IPCLST = 0
    PARNAM = ' '
    PARTYP = ' '
    INAME = ' '
    ! Pointers
    ALLOCATE(IPSUM,ICLSUM,INAMLOC,NPVAL)
    IPSUM = 0
    ICLSUM = 0
    INAMLOC = 1
    NPVAL = 0
    ALLOCATE(NSTRM)
    NSTRM = 0
    ! Derived type
    DO I=1,MAXDEPPKGS
      CALL DEPENDENT_PKG_INIT(DEPPKGS(I))
    ENDDO
    RETURN
  END SUBROUTINE INITIALIZE_MMPMOD
  !-------------------------------------------------------------------
  SUBROUTINE READ_FILENAMES(IUMAIN,OK)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IUMAIN
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=40) :: MESS
    !
    OKLOCAL = .TRUE.
    IF (NSFRFDEP+NSFRLDEP>0) THEN
      MESS = 'Error reading name of SFRFILE'
      READ(IUMAIN,*,END=900)SFRFILE
    ENDIF
    IF (NDRNFDEP+NDRNVDEP>0) THEN
      MESS = 'ERROR READING NAME OF DRN AUX VARIABLE'
      READ(IUMAIN,*,END=900)AUXVAR
    ENDIF
    RETURN
    !
    ! Error handling
    900 CONTINUE
    OKLOCAL = .FALSE.
    CALL APPEND_MESSAGE(MESS)
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_FILENAMES
  !-------------------------------------------------------------------
  SUBROUTINE MMP_ALLOC()
    ! Allocate arrays
    IMPLICIT NONE
    ALLOCATE(BUFF2D(NCOL,NROW))
    ALLOCATE(DBUFF2D(NCOL,NROW))
    ALLOCATE(BUFF3D(NCOL,NROW,NLAY))
    ALLOCATE(DBUFF3D(NCOL,NROW,NLAY))
    ALLOCATE(BUF(NCOL*NROW*NLAY))
    ALLOCATE(DBUF(NCOL*NROW*NLAY))
    ALLOCATE(IBUF(NCOL*NROW*NLAY))
    RETURN
  END SUBROUTINE MMP_ALLOC
  !-------------------------------------------------------------------
  SUBROUTINE WRITE_WEL_FILE(WT,OK)
    ! Using parameter data in module variables and non-parameter data
    ! in WELTYPE structure WT, write a WEL input file.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    USE GWM_UTLS, ONLY: OPEN_FILE
    USE MMPWELMOD, ONLY: ACTIVEPARS, IWELCB, MXPW, NAUX, NPAR, NPWEL, &
                         NWELVL, WELLP, WELOPTIONLIST
    USE PARAMMODULE, ONLY: B, INAME, IPLOC, MXPAR, PARNAM, PARTYP
    IMPLICIT NONE
    ! Arguments
    TYPE (WELTYPE), INTENT(IN) :: WT
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: COL, I, ILOC, IPL4, IST, ISTAT, IU, J, K, KP, KW, &
               L, LAY, NLST, NNPWEL, MAXTRIES, MXACTW, NP, &
               NREAD1, NREAD2, NUMINST, ROW
    REAL :: PVAL
    DOUBLE PRECISION :: WAITSECS, QFACT
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=10) :: INSTNAM, PNAM
    ! Format statements
    10 FORMAT('# Well package input file written by MMProc')
    20 FORMAT('PARAMETER',2X,I5,2X,I8)
    30 FORMAT(2I10,2X,A)
    40 FORMAT(A10,2X,A4,2X,G24.16,2X,I5,2X,'INSTANCES',2X,I5)
    50 FORMAT(A10,2X,A4,2X,G24.16,2X,I5)
    60 FORMAT(A10)
    70 FORMAT(3I10,1X,G24.16,20(:,1X,G24.16))
    80 FORMAT(2I10,2X,'Stress period: ',I5)
    !
    OKLOCAL = .TRUE.
    IF (WELFILE == '') THEN
      CALL APPEND_MESSAGE('Error: WELFILE is blank')
      OKLOCAL = .FALSE.
      GOTO 9000
    ENDIF
    MAXTRIES = 1000
    WAITSECS = 0.1
    ! Find the unit number connected to file WELFILE.
    INQUIRE(FILE=WELFILE,NUMBER=IU)
    IF (IU < 1) THEN
      ! No unit is connected to the file.  Find an unused unit number 
      ! that can be used by PLL_CLOSE_AND_DELETE to close file WELFILE.
      IU = UTL_GETUNIT(7,10000)
    ENDIF
    ISTAT = PLL_CLOSE_AND_DELETE(IU,WELFILE,MAXTRIES,WAITSECS)
    IF (ISTAT==0) THEN
      IU = UTL_GETUNIT(7,9000)
      OKLOCAL = OPEN_FILE(IU,WELFILE,'REPLACE')
      IF (.NOT. OKLOCAL) GOTO 9000
      MXACTW = WT%MXNPW
      WRITE(IU,10)                   ! Comment line
      IF (NPWEL>0) THEN
        WRITE(IU,20)NPWEL,MXPW       ! Item 1: "PARAMETER" line
        MXACTW = MXACTW + MXPW
      ENDIF
      WRITE(IU,30)MXACTW,IWELCB,TRIM(WELOPTIONLIST)    ! Item 2
      !
      IF (NPWEL>0) THEN
        ! Write items 3 and 4, which define parameters
        I = 0     ! Counter for WEL parameters (type = "Q")
        IST = 1   ! Search start point in list of all parameters
        KW = 0    ! Counter for parameter wells
        DO WHILE (I<NPWEL)
          ! Find next "Q" parameter in list of parameters.
          FINDPAR: DO NP=IST,MXPAR
            IF (PARTYP(NP)=='Q') THEN
              PNAM = PARNAM(NP)
              PVAL = B(NP)
              NLST=IPLOC(2,NP)-IPLOC(1,NP)+1
              NUMINST=IPLOC(3,NP)
              ILOC=IPLOC(4,NP)
              I = I + 1
              IST = NP + 1
              EXIT FINDPAR
            ENDIF
          ENDDO FINDPAR
          IF (NUMINST>0) THEN
            WRITE(IU,40)PNAM,'Q   ',PVAL,NLST,NUMINST  ! Item 3
            DO J=1,NUMINST
             IPL4 = IPLOC(4,NP)
             ILOC = IPL4+J-1
             INSTNAM = INAME(ILOC)
              WRITE(IU,60)INSTNAM                      ! Item 4a
              DO K=1,NLST
                KW = KW + 1
                LAY = WELLP(1,KW)
                ROW = WELLP(2,KW)
                COL = WELLP(3,KW)
                ! QFACT is WELLP(4,KW).  Aux variables, if any, follow.
                WRITE(IU,70)LAY,ROW,COL,(WELLP(L,KW),L=4,NWELVL-1) ! Item 4b
              ENDDO
            ENDDO
          ELSE
            WRITE(IU,50)PNAM,'Q   ',PVAL,NLST          ! Item 3
            DO K=1,NLST
              KW = KW + 1
              LAY = WELLP(1,KW)
              ROW = WELLP(2,KW)
              COL = WELLP(3,KW)
              ! QFACT is WELLP(4,KW).  Aux variables, if any, follow.
              WRITE(IU,70)LAY,ROW,COL,(WELLP(L,KW),L=4,NWELVL-1)  ! Item 4b
            ENDDO
          ENDIF
        ENDDO
      ENDIF 
      ! Done writing parameter definitions
      !
      NREAD2=NWELVL-1
      NREAD1=NREAD2-NAUX
      ! Write input for each stress period
      DO KP=1,NPER
        NNPWEL = WT%NNPWEL(KP)
        NP = NPAR(KP)
        WRITE(IU,80)NNPWEL,NP,KP                                  ! Item 5
        DO I=1,WT%MXNPW
          LAY = WT%WELLNP(1,I,KP)
          IF (LAY>0) THEN
            ROW = WT%WELLNP(2,I,KP)
            COL = WT%WELLNP(3,I,KP)
            ! Q is WELLNP(4,I,KP).  Aux variables, if any, follow.
            WRITE(IU,70)LAY,ROW,COL,(WT%WELLNP(L,I,KP),L=4,NREAD1) ! Item 6
          ENDIF
        ENDDO
        DO I=1,NP
          WRITE(IU,'(A)')TRIM(ACTIVEPARS(I,KP))                    ! Item 7
        ENDDO
      ENDDO
      CLOSE(IU)
    ELSE
      CALL APPEND_MESSAGE('Error: Unable to delete file "' // &
                          TRIM(WELFILE) // '"')
      OKLOCAL = .FALSE.
    ENDIF
    !
    9000 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE WRITE_WEL_FILE
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_HEADS(NCOL,NROW,DHEADS,MANAGED_WELLS,NMQMNW, &
                           MNW2WELLS,OK,HEADFILEEMPTY)
    ! Extract head data for dependents and managed-flow cells from
    ! binary head file.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_SUBS, ONLY: FUZZY_EQUALS
    USE GWM_UTLS, ONLY: OPEN_FILE
    USE PRECUTLSMOD, ONLY: HEADPRECISION, GET_BINARY_HEAD_DATASET, &
                           GET_BINARY_FILE_LENGTH
    USE MMPMNWMOD, ONLY: MNW2WELL
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NCOL, NROW, NMQMNW
    DOUBLE PRECISION, DIMENSION(NCOL,NROW), INTENT(INOUT) :: DHEADS
    TYPE (WELTYPE), INTENT(INOUT) :: MANAGED_WELLS
    TYPE (MNW2WELL), DIMENSION(NMQMNW), INTENT(INOUT) :: MNW2WELLS
    LOGICAL, INTENT(INOUT) :: OK
    LOGICAL, INTENT(OUT)   :: HEADFILEEMPTY
    ! Local variables
    INTEGER :: KH, KPERLAST, KSP, KSTPLAST, IU, M, PRECISION
    LOGICAL :: OKLOCAL
    ! Local variables
    INTEGER :: FILELEN, I, ILAY, J, K, KACTIVE, KPER, KSTP, L, &
               N, NCELLS
    CHARACTER(LEN=16) :: TEXT
    !
    OKLOCAL = .TRUE.
    HEADFILEEMPTY = .FALSE.
    KPER = 0
    KSTP = 0
    IU = UTL_GETUNIT(7,9000)
    OKLOCAL = OPEN_FILE(IU,HEADFILE,'OLD',.TRUE.)
    ! If OKLOCAL is true, HEADFILE exists
    ! Determine precision of binary head file
    PRECISION = HEADPRECISION(IU)
    IF (PRECISION<0) THEN
      CALL APPEND_MESSAGE('Error: Unable to determine precision ' // &
                          'of binary head file "' // TRIM(HEADFILE) // '"')
      OKLOCAL = .FALSE.
      ! If HEADFILE exists and precision cannot be determined, HEADFILE may
      ! be empty.  If so, it's an indication that Modflow likely failed.
      FILELEN = GET_BINARY_FILE_LENGTH(IU)
      IF (FILELEN==0) THEN
        OKLOCAL = .TRUE.
        HEADFILEEMPTY = .TRUE.
      ENDIF
    ENDIF
    IF (OKLOCAL .AND. .NOT. HEADFILEEMPTY) THEN
      KH = 0
      ! Iterate through data sets in binary head file
      ! Ned TODO: Also need to count the WEL and MNW managed wells 
      ! to ensure loop does not terminate prematurely
      ! Instead, probably just iterate through entire binary file
      DO WHILE (KH<NHEDDEP) 
        KPERLAST = KPER
        KSTPLAST = KSTP
        CALL GET_BINARY_HEAD_DATASET(IU,PRECISION,KSTP,KPER,TEXT,NCOL, &
                                 NROW,ILAY,DHEADS,OKLOCAL)
        IF (.NOT. OKLOCAL) GOTO 100
        IF (KPER.NE.KPERLAST .OR. KSTP.NE.KSTPLAST) THEN
          IF (NMQMNW > 0) THEN
            DO N=1,NMQMNW
              MNW2WELLS(N)%KDRY = 0
            ENDDO
          ENDIF
        ENDIF
        ! Consider only head values, not drawdown
        IF (TEXT=='HEAD') THEN
          ! Find and save head values for any head dependents that
          ! match stress period, time step, and layer of data set.
          DO I=1,NSV
            ! Consider only head dependents
            IF (DEPS(I)%CTYPE=='HEAD') THEN
              ! Stress period has to match
              IF (KPER==DEPS(I)%INDX(4)) THEN
                ! Consider only last time step in stress period, and
                ! see if layer matches
                IF (KSTP==NSTP(KPER) .AND. ILAY==DEPS(I)%INDX(3)) THEN
                  ! Assign simulated head value for this head DEPENDENT
                  DEPS(I)%SIMVAL = DHEADS(DEPS(I)%INDX(1),DEPS(I)%INDX(2))
                  KH = KH + 1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          ! Check head values for managed wells for HDRY
          ! Iterate through all managed wells (WEL Package)
          IF (MANAGED_WELLS%MXNPW > 0) THEN
            WELLLOOP: DO N=1,MANAGED_WELLS%MXNPW
              ! If layer of head array matches managed-well layer
              L = MANAGED_WELLS%WELLNP(1,N,KPER)
              IF (ILAY==L) THEN
                ! If stress period is one in which managed well is active
                IF (MANAGED_WELLS%WELLNP(4,N,KPER) .NE. 0.0) THEN
                  I = MANAGED_WELLS%WELLNP(2,N,KPER)
                  J = MANAGED_WELLS%WELLNP(3,N,KPER)
                  IF (FUZZY_EQUALS(DHEADS(J,I),HDRY,1.0D-7)) THEN 
                    MANAGED_WELLS%DRY(N) = .TRUE.
                    CYCLE WELLLOOP
                  ENDIF
                ENDIF
              ENDIF
            ENDDO WELLLOOP
          ENDIF
          ! Iterate through all managed wells (MNW2 Package)
          IF (NMQMNW > 0) THEN
            MNWLOOP: DO N=1,NMQMNW
              ! If stress period KPER is one in which managed well is active
              IF (MNW2WELLS(N)%IACTIVE(KPER) > 0) THEN
                NCELLS = MNW2WELLS(N)%NC
                DO M=1,NCELLS
                  ! If layer of head array matches layer of current cell
                  ! included in MNW well
                  L = MNW2WELLS(N)%INDX(M,1)   ! Layer
                  IF (ILAY==L) THEN
                    I = MNW2WELLS(N)%INDX(M,2) ! Row
                    J = MNW2WELLS(N)%INDX(M,3) ! Column
                    IF (FUZZY_EQUALS(DHEADS(J,I),HDRY,1.0D-7)) THEN
                      MNW2WELLS(N)%KDRY = MNW2WELLS(N)%KDRY + 1
                    ENDIF
                  ENDIF
                ENDDO
                IF (MNW2WELLS(N)%KDRY==NCELLS) THEN
                  MNW2WELLS(N)%DRY = .TRUE.
                  CYCLE MNWLOOP
                ENDIF
              ENDIF
            ENDDO MNWLOOP
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    !
    100 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_HEADS
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_BUD_GRID_DEPS(OK)
    ! Extract STORAGE data sets from cell-by-cell budget file, and 
    ! assign appropriate simulated values for dependents
    ! This subroutine could be extended to add support for any
    ! DEPENDENT type that uses data based on the 3D model grid
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    USE PRECUTLSMOD, ONLY: BUDGETPRECISION
    USE PRECUTLSMOD, ONLY: GET_BUDGET_3D_ARRAY
    IMPLICIT NONE
    ! Arguments
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, IPREC, IU, J, K, KPER, KSTOR, KSTP, N, NCB
    DOUBLE PRECISION :: DVAL
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=16) :: QUANTITY
    !
    OKLOCAL = .TRUE.
    IPREC = 0
    IF (FLOWBUDFILE=='') OKLOCAL = .FALSE.
    IF (OKLOCAL) THEN
      IU = UTL_GETUNIT(7,9000)
      OKLOCAL = OPEN_FILE(IU,FLOWBUDFILE,'OLD',.TRUE.)
      IF (OKLOCAL) THEN
        IPREC = BUDGETPRECISION(IU)
        ! Determine precision of binary budget file
        IF (IPREC<1) THEN
          CALL APPEND_MESSAGE('Error: Unable to determine precision ' // &
                              'of budget file "' // TRIM(FLOWBUDFILE) // '"')
          OKLOCAL = .FALSE.
        ENDIF
      ENDIF
    ENDIF
    !
    IF (OKLOCAL) THEN
      KSTOR = 0
      ! Iterate through data sets in binary budget file
      DO WHILE (KSTOR<NSTORDEP) 
        CALL GET_BUDGET_3D_ARRAY(IU,IPREC,NCOL,NROW,NLAY,QUANTITY,NCB,KPER, &
                                 KSTP,DBUFF3D,OKLOCAL)
        IF (.NOT. OKLOCAL) GOTO 100
        ! For storage dependents, only interested in data sets with QUANTITY = "STORAGE"
        IF (QUANTITY=='STORAGE') THEN 
          DO N=1,NSV
            ! Consider only storage dependents
            IF (DEPS(N)%CTYPE=='STORAGE-CHANGE' .AND. DEPS(N)%INDX(3)==-999) THEN
              ! See if stress period is in range Begin-stress-period to End-stress-period
              IF (KPER>=DEPS(N)%INDX(1) .AND. KPER<=DEPS(N)%INDX(2)) THEN
                ! Accumulate simulated values in this storage DEPENDENT
                DO K=1,NLAY
                  DO I=1,NROW
                    DO J=1,NCOL
                      DVAL = DBUFF3D(J,I,K)
                      IF (DEPS(N)%CZONE=='ZONE') THEN
                        IF (DEPS(N)%ZONE(J,I,K).LE.0) DVAL = 0.0D0
                      ENDIF
                      DEPS(N)%SIMVAL = DEPS(N)%SIMVAL + DVAL*DELT(KSTP,KPER)
                    ENDDO
                  ENDDO
                ENDDO
                IF (KPER==DEPS(N)%INDX(2) .AND. KSTP==NSTP(KPER)) THEN
                  ! Calculation of DEPS(N)%SIMVAL is complete
                  DEPS(N)%INDX(3) = 999  
                  KSTOR = KSTOR + 1
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    !
    ! Normal return
    IF (OKLOCAL) RETURN
    !
    ! Error handling
    100 CONTINUE
    CALL APPEND_MESSAGE('Error encountered in extracting storage data from file "' &
                        // TRIM(FLOWBUDFILE) // '"')
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_BUD_GRID_DEPS
  !-------------------------------------------------------------------
  SUBROUTINE SFR_IDENTIFY(OK)
    ! For SFR-FLOW and SFR-LEAK dependents, find corresponding indices
    ! in lists to be read from binary budget files
    USE GWFSFRMODULE, ONLY: NSTRM
    IMPLICIT NONE
    ! Arguments
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, IRCH, ISEG, IUSFR, J
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    !
    IF (OKLOCAL) THEN
      IF (SFRFILE .NE. '' .AND. (NSFRFDEP>0 .OR. NSFRLDEP>0)) THEN
        CALL MMP_READ_SFR_FILE(OKLOCAL)
        IF (OKLOCAL) THEN
          DO I=1,NSTRM
            ISEG = ISTRMSFR(4,I)
            IRCH = ISTRMSFR(5,I)
            DO J=1,NSV
              IF (DEPS(J)%CTYPE=='SFR-FLOW' .OR. DEPS(J)%CTYPE=='SFR-LEAK') THEN
                IF (DEPS(J)%INDX(1)==ISEG .AND. DEPS(J)%INDX(2)==IRCH) THEN
                  DEPS(J)%INDX(4) = I
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ELSE
          CALL APPEND_MESSAGE('Error reading SFR input file "' // TRIM(SFRFILE) // '"')
          OKLOCAL = .FALSE.
        ENDIF
      ENDIF
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE SFR_IDENTIFY
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_DRN_DEPS(OK)
    ! Extract DRN dependents from budget file
    USE PRECUTLSMOD, ONLY: EXTRACT_BUD_DRN_DEPS
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    CHARACTER(LEN=2000) :: DRNBFILE
    CHARACTER(LEN=200) :: MESS
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    CALL GET_FILENAME(ABS(IDRNCB),DRNBFILE)
    CALL EXTRACT_BUD_DRN_DEPS(NCOL,NROW,NLAY,NPER,NSV,MXACTD,MAXTS, &
                              NSTP,DRNBFILE,AUXVAR,DEPS,DELT,OKLOCAL)
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_DRN_DEPS
  !-------------------------------------------------------------------
  SUBROUTINE CHECK_MNW_QNET(NMQMNW,MNW2WELLS,OK)
    ! Extract Qnet values from binary output file containing MNW2 
    ! data, then, for each managed MNW2 well, compare Qnet value with 
    ! Qdes value provided in mmproc.in.  If values are not 
    ! FUZZY_EQUALS, flag managed well as dewatered.
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE PRECUTLSMOD, ONLY: EXTRACT_QNET_MNW2
    USE MMPMNWMOD, ONLY: MNW2WELL
    IMPLICIT NONE
    ! Argument-list variables
    INTEGER, INTENT(IN) :: NMQMNW
    TYPE (MNW2WELL), DIMENSION(NMQMNW), INTENT(INOUT) :: MNW2WELLS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, MAXDIM
    LOGICAL :: LEX, OKLOCAL
    CHARACTER(LEN=2000) :: MNWBUDFILE
    !
    OKLOCAL = .TRUE.
    IF(MNWMAX.GT.0)THEN
      ! ERB 2/20/2014 changed MAXDIM formula for consistency with MNW2
      MAXDIM = (MNWMAX + NMQMNW)*(3*NLAY + NLAY*NLAY) + (10*NLAY) + 25
      ! ERB 6/4/2014 increased MAXDIM again to avoid too-small value
      MAXDIM = MAXDIM * 10
    ELSE  ! Use the value of MAXDIM provided in the input file
      MAXDIM = NODTOT
    ENDIF
    CALL GET_FILENAME(ABS(IMNWCB),MNWBUDFILE)
    INQUIRE(FILE=MNWBUDFILE,EXIST=LEX)
    IF (LEX) THEN
      CALL EXTRACT_QNET_MNW2(NCOL,NROW,NLAY,NPER,NMQMNW,MAXDIM, &
                             MAXTS,NSTP,MNWBUDFILE,AUXVAR, &
                             MNW2WELLS,DELT,OKLOCAL)
    ELSE
      OKLOCAL = .FALSE.
      ! Ned TODO: Create error message here
      IF (MNWBUDFILE=='' .OR. MNWBUDFILE==' ') THEN
        AMESSAGE = 'ERROR: Binary file containing MNW2 budget data not found'
      ELSE
        AMESSAGE = 'ERROR: Binary file "'//TRIM(MNWBUDFILE)//'" does not exist'
      ENDIF
    ENDIF
    IF (OK) THEN
      OK = OKLOCAL
    ENDIF
    RETURN
  END SUBROUTINE CHECK_MNW_QNET
  !-------------------------------------------------------------------
  SUBROUTINE EXTRACT_SFR_DEPS(OK)
    ! Extract all SFR dependents from either or both budget files
    ! SFRFLOWBUDFILE and SFRLEAKBUDFILE
    USE GWFSFRMODULE, ONLY: NSTRM
    USE PRECUTLSMOD, ONLY: EXTRACT_BUD_SFR_DEPS
    IMPLICIT NONE
    ! Arguments
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: ITYP
    CHARACTER(LEN=200) :: MESS
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    IF (SFRFLOWBUDFILE.NE.'' .AND. SFRLEAKBUDFILE.NE.'') THEN
      IF (SFRFLOWBUDFILE==SFRLEAKBUDFILE) THEN
        MESS = 'Error: Binary budget files for SFR-FLOW' // &
        ' and SFR-LEAK have the same name'
        CALL APPEND_MESSAGE(MESS)
        OKLOCAL = .FALSE.
      ENDIF
    ENDIF
    IF (SFRFLOWBUDFILE.NE.'') THEN
      CALL EXTRACT_BUD_SFR_DEPS(NCOL,NROW,NLAY,NPER,NSTRM,NSV, &
                                SFRFLOWBUDFILE,1,DEPS,NSTP,DBUFF3D, &
                                OKLOCAL)
    ENDIF
    IF (SFRLEAKBUDFILE.NE.'') THEN
      CALL EXTRACT_BUD_SFR_DEPS(NCOL,NROW,NLAY,NPER,NSTRM,NSV, &
                                SFRLEAKBUDFILE,2,DEPS,NSTP,DBUFF3D, &
                                OKLOCAL)
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE EXTRACT_SFR_DEPS
  !-------------------------------------------------------------------
  SUBROUTINE WRITE_DEPENDENT_DATA(OK)
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, IU
    LOGICAL :: OKLOCAL
    ! Format statments
    10 FORMAT(A,2X,G23.16)
    !
    OKLOCAL = .TRUE.
    IU = UTL_GETUNIT(7,9000)
    OKLOCAL = OPEN_FILE(IU,SIMVALOUT,'REPLACE')
    IF (OKLOCAL) THEN
      DO I=1,NSV
        WRITE(IU,10)DEPS(I)%CID,DEPS(I)%SIMVAL
      ENDDO
      CLOSE(IU)
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE WRITE_DEPENDENT_DATA
  !-------------------------------------------------------------------
  SUBROUTINE WRITE_STATUS(MGWELLS,NMQMNW,MNW2WELLS,HEADFILEEMPTY,OK, &
                          NPER,NSTP)
    USE MGT, ONLY: MGT_FIND_STATUS, MGT_WRITESTATUS
    USE MMPMNWMOD, ONLY: MNW2WELL
    IMPLICIT NONE
    ! Arguments
    TYPE (WELTYPE), INTENT(IN) :: MGWELLS
    INTEGER, INTENT(IN) :: NMQMNW
    TYPE (MNW2WELL), DIMENSION(NMQMNW), INTENT(IN) :: MNW2WELLS
    LOGICAL, INTENT(IN) :: HEADFILEEMPTY
    LOGICAL, INTENT(INOUT) :: OK
    INTEGER, INTENT(IN) :: NPER
    INTEGER, DIMENSION(NPER), INTENT(IN) :: NSTP
    ! Local variables
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    CALL MGT_FIND_STATUS(NCOL,NROW,LISTFILE,MGWELLS,NMQMNW, &
                         MNW2WELLS,HEADFILEEMPTY,OKLOCAL,NPER,NSTP)
    CALL MGT_WRITESTATUS(OKLOCAL)
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE WRITE_STATUS
  !-------------------------------------------------------------------
  ! UTILITIES
  !-------------------------------------------------------------------
  LOGICAL FUNCTION IS_NEW_DEPPKG(DEP) 
    IMPLICIT NONE
    ! Arguments
    TYPE (DEPENDENT), INTENT(IN) :: DEP
    ! Local variables
    INTEGER :: J
    !
    IS_NEW_DEPPKG = .TRUE.
    IF (KCP > 0) THEN
      SEARCH: DO J=1,KCP
        IF (DEPPKGS(J)%PKGID == DEP%PKGID) THEN
          IS_NEW_DEPPKG = .FALSE.
          EXIT SEARCH
        ENDIF
      ENDDO SEARCH
    ENDIF
    RETURN 
  END FUNCTION IS_NEW_DEPPKG
  !-------------------------------------------------------------------
  INTEGER FUNCTION GET_BINARY_UNIT(PKGID)
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: PKGID
    ! Local variables
    INTEGER :: IU
    !
    SELECT CASE (PKGID)
    CASE ('HEAD')
      IU = IHEDUN
    CASE ('SFRL')
      IU = ISTCB2
    CASE ('SFRF')
      IU = ISTCB1
    CASE ('DRNF')
      IU = IDRNCB
    CASE ('DRNV')
      IU = IDRNCB
    CASE ('STOR')
      IU = IFLOWCB
    CASE DEFAULT
      IU = 0
    END SELECT
    !
    GET_BINARY_UNIT = ABS(IU)
    RETURN
  END FUNCTION GET_BINARY_UNIT
  !-------------------------------------------------------------------
  SUBROUTINE GET_FILENAME(UNIT,FILENAME)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: UNIT
    CHARACTER(LEN=*), INTENT(OUT) :: FILENAME
    ! Local variables
    INTEGER :: I
    !
    FILENAME = ''
    !
    IF (UNIT>0) THEN
      SEARCH: DO I=1,MAXUNITS
        IF (IUNIT(I) == UNIT) THEN
          FILENAME = FNAME(I)
          EXIT SEARCH
        ENDIF
      ENDDO SEARCH
    ENDIF
    !
    RETURN
  END SUBROUTINE GET_FILENAME
  !-------------------------------------------------------------------
  LOGICAL FUNCTION SAVE_WEL_FILE()
    USE GWM_UTLS, ONLY: COPY_FILE
    IMPLICIT NONE
    !
    WELFILECOPY = TRIM(WELFILE)//'_COPY.TXT'
    SAVE_WEL_FILE = COPY_FILE(WELFILE,WELFILECOPY)
    RETURN
  END FUNCTION SAVE_WEL_FILE
  !-------------------------------------------------------------------
  SUBROUTINE RESTORE_WEL_FILE(OK)
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE PARALLEL_PROCESSING, ONLY: PLL_CLOSE_AND_DELETE
    USE GWM_UTLS, ONLY: COPY_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    LOGICAL :: OKLOCAL
    INTEGER :: ISTAT, IU, MAXTRIES
    DOUBLE PRECISION :: WAITSECS
    !
    OKLOCAL = .TRUE.
    IU = UTL_GETUNIT(7,9000)
    MAXTRIES = 500
    WAITSECS = 0.1
    ISTAT = PLL_CLOSE_AND_DELETE(IU,WELFILE,MAXTRIES,WAITSECS)
    IF (ISTAT==0) THEN
      OKLOCAL = COPY_FILE(WELFILECOPY,WELFILE)
    ELSE
      CALL APPEND_MESSAGE( &
        'Error: Unable to restore original WEL package input file')
      OKLOCAL = .FALSE.
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE RESTORE_WEL_FILE
  !-------------------------------------------------------------------
  SUBROUTINE MMP_READ_NAME_FILE(OK)
    ! Open Modflow name file and read required information.
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_CASETRANS, UTL_RWORD, &
                         UTL_SAMENAME
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: INBAS, INDIS, IST, ISTART, ISTOP, IUNAME, IUOC, K
    DOUBLE PRECISION :: R
    LOGICAL :: LEX, OKLOCAL
    CHARACTER(LEN=200) :: LINE
    ! Format statements
    600 FORMAT('Error opening file: ',A)
    910 FORMAT('Error reading data from line:',/,A,/)
    !
    OKLOCAL = .TRUE.
    ! Initialize variables
    IUOC = 0
    !
    ! Open and read name file to get name of DIS file
    IUNAME = UTL_GETUNIT(7,1000)
    OPEN(IUNAME,FILE=NAMEFILE,STATUS='OLD',IOSTAT=IST)
    IF (IST .NE. 0) THEN
      WRITE(*,600)TRIM(NAMEFILE)
      GOTO 9100
    ENDIF
    !
    ! Read names of files from name file
    K = 0
    FIND: DO
      READ(IUNAME,'(A200)',END=100)LINE
      IF (LINE(1:1)=='#') CYCLE FIND
      K = K + 1
      READ(LINE,*,IOSTAT=IST,END=100,ERR=9000)FTYPE(K),IUNIT(K),FNAME(K)
      CALL UTL_CASETRANS(FTYPE(K),'hi')
      IF (FTYPE(K)=='DIS') THEN
        DISFILE=FNAME(K)
      ELSEIF (FTYPE(K)=='DRN') THEN
        DRNFILE=FNAME(K)
      ELSEIF (FTYPE(K)=='MNW2') THEN
        MNW2FILE=FNAME(K)
      ENDIF
    ENDDO FIND
    100 CONTINUE
    CLOSE(IUNAME)  ! Close name file
    !
    ! Open and read DIS file to populate DELT array
    IF (OKLOCAL) CALL READ_DIS_FILE(OKLOCAL)
    IF (.NOT. OKLOCAL) GOTO 9100
    !
    ! Normal return
    IF (OK) OK = OKLOCAL
    RETURN
    !
    ! Error handling
    9000 CONTINUE ! Error reading data from LINE
    WRITE(IOUT,910)TRIM(LINE)
    CALL APPEND_MESSAGE('Error reading name file')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
    9100 CONTINUE ! Error, nature of error already described
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
  END SUBROUTINE MMP_READ_NAME_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_DIS_FILE(OK)
    ! Open DIS file and read NLAY, NROW, NCOL, NPER, and NSTP for each
    ! stress period.  Calculate and store DELT for all time steps.
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_RWORD
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, INDIS, IST, K, KK, NCNFBD, ISS, ITR, N, LLOC, ISTART, &
               ISTOP
    REAL :: R, ZERO
    DOUBLE PRECISION :: ONE
    CHARACTER(LEN=2000) :: LINE
    CHARACTER*24 ANAME(5)
    LOGICAL :: OKLOCAL
    DATA ANAME(1) /'                    DELR'/
    DATA ANAME(2) /'                    DELC'/
    DATA ANAME(3) /'TOP ELEVATION OF LAYER 1'/
    DATA ANAME(4) /'  MODEL LAYER BOTTOM EL.'/
    DATA ANAME(5) /'BOT. EL. OF QUASI-3D BED'/
    ! Format statements
    600 FORMAT('Error opening file: ',A)
    910 FORMAT('Error reading data from line:',/,A,/)
    !
    ! Open DIS file
    INDIS=UTL_GETUNIT(20,1000)
    OKLOCAL = OPEN_FILE(INDIS,DISFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      ! Read DIS file to get NLAY, NROW, NCOL, NPER
      CALL URDCOM(INDIS,IOUT,LINE)
      READ(LINE,*,ERR=9000)NLAY,NROW,NCOL,NPER
      !
      !C7------ALLOCATE LAYER FLAGS.
      ALLOCATE(LBOTM(NLAY))
      ALLOCATE(LAYCBD(NLAY))
      !C
      !C8------Read confining bed information
      READ(INDIS,*) (LAYCBD(K),K=1,NLAY)
      LAYCBD(NLAY)=0
      !C
      !C9------Count confining beds, setup the pointer to each layer's
      !C9------bottom array (LBOTM), and setup LAYCBD to be the confining
      !C9------bed number for each layer.
      NCNFBD=0
      DO 100 K=1,NLAY
        LBOTM(K)=K+NCNFBD
        IF(LAYCBD(K).NE.0) THEN
          NCNFBD=NCNFBD+1
          LAYCBD(K)=NCNFBD
        END IF
      100 CONTINUE
      NBOTM=NLAY+NCNFBD
      !C
      !C10-----Allocate space for discretization arrays
      !C10-----Note that NBOTM+1 arrays are allocated for BOTM
      !C10-----because BOTM(J,I,0) contains the top elevation of layer 1.
      ALLOCATE (DELR(NCOL))
      ALLOCATE (DELC(NROW))
      ALLOCATE (BOTM(NCOL,NROW,0:NBOTM))
      ALLOCATE (PERLEN(NPER),TSMULT(NPER),ISSFLG(NPER))
      !C
      !C11-----Read the DELR and DELC arrays.
      CALL U1DREL(DELR,ANAME(1),NCOL,INDIS,IOUT)
      CALL U1DREL(DELC,ANAME(2),NROW,INDIS,IOUT)
      !C
      !C12-----Read the top elevation of layer 1.
      CALL U2DREL(BOTM(:,:,0),ANAME(3),NROW,NCOL,0,INDIS,IOUT)
      !C
      !C13-----Read the bottom elevations.
      DO 120 K=1,NLAY
        KK=K
        CALL U2DREL(BOTM(:,:,LBOTM(K)),ANAME(4),NROW,NCOL,KK,INDIS,IOUT)
        IF(LAYCBD(K).NE.0) CALL U2DREL(BOTM(:,:,LBOTM(K)+1),ANAME(5), &
                 NROW,NCOL,KK,INDIS,IOUT)
      120 CONTINUE
      !C
      !C14-----READ AND WRITE LENGTH OF STRESS PERIOD, NUMBER OF TIME STEPS,
      !C14-----TIME STEP MULTIPLIER, AND STEADY-STATE FLAG..
      ISS=0
      ITR=0
      MAXTS = 0
      DO 200 N=1,NPER
        READ(INDIS,'(A)') LINE
        LLOC=1
        ! Use UTL_RWORD instead URWORD because it supports double precision.
        CALL UTL_RWORD(IOUT,3,.FALSE.,LLOC,LINE,ISTART,ISTOP,I,PERLEN(N),INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSTP(N),R,IOUT,INDIS)
        CALL UTL_RWORD(IOUT,3,.FALSE.,LLOC,LINE,ISTART,ISTOP,I,TSMULT(N),INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,INDIS)
        IF (LINE(ISTART:ISTOP).EQ.'TR') THEN
          ISSFLG(N)=0
          ITR=1
        ELSE IF (LINE(ISTART:ISTOP).EQ.'SS') THEN
          ISSFLG(N)=1
          ISS=1
        END IF
        IF (NSTP(N)>MAXTS) MAXTS = NSTP(N)
        ZERO=0.
      200 CONTINUE
      !C
      !C16-----Assign ITRSS.
      IF(ISS.EQ.0 .AND. ITR.NE.0) THEN
        ITRSS=1
      ELSE IF(ISS.NE.0 .AND. ITR.EQ.0) THEN
        ITRSS=0
      ELSE
        ITRSS=-1
      END IF
      !
      ! Allocate and populate DELT array
      ALLOCATE(DELT(MAXTS,NPER))
      ONE=1.0D0
      DO N=1,NPER
        ! Calculate the length of the first time step.
        ! Assume time step multiplier is equal to one.
        DELT(1,N)=PERLEN(N)/REAL(NSTP(N),8)
        ! If time step multiplier is not one, then calculate first
        ! term of geometric progression.
        IF(TSMULT(N).NE.ONE) THEN
          DELT(1,N) = PERLEN(N)*(ONE-TSMULT(N))/ &
                 (ONE-TSMULT(N)**NSTP(N))
        ENDIF
        ! Populate remaining DELT values for stress period N
        IF (NSTP(N)>1) THEN
          DO I=2,NSTP(N)
            DELT(I,N) = DELT(I-1,N)*TSMULT(N)
          ENDDO
        ENDIF
      ENDDO
      !
      CLOSE(INDIS)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
    !
    ! Error handling
    9000 CONTINUE ! Error reading data from LINE
    WRITE(IOUT,910)TRIM(LINE)
    CALL APPEND_MESSAGE('Error reading name file')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
    9100 CONTINUE ! Error, nature of error already described
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
  END SUBROUTINE READ_DIS_FILE
  !-------------------------------------------------------------------
END MODULE MMPROCMOD
