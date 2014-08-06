MODULE MKMMPROC_RDNAM
  USE GWM1BAS3, ONLY: GWMOUT,SMALLEPS
  USE GWM_SUBS, ONLY: MNW2FILE
  USE GLOBAL,   ONLY: NCOL,NROW,NLAY,NPER,NSTP,PERLEN,TSMULT,BOTM, &
                      IFREFM
  USE GWFBASMODULE,ONLY: HNOFLO,HDRY
  ! Arrays and scalars used to read past unneeded data
  USE GLOBAL,   ONLY:  IBOUND, DELR, DELC, LBOTM, LAYCBD, IXSEC, NBOTM
  PRIVATE
  ! PUBLIC VARIABLES
  PUBLIC :: HCLOSE, &
            LISTFILE,DRNFILE,GWMWFILENAME,  &
            IHEDUNFILE,IFLOWCBFILE,SFRFILE,WELFILE,MNW2FILE, &
            SFRB1FILE,SFRB2FILE,DRNBFILE         
  ! PUBLIC PROCEDURES
  PUBLIC :: READ_NAME_FILE,READ_NAME_FILE_DIS   
  ! Scalars used to store data that will be written to MMProc.in.jtf
  REAL    :: HCLOSE
  CHARACTER(LEN=2000) :: LISTFILE,SFRFILE,IHEDUNFILE,GWMWFILENAME,&
                         SFRB1FILE,SFRB2FILE,DRNBFILE, &
                         WELFILE,IFLOWCBFILE,DRNFILE
  ! Local Scalars
  INTEGER :: IBCFCB, IHEDUN, IHUFCB, ILPFCB, IUPWCB, ISTCB1, &
             ISTCB2, IFLOWCB, IDRNCB, ISTORUN, ISCRCH,ICHFLG,IPRTIM
  INTEGER :: BASUNIT,DISUNIT,BCFUNIT,LPFUNIT,HUFUNIT,UPWUNIT
  CHARACTER(LEN=2000) :: BASFILE, BCFFILE, DISFILE, DE4FILE, &
                         LMGFILE, GMGFILE, HUFFILE, UPWFILE,   &
                         LPFFILE, NWTFILE, OCFILE, PCGFILE, &
                         PCGNFILE,  SIPFILE, SORFILE  
  ! Local Arrays
  INTEGER, PARAMETER :: MAXUNITS = 100
  CHARACTER(LEN=200), DIMENSION(MAXUNITS):: FNAME, FTYPE
  INTEGER, DIMENSION(MAXUNITS)           :: IMUNIT
 
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE READ_NAME_FILE(NAMFILE)
    ! Open Modflow name file and read required information.
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_CASETRANS, UTL_RWORD, &
                         UTL_SAMENAME, UTL_STOP
    USE GWM1RMS3, ONLY: HCLOSEG
    USE GWM1BAS3, ONLY: GWMWFILE
    IMPLICIT NONE
    ! Argument
    CHARACTER(LEN=2000),INTENT(IN) ::NAMFILE
    ! Local variables
    INTEGER :: INBAS, INDIS, IST, ISTART, ISTOP, IUNAME, K
    DOUBLE PRECISION :: R
    LOGICAL :: LEX, OKLOCAL
    CHARACTER(LEN=200) :: LINE
    ! Format statements
    600 FORMAT('Error opening file: ',A)
    910 FORMAT('Error reading data from line:',/,A,/)
    !
    OKLOCAL = .TRUE.
    ! Initialize variables
    CALL INITIALIZE_MMPMOD()
    GWMWFILENAME = ''
    !
    ! Open and read name file to get names of package files
    IUNAME = UTL_GETUNIT(7,1000)
    OPEN(IUNAME,FILE=NAMFILE,STATUS='OLD',IOSTAT=IST)
    IF (IST .NE. 0) THEN
      WRITE(*,600)TRIM(NAMFILE)
      GOTO 9100
    ENDIF
    !
    ! Read names of required/supported files from name file
    K = 0
    FIND: DO
      READ(IUNAME,'(A200)',END=100)LINE
      IF (LINE(1:1)=='#') CYCLE FIND
      K = K + 1
      READ(LINE,*,IOSTAT=IST,END=100,ERR=9000)FTYPE(K),IMUNIT(K),FNAME(K)
      CALL UTL_CASETRANS(FTYPE(K),'hi')
      SELECT CASE (FTYPE(K))
      CASE ('LIST')
        LISTFILE=FNAME(K)
      CASE ('DIS')
        DISFILE=FNAME(K)
        DISUNIT=IMUNIT(K)
      CASE ('BAS6')
        BASFILE=FNAME(K)
        BASUNIT=IMUNIT(K)
      CASE ('BCF6')
        BCFFILE=FNAME(K)
        BCFUNIT=IMUNIT(K)
      CASE ('LPF')
        LPFFILE=FNAME(K)
        LPFUNIT=IMUNIT(K)
      CASE ('HUF')
        HUFFILE=FNAME(K)
        HUFUNIT=IMUNIT(K)
      CASE ('UPW')
        UPWFILE=FNAME(K)
        UPWUNIT=IMUNIT(K)
      CASE ('SIP')
        SIPFILE=FNAME(K)
      CASE ('SOR')
        SORFILE=FNAME(K)
      CASE ('PCG')
        PCGFILE=FNAME(K)
      CASE ('PCGN')
        PCGNFILE=FNAME(K)
      CASE ('DE4')
        DE4FILE=FNAME(K)
      CASE ('GMG')
        GMGFILE=FNAME(K)
      CASE ('LMG')
        LMGFILE=FNAME(K)
      CASE ('NWT')
        NWTFILE=FNAME(K)
      CASE ('OC')
        OCFILE=FNAME(K)
      CASE ('SFR')
        SFRFILE=FNAME(K)
      CASE ('DRN')
        DRNFILE=FNAME(K)
      CASE ('WEL')
        WELFILE=FNAME(K)
      CASE ('MNW2')
        MNW2FILE=FNAME(K)
      CASE ('DATA')
        IF(GWMWFILE(1).EQ.IMUNIT(K))THEN          ! A GWMWFILE is called for
          GWMWFILENAME=FNAME(K)                  ! Store the file name
        ENDIF
      CASE DEFAULT
      END SELECT
    ENDDO FIND
    100 CONTINUE
    CLOSE(IUNAME)  ! Close name fil
    !
    ! Open and read files 
    !
    ! Open Scratch file for unneeded file output
    ISCRCH  = UTL_GETUNIT(7,1000)
    !OPEN(UNIT=ISCRCH,STATUS='SCRATCH',ACTION='WRITE')
    ! Writing output to scratch file is a problem if input contains
    ! errors, because error messages are then lost!
    OPEN(UNIT=ISCRCH,FILE='messages_from_gwmvi.txt',STATUS='REPLACE', &
         ACTION='WRITE')
    ! DIS file - read NLAY, NROW, NCOL, NPER, and NSTP 
    !
    ! it has already been read IF (OKLOCAL) CALL READ_DIS_FILE(OKLOCAL)
    ! BAS package - read HNOFLO and IFREFM
    IF (OKLOCAL) CALL READ_BAS_FILE(OKLOCAL)
    IF (.NOT. OKLOCAL) GOTO 9100
    ! Flow package - read HDRY and unit number for flow save (IFLOWCB)
    HDRY = 0.0
    IFLOWCB = 999999
    IF (OKLOCAL) THEN
      IF (BCFFILE .NE. '') THEN
        CALL READ_BCF_FILE(OKLOCAL)
      ELSEIF (LPFFILE .NE. '') THEN
        CALL READ_LPF_FILE(OKLOCAL)
      ELSEIF (HUFFILE .NE. '') THEN
        CALL READ_HUF_FILE(OKLOCAL)
      ELSEIF (UPWFILE .NE. '') THEN
        CALL READ_UPW_FILE(OKLOCAL)
      ELSE       ! No flow package file found; terminate run
        WRITE(GWMOUT,9101)        
        GOTO 9100
      ENDIF
    ENDIF
    ! Get the name of the flow save file - if it exists
    CALL GET_FILENAME(IFLOWCB,IFLOWCBFILE)
    !
    ! Solver package - read HCLOSE
    HCLOSE = SMALLEPS  ! Initialize HCLOSE
    IF (OKLOCAL) THEN
      IF (SIPFILE .NE. '') THEN
        CALL READ_SIP_FILE(OKLOCAL)
      ELSEIF (PCGFILE .NE. '') THEN
        CALL READ_PCG_FILE(OKLOCAL)
      ELSEIF (PCGNFILE .NE. '') THEN
        CALL READ_PCGN_FILE(OKLOCAL)
      ELSEIF (DE4FILE .NE. '') THEN
        CALL READ_DE4_FILE(OKLOCAL)
      ELSEIF (GMGFILE .NE. '') THEN
        CALL READ_GMG_FILE(OKLOCAL)
      ELSEIF (LMGFILE .NE. '') THEN
        CALL READ_LMG_FILE(OKLOCAL)
      ELSEIF (NWTFILE .NE. '') THEN
        CALL READ_NWT_FILE(OKLOCAL)
      ELSEIF (SORFILE .NE. '') THEN
        CALL READ_SOR_FILE(OKLOCAL)
      ENDIF
      HCLOSEG = HCLOSE
    ENDIF
    ! OC package - read unit number for streamflow save (IHEDUN)
    IF (OKLOCAL) THEN
      IF (OCFILE .NE. '')CALL READ_OC_FILE(OKLOCAL)
    ENDIF
    ! Get the name of the head save file - if it exists
    CALL GET_FILENAME(IHEDUN,IHEDUNFILE)
    ! SFR package - read unit number for streamflow save (ISTCB2)
    IF (OKLOCAL) THEN
      IF (SFRFILE .NE. '')CALL READ_SFR_FILE(OKLOCAL)
    ENDIF
    ! Get the name of the streamflow leakage save file - if it exists
    CALL GET_FILENAME(ISTCB1,SFRB1FILE) ! ISTCB1 is positive to make binary file
    ! Get the name of the streamflow flow save file - if it exists
    CALL GET_FILENAME(-ISTCB2,SFRB2FILE) ! ISTCB2 is negative to make binary file
    !
    ! DRN package
    IF (OKLOCAL) THEN
      IF (DRNFILE .NE. '' ) CALL READ_DRN_FILE(OKLOCAL)
    ENDIF
    ! Get the name of the drain save file - if it exists
    CALL GET_FILENAME(IDRNCB,DRNBFILE) ! IDRNCB is positive to make binary file
    IF (.NOT. OKLOCAL) GOTO 9100
    !
    ! Normal return
    ! Close and erase scratch file
    CLOSE(ISCRCH)
    RETURN
    !
    ! Error handling
9000 CONTINUE ! Error reading data from LINE
    WRITE(GWMOUT,910)TRIM(LINE)
    WRITE(GWMOUT,9990)TRIM(NAMFILE)
    CALL UTL_STOP()
9101 FORMAT(/,' STOP EXECUTION - GWM-VI CAN NOT FIND FLOW PACKAGE OF ',&
     /,'   TYPE BCF, LPF, HUF or UPW')    
9990 FORMAT(/,1X,'*** ERROR READING LINE OF FILE "',A,'"',/,&
      2X,'-- STOP EXECUTION (READ_NAME_FILE)')!
    !
9100 CONTINUE ! Error, nature of error already described
    CALL UTL_STOP()
    !
  END SUBROUTINE READ_NAME_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_NAME_FILE_DIS(NAMFILE)
    ! Open Modflow name file and read DIS file information only
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_CASETRANS, UTL_STOP
    IMPLICIT NONE
    ! Argument
    CHARACTER(LEN=2000),INTENT(IN) ::NAMFILE
    ! Local variables
    INTEGER :: INDIS, IST, IUNAME, K
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=200) :: LINE
    ! Format statements
    600 FORMAT('Error opening file: ',A)
    910 FORMAT('Error reading data from line:',/,A,/)
    !
    OKLOCAL = .TRUE.
    ! Initialize variables
    CALL INITIALIZE_MMPMOD()
    !
    ! Open and read name file to get names of package files
    IUNAME = UTL_GETUNIT(7,1000)
    OPEN(IUNAME,FILE=NAMFILE,STATUS='OLD',IOSTAT=IST)
    IF (IST .NE. 0) THEN
      WRITE(*,600)TRIM(NAMFILE)
      GOTO 9100
    ENDIF
    !
    ! Read names of required/supported files from name file
    K = 0
    FIND: DO
      READ(IUNAME,'(A200)',END=100)LINE
      IF (LINE(1:1)=='#') CYCLE FIND
      K = K + 1
      READ(LINE,*,IOSTAT=IST,END=100,ERR=9000)FTYPE(K),IMUNIT(K),FNAME(K)
      CALL UTL_CASETRANS(FTYPE(K),'hi')
      SELECT CASE (FTYPE(K))
      CASE ('DIS')
        DISFILE=FNAME(K)
        DISUNIT=IMUNIT(K)
      CASE DEFAULT
      END SELECT
    ENDDO FIND
    100 CONTINUE
    CLOSE(IUNAME)  ! Close name fil
    ! Open Scratch file for unneeded file output
    ISCRCH  = UTL_GETUNIT(7,1000)
    OPEN(UNIT=ISCRCH,STATUS='SCRATCH',ACTION='WRITE')
    ! DIS file - read NLAY, NROW, NCOL, NPER, and NSTP 
    IF (OKLOCAL) CALL READ_DIS_FILE(OKLOCAL)
    ! Normal return
    ! Close and erase scratch file
    CLOSE(ISCRCH)
    RETURN
    !
    ! Error handling
9000 CONTINUE ! Error reading data from LINE
    WRITE(GWMOUT,910)TRIM(LINE)
    WRITE(GWMOUT,9990)TRIM(NAMFILE)
    CALL UTL_STOP()
9990 FORMAT(/,1X,'*** ERROR READING LINE OF FILE "',A,'"',/,&
      2X,'-- STOP EXECUTION (READ_NAME_FILE)')!
    !
9100 CONTINUE ! Error, nature of error already described
    CALL UTL_STOP()
    !
  END SUBROUTINE READ_NAME_FILE_DIS
  !-------------------------------------------------------------------
  SUBROUTINE READ_BAS_FILE(OK)
    ! Open BAS6 file and read HNOFLO 
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_RWORD
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: ICOL, INBAS, IST, ISTART, ISTOP, K, KK, N
    DOUBLE PRECISION :: R
    CHARACTER(LEN=200) :: LINE
    CHARACTER*24 ANAME(2)
    LOGICAL :: OKLOCAL
    DATA ANAME(1) /'          BOUNDARY ARRAY'/
    DATA ANAME(2) /'            INITIAL HEAD'/
    ! Format statements
    600 FORMAT('Error opening file: ',A)
    910 FORMAT('Error reading data from line:',/,A,/)
    !
    OKLOCAL = .TRUE.
    INBAS=BASUNIT   ! Use unit number specified in NAM file to match any array-reading unit numbers
    OKLOCAL = OPEN_FILE(INBAS,BASFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      ALLOCATE (IBOUND(NCOL,NROW,NLAY))
      !
      ! Read item 1: Options
      CALL URDCOM(INBAS,ISCRCH,LINE)
      !
      ! Look for OPTIONS in the first item after the heading.
      IXSEC=0
      ICHFLG=0
      IFREFM=0
      IPRTIM=0
      ICOL=1
      25 CONTINUE
      CALL UTL_RWORD(ISCRCH,1,.FALSE.,ICOL,LINE,ISTART,ISTOP,N,R,INBAS)
      IF(LINE(ISTART:ISTOP).EQ.'XSECTION') THEN
        IXSEC=1
      ELSE IF(LINE(ISTART:ISTOP).EQ.'CHTOCH') THEN
        ICHFLG=1
      ELSE IF(LINE(ISTART:ISTOP).EQ.'FREE') THEN
        IFREFM=1
      ELSEIF(LINE(ISTART:ISTOP).EQ.'PRINTTIME') THEN
        IPRTIM=1
      END IF
      IF(ICOL.LT.200) GO TO 25
      !
      ! IFREFM must = 1 so that stresses can be read with maximum precision
      IF (IFREFM .NE. 1) THEN
        OKLOCAL = .FALSE.
        CALL APPEND_MESSAGE('Error: BASIC Package input file does not specify FREE option')
        GOTO 9000
      ENDIF
      !
      ! Read IBOUND arrays
      IF(IXSEC.EQ.0) THEN
         DO K=1,NLAY
           KK=K
           CALL U2DINT(IBOUND(:,:,KK),ANAME(1),NROW,NCOL,KK,INBAS,ISCRCH)
         ENDDO
      ELSE
         CALL U2DINT(IBOUND(:,:,1),ANAME(1),NLAY,NCOL,-1,INBAS,ISCRCH)
      END IF
      !
      ! Read HNOFLO
      IF(IFREFM.EQ.0) THEN
        READ(INBAS,'(F10.0)') HNOFLO
      ELSE
        READ(INBAS,*) HNOFLO
      END IF
      !
      CLOSE(INBAS)
    ENDIF
    OKLOCAL = .TRUE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
    ! Error handling
    9000 CONTINUE ! Error encountered
    CALL APPEND_MESSAGE('Error encountered in BASIC Package input file')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_BAS_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_DIS_FILE(OK)
    ! Open DIS file and read NLAY, NROW, NCOL, NPER, and NSTP for each
    ! stress period.  
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_RWORD
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, INDIS, IST, K, KK, NCNFBD, ISS, ITR, N, LLOC, ISTART, &
               ISTOP, MAXTS
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
    INDIS=DISUNIT   ! Use unit number specified in NAM file to match any array-reading unit numbers
    OKLOCAL = OPEN_FILE(INDIS,DISFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      ! Read DIS file to get NLAY, NROW, NCOL, NPER
      CALL URDCOM(INDIS,ISCRCH,LINE)
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
      ALLOCATE (PERLEN(NPER),NSTP(NPER),TSMULT(NPER))
      !C
      !C11-----Read the DELR and DELC arrays.
      CALL U1DREL(DELR,ANAME(1),NCOL,INDIS,ISCRCH)
      CALL U1DREL(DELC,ANAME(2),NROW,INDIS,ISCRCH)
      !C
      !C12-----Read the top elevation of layer 1.
      CALL U2DREL(BOTM(:,:,0),ANAME(3),NROW,NCOL,0,INDIS,ISCRCH)
      !C
      !C13-----Read the bottom elevations.
      DO 120 K=1,NLAY
        KK=K
        CALL U2DREL(BOTM(:,:,LBOTM(K)),ANAME(4),NROW,NCOL,KK,INDIS,ISCRCH)
        IF(LAYCBD(K).NE.0) CALL U2DREL(BOTM(:,:,LBOTM(K)+1),ANAME(5), &
                 NROW,NCOL,KK,INDIS,ISCRCH)
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
        ! Use UTL_RWORD because it supports double precision.
        CALL UTL_RWORD(ISCRCH,3,.FALSE.,LLOC,LINE,ISTART,ISTOP,I,PERLEN(N),INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSTP(N),R,ISCRCH,INDIS)
        CALL UTL_RWORD(ISCRCH,3,.FALSE.,LLOC,LINE,ISTART,ISTOP,I,TSMULT(N),INDIS)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,ISCRCH,INDIS)
      200 CONTINUE
      !
 888     CLOSE(INDIS)
    ENDIF
    IF (OK) OK = OKLOCAL
    CALL SGWF2BAS7PSV(1)  ! Save read data to GLOBALDAT
    RETURN
    !
    ! Error handling
    9000 CONTINUE ! Error reading data from LINE
    WRITE(GWMOUT,910)TRIM(LINE)
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
  SUBROUTINE READ_BCF_FILE(OK)
    ! Open file and read HDRY and unit number used for budget data.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU
    LOGICAL :: OKLOCAL
    !
    IU = BCFUNIT ! Use unit number specified in NAM file to match any array-reading unit numbers
    OKLOCAL = OPEN_FILE(IU,BCFFILE,'OLD')
    IF (OKLOCAL) THEN
      IF(IFREFM.EQ.0) THEN
        READ(IU,'(I10,F10.0)') IBCFCB,HDRY
      ELSE
        READ(IU,*) IBCFCB,HDRY
      END IF
      IFLOWCB = IBCFCB
      !
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_BCF_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_LPF_FILE(OK)
    ! Open file and read HDRY and unit number used for budget data.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ISTART, ISTOP, IU, LLOC
    REAL :: R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = LPFUNIT ! Use unit number specified in NAM file to match any array-reading unit numbers
    OKLOCAL = OPEN_FILE(IU,LPFFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ILPFCB,R,ISCRCH,IU)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,ISCRCH,IU)
      IFLOWCB = ILPFCB
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_LPF_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_HUF_FILE(OK)
    ! Open file and read HDRY and unit number used for budget data.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ISTART, ISTOP, IU, LLOC
    REAL :: R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = HUFUNIT ! Use unit number specified in NAM file to match any array-reading unit numbers
    OKLOCAL = OPEN_FILE(IU,HUFFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHUFCB,R,ISCRCH,IU)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,ISCRCH,IU)
      IFLOWCB = IHUFCB
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_HUF_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_UPW_FILE(OK)
    ! Open file and read HDRY and unit number used for budget data.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ISTART, ISTOP, IU, LLOC
    REAL :: R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UPWUNIT ! Use unit number specified in NAM file to match any array-reading unit numbers
    OKLOCAL = OPEN_FILE(IU,UPWFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUPWCB,R,ISCRCH,IU)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,ISCRCH,IU)
      IFLOWCB = IUPWCB
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_UPW_FILE
  !-------------------------------------------------------------------
  ! READERS FOR SOLVER PACKAGES
  !-------------------------------------------------------------------
  SUBROUTINE READ_SIP_FILE(OK)
    ! Open file and read HCLOSE
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ISTART, ISTOP, IU, LLOC
    DOUBLE PRECISION :: ACCL, R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,SIPFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      READ(IU,*) ACCL,HCLOSE
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_SIP_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_PCG_FILE(OK)
    ! Open file and read HCLOSE
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,PCGFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      READ (IU,*) HCLOSE
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_PCG_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_PCGN_FILE(OK)
    ! Open file and read HCLOSE
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU, ITER_MO, ITER_MI
    REAL :: CLOSE_R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,PCGNFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      READ(LINE,*,ERR=50)ITER_MO, ITER_MI, CLOSE_R, HCLOSE
      GOTO 60
      50 CONTINUE
      OKLOCAL = .FALSE.
      60 CONTINUE
      CLOSE(IU)
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_PCGN_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_SOR_FILE(OK)
    ! Open file and read HCLOSE
    ! MODFLOW-2005 does not support SSOR; however, this subroutine
    ! is included to provide support for SSOR in versions of MODFLOW
    ! other than MODFLOW-2005, just in case.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU, MXITER
    REAL :: ACCL
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,SORFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      READ(LINE,*,ERR=50)MXITER
      READ(IU,*,ERR=50)ACCL,HCLOSE
      GOTO 60
      50 CONTINUE
      OKLOCAL = .FALSE.
      60 CONTINUE
      CLOSE(IU)
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_SOR_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_DE4_FILE(OK)
    ! Open file and read HCLOSE
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, ISTART, ISTOP, IU, LLOC
    INTEGER :: IFREQ, ITMX, MUTD4, MXUP, MXLOW, MXBW
    REAL :: ACCLDE4, R
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,DE4FILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITMX,R,ISCRCH,IU)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXUP,R,ISCRCH,IU)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXLOW,R,ISCRCH,IU)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXBW,R,ISCRCH,IU)
      READ(IU,*) IFREQ,MUTD4,ACCLDE4,HCLOSE
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_DE4_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_GMG_FILE(OK)
    ! Open file and read HCLOSE
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IITER, IU
    DOUBLE PRECISION :: RCLOSEGMG
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,GMGFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      READ(LINE,*) RCLOSEGMG,IITER,HCLOSE
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_GMG_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_LMG_FILE(OK)
    ! Open file and read ... WHAT?
    ! LMG does not use HCLOSE or equivalent, 
    ! so just return with OK unchanged
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    !    
    RETURN
    !
  END SUBROUTINE READ_LMG_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_NWT_FILE(OK)
    ! Open file and read 
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,NWTFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      READ(LINE,*,ERR=50)HCLOSE ! HEADTOL in NWT documentation
      GOTO 60
      50 CONTINUE
      OKLOCAL = .FALSE.
      60 CONTINUE
      CLOSE(IU)
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_NWT_FILE
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  SUBROUTINE READ_SFR_FILE(OK)
    ! Open file and read ISTCB1, ISTCB2, and list of reaches
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IN,IU,LLOC,ISFROPT,IUZT,IRTFLG,NUMTIM,&
               ISTART, ISTOP,I,NSS,NSFRPAR,NPARSEG,NSTRM
    REAL    :: R,CONST,DLEAK
    CHARACTER(LEN=200) :: LINE
    LOGICAL :: OKLOCAL
    !
    IU = UTL_GETUNIT(7,1000)
    IN = IU
    OKLOCAL = OPEN_FILE(IU,SFRFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU, ISCRCH, LINE)
      LLOC = 1
      ISFROPT = 0
      IUZT = 0
      IRTFLG = 0
      NUMTIM = 1
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NSTRM, R, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NSS, R, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NSFRPAR, R, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, NPARSEG, R, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, CONST, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 3, I, DLEAK, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, ISTCB1, R, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, ISTCB2, R, ISCRCH, IU)
      !
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_SFR_FILE
  !-------------------------------------------------------------------
  SUBROUTINE READ_DRN_FILE(OK)
    ! Open file and read DRN file to get IDRNCB
    USE UTILITIES, ONLY: UTL_GETUNIT, UTL_SAMENAME
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU,MXACTD,LLOC,ISTART,ISTOP
    LOGICAL :: OKLOCAL
    REAL    :: R
    CHARACTER(LEN=200) :: LINE, WORD
   !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,DRNFILE,'OLD')
    !
   IF (OKLOCAL) THEN
      CALL URDCOM(IU, ISCRCH, LINE)
      READ(LINE,*)WORD
      IF (UTL_SAMENAME(WORD,'PARAMETER')) CALL URDCOM(IU, ISCRCH, LINE)
      LLOC = 1
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, MXACTD, R, ISCRCH, IU)
      CALL URWORD(LINE, LLOC, ISTART, ISTOP, 2, IDRNCB, R, ISCRCH, IU)
      CLOSE(IU)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_DRN_FILE
  !-------------------------------------------------------------------
  SUBROUTINE INITIALIZE_MMPMOD()
    USE PARAMMODULE, ONLY: MXPAR, MXCLST, MXINST, B, IPLOC, IPCLST, &
                           PARNAM, PARTYP, INAME, IACTIVE, IPSUM, &
                           ICLSUM, INAMLOC, NPVAL
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
    ISTORUN = 0
    IDRNCB = 0
    IMUNIT = 0
    DISUNIT = 0
    BASUNIT = 0
    BCFUNIT = 0
    LPFUNIT = 0
    HUFUNIT = 0
    UPWUNIT = 0
    ! Character
    FTYPE = ''
    FNAME = ''
    BASFILE = ''
    BCFFILE = ''
    DISFILE = ''
    DE4FILE = ''
    DRNFILE = ''
    DRNBFILE = ''
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
    SIPFILE = ''
    SORFILE = ''
    UPWFILE = ''
    WELFILE = ''
    SFRFILE = ''
    SFRB1FILE = ''
    SFRB2FILE = ''
    IHEDUNFILE = ''
    IFLOWCBFILE = ''
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
    RETURN
  END SUBROUTINE INITIALIZE_MMPMOD
  !-------------------------------------------------------------------
  SUBROUTINE READ_OC_FILE(OK)
    ! Open file and read unit number (IHEDUN) for binary file 
    ! containing head arrays.
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Arguments
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IU,LLOC,ISTART,ISTOP,N,IHEDFM,IDDNFM,IDDNUN
    CHARACTER(LEN=200) :: LINE
    REAL :: R
    LOGICAL :: OKLOCAL
    !
    ! Assign default values
    IHEDFM=0
    IDDNFM=0
    IHEDUN=0
    IDDNUN=0
    !C
    !C3------OUTPUT CONTROL IS ACTIVE.  READ FIRST RECORD AND DECODE FIRST
    !C3------WORD.  MUST USE URWORD IN CASE FIRST WORD IS ALPHABETIC.
    !
    IU = UTL_GETUNIT(7,1000)
    OKLOCAL = OPEN_FILE(IU,OCFILE,'OLD')
    !
    IF (OKLOCAL) THEN
      CALL URDCOM(IU,ISCRCH,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,IU)
      !C
      !C4------TEST FOR NUMERIC OUTPUT CONTROL.  FIRST WORD WILL NOT BE
      !C4------"PERIOD", "HEAD", "DRAWDOWN", OR "COMPACT".
      IF(   LINE(ISTART:ISTOP).NE.'PERIOD' .AND. &
            LINE(ISTART:ISTOP).NE.'HEAD' .AND.   &
            LINE(ISTART:ISTOP).NE.'DRAWDOWN' .AND. &
            LINE(ISTART:ISTOP).NE.'COMPACT' .AND. &
            LINE(ISTART:ISTOP).NE.'IBOUND') THEN
      !C4A-----NUMERIC OUTPUT CONTROL.  DECODE THE INITIAL RECORD ACCORDINGLY.
        IF(IFREFM.EQ.0) THEN
          READ(LINE,'(4I10)') IHEDFM,IDDNFM,IHEDUN,IDDNUN
        ELSE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDFM,R,ISCRCH,IU)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNFM,R,ISCRCH,IU)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDUN,R,ISCRCH,IU)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNUN,R,ISCRCH,IU)
        END IF
      ELSE
      !C4B-----ALPHABETIC OUTPUT CONTROL.  CALL MODULE TO READ INITIAL RECORDS.
        CALL MMPSGWF2BAS7J(IU,LINE,LLOC,ISTART,ISTOP)
      END IF
      !
      CLOSE(IU)
    ENDIF
    1000 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_OC_FILE
  !-------------------------------------------------------------------
  SUBROUTINE MMPSGWF2BAS7J(INOC,LINE,LLOC,ISTART,ISTOP)
    !C     ******************************************************************
    !C     READ INITIAL ALPHABETIC OUTPUT CONTROL RECORDS.
    !C     ******************************************************************
    !C
    !C        SPECIFICATIONS:
    !C     ------------------------------------------------------------------
    USE UTILITIES, ONLY: UTL_STOP
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: INOC
    INTEGER, INTENT(INOUT) :: ISTART, ISTOP, LLOC 
    INTEGER :: IHEDFM, IDDNFM, IDDNUN, IBDOPT, LBHDSV, &
               LBDDSV, IAUXSV, IDDREFNEW, N
    REAL :: R
    CHARACTER(LEN=20) :: CHEDFM, CDDNFM
    !C
    CHARACTER*200 LINE
    !C     ------------------------------------------------------------------
    !C
    !C1------ALPHABETIC OUTPUT CONTROL.  
    !C
    !C2B-----LOOK FOR "HEAD PRINT ..." AND "HEAD SAVE ...".  IF
    !C2B-----FOUND, SET APPROPRIATE FLAGS.
100 IF(LINE(ISTART:ISTOP).EQ.'PERIOD') THEN
      GOTO 1000  ! END OF DATA SET ONE HAS BEEN REACHED
    ELSEIF(LINE(ISTART:ISTOP).EQ.'HEAD') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
      IF(LINE(ISTART:ISTOP).EQ.'PRINT') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
        IF(LINE(ISTART:ISTOP).NE.'FORMAT') GO TO 2000
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDFM,R,ISCRCH,INOC)
      ELSE IF(LINE(ISTART:ISTOP).EQ.'SAVE') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
        IF(LINE(ISTART:ISTOP).EQ.'UNIT') THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDUN,R,ISCRCH, &
                 INOC)
        ELSE IF(LINE(ISTART:ISTOP).EQ.'FORMAT') THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,ISCRCH,INOC)
          CHEDFM=LINE(ISTART:ISTOP)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
          IF(LINE(ISTART:ISTOP).EQ.'LABEL') THEN
            LBHDSV=1
          END IF
        ELSE
          GO TO 2000
        END IF
      ELSE
        GO TO 2000
      END IF
      !C
      !C2C-----LOOK FOR "DRAWDOWN PRINT ..." AND "DRAWDOWN SAVE ...".
      !C2C-----IF FOUND, SET APPROPRIATE FLAGS
    ELSE IF(LINE(ISTART:ISTOP).EQ.'DRAWDOWN') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
      IF(LINE(ISTART:ISTOP).EQ.'PRINT') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
        IF(LINE(ISTART:ISTOP).NE.'FORMAT') GO TO 2000
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNFM,R,ISCRCH,INOC)
      ELSE IF(LINE(ISTART:ISTOP).EQ.'SAVE') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
        IF(LINE(ISTART:ISTOP).EQ.'UNIT') THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNUN,R,ISCRCH, &
                        INOC)
        ELSE IF(LINE(ISTART:ISTOP).EQ.'FORMAT') THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,ISCRCH,INOC)
          CDDNFM=LINE(ISTART:ISTOP)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
          IF(LINE(ISTART:ISTOP).EQ.'LABEL') THEN
            LBDDSV=1
          END IF
        ELSE
          GO TO 2000
        END IF
      ELSE
        GO TO 2000
      END IF
      !C
      !C2D-----LOOK FOR "COMPACT BUDGET FILES" -- "COMPACT" IS SUFFICIENT.
      !C2D-----IF FOUND, SET APPROPRIATE FLAG.
    ELSE IF(LINE(ISTART:ISTOP).EQ.'COMPACT') THEN
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
      IF(LINE(ISTART:ISTOP).EQ.'BUDGET') THEN
        IBDOPT=2
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
        IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR. &
            LINE(ISTART:ISTOP).EQ.'AUX') THEN
          IAUXSV=1
        END IF
      ELSE
        GO TO 2000
      END IF
      !C
      !C2E-----LOOK FOR  "IBOUND SAVE ...".  
    ELSE IF(LINE(ISTART:ISTOP).EQ.'IBOUND') THEN
      !C
      !C2F-----ERROR IF UNRECOGNIZED WORD.
    ELSE
      GO TO 2000
    END IF
    !C
    !C3------FINISHED READING A RECORD.  READ NEXT RECORD, IGNORING BLANK
    !C3------LINES.  GO BACK AND DECODE IT.
    110 CONTINUE
    READ(INOC,'(A)',END=1000) LINE
    IF(LINE.EQ.' ') GO TO 110
    LLOC=1
    CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,ISCRCH,INOC)
    GO TO 100
    !C
    !C4------RETURN.
    1000 CONTINUE
    RETURN
    !C
    !C5------ERROR DECODING INPUT DATA.
    2000 CONTINUE
    WRITE(GWMOUT,2001) LINE
    2001 FORMAT(1X,/1X,'ERROR READING OUTPUT CONTROL INPUT DATA:'/1X,A80)
    CALL UTL_STOP()
  END SUBROUTINE MMPSGWF2BAS7J
  !-------------------------------------------------------------------
  ! UTILITIES
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
        IF (IMUNIT(I) == UNIT) THEN
          FILENAME = FNAME(I)
          EXIT SEARCH
        ENDIF
      ENDDO SEARCH
    ENDIF
    !
    RETURN
  END SUBROUTINE GET_FILENAME
  !-------------------------------------------------------------------
END MODULE MKMMPROC_RDNAM
