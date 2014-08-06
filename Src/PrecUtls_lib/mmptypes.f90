MODULE MMPTYPES
  USE GWM_UTLS, ONLY: APPEND_MESSAGE
  IMPLICIT NONE
  PRIVATE
  ! Public data type
  PUBLIC :: DEPENDENT, DEPENDENT_PKG
  ! Public procedures
  PUBLIC :: DEPENDENT_READ, DEPENDENT_PKG_INIT
  !
  TYPE :: DEPENDENT
    CHARACTER(LEN=15) :: CTYPE
    CHARACTER(LEN=20) :: CID
    CHARACTER(LEN=4)  :: CZONE
    CHARACTER(LEN=4)  :: PKGID
    INTEGER :: NCOL, NROW, NLAY
    DOUBLE PRECISION :: SIMVAL
    INTEGER, DIMENSION(7) :: INDX  ! COL, ROW, LAY, SP, except:
    ! For STORAGE-CHANGE dependents:
    !   INDX(1) = Begin-Stress-Period 
    !   INDX(2) = End-Stress-Period
    !   INDX(3) = 999 is used to indicate that computation of
    !                 storage DEPENDENT is complete.
    ! For SFR-FLOW and SFR-LEAK dependents:
    !   INDX(1) = SEG
    !   INDX(2) = REACH
    !   INDX(3) = SP
    !   INDX(4) = SEQNUM (index in list written to binary budget file)
    ! For DRN-FLOW dependents:
    !   INDX(1) = COL
    !   INDX(2) = ROW
    !   INDX(3) = LAY
    !   INDX(4) = AUXVAL
    !   INDX(5) = Stress-Period
    ! For DRN-VOLUME dependents:
    !   INDX(1) = COL
    !   INDX(2) = ROW
    !   INDX(3) = LAY
    !   INDX(4) = AUXVAL
    !   INDX(5) = Begin-Stress-Period
    !   INDX(6) = End-Stress-Period
    !   INDX(7) = 999 is used to indicate that computation of
    !             DRN-VOLUME DEPENDENT is complete
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ZONE
  END TYPE DEPENDENT
  !
  INTERFACE DEPENDENT_CONSTRUCT
    ! A generic interface for allocating one DEPENDENT of any type
    MODULE PROCEDURE DEPENDENT_HD_CONSTRUCT
    MODULE PROCEDURE DEPENDENT_SFR_CONSTRUCT
    MODULE PROCEDURE DEPENDENT_STOR_CONSTRUCT
    MODULE PROCEDURE DEPENDENT_DRNF_CONSTRUCT
    MODULE PROCEDURE DEPENDENT_DRNV_CONSTRUCT
  END INTERFACE
  !
  TYPE :: DEPENDENT_PKG
    CHARACTER(LEN=4)   :: PKGID
    INTEGER            :: UNIT
    CHARACTER(LEN=200) :: FILENAME
  END TYPE DEPENDENT_PKG
  !
CONTAINS
  !
  !-------------------------------------------------------------------
  ! Specific subroutines that implement interface DEPENDENT_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_HD_CONSTRUCT(ID,DEPTYPE,COL,ROW,LAY,SP,DEP)
    ! Construct a DEPENDENT of type HEAD
    USE UTILITIES, ONLY: UTL_CASE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: ID
    CHARACTER(LEN=*), INTENT(IN) :: DEPTYPE
    INTEGER, INTENT(IN) :: COL, ROW, LAY, SP
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    ! Local variables
    CHARACTER(LEN=15) :: TYPELOCAL
    !
    CALL UTL_CASE(DEPTYPE,TYPELOCAL,1)
    DEP%CTYPE = TYPELOCAL
    DEP%CID = ID
    DEP%CZONE = ''
    IF (TYPELOCAL == 'HEAD') THEN
      DEP%PKGID = 'HEAD'
    ENDIF
    ! For a HEAD DEPENDENT, INDX contains
    ! (Column, Row, Layer, StressPeriod)
    DEP%INDX(1) = COL
    DEP%INDX(2) = ROW
    DEP%INDX(3) = LAY
    DEP%INDX(4) = SP
    DEP%INDX(5) = -999
    DEP%INDX(6) = -999
    DEP%INDX(7) = -999
    DEP%NCOL = -999
    DEP%NROW = -999
    DEP%NLAY = -999
    DEP%SIMVAL = -99999.0D0
    RETURN
  END SUBROUTINE DEPENDENT_HD_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_SFR_CONSTRUCT(ID,DEPTYPE,SEG,REACH,SP,DEP)
    ! Construct a DEPENDENT of type SFR-FLOW or SFR-LEAK
    USE UTILITIES, ONLY: UTL_CASE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: ID
    CHARACTER(LEN=*), INTENT(IN) :: DEPTYPE
    INTEGER, INTENT(IN) :: SEG, REACH, SP
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    ! Local variables
    CHARACTER(LEN=15) :: TYPELOCAL
    !
    CALL UTL_CASE(DEPTYPE,TYPELOCAL,1)
    DEP%CTYPE = TYPELOCAL
    DEP%CID = ID
    DEP%CZONE = ''
    IF (TYPELOCAL == 'SFR-FLOW') THEN
      DEP%PKGID = 'SFRF'
    ELSEIF (TYPELOCAL == 'SFR-LEAK') THEN
      DEP%PKGID = 'SFRL'
    ENDIF
    ! For an SFR DEPENDENT, INDX contains (Segment,Reach,StressPeriod)
    DEP%INDX(1) = SEG
    DEP%INDX(2) = REACH
    DEP%INDX(3) = SP     ! Stress period
    DEP%INDX(4) = -999
    DEP%INDX(5) = -999
    DEP%INDX(6) = -999
    DEP%INDX(7) = -999
    DEP%NCOL = -999
    DEP%NROW = -999
    DEP%NLAY = -999
    DEP%SIMVAL = -99999.0D0
    RETURN
  END SUBROUTINE DEPENDENT_SFR_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_STOR_CONSTRUCT(ID,NCOL,NROW,NLAY,BEGSP, &
                                       ENDSP,CZONE,DEP)
    ! Construct a DEPENDENT of type STORAGE-CHANGE
    USE UTILITIES, ONLY: UTL_CASE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: ID
    INTEGER, INTENT(IN) :: NCOL, NROW, NLAY, BEGSP, ENDSP
    CHARACTER(LEN=*), INTENT(IN) :: CZONE
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    ! Local variables
    CHARACTER(LEN=4) :: CZONELOCAL
    !
    DEP%CTYPE = 'STORAGE-CHANGE'
    DEP%CID = ID
    CALL UTL_CASE(CZONE,CZONELOCAL,1)
    DEP%CZONE = CZONELOCAL
    DEP%PKGID = 'STOR'
    ! For a STORAGE-CHANGE DEPENDENT, INDX contains 
    ! (Begin-Stress-Period, End-Stress-Period)
    DEP%INDX(1) = BEGSP
    DEP%INDX(2) = ENDSP
    DEP%INDX(3) = -999 ! Change to 999 when calculation of SIMVAL is complete
    DEP%INDX(4) = -999
    DEP%INDX(5) = -999
    DEP%INDX(6) = -999
    DEP%INDX(7) = -999
    DEP%NCOL = NCOL
    DEP%NROW = NROW
    DEP%NLAY = NLAY
    DEP%SIMVAL = 0.0D0
    IF (CZONELOCAL=='ZONE') THEN
      ALLOCATE(DEP%ZONE(NCOL,NROW,NLAY))
    ENDIF
    RETURN
  END SUBROUTINE DEPENDENT_STOR_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_DRNF_CONSTRUCT(ID,COL,ROW,LAY,AUXVAL, &
                                      SP,DEP)
    ! Construct a DEPENDENT of type DRN-FLOW
    USE UTILITIES, ONLY: UTL_CASE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: ID
    INTEGER, INTENT(IN) ::  COL, ROW, LAY, AUXVAL, SP
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    ! Local variables
    !
    DEP%CTYPE = 'DRN-FLOW'
    DEP%CID = ID
    DEP%PKGID = 'DRNF'
    DEP%INDX(1) = COL
    DEP%INDX(2) = ROW
    DEP%INDX(3) = LAY
    DEP%INDX(4) = AUXVAL
    DEP%INDX(5) = SP
    DEP%INDX(6) = -999
    DEP%INDX(7) = -999
    DEP%NCOL = -999
    DEP%NROW = -999
    DEP%NLAY = -999
    DEP%SIMVAL = 0.0D0
    RETURN
  END SUBROUTINE DEPENDENT_DRNF_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_DRNV_CONSTRUCT(ID,COL,ROW,LAY,AUXVAL, &
                                      BEGSP,ENDSP,DEP)
    ! Construct a DEPENDENT of type DRN-VOLUME
    USE UTILITIES, ONLY: UTL_CASE
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: ID
    INTEGER, INTENT(IN) :: COL, ROW, LAY, AUXVAL, BEGSP, ENDSP
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    ! Local variables
    !
    DEP%CTYPE = 'DRN-VOLUME'
    DEP%CID = ID
    DEP%PKGID = 'DRNV'
    DEP%INDX(1) = COL
    DEP%INDX(2) = ROW
    DEP%INDX(3) = LAY
    DEP%INDX(4) = AUXVAL
    DEP%INDX(5) = BEGSP
    DEP%INDX(6) = ENDSP
    DEP%INDX(7) = -999 ! Change to 999 when calculation of SIMVAL is complete
    DEP%NCOL = -999
    DEP%NROW = -999
    DEP%NLAY = -999
    DEP%SIMVAL = 0.0D0
    RETURN
  END SUBROUTINE DEPENDENT_DRNV_CONSTRUCT
  !-------------------------------------------------------------------
  ! End of subroutines that implement interface DEPENDENT_CONSTRUCT
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_READ(IU,NCOL,NROW,NLAY,DEP,OK)
    ! Read one DEPENDENT from unit IU; store data in DEPENDENT DEP.
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_RWORD
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, NCOL, NROW, NLAY
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: N
    DOUBLE PRECISION :: R
    CHARACTER(LEN=200) :: LINE
    CHARACTER(LEN=20)  :: CID
    CHARACTER(LEN=15)  :: CTYPE
    CHARACTER(LEN=4)   :: CZONE
    INTEGER :: AUXVAL, ICOL, ISTART, ISTOP, COL, ROW, LAY, SEG, &
               REACH, SP1, SP2
    LOGICAL :: OKLOCAL
    !
    ! Read item A
    OKLOCAL = .TRUE.
    READ(IU,'(A200)',END=900)LINE
    ICOL = 1
    CALL UTL_RWORD(0,0,.TRUE.,ICOL,LINE,ISTART,ISTOP,N,R)
    CID = LINE(ISTART:ISTOP)
    CALL UTL_RWORD(0,1,.TRUE.,ICOL,LINE,ISTART,ISTOP,N,R)
    CTYPE = LINE(ISTART:ISTOP)
    !
    ! Read items B and (for type STORAGE-CHANGE) C
    SELECT CASE (CTYPE)
    CASE ('HEAD')
      READ(IU,*,ERR=910)LAY,ROW,COL,SP1
      CALL DEPENDENT_CONSTRUCT(CID,CTYPE,COL,ROW,LAY,SP1,DEP)
    CASE ('SFR-FLOW', 'SFR-LEAK')
      READ(IU,*,ERR=910)SEG,REACH,SP1
      CALL DEPENDENT_CONSTRUCT(CID,CTYPE,SEG,REACH,SP1,DEP)
    CASE ('STORAGE-CHANGE')
      READ(IU,*,ERR=910)SP1,SP2,CZONE
      CALL DEPENDENT_CONSTRUCT(CID,NCOL,NROW,NLAY,SP1,SP2,CZONE,DEP)
      ! Read item C
      CALL DEPENDENT_READ_ZONES(IU,DEP,OK)
      IF (.NOT. OK) GOTO 905
    CASE ('DRN-FLOW')
      READ(IU,*,ERR=910)LAY,ROW,COL,AUXVAL,SP1
      CALL DEPENDENT_CONSTRUCT(CID,COL,ROW,LAY,AUXVAL,SP1,DEP)
    CASE ('DRN-VOLUME')
      READ(IU,*,ERR=910)LAY,ROW,COL,AUXVAL,SP1,SP2
      CALL DEPENDENT_CONSTRUCT(CID,COL,ROW,LAY,AUXVAL,SP1,SP2,DEP)
    CASE DEFAULT
      CALL APPEND_MESSAGE('Error: Unrecognized DEPENDENT type "' // &
                 TRIM(CTYPE) // '"')
      GOTO 905
    END SELECT
    !
    RETURN
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error encountered reading DEPENDENT')
    905 CONTINUE
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    !
    910 CONTINUE
    CALL APPEND_MESSAGE('Error reading DEPENDENT "' // TRIM(CID) // '" of type ' &
               // TRIM(CTYPE))
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE DEPENDENT_READ
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_READ_ZONES(IU,DEP,OK)
    ! Read and store a 3-D Zone array
    USE GLOBAL_DATA, ONLY: AMESSAGE
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU
    TYPE (DEPENDENT), INTENT(INOUT) :: DEP
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, J, K
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    IF (DEP%CZONE == 'ZONE') THEN
      DO K=1,DEP%NLAY
        DO I=1,DEP%NROW
          READ(IU,*,ERR=900)(DEP%ZONE(J,I,K),J=1,DEP%NCOL)
        ENDDO
      ENDDO
    ENDIF
    RETURN
    !
    ! Error handling
    900 CONTINUE
    CALL APPEND_MESSAGE('Error reading ZONE array for DEPENDENT "' // &
               TRIM(DEP%CID) // '"')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE DEPENDENT_READ_ZONES
  !-------------------------------------------------------------------
  SUBROUTINE DEPENDENT_PKG_INIT(CP)
    IMPLICIT NONE
    ! Argument
    TYPE (DEPENDENT_PKG), INTENT(INOUT) :: CP
    CP%UNIT = 0
    CP%PKGID = ''
    CP%FILENAME = ''
    RETURN
  END SUBROUTINE DEPENDENT_PKG_INIT
  !-------------------------------------------------------------------
END MODULE MMPTYPES