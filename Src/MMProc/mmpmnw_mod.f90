MODULE MMPMNWMOD
  PRIVATE
  ! Public variables
  ! Public type and program units
  PUBLIC :: MNW2WELL, MNW2WELL_ALLOCATE, MNW2WELL_DEALLOC, &
            READ_MNW2_DATA
  !
  TYPE MNW2WELL
    ! Holds data related to one MNW2 well for which heads in cells
    ! included in the well need to be evaluated.
    CHARACTER(LEN=20) :: ID
    INTEGER :: NC  ! Number of cells included in this MNW2 well
    ! INDX contains layer, row, column indices of each cell
    ! INDX(i,1) = Layer
    ! INDX(i,2) = Row
    ! INDX(i,3) = Column
    INTEGER :: NPER ! Number of stress periods, and dimension of IACTIVE
    INTEGER :: KDRY ! Number of dry cells along MNW well in current time step
    INTEGER :: AUXVALUE ! Unique integer identifier for this MNW2 well
    LOGICAL :: DRY
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDX  
    ! IACTIVE includes an element for each stress period
    ! If MNW2 well is not active, IACTIVE(i) = 0
    ! If MNW2 well is active, IACTIVE(i) = 1
    INTEGER, ALLOCATABLE, DIMENSION(:) :: IACTIVE
    DOUBLE PRECISION :: QDES
  END TYPE MNW2WELL
  !
CONTAINS
  !---------------------------------------------------------------------
  SUBROUTINE MNW2WELL_ALLOCATE(WELL,NCELLS,NPER,WELLID,QDES,AUXVAL)
    IMPLICIT NONE
    ! Arguments
    TYPE (MNW2WELL), INTENT(INOUT) :: WELL
    INTEGER, INTENT(IN) :: NCELLS, NPER, AUXVAL
    CHARACTER(LEN=*), INTENT(IN) :: WELLID
    DOUBLE PRECISION, INTENT(IN) :: QDES
    !
    IF (ALLOCATED(WELL%INDX)) DEALLOCATE(WELL%INDX)
    IF (ALLOCATED(WELL%IACTIVE)) DEALLOCATE(WELL%IACTIVE)
    !
    WELL%ID = WELLID
    WELL%NC = NCELLS
    WELL%NPER = NPER
    WELL%KDRY = 0
    WELL%AUXVALUE = AUXVAL
    WELL%DRY = .FALSE.
    WELL%QDES = QDES
    IF (NCELLS>0 .AND. NPER>0) THEN
      ALLOCATE(WELL%INDX(NCELLS,3))
      WELL%INDX = 0
      ALLOCATE(WELL%IACTIVE(NPER))
      WELL%IACTIVE = 0
    ENDIF
    !
    RETURN
  END SUBROUTINE MNW2WELL_ALLOCATE
  !---------------------------------------------------------------------
  SUBROUTINE MNW2WELL_DEALLOC(WELL)
    IMPLICIT NONE
    ! Argument
    TYPE (MNW2WELL), INTENT(INOUT) :: WELL
    IF (ALLOCATED(WELL%INDX)) DEALLOCATE(WELL%INDX)
    IF (ALLOCATED(WELL%IACTIVE)) DEALLOCATE(WELL%IACTIVE)
    WELL%NC = 0
    WELL%NPER = 0
    WELL%KDRY = 0
    WELL%AUXVALUE = -999
    WELL%ID = ''
    WELL%QDES = 0.0D0
    RETURN
  END SUBROUTINE MNW2WELL_DEALLOC
  !-------------------------------------------------------------------
  SUBROUTINE READ_MNW2_DATA(IU,NMQMNW,NPER,WELLS,OK)
    USE GWM_UTLS, ONLY: APPEND_MESSAGE
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, NMQMNW, NPER
    TYPE (MNW2WELL), DIMENSION(NMQMNW) :: WELLS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: AUXVAL, LAY, ROW, COL, M, N, NCELLS
    DOUBLE PRECISION :: QDES
    CHARACTER(LEN=20) :: LASTID, WELLID
    CHARACTER(LEN=80) :: ERRMSG
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    ERRMSG = ' '
    LASTID = ' '
    DO N=1,NMQMNW
      ! Read item 15a
      WELLID = ' '
      READ(IU,*,ERR=100)WELLID,NCELLS,QDES,AUXVAL
      CALL MNW2WELL_ALLOCATE(WELLS(N),NCELLS,NPER,WELLID,QDES,AUXVAL)
      DO M=1,NCELLS
        ! Read and store item 15b
        READ(IU,*,ERR=100)LAY,ROW,COL
        WELLS(N)%INDX(M,1) = LAY
        WELLS(N)%INDX(M,2) = ROW
        WELLS(N)%INDX(M,3) = COL
      ENDDO
      ! Read and store item 15c
      READ(IU,*,ERR=100)(WELLS(N)%IACTIVE(M),M=1,NPER)
      LASTID = WELLID
    ENDDO
    GOTO 200
    !
    ! Error processing
    100 CONTINUE
    OKLOCAL = .FALSE.
    IF (WELLID .NE. ' ') THEN
      ERRMSG = 'Error reading MNW data for well: ' // TRIM(WELLID)
    ELSE
      ERRMSG = 'Error reading MNW data. Last well successfully read: ' &
               // TRIM(LASTID)
    ENDIF
    CALL APPEND_MESSAGE(ERRMSG)
    !
    200 CONTINUE
    IF (OK) OK = OKLOCAL
    !
    RETURN
  END SUBROUTINE READ_MNW2_DATA
  !---------------------------------------------------------------------
END MODULE MMPMNWMOD
