MODULE MGT
  PUBLIC MGT_FIND_STATUS, MGT_INI, MGT_SETSTATUS, MGT_WRITESTATUS, &
         HCLOSEM, DEWATERQ, PLL_ACTIVE
  PRIVATE
    INTEGER :: MFSTATUS
    REAL :: HCLOSEM, HDRYM, HNOFLOM
    LOGICAL ::DEWATERQ, PLL_ACTIVE   
    CHARACTER(LEN=40), DIMENSION(-1:2) :: STATUSTEXT
    CHARACTER(LEN=50), DIMENSION(0:1)  :: DEWATERTEXT
    CHARACTER(LEN=2000) :: MODFLOWSTATUS = 'modflow.status'
    DATA STATUSTEXT/ &
        'Simulation started','Simulation successfully completed', &
        'Error: Failure to converge', &
        'Error: Other than failure to converge'/
    DATA DEWATERTEXT/'No managed-flow cells converted to dry',  &
         'At least one managed-flow cell converted to dry'/
CONTAINS
  !-------------------------------------------------------------------------
  SUBROUTINE MGT_INI()
    IMPLICIT NONE
    ! Local variables
    LOGICAL :: OKLOCAL
    !
    MFSTATUS = -1
    DEWATERQ = .FALSE.
    CALL MGT_WRITESTATUS(OKLOCAL)
    RETURN
  END SUBROUTINE MGT_INI
  !-------------------------------------------------------------------------
  SUBROUTINE MGT_FIND_STATUS(NCOL,NROW,LISTFILE,MANAGED_WELLS,NMQMNW, &
                             MNW2WELLS,HEADFILEEMPTY,OK,NPER,NSTP)
    ! Assign MFSTATUS and DEWATERQ by evaluating Modflow output
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE PARALLEL_PROCESSING, ONLY: PLL_WAIT
    USE GWM_UTLS, ONLY: APPEND_MESSAGE, OPEN_FILE
    USE MMPWELMOD, ONLY: WELTYPE
    USE MMPMNWMOD, ONLY: MNW2WELL
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: NCOL, NROW, NMQMNW
    CHARACTER(LEN=*), INTENT(IN) :: LISTFILE
    TYPE (WELTYPE), INTENT(IN) :: MANAGED_WELLS
    TYPE (MNW2WELL), DIMENSION(NMQMNW), INTENT(IN) :: MNW2WELLS
    LOGICAL, INTENT(IN) :: HEADFILEEMPTY
    LOGICAL, INTENT(INOUT) :: OK
    INTEGER, INTENT(IN) :: NPER
    INTEGER, DIMENSION(NPER), INTENT(IN) :: NSTP
    ! Local variables
    INTEGER :: I, IU, IFIND, ISTAT, KTRY
    CHARACTER(LEN=200) :: LINE, TIMESUMMARY
    LOGICAL :: COMPLETE, CONVERGED, FINALTIMESUMMARY, LOPEN, MGTCELLDRY, OKLOCAL, RUNEND
    !
    MFSTATUS = 2 ! Default: Error: Other than failure to converge
    COMPLETE = .FALSE.
    CONVERGED = .TRUE.
    MGTCELLDRY = .FALSE.
    RUNEND = .FALSE.
    FINALTIMESUMMARY = .FALSE.
    WRITE(TIMESUMMARY,15)NSTP(NPER),NPER
    15 FORMAT('TIME SUMMARY AT END OF TIME STEP',I4,' IN STRESS PERIOD',I5)
    IU = UTL_GETUNIT(7,9000)
    OKLOCAL = OPEN_FILE(IU,LISTFILE,'OLD')
    IF (OKLOCAL) THEN
      DO
        READ(IU,'(A200)',END=10)LINE
        IFIND = INDEX(LINE,'Normal termination')
        IF (IFIND>0) COMPLETE = .TRUE.
        IFIND = INDEX(LINE,'Run end date and time')
        IF (IFIND>0) RUNEND = .TRUE.
        IFIND = INDEX(LINE,'FAILED TO MEET SOLVER CONVERGENCE')
        IF (IFIND>0) CONVERGED = .FALSE.
        IFIND = INDEX(LINE,'Failure to converge')
        IF (IFIND>0) CONVERGED = .FALSE.
        IFIND = INDEX(LINE,'FAILED TO CONVERGE') ! This is what SEAWAT writes
        IF (IFIND>0) CONVERGED = .FALSE.
        IFIND = INDEX(LINE,TRIM(TIMESUMMARY))
        IF (IFIND>0) FINALTIMESUMMARY = .TRUE.
      ENDDO
      10 CONTINUE
      !
      ! SEAWAT does not write 'Normal termination' or
      ! 'Run end date and time' to list file.  Check for this case.
      !
      ! The potential problem with this logic is that if a version of
      ! Modflow uses a convergence-failure messages other than one of
      ! those that appear in the above DO loop, CONVERGED could be true 
      ! even if the Modflow run did not converge.
      !
      IF (CONVERGED .AND. .NOT. COMPLETE .AND. .NOT. RUNEND) THEN
        IF (FINALTIMESUMMARY) COMPLETE = .TRUE.
      ENDIF
      !
      KTRY = 0
      20 CONTINUE
      CLOSE(IU)
      INQUIRE(UNIT=IU,OPENED=LOPEN)
      IF (LOPEN) THEN
        KTRY = KTRY + 1
        IF (KTRY < 200) THEN
          CALL PLL_WAIT()
          GOTO 20
        ENDIF
      ENDIF
      !
      IF (.NOT. CONVERGED) THEN
        ISTAT = 1   ! Error: Failure to converge
      ELSE
        IF (COMPLETE .OR. RUNEND) THEN
          ISTAT = 0 ! Simulation successfully completed
        ELSEIF (HEADFILEEMPTY) THEN
          ISTAT = 2 ! Error: Other than failure to converge
        ELSE
          ISTAT = 2 ! Error: Other than failure to converge
        ENDIF
      ENDIF
      MFSTATUS = ISTAT
      !
!      IF (ISTAT==0) THEN
        ! Assign DEWATERQ
        IF (MANAGED_WELLS%MXNPW > 0) THEN
          DO I=1,MANAGED_WELLS%MXNPW
            IF (MANAGED_WELLS%DRY(I)) DEWATERQ = .TRUE.
          ENDDO
        ENDIF
        IF (NMQMNW > 0) THEN
          DO I=1,NMQMNW
            IF (MNW2WELLS(I)%DRY) DEWATERQ = .TRUE.
          ENDDO
        ENDIF
!      ENDIF
      !
    ELSE
      CALL APPEND_MESSAGE('Error: Unable to open MODFLOW list file "' &
                          // TRIM(LISTFILE) // '"')
    ENDIF
    !
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE MGT_FIND_STATUS
  !-------------------------------------------------------------------------
  SUBROUTINE MGT_SETSTATUS(I,HDRY,HNOFLO)
    !  I=-1 -- Simulation started
    !  I=0  -- Simulation successfully completed
    !  I=1  -- Error: Failure to converge
    !  I=2  -- Error: Other than failure to converge
    IMPLICIT NONE
    !   Argument-list variables
    INTEGER, INTENT(IN) :: I
    REAL, INTENT(IN) :: HDRY, HNOFLO
    !
    IF (I.GE.-1 .AND. I.LE.2) THEN
      MFSTATUS = I
    ENDIF
    IF (I>-1 .AND. I<2) THEN
      HDRYM = HDRY
      HNOFLOM = HNOFLO
    ENDIF
    RETURN
  END SUBROUTINE MGT_SETSTATUS
  !-------------------------------------------------------------------------
  SUBROUTINE MGT_WRITESTATUS(OK)
    USE UTILITIES, ONLY: UTL_GETUNIT
    USE GWM_UTLS, ONLY: OPEN_FILE
    IMPLICIT NONE
    ! Argument
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: N, IDEWATER
    LOGICAL :: OKLOCAL
    !   Formats
    10 FORMAT(1X,'Modflow status: ',I2,2X,'(',A,')')
    20 FORMAT(1X,'Dewater status: ',I2,2X,'(',A,')')
    !
    N = UTL_GETUNIT(7,2000)
    OKLOCAL = OPEN_FILE(N,MODFLOWSTATUS,'REPLACE')
    IF (OKLOCAL) THEN
      WRITE(N,10)MFSTATUS,TRIM(STATUSTEXT(MFSTATUS))
      IF (DEWATERQ) THEN
        IDEWATER = 1
      ELSE
        IDEWATER = 0
      ENDIF
      WRITE(N,20)IDEWATER,TRIM(DEWATERTEXT(IDEWATER))
      CLOSE(N)
    ENDIF
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE MGT_WRITESTATUS
  !-------------------------------------------------------------------------
END MODULE MGT