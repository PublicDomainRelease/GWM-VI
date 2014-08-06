MODULE GWM_STOP
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE GSTOP(STOPMESS,IOUT)
  ! A fatal error has occurred during execution. Terminate the run.
    ! Version called from GWM-VI
    USE GWM_UTLS, ONLY: APPEND_MESSAGE
    USE GWM_SUBS, ONLY: CLEAN_UP
    USE PLLM, ONLY: PLLM_USTOP
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: STOPMESS
    INTEGER, OPTIONAL, INTENT(IN) :: IOUT
    ! Format statements
    10 FORMAT(/,1X,A)
    !
    ! Check to see if it necessary to restore original MNW2 input file
    CALL CLEAN_UP()
    !
    IF (PRESENT(STOPMESS)) THEN
      IF (STOPMESS .NE. ' ') THEN
        CALL APPEND_MESSAGE(STOPMESS)
      ENDIF
    ENDIF
    CALL PLLM_USTOP()
    !
  END SUBROUTINE GSTOP
  !-----------------------------------------------------------------------------
END MODULE GWM_STOP