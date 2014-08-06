MODULE GWM_UTLS
  PUBLIC :: APPEND_MESSAGE, COPY_FILE, FINDCELL, &
            OPEN_FILE, SAVE_FILE, TEXT_FILE_COPY
CONTAINS
  !-------------------------------------------------------------------
  SUBROUTINE APPEND_MESSAGE(ERRMSG)
    USE GLOBAL_DATA, ONLY: AMESSAGE
    IMPLICIT NONE
    ! Argument
    CHARACTER(LEN=*), INTENT(IN) :: ERRMSG
    ! Local variables
    !
    IF (AMESSAGE == '') THEN
      AMESSAGE = ERRMSG
    ELSE
      AMESSAGE = TRIM(AMESSAGE) // ' :: ' // ERRMSG
    ENDIF
    RETURN
  END SUBROUTINE APPEND_MESSAGE
  !-------------------------------------------------------------------
  LOGICAL FUNCTION COPY_FILE(FROMFILE,TOFILE)
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_GETUNIT
    IMPLICIT NONE
    ! Arguments
    CHARACTER(LEN=*), INTENT(IN) :: FROMFILE, TOFILE
    ! Local variables
    INTEGER, PARAMETER :: MAXLEN=1000
    INTEGER :: IUFROM, IUTO
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=MAXLEN) :: LINE
    CHARACTER(LEN=200) :: ERRMSG
    !
    IUFROM = UTL_GETUNIT(7,9000)
    OKLOCAL = OPEN_FILE(IUFROM,FROMFILE,'OLD')
    IF (OKLOCAL) THEN
      IUTO = UTL_GETUNIT(7,9000)
      OKLOCAL = OPEN_FILE(IUTO,TOFILE,'REPLACE')
      IF (OKlOCAL) THEN
        DO
          READ(IUFROM,'(A1000)',END=20)LINE
          WRITE(IUTO,'(A)',ERR=10)TRIM(LINE)
        ENDDO
        10 CONTINUE
        OKLOCAL = .FALSE.
        20 CONTINUE
        CLOSE(IUTO)  
      ENDIF
      CLOSE(IUFROM)
    ENDIF
    IF (.NOT. OKLOCAL) THEN
      ERRMSG = 'Error copying file "' // TRIM(FROMFILE) // &
                '" to "' // TRIM(TOFILE) // '"'
      CALL APPEND_MESSAGE(ERRMSG)
    ENDIF
    COPY_FILE = OKLOCAL
    RETURN
  END FUNCTION COPY_FILE
  !-------------------------------------------------------------------
  SUBROUTINE FINDCELL(NODE,NROW,NCOL,NLAY,I,J,K)
    ! Return row (I), column (J), and layer (K) that correspond to a
    ! node number
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: NODE, NROW, NCOL, NLAY
    INTEGER, INTENT(OUT) :: I, J, K
    INTEGER :: N, NODES, NRC
    !
    NRC=NROW*NCOL
    NODES=NRC*NLAY
    !
    IF (NODE.LE.0 .OR. NODE.GT.NODES) THEN
      K = 0
      I = 0
      J = 0
    ELSE
      K=(NODE-1)/NRC+1  ! Layer index
      N=NODE-(K-1)*NRC
      I=(N-1)/NCOL+1    ! Row index
      J=N-(I-1)*NCOL    ! Column index
    ENDIF
    !
    RETURN
  END SUBROUTINE FINDCELL
  !-------------------------------------------------------------------
  LOGICAL FUNCTION OPEN_FILE(IU,FNAME,STATUS,BINARY)
    ! Open file FILENAME on unit IU
    USE GLOBAL_DATA, ONLY: AMESSAGE
    USE UTILITIES, ONLY: UTL_CASE
    USE PARALLEL_PROCESSING, ONLY: PLL_WAIT
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU
    CHARACTER(LEN=*), INTENT(IN) :: FNAME, STATUS
    LOGICAL, INTENT(IN), OPTIONAL :: BINARY
    ! Local variables
    INTEGER :: I, ISTAT, KTRY
    LOGICAL :: BINLOCAL, LEX, LOPEN
    CHARACTER(LEN=10) :: STATUSLOCAL
    !
    INCLUDE 'openspec_f90.inc'
    OPEN_FILE = .FALSE.
    BINLOCAL = .FALSE.
    IF (PRESENT(BINARY)) BINLOCAL = BINARY
    CALL UTL_CASE(STATUS,STATUSLOCAL,1)
    IF (STATUSLOCAL == 'OLD') THEN
      INQUIRE(FILE=FNAME,EXIST=LEX)
      IF (.NOT. LEX) THEN
        CALL APPEND_MESSAGE('Error: File "' // TRIM(FNAME) // &
                            '" does not exist.')
      ENDIF
    ENDIF
    !
    KTRY = 0
    10 CONTINUE
    IF (BINLOCAL) THEN
      OPEN(IU,FILE=FNAME,STATUS=STATUSLOCAL,FORM=FORM,IOSTAT=ISTAT)
    ELSE
      OPEN(IU,FILE=FNAME,STATUS=STATUSLOCAL,IOSTAT=ISTAT)
    ENDIF
    IF (ISTAT==0) THEN
      OPEN_FILE = .TRUE.
    ELSE
      KTRY = KTRY + 1
      IF (KTRY < 20) THEN
        CALL PLL_WAIT()
        GOTO 10
      ELSE
        CALL APPEND_MESSAGE('Error opening file "' // TRIM(FNAME) &
                            // '"')
      ENDIF
    ENDIF
    RETURN
  END FUNCTION OPEN_FILE
  !-----------------------------------------------------------------------------
  SUBROUTINE SAVE_FILE(IU,FNAME)
    USE UTILITIES
    IMPLICIT NONE
    !   Argument-list variables
    INTEGER,          INTENT(IN) :: IU
    CHARACTER(LEN=*), INTENT(IN) :: FNAME
    !   Local variables
    CHARACTER(LEN=2000) :: LINE
    INTEGER :: ISTAT, IO
    !   Format statements
    10 FORMAT(A2000)
    !
    IO = UTL_GETUNIT(7,2000)
    OPEN(IO,FILE=FNAME,STATUS='REPLACE')
    REWIND(IU)
    !
    COPY: DO
      READ(IU,10,IOSTAT=ISTAT)LINE
      IF (ISTAT .NE. 0) EXIT COPY
      WRITE(IO,*)TRIM(LINE)
    ENDDO COPY
    CLOSE(IO)
    !
    RETURN
  END SUBROUTINE SAVE_FILE
  !-------------------------------------------------------------------
  LOGICAL FUNCTION TEXT_FILE_COPY(TXT_FILE,TXT_FILE_COPY)
    IMPLICIT NONE
    ! Argument
    CHARACTER(LEN=*), INTENT(IN) ::  TXT_FILE
    CHARACTER(LEN=*), INTENT(OUT) :: TXT_FILE_COPY
    !
    TXT_FILE_COPY = TRIM(TXT_FILE)//'_COPY.TXT'
    TEXT_FILE_COPY = COPY_FILE(TXT_FILE,TXT_FILE_COPY)
    RETURN
  END FUNCTION TEXT_FILE_COPY
  !-------------------------------------------------------------------
END MODULE GWM_UTLS

