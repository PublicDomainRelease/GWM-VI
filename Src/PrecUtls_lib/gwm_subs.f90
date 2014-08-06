MODULE GWM_SUBS
  ! This module contains some utility subprograms that do
  ! not require modules of the JUPITER API.
  PRIVATE
  ! Public data
  PUBLIC :: EPSQNET, MNWCOPYOK, MNW2FILE, MNW2FILEORIG
  ! Public program units
  PUBLIC :: BASE36_I, CLEAN_UP, FUZZY_EQUALS, I_BASE36, IGETUNIT, &
            TO_LONGNAME, URWORD2
  !
  INTERFACE URWORD2
    MODULE PROCEDURE URWORD2_D
  END INTERFACE
  !
  CHARACTER(LEN=2000) :: MNW2FILE, MNW2FILEORIG
  LOGICAL :: MNWCOPYOK
  DOUBLEPRECISION :: EPSQNET = 1.0D-5
CONTAINS
  !-----------------------------------------------------------------------------
  INTEGER FUNCTION BASE36_I(CH2) RESULT(IRETURN)
    !   Convert 2-character string representing a base-36 integer
    !   to a decimal integer.  "Digits" in CH2 may be any of 
    !   10 digits 0-9 and 26 letters A-Z (lowercase letters also
    !   are acceptable).  Error is indicated by a return value < 0.
    IMPLICIT NONE
    !   Argument-list variable
    CHARACTER(LEN=2), INTENT(IN) :: CH2
    !   Local variables
    INTEGER :: I, J, KFOUND
    INTEGER, DIMENSION(2) :: IA
    CHARACTER(LEN=1), DIMENSION(0:35) :: DIGIT36 
    CHARACTER(LEN=1) :: CH1
    DATA DIGIT36/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E', &
        'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V', &
        'W','X','Y','Z'/
    !
    IRETURN = -1
    IA = -1
    KFOUND = 0
    EACHCHAR: DO I=1,2
      CH1 = CH2(I:I)
      EACHDIG: DO J=0,35
        IF (SAMENAME(DIGIT36(J),CH1)) THEN
          IA(I) = J
          KFOUND = KFOUND+1
          EXIT EACHDIG
        ENDIF
      ENDDO EACHDIG
    ENDDO EACHCHAR
    !
    IF (KFOUND==2) THEN
      IRETURN = 0
      DO I=1,2
        IRETURN = IRETURN+IA(I)*36**(I-1)
      ENDDO
    ENDIF
    RETURN
  END FUNCTION BASE36_I
  !-----------------------------------------------------------------------------
  SUBROUTINE CHGCASE(WORDIN,WORDOUT,ICASE)
    !   CONVERT A CHARACTER STRING TO ALL UPPER (ICASE > 0) OR ALL
    !   LOWER (ICASE < 0) CASE.  IF ICASE = 0, NO CASE CONVERSION IS DONE.
    IMPLICIT NONE
    !
    !   Argument-list variables
    CHARACTER(LEN=*),  INTENT(IN)  :: WORDIN
    CHARACTER(LEN=*),  INTENT(OUT) :: WORDOUT
    INTEGER, OPTIONAL, INTENT(IN)  :: ICASE   ! Default = 1
    !
    !   Local variables
    INTEGER :: ICAS, IDIFF, K, LENGNB, LENGOUT
    !
    ICAS = 1
    IF (PRESENT(ICASE)) ICAS = ICASE
    !   DETERMINE LENGTH OF STRING VARIABLES
    LENGOUT = LEN(WORDOUT)

    !   DETERMINE IF WORDOUT IS LONG ENOUGH TO CONTAIN NON-BLANK LENGTH
    !   OF WORDIN
    LENGNB = LEN_TRIM(WORDIN)
    IF (LENGNB.GT.LENGOUT) THEN
      WORDOUT = 'STRING-LENGTH ERROR IN CHGCASE'
      RETURN
    ENDIF
    !
    WORDOUT = WORDIN
    IDIFF=ICHAR('a')-ICHAR('A')
    IF (ICAS.GT.0) THEN
    !     CONVERT STRING TO UPPER CASE
      DO 10 K=1,LENGNB
        IF(WORDIN(K:K).GE.'a' .AND. WORDIN(K:K).LE.'z')   &
            WORDOUT(K:K)=CHAR(ICHAR(WORDIN(K:K))-IDIFF)
10    CONTINUE
    !
    ELSEIF (ICAS.LT.0) THEN
    !     CONVERT STRING TO LOWER CASE
      DO 20 K=1,LENGNB
        IF(WORDIN(K:K).GE.'A' .AND. WORDIN(K:K).LE.'Z')   &
            WORDOUT(K:K)=CHAR(ICHAR(WORDIN(K:K))+IDIFF)
20    CONTINUE
    !
    ENDIF
    !
    RETURN
  END SUBROUTINE CHGCASE
  !-------------------------------------------------------------------
  SUBROUTINE CLEAN_UP()
    IMPLICIT NONE
    ! Restore original MNW2 input file
    IF (MNWCOPYOK) THEN
      MNWCOPYOK = COPY_FILE(MNW2FILEORIG,MNW2FILE)
    ENDIF
    RETURN
  END SUBROUTINE CLEAN_UP
  !-------------------------------------------------------------------
  LOGICAL FUNCTION COPY_FILE(FROMFILE,TOFILE)
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
    IUFROM = IGETUNIT(7,9000)
    OKLOCAL = OPEN_FILE(IUFROM,FROMFILE,'OLD')
    IF (OKLOCAL) THEN
      IUTO = IGETUNIT(7,9000)
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
    ENDIF
    COPY_FILE = OKLOCAL
    RETURN
  END FUNCTION COPY_FILE
  !-----------------------------------------------------------------------------
  LOGICAL FUNCTION FUZZY_EQUALS(A,B,EPS) RESULT(C)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: A, B, EPS
    DOUBLE PRECISION :: DZERO = 0.0D0, EPSLOCAL
    !
    EPSLOCAL = ABS(EPS)
    IF (B .NE. DZERO) THEN
      IF (ABS((A-B)/B) .LT. EPSLOCAL) THEN
        C = .TRUE.
      ELSE
        C = .FALSE.
      ENDIF
    ELSE 
      ! B==0
      IF (ABS(A) .LT. EPSLOCAL) THEN
        C = .TRUE.
      ELSE
        C = .FALSE.
      ENDIF
    ENDIF
    !
    RETURN
  END FUNCTION FUZZY_EQUALS
  !-----------------------------------------------------------------------------
  FUNCTION I_BASE36(I) RESULT(CH2)
    !   Return a LEN=2 character string based on integer KSP
    !   The string may contain any of 10 digits 0-9 and 26 letters A-Z
    !   The result is a two-digit, base-36 representation of integer KSP
    !   The largest integer that can be represented this way
    !   is 36*36-1 = 1295.  I may be in range 0-1295.  If I is outside
    !   this range, the function returns '**'
    IMPLICIT NONE
    !   Argument-list and result variables
    INTEGER, INTENT(IN) :: I
    CHARACTER(LEN=2) :: CH2
    !   Local variables
    INTEGER :: I1, I2
    CHARACTER(LEN=1), DIMENSION(0:35) :: DIGIT36 
    DATA DIGIT36/'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E', &
        'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V', &
        'W','X','Y','Z'/
    !
    CH2 = '**'   ! Indicates error
    IF (I.GE.0 .AND. I.LE.1295) THEN
      I1 = I/36             ! Find the 36's digit, in decimal 
      I2 = I-I1*36          ! Find the 1's digit, in decimal
      CH2(1:1) = DIGIT36(I1)  ! Assign 36's digit, in base-36
      CH2(2:2) = DIGIT36(I2)  ! Assign 1's digit, in base-36
    ENDIF
    RETURN
  END FUNCTION I_BASE36
  !-----------------------------------------------------------------------------
  FUNCTION IGETUNIT(IFIRST,MAXUNIT) RESULT (KUNIT)
    !   Find first unused file unit number between IFIRST and MAXUNIT
    IMPLICIT NONE
    !
    !   Argument-list and result variables
    INTEGER, INTENT(IN) :: IFIRST
    INTEGER, INTENT(IN) :: MAXUNIT
    INTEGER             :: KUNIT
    !
    !   Local variables
    INTEGER I, IOST
    LOGICAL LOP
    !
    !   Loop through range provided to find first unused unit number
    DO I=IFIRST,MAXUNIT
      INQUIRE(UNIT=I,IOSTAT=IOST,OPENED=LOP,ERR=5)
      IF (IOST==0) THEN
        IF (.NOT. LOP) THEN
          KUNIT = I
          RETURN
        ENDIF
      ENDIF
 5    CONTINUE
    ENDDO
    !
    !   If there are no unused unit numbers in range provided, return
    !   a value indicating an error
    KUNIT = -1
    !
    RETURN
  END FUNCTION IGETUNIT
  !-------------------------------------------------------------------
  LOGICAL FUNCTION OPEN_FILE(IU,FNAME,STATUS,BINARY)
    ! Open file FILENAME on unit IU
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
    CALL CHGCASE(STATUS,STATUSLOCAL,1)
    IF (STATUSLOCAL == 'OLD') THEN
      INQUIRE(FILE=FNAME,EXIST=LEX)
      IF (.NOT. LEX) THEN
!        CALL APPEND_MESSAGE('Error: File "' // TRIM(FNAME) // &
!                            '" does not exist.')
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
        !CALL PLL_WAIT()
        GOTO 10
      ELSE
!        CALL APPEND_MESSAGE('Error opening file "' // TRIM(FNAME) &
!                            // '"')
      ENDIF
    ENDIF
    RETURN
  END FUNCTION OPEN_FILE
  !-----------------------------------------------------------------------------
  FUNCTION SAMENAME(NAME1,NAME2) RESULT(SAME)
    !   Make a case-insensitive comparison of two names
    IMPLICIT NONE
    !
    !   Argument-list and result variables
    CHARACTER(LEN=*), INTENT(IN) :: NAME1 ! Name to be compared
    CHARACTER(LEN=*), INTENT(IN) :: NAME2 ! Name to be compared
    LOGICAL :: SAME
    !
    !   Local variables
    CHARACTER(LEN=2000) :: NAME1UP, NAME2UP
    !
    !   Convert both names to upper case, and make comparison
    CALL CHGCASE(NAME1,NAME1UP,1)
    CALL CHGCASE(NAME2,NAME2UP,1)
    IF (NAME1UP==NAME2UP) THEN
      SAME = .TRUE.
    ELSE
      SAME = .FALSE.
    ENDIF
    !
    RETURN
  END FUNCTION SAMENAME
  !-----------------------------------------------------------------------------
  SUBROUTINE TO_LONGNAME(SHORTNAME,LONGNAME,NUM)
    !   Define LONGNAME by concatenating base-36 representation
    !   of NUM to SHORTNAME, padded with _ characters between
    IMPLICIT NONE
    !   Argument-list variables
    CHARACTER(LEN=*), INTENT(IN)  :: SHORTNAME
    CHARACTER(LEN=*), INTENT(OUT) :: LONGNAME
    INTEGER,          INTENT(IN)  :: NUM
    !   Local variables
    INTEGER :: LENL, LENS, LENSP1, LENSTRIM, LENOUT
    CHARACTER(LEN=40) :: TEMPNAM, UNDERSCORES
    CHARACTER(LEN=2)  :: CH2
    DATA UNDERSCORES/'________________________________________'/
    !
    LENS = LEN(SHORTNAME)
    LENL = LEN(LONGNAME)
    LENSP1 = LENS+1
    LENOUT = LENS+2
    CH2 = I_BASE36(NUM)
    !
    !   Error checks
    IF (LENOUT .GT. 40) CALL USTOP('SHORTNAME TOO LONG IN TO_LONGNAME')
    IF (LENL .LT. LENOUT) CALL USTOP('PROGRAMMER ERROR IN TO_LONGNAME')
    IF (CH2 == '**') CALL USTOP('ILLEGAL VALUE OF NUM IN TO_LONGNAME')
    !
    LENSTRIM = LEN_TRIM(SHORTNAME)
    TEMPNAM = UNDERSCORES
    TEMPNAM(1:LENSTRIM) = TRIM(SHORTNAME)
    TEMPNAM(LENSP1:LENOUT) = CH2
    TEMPNAM(LENOUT+1:LENL) = ' '
    LONGNAME = TEMPNAM
    RETURN
  END SUBROUTINE TO_LONGNAME
  !-----------------------------------------------------------------------------
  SUBROUTINE URWORD2_D(LINE,ICOL,ISTART,ISTOP,NCODE,N,R,IOUT,IN)
!C     ******************************************************************
!C     ROUTINE TO EXTRACT A WORD FROM A LINE OF TEXT, AND OPTIONALLY
!C     CONVERT THE WORD TO A NUMBER.
!C        ISTART AND ISTOP WILL BE RETURNED WITH THE STARTING AND
!C          ENDING CHARACTER POSITIONS OF THE WORD.
!C        THE LAST CHARACTER IN THE LINE IS SET TO BLANK SO THAT IF ANY
!C          PROBLEMS OCCUR WITH FINDING A WORD, ISTART AND ISTOP WILL
!C          POINT TO THIS BLANK CHARACTER.  THUS, A WORD WILL ALWAYS BE
!C          RETURNED UNLESS THERE IS A NUMERIC CONVERSION ERROR.  BE SURE
!C          THAT THE LAST CHARACTER IN LINE IS NOT AN IMPORTANT CHARACTER
!C          BECAUSE IT WILL ALWAYS BE SET TO BLANK.
!C        A WORD STARTS WITH THE FIRST CHARACTER THAT IS NOT A SPACE OR
!C          COMMA, AND ENDS WHEN A SUBSEQUENT CHARACTER THAT IS A SPACE
!C          OR COMMA.  NOTE THAT THESE PARSING RULES DO NOT TREAT TWO
!C          COMMAS SEPARATED BY ONE OR MORE SPACES AS A NULL WORD.
!C        FOR A WORD THAT BEGINS WITH "'", THE WORD STARTS WITH THE
!C          CHARACTER AFTER THE QUOTE AND ENDS WITH THE CHARACTER
!C          PRECEDING A SUBSEQUENT QUOTE.  THUS, A QUOTED WORD CAN
!C          INCLUDE SPACES AND COMMAS.  THE QUOTED WORD CANNOT CONTAIN
!C          A QUOTE CHARACTER.
!C        IF NCODE IS 1, THE WORD IS CONVERTED TO UPPER CASE.
!C        IF NCODE IS 2, THE WORD IS CONVERTED TO AN INTEGER.
!C        IF NCODE IS 3, THE WORD IS CONVERTED TO A REAL NUMBER.
!C        NUMBER CONVERSION ERROR IS WRITTEN TO UNIT IOUT IF IOUT IS
!C          POSITIVE; ERROR IS WRITTEN TO DEFAULT OUTPUT IF IOUT IS 0;
!C          NO ERROR MESSAGE IS WRITTEN IF IOUT IS NEGATIVE.
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      CHARACTER*(*) LINE
      CHARACTER*20 STRING
      CHARACTER*30 RW
      CHARACTER*1 TAB
      DOUBLE PRECISION :: R
!C     ------------------------------------------------------------------
      TAB=CHAR(9)
!C
!C1------Set last char in LINE to blank and set ISTART and ISTOP to point
!C1------to this blank as a default situation when no word is found.  If
!C1------starting location in LINE is out of bounds, do not look for a
!C1------word.
      LINLEN=LEN(LINE)
      LINE(LINLEN:LINLEN)=' '
      ISTART=LINLEN
      ISTOP=LINLEN
      LINLEN=LINLEN-1
      IF(ICOL.LT.1 .OR. ICOL.GT.LINLEN) GO TO 100
!C
!C2------Find start of word, which is indicated by first character that
!C2------is not a blank, a comma, or a tab.
      DO 10 I=ICOL,LINLEN
      IF(LINE(I:I).NE.' ' .AND. LINE(I:I).NE.',' &
          .AND. LINE(I:I).NE.TAB) GO TO 20
10    CONTINUE
      ICOL=LINLEN+1
      GO TO 100
!C
!C3------Found start of word.  Look for end.
!C3A-----When word is quoted, only a quote can terminate it.
20    IF(LINE(I:I).EQ.'''') THEN
         I=I+1
         IF(I.LE.LINLEN) THEN
            DO 25 J=I,LINLEN
            IF(LINE(J:J).EQ.'''') GO TO 40
25          CONTINUE
         END IF
!C
!C3B-----When word is not quoted, space, comma, or tab will terminate.
      ELSE
         DO 30 J=I,LINLEN
         IF(LINE(J:J).EQ.' ' .OR. LINE(J:J).EQ.',' &
          .OR. LINE(J:J).EQ.TAB) GO TO 40
30       CONTINUE
      END IF
!C
!C3C-----End of line without finding end of word; set end of word to
!C3C-----end of line.
      J=LINLEN+1
!C
!C4------Found end of word; set J to point to last character in WORD and
!C-------set ICOL to point to location for scanning for another word.
40    ICOL=J+1
      J=J-1
      IF(J.LT.I) GO TO 100
      ISTART=I
      ISTOP=J
!C
!C5------Convert word to upper case and RETURN if NCODE is 1.
      IF(NCODE.EQ.1) THEN
         IDIFF=ICHAR('a')-ICHAR('A')
         DO 50 K=ISTART,ISTOP
            IF(LINE(K:K).GE.'a' .AND. LINE(K:K).LE.'z') &
                   LINE(K:K)=CHAR(ICHAR(LINE(K:K))-IDIFF)
50       CONTINUE
         RETURN
      END IF
!C
!C6------Convert word to a number if requested.
100   IF(NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
         RW=' '
         L=30-ISTOP+ISTART
         IF(L.LT.1) GO TO 200
         RW(L:30)=LINE(ISTART:ISTOP)
         IF(NCODE.EQ.2) READ(RW,'(I30)',ERR=200) N
         IF(NCODE.EQ.3) READ(RW,'(F30.0)',ERR=200) R
      END IF
      RETURN
!C
!C7------Number conversion error.
200   IF(NCODE.EQ.3) THEN
         STRING= 'A REAL NUMBER'
         L=13
      ELSE
         STRING= 'AN INTEGER'
         L=10
      END IF
!C
!C7A-----If output unit is negative, set last character of string to 'E'.
      IF(IOUT.LT.0) THEN
         N=0
         R=0.
         LINE(LINLEN+1:LINLEN+1)='E'
         RETURN
!C
!C7B-----If output unit is positive; Write a message to output unit.
      ELSE IF(IOUT.GT.0) THEN
         IF(IN.GT.0) THEN
            WRITE(IOUT,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
         ELSE
            WRITE(IOUT,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
         END IF
201      FORMAT(1X,/1X,'FILE UNIT ',I4,' : ERROR CONVERTING "',A, &
             '" TO ',A,' IN LINE:',/1X,A)
202      FORMAT(1X,/1X,'KEYBOARD INPUT : ERROR CONVERTING "',A, &
             '" TO ',A,' IN LINE:',/1X,A)
!C
!C7C-----If output unit is 0; Write a message to default output.
      ELSE
         IF(IN.GT.0) THEN
            WRITE(*,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
         ELSE
            WRITE(*,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
         END IF
      END IF
!C
!C7D-----STOP after writing message.
      CALL USTOP(' ')
      END SUBROUTINE URWORD2_D
  !-------------------------------------------------------------------
END MODULE GWM_SUBS
