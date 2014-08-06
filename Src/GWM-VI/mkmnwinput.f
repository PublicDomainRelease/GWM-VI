!
!***********************************************************************
      SUBROUTINE MKMNWINPUT(MNW2JTF_FILE)
!***********************************************************************
!     Open the MNW2 input file and call GWM1DCV3MNW2 which will read
!       the MNW2 input file, close it, and write a Jupiter Template File
!       named mnw2_input.jtf which is a template for an MNW2 input file
!       that contains both managed and unmanaged MNW2 wells.
!-----------------------------------------------------------------------
      USE GWM_SUBS, ONLY: IGETUNIT
      USE GWM1DCV3MNW, ONLY: GWM1DCV3MNW2, NAUX_EXTERNAL
      USE GWM1DCV3, ONLY: MNW2HERE
      USE MKMMPROC_RDNAM, ONLY: MNW2FILE
      IMPLICIT NONE
      CHARACTER(LEN=2000),INTENT(IN)::MNW2JTF_FILE
      ! Local variables
      LOGICAL FIRSTSIM
      INTEGER MNW2UNIT,LOCAT
      CHARACTER(LEN=2000) :: FNAME
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF(.NOT.MNW2HERE)RETURN  ! NO MNW-type wells in flow variables
        IF(MNW2FILE.EQ.'')THEN   ! NO MNW file
          CALL USTOP('MNW-type in DECVAR but No MNW2 input file found')
        ENDIF
        ! Open the original MNW2 input file on an arbitrary unit number
        MNW2UNIT = IGETUNIT(7,1000)  
        FNAME = MNW2FILE 
        LOCAT = MNW2UNIT
        OPEN(UNIT=MNW2UNIT,FILE=FNAME,STATUS='OLD',ERR=999)
        READ(MNW2UNIT,'(A)')FNAME  ! Read first line of file into FNAME
        REWIND (MNW2UNIT)          ! Rewind file
        IF(TRIM(FNAME).EQ.'# MNW2 input file written by GWM')THEN
          ! User must place original MNW2 file content into this file 
          CALL USTOP('Combined MNW found instead of original MNW2 file')
        ENDIF
        ! Read the original MNW2 input file and write mnw2_input.jtf
        FIRSTSIM = .TRUE.
        ! GWM-VI uses an auxiliary variable in MNW2 input to identify wells        
        NAUX_EXTERNAL = 1
        ! The first four zero arguments define the solver used
        !   which is used to determine SMALL within MNW2.  However,
        !   SMALL is irrelevant for reading the data.
        CALL GWM1DCV3MNW2(MNW2UNIT,FIRSTSIM,0,0,0,0,0,MNW2JTF_FILE)
        ! Upon return, MNW2UNIT is the unit number for mnw2_input.jtf
        CLOSE(UNIT=MNW2UNIT)
      RETURN
  999 CONTINUE ! FILE-OPENING ERROR
      WRITE(*,9990)TRIM(FNAME),LOCAT
 9990 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (MKMNWINPUT)')
      CALL USTOP(' ')
      END SUBROUTINE MKMNWINPUT 
