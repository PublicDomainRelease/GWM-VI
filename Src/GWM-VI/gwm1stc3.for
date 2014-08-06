      MODULE GWM1STC3
C     VERSION: 27AUG2013 
      USE GWM_SUBS, ONLY: IGETUNIT
      USE GWM_STOP, ONLY:   GSTOP
      USE MF2005_UTLS, ONLY: URDCOM, URWORD

      IMPLICIT NONE
      PRIVATE
      PUBLIC::STCNUM,STCSTATE,STCSTATE0,NSF,NSD,NLK,STCNAME,
     &        STCSLOC,STCRLOC
      PUBLIC::GWM1STC3AR,GWM1STC3FP,GWM1STC3FPR,GWM1STC3FM,
     &        GWM1STC3OT,GWM1STC3OS,GWM1STC3OS2
C-----MAKE PUBLIC FOR GWM-VI
      PUBLIC ::STCSP
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
C
C-----VARIABLES FOR STREAM CONSTRAINTS
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::STCNAME(:)
      REAL(DP),         SAVE,ALLOCATABLE::STCSTATE(:),STCSTATE0(:)
      REAL(DP),         SAVE,ALLOCATABLE::STCRHS(:)
      INTEGER(I4B),     SAVE,ALLOCATABLE::STCSLOC(:),STCRLOC(:)
      INTEGER(I4B),     SAVE,ALLOCATABLE::STCDIR(:),STCSP(:)
      INTEGER(I4B),     SAVE,ALLOCATABLE::GRDLOCSTC(:)
      INTEGER(I4B),     SAVE            ::STCNUM,NSF,NSD,NLK
C
C      STCNAME  -name of constraint 
C      STCSTATE -array of constraint left hand side values resulting from most 
C                 recent MODFLOW simulation
C      STCSTATE0-array used to store base constraint values during pertubation
C                 for response matrix generation
C      STCRHS   -value of the right hand side for this constraint 
C      STCSLOC  -stream segment at which constraint is applied
C      STCRLOC  -stream reach at which constraint is applied
C      STCDIR   -flag to determine direction of inequality for this constraint
C                 1 => left hand side < right hand side
C                 2 => left hand side > right hand side
C      STCSP    -stress period during which this constraint is active 
C      NSF      -number of streamflow constraints
C      NSD      -number of stream depletion constraints
C      NLK      -number of stream leakage constraints
C      STCNUM   -total number of stream constraints
C      GRDLOCSTC-grid number on which the constraint is located
C
C-----FOR ERROR HANDLING
      INTEGER(I2B)::ISTAT  
      CHARACTER(LEN=200)::FLNM
      CHARACTER(LEN=20)::FILACT,FMTARG,ACCARG 
      INTEGER(I4B)::NDUM
      REAL(SP)::RDUM
C
      CONTAINS
C***********************************************************************
      SUBROUTINE GWM1STC3AR(FNAMEN,IOUT,NV,NCON,NRMC,NGRIDS)
C***********************************************************************
C   PURPOSE: READ INPUT FROM THE STREAM CONSTRAINTS FILE
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : CUTCOM
      USE GWM1RMS3, ONLY : IREF
      USE GLOBAL  , ONLY : IUNIT
      CHARACTER(LEN=200),INTENT(IN),DIMENSION(NGRIDS)::FNAMEN
      INTEGER(I4B),INTENT(INOUT)::NCON,NRMC,NV
      INTEGER(I4B),INTENT(IN)::IOUT,NGRIDS
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,II,ITYPES,ITYPEF,LLOC,INMS,INMF,BYTES
      INTEGER(I4B)::LOCAT,ISTART,ISTOP,IPRN,IPRNG,NSP
      INTEGER(I4B)::NSEG,NRCH
      REAL(SP)::TBND
      CHARACTER(LEN=10)::TCTNM
      CHARACTER(LEN=2)::FG2
      CHARACTER(LEN=200)::LINE
      INTEGER(I4B),DIMENSION(NGRIDS)::NUNOPN
      INTEGER(I4B),DIMENSION(NGRIDS)::TNSF,TNSD,TNLK,TNS
      INTEGER(I4B)::G,JROW,NS
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      TNSF = 0
      TNSD = 0
      TNLK = 0 
      IPRN = -1
C
C-----LOOP OVER ALL GRIDS TO COUNT CONSTRAINTS
C
      DO 100 G=1,NGRIDS
        IF(FNAMEN(G).NE.' ')THEN
C---------CHECK THAT SFR AND STR PACKAGES ARE NOT BOTH ACTIVE ON THIS GRID           
          CALL SGWF2BAS7PNT(G)                   ! POINT TO CORRECT IUNIT
          IF(IUNIT(18).GT.0 .AND. IUNIT(44).GT.0)THEN ! STR AND SFR ACTIVE
            WRITE(IOUT,2000,ERR=990)G            ! GWM DOES NOT ALLOW THIS
            CALL GSTOP(' ')                      ! SO STOP THE RUN
          ENDIF
C---------OPEN FILE
          NUNOPN(G)=IGETUNIT(10,200)
          LOCAT=NUNOPN(G)
          WRITE(IOUT,1000,ERR=990) LOCAT,FNAMEN(G)
          OPEN(UNIT=LOCAT,FILE=FNAMEN(G),ACTION='READ',ERR=999)
C
C---------CHECK FOR COMMENT LINES
          CALL URDCOM(LOCAT,IOUT,LINE)
C
C---------READ IPRN
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPRNG,RDUM,IOUT,LOCAT)
          IF(IPRNG.EQ.0 .OR. IPRNG.EQ.1)THEN
            IPRN = MAX(IPRN,IPRNG)               ! USE MOST DETAILED ECHO
          ELSE
            WRITE(IOUT,2010,ERR=990)IPRNG        ! INVALID VALUE OF IPRN
            CALL GSTOP(' ')
          ENDIF
C
C---------READ TOTAL NUMBER OF STREAM-CONSTRAINT TYPES
          READ(LOCAT,'(A)',ERR=991)LINE
          CALL CUTCOM(LINE,200)                  ! REMOVE COMMENTS FROM LINE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,TNSF(G),RDUM,IOUT,LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,TNSD(G),RDUM,IOUT,LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,TNLK(G),RDUM,IOUT,LOCAT)
          IF(TNSF(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NSF IS NEGATIVE')
          IF(TNSD(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NSD IS NEGATIVE')
          IF(TNLK(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NLK IS NEGATIVE')
          WRITE(IOUT,3000,ERR=990)TNSF(G),TNSD(G)
          WRITE(IOUT,3010,ERR=990)TNLK(G)
        ENDIF
  100 ENDDO 
C
C-----ADD CONSTRAINTS FROM ALL GRIDS
      NSF=SUM(TNSF)
      NSD=SUM(TNSD)
      NLK=SUM(TNLK)
      IF(NSD.GT.0)IREF=1                         ! REFERENCE SIMULATION NEEDED
      STCNUM= NSF + NSD + NLK
      NCON = NCON + STCNUM                       ! INCREMENT CONSTRAINT NUMBER
      NRMC = NRMC + STCNUM                       ! INCREMENT RESPONSE CONSTRAINTS
      NV   = NV   + STCNUM                       ! ADD CONSTRAINT SLACKS TO NV
C
      IF(STCNUM.GT.0)THEN
        ALLOCATE (STCNAME(STCNUM),               ! ALLOCATE SPACE FOR STREAM 
     &            STCSTATE(STCNUM),STCSTATE0(STCNUM),STCRHS(STCNUM),
     &            STCSLOC(STCNUM),STCRLOC(STCNUM),  ! CONSTRAINT INFORMATION
     &            STCDIR(STCNUM),STCSP(STCNUM),GRDLOCSTC(STCNUM), 
     &            STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992
        BYTES = 10*STCNUM + 8*3*STCNUM + 4*5*STCNUM
      ELSE
        RETURN                                   ! NO STREAM CONSTRAINTS 
      ENDIF
C
C-----READ CONSTRAINT INFORMATION
C
      JROW=0                                     ! INITIALIZE MASTER ROW COUNTER
      DO 230 NS =1,3
      IF(NS.EQ.1) TNS = TNSF                     ! READ STREAMFLOW CONSTRAINTS
      IF(NS.EQ.2) TNS = TNSD                     ! READ STREAM DEPLETION CONST.
      IF(NS.EQ.3) TNS = TNLK                     ! READ STREAM LEAKAGE CONST.
C-----LOOP OVER ALL GRIDS TO READ EACH STREAMFLOW CONSTRAINT SET
      DO 220 G=1,NGRIDS
        IF(TNS(G).GT.0)THEN     
          LOCAT=NUNOPN(G) 
          DO 210 I=1+JROW,TNS(G)+JROW
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSEG,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRCH,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ITYPES,ITYPEF,1,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,NDUM,TBND,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSP,RDUM,IOUT,LOCAT)
C-----------PROCESS THE CONSTRAINT NAME
            TCTNM=LINE(INMS:INMF)
C-----------CHECK THAT CONSTRAINT NAME HAS NOT BEEN USED
            DO 200 II=1,I-1
              IF(STCNAME(II).EQ.TCTNM)THEN
                WRITE(IOUT,4000,ERR=990)TCTNM
                CALL GSTOP(' ')
              ENDIF
  200       CONTINUE
            STCNAME(I)=TCTNM                     ! STORE CONSTRAINT NAME
            GRDLOCSTC(I)=G                       ! STORE GRID NUMBER
C-----------PROCESS SEGMENT AND REACH NUMBER 
            IF(NSEG.LT.0)CALL GSTOP('PROGRAM STOPPED: NSEG IS NEGATIVE')
            STCSLOC(I)=NSEG                      ! STORE SEGMENT NUMBER 
            IF(NRCH.LT.0)CALL GSTOP('PROGRAM STOPPED: NRCH IS NEGATIVE')
            STCRLOC(I)=NRCH                      ! STORE REACH NUMBER 
C-----------PROCESS THE TYPE OF CONSTRAINT
            FG2=LINE(ITYPES:ITYPEF)
            IF(FG2.EQ.'LE') THEN
              STCDIR(I) = 1
            ELSEIF(FG2.EQ.'GE') THEN
              STCDIR(I) = 2
            ELSE
              CALL GSTOP('STREAMFLOW CONSTRAINT NOT LE OR GE  ')
            ENDIF
C-----------PROCESS BOUND AND STRESS PERIOD TO WHICH CONSTRAINT APPLIES 
            STCRHS(I)=REAL(TBND,DP)              ! STORE CONSTRAINT BOUND 
            STCSP(I)=NSP                         ! STORE STRESS PERIOD 
  210     ENDDO
        ENDIF
        JROW=JROW+TNS(G)
  220 ENDDO
  230 ENDDO
C
C-----WRITE THE INFORMATION TO THE OUTPUT FILE
      IF(IPRN.EQ.1)THEN
        IF(NSF.GT.0)THEN
          WRITE(IOUT,5000,ERR=990)
          DO 300 I=1,NSF
            IF(STCDIR(I).EQ.1)THEN
               WRITE(IOUT,5010,ERR=990)I,STCNAME(I),STCSLOC(I),
     1                                STCRLOC(I),'<',STCRHS(I),STCSP(I)
            ENDIF
            IF(STCDIR(I).EQ.2)THEN
              WRITE(IOUT,5010,ERR=990)I,STCNAME(I),STCSLOC(I),
     1                                STCRLOC(I),'>',STCRHS(I),STCSP(I)
            ENDIF
  300     ENDDO
        ENDIF
        IF(NSD.GT.0)THEN
          WRITE(IOUT,5020,ERR=990)
          DO 400 I=NSF+1,NSF+NSD
            IF(STCDIR(I).EQ.1)THEN
              WRITE(IOUT,5010,ERR=990)I,STCNAME(I),STCSLOC(I),
     1                                STCRLOC(I),'<',STCRHS(I),STCSP(I)
            ENDIF
            IF(STCDIR(I).EQ.2)THEN
              WRITE(IOUT,5010,ERR=990)I,STCNAME(I),STCSLOC(I),
     1                                STCRLOC(I),'>',STCRHS(I),STCSP(I)
            ENDIF
  400     ENDDO
        ENDIF
        IF(NLK.GT.0)THEN
          WRITE(IOUT,5040,ERR=990)
          DO 500 I=NSF+NSD+1,NSF+NSD+NLK
            IF(STCDIR(I).EQ.1)THEN
              WRITE(IOUT,5010,ERR=990)I,STCNAME(I),STCSLOC(I),
     1                                STCRLOC(I),'<',STCRHS(I),STCSP(I)
            ENDIF
            IF(STCDIR(I).EQ.2)THEN
              WRITE(IOUT,5010,ERR=990)I,STCNAME(I),STCSLOC(I),
     1                                STCRLOC(I),'>',STCRHS(I),STCSP(I)
            ENDIF
  500     ENDDO
        ENDIF
      ENDIF
C      
 600  WRITE(IOUT,6000,ERR=990)BYTES
      WRITE(IOUT,7000,ERR=990)
C
C-----CLOSE FILES
      DO I=1,NGRIDS
        CLOSE(UNIT = NUNOPN(I))
      ENDDO
C
 1000 FORMAT(1X,/1X,'OPENING STREAMFLOW CONSTRAINTS FILE',/,
     1  ' ON UNIT ',I4,':',/1X,A200)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. BOTH SFR AND ',
     1  'STR PACKAGE ACTIVE ON GRID',I4)
 2010 FORMAT(1X,/1X,'PROGRAM STOPPED. IPRN MUST BE EQUAL TO 0 OR 1: ',
     1  I4)
 3000 FORMAT(/,1X,'NUMBER OF STREAMFLOW (NSF) AND STREAMFLOW-',
     1  'DEPLETION',/,' CONSTRAINTS (NSD) ARE ',I5,' AND ',I5,' ,',
     2  ' RESPECTIVELY.')
 3010 FORMAT(/,1X,'NUMBER OF STREAM LEAKAGE (NLK) CONSTRAINTS IS',I5)
 4000 FORMAT(1X,/1X,'PROGRAM STOPPED. CONSTRAINT NAME ',A10,' HAS',
     1  ' ALREADY BEEN USED.') 
 5000 FORMAT(1X,/1X,'STREAMFLOW CONSTRAINTS:',/,T43,'RIGHT-HAND',
     1  T56,'STRESS',/,T2,'NUMBER',T10,'NAME',T20,'SEGMENT',T29,
     2  'REACH',T36,'TYPE',T46,'SIDE',T56,'PERIOD',/,' --------',
     3  '-------------------------------------------------------')
 5010 FORMAT(1X,I5,T11,A10,T21,I3,T30,I4,T37,A1,T42,ES11.4,T56,I5)
 5020 FORMAT(1X,/1X,'STREAMFLOW-DEPLETION CONSTRAINTS:',/,T43,
     1  'RIGHT-HAND',T56,'STRESS',/,T2,'NUMBER',T10,'NAME',T20,
     2  'SEGMENT',T29,'REACH',T36,'TYPE',T46,'SIDE',T56,'PERIOD',
     3  /,' ---------------------------------------------------',
     4  '------------')
 5040 FORMAT(1X,/1X,'STREAM LEAKAGE CONSTRAINTS:',/,T43,
     1  'RIGHT-HAND',T56,'STRESS',/,T2,'NUMBER',T10,'NAME',T20,
     2  'SEGMENT',T29,'REACH',T36,'TYPE',T46,'SIDE',T56,'PERIOD',
     3  /,' ---------------------------------------------------',
     4  '------------')
 6000 FORMAT(/,1X,I8,' BYTES OF MEMORY ALLOCATED TO STORE DATA FOR',
     1               ' STREAMFLOW CONSTRAINTS')
 7000 FORMAT(/,1X,'CLOSING STREAMFLOW CONSTRAINTS FILE',/)
C
      RETURN
C

C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),IOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1STC3AR)')
      CALL GSTOP(' ')
C
  991 CONTINUE
C-----FILE-READING ERROR
      INQUIRE(LOCAT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9910)TRIM(FLNM),LOCAT,FMTARG,ACCARG,FILACT
      WRITE(IOUT,9910)TRIM(FLNM),LOCAT,FMTARG,ACCARG,FILACT
 9910 FORMAT(/,1X,'*** ERROR READING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1STC3AR)')
      CALL GSTOP(' ')
C
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1STC3AR)')
      CALL GSTOP(' ')
C
  999 CONTINUE
C-----FILE-OPENING ERROR
      INQUIRE(LOCAT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9990)TRIM(FLNM),LOCAT,'OLD',FMTARG,ACCARG,FILACT
      WRITE(IOUT,9990)TRIM(FLNM),LOCAT,'OLD',FMTARG,ACCARG,
     &                 FILACT
 9990 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE STATUS: ',A,/
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1STC3AR)')
      CALL GSTOP(' ')
C
      END SUBROUTINE GWM1STC3AR
C
C
C***********************************************************************
      SUBROUTINE GWM1STC3OS(IGRID,KPER,IC1,NDEP,DEPVALS)
C***********************************************************************
C  VERSION: 16JULY2009
C  PURPOSE: ASSIGN COMPUTED STREAMFLOW VALUES TO STATE ARRAY FROM STR PACKAGE
C-----------------------------------------------------------------------
      USE GWFSTRMODULE, ONLY: NSTREM, STRM, ISTRM
      INTEGER(I4B),INTENT(IN)::IGRID,KPER,IC1,NDEP
      DOUBLE PRECISION, DIMENSION(NDEP), INTENT(IN), OPTIONAL :: DEPVALS
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,ISEG,IREC,J,LOC,IC2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(NDEP.EQ.0)THEN                          ! CALL NOT FROM GWM-VI
C
      DO 100 I= 1,NSF+NSD+NLK 
        IF(IGRID.EQ.GRDLOCSTC(I))THEN            ! CONSTRAINT IS ON THIS GRID
          IF(KPER.EQ.STCSP(I)) THEN              ! CONSTRAINT ACTIVE THIS SP
            CALL SGWF2STR7PNT(IGRID)             ! POINT TO CORRECT GRID
            ISEG = STCSLOC(I)                    ! STREAM SEGMENT FOR CONSTRAINT
            IREC = STCRLOC(I)                    ! STREAM REACH FOR CONSTRAINT
            DO 80 J = 1,NSTREM
              IF(ISTRM(4,J).EQ.ISEG .AND. ISTRM(5,J).EQ.IREC)THEN
                LOC = J                          ! THIS CONSTRAINT AT LOC
                GOTO 90
              ENDIF
   80       ENDDO
            GOTO 100
   90       CONTINUE
            IF(I.LE.NSF+NSD)THEN                 ! THIS IS A STREAMFLOW CONSTR.
C             STCSTATE(I) = REAL(STRM(9,LOC),DP) ! IF STRM IS SINGLE PRECISION NEED THIS
              STCSTATE(I) = STRM(9,LOC)          ! LOAD STREAMFLOW VALUE
            ELSE                                 ! THIS IS A STREAM LEAK CONSTR.
              STCSTATE(I) = STRM(11,LOC)         ! LOAD STREAM LEAKAGE VALUE
            ENDIF
          ENDIF
        ENDIF
  100 ENDDO
C
      ELSEIF(NDEP.GT.0 .AND. PRESENT(DEPVALS))THEN ! CALL IS FROM GWM-VI
        IC2 = IC1+NSF+NSD+NLK-1
        J = 0
        DO 200 I=IC1,IC2 
          J = J+1
          STCSTATE(J) = DEPVALS(I)
  200   ENDDO
C
      ENDIF
C
      RETURN
      END SUBROUTINE GWM1STC3OS
C
C***********************************************************************
      SUBROUTINE GWM1STC3OS2(IGRID,KPER,IC1,NDEP,DEPVALS)
C***********************************************************************
C  VERSION: 16JULY2009
C  PURPOSE: ASSIGN COMPUTED STREAMFLOW VALUES TO STATE ARRAY FROM SFR PACKAGE
C-----------------------------------------------------------------------
      USE GWFSFRMODULE, ONLY: NSTRM, STRM, ISTRM
      INTEGER(I4B),INTENT(IN)::IGRID,KPER,IC1,NDEP
      DOUBLE PRECISION, DIMENSION(NDEP), INTENT(IN), OPTIONAL :: DEPVALS
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,ISEG,IREC,J,LOC,IC2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(NDEP.EQ.0)THEN                          ! CALL NOT FROM GWM-VI
C
      DO 100 I= 1,NSF+NSD+NLK 
        IF(IGRID.EQ.GRDLOCSTC(I))THEN            ! CONSTRAINT IS ON THIS GRID
          IF(KPER.EQ.STCSP(I)) THEN              ! CONSTRAINT ACTIVE THIS SP
            CALL SGWF2SFR7PNT(IGRID)             ! POINT TO CORRECT GRID
            ISEG = STCSLOC(I)                    ! STREAM SEGMENT FOR CONSTRAINT
            IREC = STCRLOC(I)                    ! STREAM REACH FOR CONSTRAINT
            DO 80 J = 1,NSTRM
              IF(ISTRM(4,J).EQ.ISEG .AND. ISTRM(5,J).EQ.IREC)THEN
                LOC = J                          ! THIS CONSTRAINT AT LOC
                GOTO 90
              ENDIF
   80       ENDDO
            GOTO 100
   90       CONTINUE
            IF(I.LE.NSF+NSD)THEN                 ! THIS IS A STREAMFLOW CONSTR.
C             STCSTATE(I) = REAL(STRM(9,LOC),DP) ! IF STRM IS SINGLE PRECISION NEED THIS
              STCSTATE(I) = STRM(9,LOC)          ! LOAD STREAMFLOW VALUE
            ELSE                                 ! THIS IS A STREAM LEAK CONSTR.
              STCSTATE(I) = STRM(11,LOC)         ! LOAD STREAM LEAKAGE VALUE
            ENDIF
          ENDIF
        ENDIF
  100 ENDDO
C
      ELSEIF(NDEP.GT.0 .AND. PRESENT(DEPVALS))THEN ! CALL IS FROM GWM-VI
        IC2 = IC1+NSF+NSD+NLK-1
        J = 0
        DO 200 I=IC1,IC2 
          J = J+1
          STCSTATE(J) = DEPVALS(I)               
  200   ENDDO
C
      ENDIF
C
      RETURN
      END SUBROUTINE GWM1STC3OS2
C
C***********************************************************************
      SUBROUTINE GWM1STC3FP(RSTRT,IPERT)
C***********************************************************************
C  VERSION: 19MAR2008
C  PURPOSE - USE SIMULATION RESULTS TO COMPUTE RESPONSE MATRIX AND 
C            AUGMENTED RIGHT HAND SIDE
C-----------------------------------------------------------------------
      USE GWM1OBJ3, ONLY : SOLNTYP
      USE GWM1DCV3, ONLY : FVBASE
      USE GWM1BAS3, ONLY : RMFILE
      USE GWM1RMS3, ONLY : SLPITCNT,IBASE,IREF,IRM,DELINC,
     &                     RHSIN,RHSINF,RANGENAME,RANGENAMEF,CONTYP,NDV
C-----AMAT HAS LOCAL NAME RESMAT
      USE GWM1RMS3, ONLY : RHS,RESMAT=>AMAT
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::IPERT
C-----LOCAL VARIABLES
      REAL(DP)::STATEN
      INTEGER(I4B)::ROW,ISTC,NPGNA,I
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IPERT.EQ.-1)THEN
C-------THIS IS A REFERENCE SIMULATION
        CALL SGWM1STC3FP                         ! MODIFY CONSTRAINTS
      ELSEIF(IPERT.EQ.0)THEN
C-------THIS IS A BASE SIMULATION
        IF(SLPITCNT.GT.0)THEN
C---------THIS IS A BASE SIMULATION FOR A LATER SLP ITERATION
        ELSEIF(.NOT.(IREF.EQ.1.AND.IBASE.EQ.1))THEN 
C---------THIS IS EITHER THE FIRST SLP OR BASE LP SIMULATION    
          CALL SGWM1STC3FP                       ! MODIFY CONSTRAINTS                   
        ENDIF
C
        ISTC = 0
        DO 100 ROW=RSTRT,RSTRT+STCNUM-1
          ISTC = ISTC+1
C
C---------RETRIEVE THE INITIAL STATE FROM MODFLOW OUTPUT
          STCSTATE0(ISTC) = STCSTATE(ISTC)
C
C---------IF SAVING RESPONSE MATRIX WRITE THE INITIAL STATE TO FILE
          IF(IRM.EQ.1 .OR. IRM.EQ.4)WRITE(RMFILE)STCSTATE0(ISTC)
C
C---------COMPUTE THE DIFFERENCE BETWEEN INITIAL STATE AND RHS
          RHS(ROW) = STCRHS(ISTC) - STCSTATE0(ISTC)
100     ENDDO
C
      ELSEIF(IPERT.GT.0)THEN
C-----  THIS IS A PERTURBATION SIMULATION
C
C-------COMPUTE THE RESPONSE MATRIX VALUE FOR EACH CONSTRAINT
        ISTC = 0
        DO 200 ROW=RSTRT,RSTRT+STCNUM-1
          ISTC = ISTC+1
C
C---------COMPUTE STATE DIFFERENCE (PERTURBED STATE MINUS INITIAL STATE)
          STATEN = STCSTATE(ISTC) - STCSTATE0(ISTC)
C
C---------DIVIDE STATE DIFFERENCE BY THE PERTURBATION SIZE 
          RESMAT(ROW,IPERT) = STATEN/DELINC(IPERT)
C
C---------AUGMENT RIGHT HAND SIDE WITH KNOWN PERTURBATION RESPONSE TERM 
          RHS(ROW) = RHS(ROW) + RESMAT(ROW,IPERT)*FVBASE(IPERT)
C
C---------WRITE RESPONSE MATRIX ELEMENT TO DISK
          IF(IRM.EQ.1 .OR. IRM.EQ.4) WRITE(RMFILE)RESMAT(ROW,IPERT)
  200   ENDDO
C
      ENDIF
C
C-----SET NEXT STARTING LOCATION
      RSTRT = RSTRT+STCNUM
C
      RETURN
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1STC3FP
C***********************************************************************
C
C  PURPOSE: USE RESULTS OF THE FIRST SIMULATION TO MODIFY CONSTRAINTS
C-----------------------------------------------------------------------
C
C-----STORE THE INITIAL RHS FOR USE IN OUTPUT 
      ISTC = 0
      DO 100 ROW=RSTRT,RSTRT+STCNUM-1
        ISTC = ISTC+1
        RHSIN(ROW)= STCRHS(ISTC)                  ! STORE THE INPUT RHS
        RANGENAME(ROW+NDV-1)=STCNAME(ISTC)        ! LOAD FOR RANGE ANALYSIS OUTPUT
        RHSINF(ROW)= STCRHS(ISTC)                 ! STORE THE INPUT RHS
        RANGENAMEF(ROW+NDV-1)=STCNAME(ISTC)       ! LOAD FOR RANGE ANALYSIS OUTPUT
        IF(ISTC.LE.NSF)THEN                       ! STREAMFLOW CONSTRAINT
          CONTYP(ROW)=1                           ! SET FOR RANGE OUTPUT
        ELSE                                      ! STREAM DEPLETION CONSTRAINT
          CONTYP(ROW)=2                           ! THIS IS TRANSFORMED CONSTRAINT
        ENDIF
  100 ENDDO
C
C-----FOR STREAM CONSTRAINTS, CONVERT STREAM DEPLETION TO STREAM FLOW
      DO 200 I=NSF+1,NSF+NSD
        IF(IRM.EQ.1 .OR. IRM.EQ.4)WRITE(RMFILE)STCSTATE(I)! WRITE REFERENCE STATE
        STCRHS(I) = STCSTATE(I) - STCRHS(I)       ! STCRHS IS NOW MINIMUM STREAMFLOW
C------SWITCH DIRECTION OF INEQUALITY: UPPER BOUND ON DEPLETION
C      BECOMES LOWER BOUND ON STREAMFLOW
        IF(STCDIR(I).EQ.1)THEN
          STCDIR(I) = 2
        ELSEIF(STCDIR(I).EQ.2)THEN
          STCDIR(I) = 1
        ENDIF
 200  ENDDO
C
      RETURN
      END SUBROUTINE SGWM1STC3FP
C
      END SUBROUTINE GWM1STC3FP
C
C
C***********************************************************************
      SUBROUTINE GWM1STC3FPR(RSTRT,IREAD,COL)
C***********************************************************************
C  VERSION: 19MAR2008
C  PURPOSE - READ RESPONSE MATRIX AND AUGMENTED RIGHT HAND SIDE
C-----------------------------------------------------------------------
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,FVBASE
      USE GWM1BAS3, ONLY : RMFILE
      USE GWM1RMS3, ONLY : RHSIN,RHSINF,RANGENAME,RANGENAMEF,CONTYP,NDV
C-----AMAT HAS LOCAL NAME RESMAT 
      USE GWM1RMS3, ONLY : RESMAT => AMAT
      USE GWM1RMS3, ONLY : RHS
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::IREAD
      INTEGER(I4B),INTENT(IN)::COL
C-----LOCAL VARIABLES
      INTEGER(I4B)::ISTC,ROW,I
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IREAD.EQ.0)THEN    
        CALL SGWM1STC3FPR                        ! READ THE REFERENCE STATE                     
        ISTC = 0
        DO 100 ROW=RSTRT,RSTRT+STCNUM-1 
          ISTC = ISTC+1
          READ(RMFILE)STCSTATE0(ISTC)            ! READ THE BASE STATE FROM A FILE
          RHS(ROW) = STCRHS(ISTC) - STCSTATE0(ISTC)
          STCSTATE(ISTC) = STCSTATE0(ISTC)       ! ASSIGN BASE STATE FOR OUTPUT
  100   ENDDO
      ELSEIF(IREAD.EQ.1)THEN                     ! READ RESPONSE MATRIX
        DO 200 ROW=RSTRT,RSTRT+STCNUM-1
          READ(RMFILE)RESMAT(ROW,COL)            ! READ EACH RESPONSE COEFFICIENT
          RHS(ROW) = RHS(ROW) + RESMAT(ROW,COL)*FVBASE(COL)
  200   ENDDO
      ENDIF
      RSTRT = RSTRT+STCNUM                       ! SET NEXT STARTING LOCATION
C
      RETURN
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1STC3FPR
C***********************************************************************
C
C  PURPOSE: READ RESULTS OF THE FIRST SIMULATION AND MODIFY CONSTRAINTS
C-----------------------------------------------------------------------
C
C-----STORE THE INITIAL RHS FOR USE IN OUTPUT 
      ISTC = 0
      DO 100 ROW=RSTRT,RSTRT+STCNUM-1
        ISTC = ISTC+1
        RHSIN(ROW)= STCRHS(ISTC)                  ! STORE THE INPUT RHS
        RANGENAME(ROW+NDV-1)=STCNAME(ISTC)        ! LOAD FOR RANGE ANALYSIS OUTPUT
        RHSINF(ROW)= STCRHS(ISTC)                 ! STORE THE INPUT RHS
        RANGENAMEF(ROW+NDV-1)=STCNAME(ISTC)       ! LOAD FOR RANGE ANALYSIS OUTPUT
        IF(ISTC.LE.NSF)THEN                       ! STREAMFLOW CONSTRAINT
          CONTYP(ROW)=1                           ! SET FOR RANGE OUTPUT
        ELSE                                      ! STREAM DEPLETION CONSTRAINT
          CONTYP(ROW)=2                           ! THIS IS TRANSFORMED CONSTRAINT
        ENDIF
  100 ENDDO
C
C-----TRANSFORM STREAM DEPLETION CONSTRAINTS
      DO 200 I=NSF+1,NSF+NSD
        READ(RMFILE)STCSTATE(I)                   ! READ THE REFERENCE STATE 
        STCRHS(I) = STCSTATE(I) - STCRHS(I)       ! STCRHS IS NOW MINIMUM HEAD
        IF(STCDIR(I).EQ.1)THEN                    ! SWITCH DIRECTION OF INEQUALITY
          STCDIR(I) = 2
        ELSEIF(STCDIR(I).EQ.2)THEN
          STCDIR(I) = 1
        ENDIF
  200 ENDDO
C
      RETURN
      END SUBROUTINE SGWM1STC3FPR
C
      END SUBROUTINE GWM1STC3FPR
C
C
C***********************************************************************
      SUBROUTINE GWM1STC3FM(RSTRT,NSLK)
C***********************************************************************
C  VERSION: 19MAR2008
C  PURPOSE - PLACE SLACK COEFFICIENTS IN LP MATRIX STARTING IN ROW RSTRT
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE
      USE GWM1RMS3, ONLY : AMAT
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::NSLK
C-----LOCAL VARIABLES
      INTEGER(I4B) ROW,ISTC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----SET STREAM FLOW CONSTRAINT SLACKS
      ISTC = 0
      DO 100 ROW=RSTRT,RSTRT+STCNUM-1 
        ISTC = ISTC+1
        IF(STCDIR(ISTC).EQ.1) THEN
          AMAT(ROW,NSLK) =  ONE
        ELSE
          AMAT(ROW,NSLK) = -ONE
        ENDIF
  100 ENDDO
C
C-----SET NEXT STARTING LOCATION
      RSTRT = RSTRT+STCNUM
C
      RETURN
      END SUBROUTINE GWM1STC3FM
C
C
C***********************************************************************
      SUBROUTINE GWM1STC3OT(RSTRT,IFLG)
C*************************************************************************
C  VERSION: 25JAN2011
C  PURPOSE - WRITE STATUS OF CONSTRAINT
C---------------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,ONE
      USE GWM1RMS3, ONLY : RHS,RANGENAME,NDV,HCLOSEG
      USE GWM1BAS3, ONLY : GWM1BAS3CS
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::IFLG
C-----LOCAL VARIABLES
      CHARACTER(LEN=25)::CTYPE
      INTEGER ROW,ISTC,DIRR
      DOUBLE PRECISION DIFF
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----WRITE STATUS FOR A FLOW PROCESS SIMULATION
      IF(IFLG.EQ.1.OR.IFLG.EQ.4)THEN            
        DO 100 ISTC=1,STCNUM
          CALL SGWM1STC3OT
          DIFF = STCSTATE(ISTC) - STCRHS(ISTC)
          DIRR = STCDIR(ISTC)           
          IF(IFLG.EQ.1)THEN
            CALL GWM1BAS3CS(CTYPE,STCNAME(ISTC),0.0,0.0,DIFF,DIRR,1)
          ELSEIF(IFLG.EQ.4)THEN
            CALL GWM1BAS3CS(CTYPE,STCNAME(ISTC),
     &                      STCSTATE(ISTC),STCRHS(ISTC),DIFF,DIRR,2)
          ENDIF
100     ENDDO
C
C-----WRITE STATUS FOR LINEAR PROGRAM OUTPUT
      ELSEIF(IFLG.EQ.2)THEN
        ISTC = 0
        DO 200 ROW=RSTRT,RSTRT+STCNUM-1
        ISTC = ISTC + 1
          IF(ABS(RHS(ROW)).GT.ZERO)THEN          ! DUAL VARIABLE IS NON-ZERO 
                                                 ! CONSTRAINT IS BINDING         
            CALL SGWM1STC3OT
            CALL GWM1BAS3CS(CTYPE,STCNAME(ISTC),0.0,0.0,RHS(ROW),0,0) 
          ENDIF
  200   ENDDO
        RSTRT = RSTRT+STCNUM                     ! SET NEXT STARTING LOCATION
      ENDIF
C
      RETURN
      CONTAINS
C*************************************************************************
      SUBROUTINE SGWM1STC3OT
C*************************************************************************
C
C  PURPOSE - ASSIGN CONSTRAINT TYPE
C---------------------------------------------------------------------------
C
      IF(ISTC.LE.NSF)THEN
        CTYPE = 'Streamflow'
      ELSEIF(ISTC.LE.NSF+NSD)THEN
        CTYPE = 'Stream Depletion'
      ELSEIF(ISTC.LE.NSF+NSD+NLK)THEN
        CTYPE = 'Stream Leakage'
      ENDIF
C
      RETURN
      END SUBROUTINE SGWM1STC3OT
C
      END SUBROUTINE GWM1STC3OT
C
C
      END MODULE GWM1STC3