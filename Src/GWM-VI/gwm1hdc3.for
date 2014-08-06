      MODULE GWM1HDC3
C     VERSION: 27AUG2013
      USE GWM_SUBS, ONLY: FUZZY_EQUALS, IGETUNIT
      USE GWM_STOP, ONLY:   GSTOP
      USE MF2005_UTLS, ONLY: URDCOM, URWORD
      IMPLICIT NONE
      PRIVATE
      PUBLIC::HDCNAME,HDCSTATE,HDCSTATE0,HDCRHS,HDCDIR,HDCNUM,
     &        HDCILOC,HDCJLOC,HDCKLOC,HDCSP,NHB,NDD,NDF,NGD
      PUBLIC::GWM1HDC3AR,GWM1HDC3FP,GWM1HDC3FPR,GWM1HDC3FM,
     &        GWM1HDC3OT,GWM1HDC3OS
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
C
C-----VARIABLES FOR HEAD CONSTRAINTS
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::HDCNAME(:)
      REAL(DP),         SAVE,ALLOCATABLE::HDCSTATE(:),HDCSTATE0(:)
      REAL(DP),         SAVE,ALLOCATABLE::HDCRHS(:)
      INTEGER(I4B),     SAVE,ALLOCATABLE::HDCILOC(:,:),HDCJLOC(:,:),
     &                                    HDCKLOC(:,:)
      INTEGER(I4B),     SAVE,ALLOCATABLE::HDCDIR(:),HDCSP(:)
      INTEGER(I4B),     SAVE,ALLOCATABLE::GRDLOCHDC(:)
      INTEGER(I4B),     SAVE            ::HDCNUM,NHB,NDD,NDF,NGD
C
C      HDCNAME  -name of constraint 
C      HDCSTATE -array of constraint left hand side values resulting from most 
C                 recent MODFLOW simulation
C      HDCSTATE0-array used to store base constraint values during pertubation
C                 for response matrix generation
C      HDCRHS   -value of the right hand side for this constraint 
C      HDCILOC  -constraint observation point, first index
C      HDCJLOC  -constraint observation point, second index
C      HDCKLOC  -constraint observation point, third index
C                 HDCILOC, HDCJLOC, and HDCKLOC contain two entries for each 
C                 constraint.  
C                 If a head bound or drawdown constraint then the first entries are the 
C                 row, column,layer location in MODFLOW grid of the head 
C                 observation point and the second entry is ignored.  
C                 If a head difference or gradient constraint then the first and second
C                 entries are the row, column,layer location in MODFLOW grid of constraint 
C                 the first and second head observation points.  
C      HDCDIR   -flag to determine direction of inequality for this constraint
C                 1 => left hand side < right hand side
C                 2 => left hand side > right hand side
C      HDCSP    -stress period during which this constraint is active 
C      NHB      -number of head bound constraints
C      NDD      -number of drawdown constraints
C      NDF      -number of head difference constraints
C      NGD      -number of gradient constraints
C      HDCNUM   -total number of head type constraints
C      GRDLOCHDC-grid number on which the constraint is located
C
C-----FOR ERROR HANDLING
      INTEGER(I2B)::ISTAT  
      CHARACTER(LEN=200)::FLNM
      CHARACTER(LEN=20)::FILACT,FMTARG,ACCARG 
      INTEGER(I4B)::NDUM
      REAL(SP)::RDUM
C
      CONTAINS
C*************************************************************************
      SUBROUTINE GWM1HDC3AR(FNAMEN,IOUT,NPER,NV,NCON,NRMC,NGRIDS) 
C*************************************************************************
C     VERSION: 19MAR2008
C     PURPOSE: READ INPUT FROM THE HEAD CONSTRAINTS FILE 
C---------------------------------------------------------------------------
      USE GWM1BAS3, ONLY : BIGINF
      USE GWM1RMS3, ONLY : IREF
      USE GWM1DCV3, ONLY : NFVAR,FVNCELL,FVMNW,GRDLOCDCV
      USE GLOBAL,   ONLY : NCOL,NROW,NLAY
      INTEGER(I4B),INTENT(INOUT)::NCON,NV,NRMC
      INTEGER(I4B),INTENT(IN)::IOUT,NPER,NGRIDS
      CHARACTER(LEN=200),INTENT(IN),DIMENSION(NGRIDS)::FNAMEN
C-----LOCAL VARIABLES
      INTEGER(I4B)::NSP,MNWCNT
      CHARACTER(LEN=10)::TCTNM
      CHARACTER(LEN=20)::TCMNW
      REAL(SP)::TBND,TLEN,TGRAD
      CHARACTER(LEN=2)::FG2
      CHARACTER(LEN=1)::DIRPRT
      INTEGER(I4B)::I,J,K,II,ITYPES,ITYPEF,IPRN,IR1,IC1,IL1,IR2,IC2,IL2
      INTEGER(I4B)::G,JROW,NH,IPRNG
      INTEGER(I4B)::BYTES,LOCAT,ISTART,ISTOP,LLOC,INMS,INMF
      CHARACTER(LEN=200)::LINE
      REAL(DP), ALLOCATABLE::LEN(:),GRAD(:)
      INTEGER(I4B),DIMENSION(NGRIDS)::NUNOPN
      INTEGER(I4B),DIMENSION(NGRIDS)::TNHD,TNHB,TNDD,TNDF,TNGD
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      TNHB = 0
      TNDD = 0
      TNDF = 0
      TNGD = 0
      IPRN = -1
C
C-----LOOP OVER ALL GRIDS TO COUNT CONSTRAINTS
C
      DO 100 G=1,NGRIDS
        IF(FNAMEN(G).NE.' ')THEN
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
            WRITE(IOUT,2000,ERR=990)IPRNG        ! INVALID VALUE OF IPRN
            CALL GSTOP(' ')
          ENDIF
C
C---------READ TOTAL NUMBER OF HEAD-TYPE CONSTRAINTS
          READ(LOCAT,*,ERR=991)TNHB(G),TNDD(G),TNDF(G),TNGD(G)
          IF(TNHB(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NHB IS NEGATIVE')
          IF(TNDD(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NDD IS NEGATIVE')
          IF(TNDF(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NDF IS NEGATIVE')
          IF(TNGD(G).LT.0)CALL GSTOP('PROGRAM STOPPED: NGD IS NEGATIVE')
        ENDIF
  100 ENDDO   
C-----ADD CONSTRAINTS FROM ALL GRIDS
      NHB=SUM(TNHB)
      NDD=SUM(TNDD)
      NDF=SUM(TNDF)
      NGD=SUM(TNGD)
      IF(NDD.GT.0)IREF=1                         ! REFERENCE SIMULATION NEEDED
      HDCNUM = NHB + NDD + NDF + NGD
      NCON = NCON + HDCNUM                       ! INCREMENT CONSTRAINT NUMBER
      NRMC = NRMC + HDCNUM                       ! INCREMENT RESPONSE CONSTRAINTS
      NV   = NV   + HDCNUM                       ! ADD CONSTRAINT SLACKS TO NV
C
      IF(HDCNUM.GT.0)THEN
C-----ALLOCATE SPACE FOR CONSTRAINT INFORMATION
        ALLOCATE (HDCNAME(HDCNUM),               ! ALLOCATE SPACE FOR HEAD 
     1            HDCSTATE(HDCNUM),HDCSTATE0(HDCNUM),HDCRHS(HDCNUM),
     2            HDCILOC(HDCNUM,2),HDCJLOC(HDCNUM,2),HDCKLOC(HDCNUM,2),
     3            HDCDIR(HDCNUM),HDCSP(HDCNUM),      ! CONSTRAINTS
     4            LEN(HDCNUM),GRAD(HDCNUM),GRDLOCHDC(HDCNUM),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992
        BYTES = 10*HDCNUM + 8*3*HDCNUM + 4*9*HDCNUM  
      ELSE
        RETURN                                   ! NO HEAD CONSTRAINTS 
      ENDIF
C
C-----READ CONSTRAINT INFORMATION
C
      JROW=0                                     ! INITIALIZE MASTER ROW COUNTER
      DO 230 NH =1,2
      IF(NH.EQ.1) TNHD = TNHB                    ! READ HEAD BOUND CONSTRAINTS
      IF(NH.EQ.2) TNHD = TNDD                    ! READ DRAWDOWN CONSTRAINTS
C-----LOOP OVER ALL GRIDS TO READ HEAD-BOUND AND DRAWDOWN CONSTRAINTS
      DO 220 G=1,NGRIDS
        IF(TNHD(G).GT.0)THEN      
          LOCAT=NUNOPN(G)  
          DO 210 I=1+JROW,TNHD(G)+JROW
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            TCTNM=LINE(INMS:INMF) ! SAVE THE CONSTRAINT NAME
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL1,RDUM,IOUT,LOCAT)
            IF(IL1.GT.0)THEN ! THIS IS A CELL CONSTRAINT; READ ROW,COL
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR1,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC1,RDUM,IOUT,LOCAT)
            ELSE             ! THIS IS AN MNW WELL CONSTRAINT;READ WELLID
              CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            ENDIF
            CALL URWORD(LINE,LLOC,ITYPES,ITYPEF,1,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,NDUM,TBND,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSP,RDUM,IOUT,LOCAT)
            
C
C-----------PROCESS THE CONSTRAINT NAME; CHECK THAT IT HAS NOT BEEN USED
            DO 200 II=1,I-1
              IF(HDCNAME(II).EQ.TCTNM)THEN
                WRITE(IOUT,3000,ERR=990)TCTNM
                CALL GSTOP(' ')
              ENDIF
  200       ENDDO
            HDCNAME(I)=TCTNM                     ! STORE CONSTRAINT NAME
            GRDLOCHDC(I)=G                       ! STORE GRID NUMBER
C-----------PROCESS ROW, COLUMN AND LAYER NUMBER OR MNW WELLID
            CALL SGWF2BAS7PNT(G)                 ! CHANGE POINTERS TO THIS GRID
            IF(IL1.EQ.0)THEN  ! THIS IS AN MNW WELL CONSTRAINT
              TCMNW=LINE(INMS:INMF) ! THIS IS THE WELLID NAME
              CALL UPCASE(TCMNW)    ! convert to uppercase
              MNWCNT = 0            ! COUNT THE NUMBER OF MNW-TYPE VARIABLES
              IC1=-1
              LOOP_FV : DO J=1,NFVAR      ! LOOP OVER GWM FLOW VARIABLES
                IF(FVNCELL(J).LT.0)THEN   ! THIS IS MNW-TYPE FLOW VARIABLE
                  DO K=1,ABS(FVNCELL(J))  ! LOOP OVER CELLS FOR THIS VARIABLE
                    MNWCNT = MNWCNT + 1
                    IF(FVMNW(J)%FVWELLID(K).EQ.TCMNW)THEN ! MATCH FOUND 
                      IF(G.NE.GRDLOCDCV(J))THEN ! FLOW VARIABLE NOT ON GRID
                        WRITE(IOUT,3010,ERR=990)TRIM(TCTNM),G  
                        CALL GSTOP(' ')
                      ENDIF
                      IR1 = MNWCNT ! STORE THE INDEX IN IR1
                      IC1 = IR1    ! THIS INDEX GETS CHANGED IN GWM1DCV3ARW
                      EXIT LOOP_FV
                    ENDIF               
                  ENDDO
                ENDIF
             ENDDO LOOP_FV
             IF(IC1.EQ.-1)WRITE(IOUT,3020,ERR=990)TRIM(TCTNM)  
             IF(IC1.EQ.-1)CALL GSTOP(' ')
            ELSE
              IF(IR1.LT.1 .OR. IR1.GT.NROW)THEN
                WRITE(IOUT,3100,ERR=990)IR1      ! NOT A VALID ROW NUMBER 
                CALL GSTOP(' ')
              ENDIF
              IF(IC1.LT.1 .OR. IC1.GT.NCOL)THEN
                WRITE(IOUT,3200,ERR=990)IC1      ! NOT A VALID COLUMN NUMBER 
                CALL GSTOP(' ')
              ENDIF
              IF(IL1.LT.1 .OR. IL1.GT.NLAY)THEN
                WRITE(IOUT,3300,ERR=990)IL1      ! NOT A VALID LAYER NUMBER 
                CALL GSTOP(' ')
              ENDIF
            ENDIF
            HDCILOC(I,1)=IR1                     ! STORE ROW NUMBER 
            HDCJLOC(I,1)=IC1                     ! STORE COLUMN NUMBER 
            HDCKLOC(I,1)=IL1                     ! STORE LAYER NUMBER 
C
C-----------PROCESS THE TYPE OF CONSTRAINT
            FG2=LINE(ITYPES:ITYPEF)
            IF(FG2.EQ.'LE') THEN
              HDCDIR(I) = 1
            ELSEIF(FG2.EQ.'GE') THEN
              HDCDIR(I) = 2
            ELSE
              CALL GSTOP('HEAD CONSTRAINT NOT LE OR GE  ')
            ENDIF
C
C-----------PROCESS BOUND AND STRESS PERIOD TO WHICH CONSTRAINT APPLIES 
            HDCRHS(I)=REAL(TBND,DP)                ! STORE RIGHT HAND SIDE
            IF(HDCRHS(I).GT.BIGINF)THEN
              WRITE(IOUT,3400,ERR=990)TCTNM,BIGINF ! RHS VALUE NOT WITHIN BOUNDS
              CALL GSTOP(' ')
            ENDIF
            HDCSP(I)=NSP                           ! STORE STRESS PERIOD
  210     ENDDO
C
C-------END OF HEAD-BOUNDS LOOP
        ENDIF
        JROW=JROW+TNHD(G)
  220 ENDDO
  230 ENDDO
C
C-----READ DIFFERENCE CONSTRAINTS
      DO 320 G=1,NGRIDS
        IF(TNDF(G).GT.0)THEN
          LOCAT=NUNOPN(G)  
          DO 310 I=1+JROW,TNDF(G)+JROW
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL1,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR1,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC1,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL2,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR2,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC2,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,NDUM,TBND,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSP,RDUM,IOUT,LOCAT)
C
C-----------PROCESS THE CONSTRAINT NAME
            TCTNM=LINE(INMS:INMF)
C-----------CHECK THAT CONSTRAINT NAME HAS NOT BEEN USED
            DO 300 II=1,I-1
              IF(HDCNAME(II).EQ.TCTNM)THEN
                WRITE(IOUT,3000,ERR=990)TCTNM
                CALL GSTOP(' ')
              ENDIF
  300       ENDDO
            HDCNAME(I)=TCTNM                     ! STORE CONSTRAINT NAME
            GRDLOCHDC(I)=G                       ! STORE GRID NUMBER
C-----------PROCESS ROW, COLUMN AND LAYER NUMBER 
            CALL SGWF2BAS7PNT(G)                 ! CHANGE POINTERS TO THIS GRID
            IF(IR1.LT.1 .OR. IR1.GT.NROW)THEN
              WRITE(IOUT,3100,ERR=990)IR1        ! NOT A VALID ROW NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IC1.LT.1 .OR. IC1.GT.NCOL)THEN
              WRITE(IOUT,3200,ERR=990)IC1        ! NOT A VALID COLUMN NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IL1.LT.1 .OR. IL1.GT.NLAY)THEN
              WRITE(IOUT,3300,ERR=990)IL1        ! NOT A VALID LAYER NUMBER
              CALL GSTOP(' ')
            ENDIF
            IF(IR2.LT.1 .OR. IR2.GT.NROW)THEN
              WRITE(IOUT,3100,ERR=990)IR2        ! NOT A VALID ROW NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IC2.LT.1 .OR. IC2.GT.NCOL)THEN
              WRITE(IOUT,3200,ERR=990)IC2        ! NOT A VALID COLUMN NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IL2.LT.1 .OR. IL2.GT.NLAY)THEN
              WRITE(IOUT,3300,ERR=990)IL2        ! NOT A VALID LAYER NUMBER
              CALL GSTOP(' ')
            ENDIF
            HDCILOC(I,1)=IR1                     ! STORE ROW NUMBER 
            HDCJLOC(I,1)=IC1                     ! STORE COLUMN NUMBER 
            HDCKLOC(I,1)=IL1                     ! STORE LAYER NUMBER 
            HDCILOC(I,2)=IR2                     ! STORE ROW NUMBER 
            HDCJLOC(I,2)=IC2                     ! STORE COLUMN NUMBER 
            HDCKLOC(I,2)=IL2                     ! STORE LAYER NUMBER 
C
C-----------PROCESS BOUND AND STRESS PERIOD TO WHICH CONSTRAINT APPLIES 
            HDCRHS(I)=REAL(TBND,DP)              ! STORE RIGHT HAND SIDE
            IF(HDCRHS(I).GT.BIGINF)THEN
              WRITE(IOUT,3400,ERR=990)TCTNM,BIGINF ! RHS VALUE NOT WITHIN BOUNDS
              CALL GSTOP(' ')
            ENDIF
            HDCSP(I)=NSP                         ! STORE STRESS PERIOD
            HDCDIR(I)=2                          ! ONLY GE TYPE USED
  310     ENDDO
        ENDIF
        JROW=JROW+TNDF(G)
  320 ENDDO
C
C-----READ GRADIENT CONSTRAINTS
      DO 420 G=1,NGRIDS
        IF(TNGD(G).GT.0)THEN
          LOCAT=NUNOPN(G)  
          DO 410 I=1+JROW,TNGD(G)+JROW
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL1,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR1,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC1,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL2,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR2,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC2,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,NDUM,TLEN,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,NDUM,TGRAD,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSP,RDUM,IOUT,LOCAT)
C 
C-----------PROCESS THE CONSTRAINT NAME
            TCTNM=LINE(INMS:INMF)
C-----------CHECK THAT CONSTRAINT NAME HAS NOT BEEN USED
            DO 400 II=1,I-1
              IF(HDCNAME(II).EQ.TCTNM)THEN
                WRITE(IOUT,3000,ERR=990)TCTNM
                CALL GSTOP(' ')
              ENDIF
  400       ENDDO
            HDCNAME(I)=TCTNM                     ! STORE CONSTRAINT NAME
            GRDLOCHDC(I)=G                       ! STORE GRID NUMBER
C-----------PROCESS ROW, COLUMN AND LAYER NUMBER 
            CALL SGWF2BAS7PNT(G)                 ! CHANGE POINTERS TO THIS GRID
            IF(IR1.LT.1 .OR. IR1.GT.NROW)THEN
              WRITE(IOUT,3100,ERR=990)IR1        ! NOT A VALID ROW NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IC1.LT.1 .OR. IC1.GT.NCOL)THEN
              WRITE(IOUT,3200,ERR=990)IC1        ! NOT A VALID COLUMN NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IL1.LT.1 .OR. IL1.GT.NLAY)THEN
              WRITE(IOUT,3300,ERR=990)IL1        ! NOT A VALID LAYER NUMBER
              CALL GSTOP(' ')
            ENDIF
            IF(IR2.LT.1 .OR. IR2.GT.NROW)THEN
              WRITE(IOUT,3100,ERR=990)IR2        ! NOT A VALID ROW NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IC2.LT.1 .OR. IC2.GT.NCOL)THEN
              WRITE(IOUT,3200,ERR=990)IC2        ! NOT A VALID COLUMN NUMBER 
              CALL GSTOP(' ')
            ENDIF
            IF(IL2.LT.1 .OR. IL2.GT.NLAY)THEN
              WRITE(IOUT,3300,ERR=990)IL2        ! NOT A VALID LAYER NUMBER
              CALL GSTOP(' ')
            ENDIF
            HDCILOC(I,1)=IR1                     ! STORE ROW NUMBER 
            HDCJLOC(I,1)=IC1                     ! STORE COLUMN NUMBER 
            HDCKLOC(I,1)=IL1                     ! STORE LAYER NUMBER 
            HDCILOC(I,2)=IR2                     ! STORE ROW NUMBER 
            HDCJLOC(I,2)=IC2                     ! STORE COLUMN NUMBER 
            HDCKLOC(I,2)=IL2                     ! STORE LAYER NUMBER 
C
C-----------PROCESS BOUND AND STRESS PERIOD TO WHICH CONSTRAINT APPLIES 
            LEN(I) =REAL(TLEN,DP)                ! LOCAL STORAGE FOR OUTPUT
            GRAD(I)=REAL(TGRAD,DP)               ! LOCAL STORAGE FOR OUTPUT
            TBND=TLEN*TGRAD                      ! CONVERT LEN * GRAD TO RHS
            HDCRHS(I)=REAL(TBND,DP)              ! STORE RIGHT HAND SIDE
            IF(HDCRHS(I).GT.BIGINF)THEN
              WRITE(IOUT,3400,ERR=990)TCTNM,BIGINF ! RHS VALUE NOT WITHIN BOUNDS
              CALL GSTOP(' ')
            ENDIF
            HDCSP(I)=NSP                           ! STORE STRESS PERIOD
            HDCDIR(I)=2                            ! ONLY GE TYPE USED
  410     ENDDO
C
C----END OF GRADIENT LOOP
        ENDIF
        JROW=JROW+TNGD(G)
  420 ENDDO
C
C----WRITE THE INFORMATION TO THE OUTPUT FILE
      IF(IPRN.EQ.1)THEN
        IF(NHB.GT.0)THEN
          WRITE(IOUT,4000,ERR=990)
          DO 500 I=1,NHB
            IF(NGRIDS.GT.1)THEN                   ! WRITE THE GRID NAME
              IF(I.EQ.1)THEN                     ! FIRST GRID
                WRITE(IOUT,4020,ERR=990)FNAMEN(GRDLOCHDC(I))
              ELSEIF(GRDLOCHDC(I).NE.GRDLOCHDC(I-1))THEN ! NEW GRID
                WRITE(IOUT,4020,ERR=990)FNAMEN(GRDLOCHDC(I))
              ENDIF
            ENDIF
            IF(HDCDIR(I).EQ.1)THEN
              DIRPRT = '<'
            ELSEIF(HDCDIR(I).EQ.2)THEN
              DIRPRT = '>'
            ENDIF
            IF(HDCKLOC(I,1).EQ.0)THEN ! THIS IS A FLAG; ITS AN MNW WELL HEAD
              MNWCNT = 0            ! COUNT THE NUMBER OF MNW-TYPE VARIABLES
              LOOP_FVW : DO J=1,NFVAR      ! LOOP OVER GWM FLOW VARIABLES
                IF(FVNCELL(J).LT.0)THEN   ! THIS IS MNW-TYPE FLOW VARIABLE
                  DO K=1,ABS(FVNCELL(J))  ! LOOP OVER CELLS FOR THIS VARIABLE
                    MNWCNT = MNWCNT + 1
                    IF(MNWCNT.EQ.HDCILOC(I,1))THEN ! MATCH FOUND 
                      TCMNW=FVMNW(J)%FVWELLID(K)
                      EXIT LOOP_FVW
                    ENDIF               
                  ENDDO
                ENDIF
              ENDDO LOOP_FVW
              WRITE(IOUT,4015,ERR=990)I,HDCNAME(I),TCMNW,
     1             DIRPRT,HDCRHS(I),HDCSP(I)
            ELSE
              WRITE(IOUT,4010,ERR=990)I,HDCNAME(I),HDCKLOC(I,1),
     1             HDCILOC(I,1),HDCJLOC(I,1),DIRPRT,HDCRHS(I),HDCSP(I)
            ENDIF
 500     ENDDO
        ENDIF
        IF(NDD.GT.0)THEN
          WRITE(IOUT,4100,ERR=990)
          DO 510 I=NHB+1,NHB+NDD
            IF(NGRIDS.GT.1)THEN                   ! WRITE THE GRID NAME
              IF(I.EQ.NHB+1)THEN                 ! FIRST GRID
                WRITE(IOUT,4110,ERR=990)FNAMEN(GRDLOCHDC(I))
              ELSEIF(GRDLOCHDC(I).NE.GRDLOCHDC(I-1))THEN ! NEW GRID
                WRITE(IOUT,4110,ERR=990)FNAMEN(GRDLOCHDC(I))
              ENDIF
            ENDIF
            IF(HDCDIR(I).EQ.1)THEN
              WRITE(IOUT,4010,ERR=990)I,HDCNAME(I),HDCKLOC(I,1),
     1                  HDCILOC(I,1),HDCJLOC(I,1),'<',HDCRHS(I),HDCSP(I)
            ENDIF
            IF(HDCDIR(I).EQ.2)THEN
              WRITE(IOUT,4010,ERR=990)I,HDCNAME(I),HDCKLOC(I,1),
     1                  HDCILOC(I,1),HDCJLOC(I,1),'>',HDCRHS(I),HDCSP(I)
            ENDIF
  510     ENDDO
        ENDIF
        IF(NDF.GT.0)THEN
          WRITE(IOUT,4200,ERR=990)
          DO 520 I=(NHB+NDD)+1,(NHB+NDD)+NDF
            IF(NGRIDS.GT.1)THEN                   ! WRITE THE GRID NAME
              IF(I.EQ.(NHB+NDD)+1)THEN           ! FIRST GRID
                WRITE(IOUT,4220,ERR=990)FNAMEN(GRDLOCHDC(I))
              ELSEIF(GRDLOCHDC(I).NE.GRDLOCHDC(I-1))THEN  ! NEW GRID
                WRITE(IOUT,4220,ERR=990)FNAMEN(GRDLOCHDC(I))
              ENDIF
            ENDIF
            WRITE(IOUT,4210,ERR=990)I,HDCNAME(I),HDCKLOC(I,1),
     1            HDCILOC(I,1),HDCJLOC(I,1),HDCKLOC(I,2),HDCILOC(I,2),
     2            HDCJLOC(I,2),HDCRHS(I),HDCSP(I)

  520     ENDDO
        ENDIF
        IF(NGD.GT.0)THEN
          WRITE(IOUT,4300,ERR=990)
          DO 530 I=(NHB+NDD+NDF)+1,(NHB+NDD+NDF)+NGD
            IF(NGRIDS.GT.1)THEN                   ! WRITE THE GRID NAME
              IF(I.EQ.(NHB+NDD+NDF)+1)THEN       ! FIRST GRID
                WRITE(IOUT,4320,ERR=990)FNAMEN(GRDLOCHDC(I))
              ELSEIF(GRDLOCHDC(I).NE.GRDLOCHDC(I-1))THEN ! NEW GRID
                WRITE(IOUT,4320,ERR=990)FNAMEN(GRDLOCHDC(I))
              ENDIF
            ENDIF
            WRITE(IOUT,4310,ERR=990)I,HDCNAME(I),HDCKLOC(I,1),
     1            HDCILOC(I,1),HDCJLOC(I,1),HDCKLOC(I,2),HDCILOC(I,2),
     2            HDCJLOC(I,2),LEN(I),GRAD(I),HDCSP(I)

  530     ENDDO
        ENDIF
      ENDIF
C
C-----DEALLOCATE LOCAL ARRAYS
      DEALLOCATE (LEN,GRAD,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 993 
C
C-----CLOSE FILES
      WRITE(IOUT,5000)BYTES
      WRITE(IOUT,6000)
      DO I=1,NGRIDS
        CLOSE(UNIT = NUNOPN(I))
      ENDDO
C
 1000 FORMAT(1X,/1X,'OPENING HEAD CONSTRAINTS FILE',
     1  ' ON UNIT ',I4,':',/1X,A200)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. IPRN MUST BE EQUAL TO 0 OR 1: ',
     1  I4)
 3000 FORMAT(1X,/1X,'PROGRAM STOPPED. CONSTRAINT NAME ',A10,' HAS',
     1  ' ALREADY BEEN USED.') 
 3010 FORMAT(/1X,'PROGRAM STOPPED READING HEDCON. MNW WELLID ',A,
     1  ' NOT FOUND ON GRID:',I5)
 3020 FORMAT(/1X,'PROGRAM STOPPED READING HEDCON. MNW WELLID ',A,
     1  ' NOT FOUND')
 3100 FORMAT(1X,/1X,'PROGRAM STOPPED. ROW NUMBER FOR CELL IS OUT OF',
     1  ' BOUNDS: ',I5)
 3200 FORMAT(1X,/1X,'PROGRAM STOPPED. COLUMN NUMBER FOR CELL IS OUT',
     1  ' OF BOUNDS: ',I5)
 3300 FORMAT(1X,/1X,'PROGRAM STOPPED. LAYER NUMBER FOR CELL IS OUT',
     1  ' OF BOUNDS: ',I5)
 3400 FORMAT(1X,/1X,'PROGRAM STOPPED. CONSTRAINT ',A10,' HAS A RIGHT-',
     1  'HAND SIDE',/,' VALUE (RHS) THAT IS GREATER THAN THE',
     2  ' VALUE',/,' ASSIGNED TO MACHINE INFINITY IN GWM: ',D5.2)
 4000 FORMAT(1X,/1X,'HEAD CONSTRAINTS:',/,
     1  T22,'MNW WELLID or',T44,'RIGHT-HAND',T58,'STRESS',/,
     1  T2,'NUMBER',T10,'NAME',T22,'LAY',T27,
     2  'ROW',T32,'COL',T37,'TYPE',T47,'SIDE',T58,'PERIOD',/,' --------'
     3  ,'-------------------------------------------------------')
 4010 FORMAT(1X,I5,T10,A10,T21,I4,T26,I4,T31,I4,T38,
     1  A1,T42,ES11.4,T58,I5)
 4015 FORMAT(1X,I5,T10,A10,T21,A,T41,
     1  A1,T42,ES11.4,T58,I5)
 4020 FORMAT('  HEAD CONSTRAINTS READ FROM FILE: ',A120)
 4100 FORMAT(1X,/1X,'DRAWDOWN CONSTRAINTS:',/,T44,'RIGHT-HAND',
     1  T58,'STRESS',/,T2,'NUMBER',T10,'NAME',T22,'LAY',T27,
     2  'ROW',T32,'COL',T37,'TYPE',T47,'SIDE',T58,'PERIOD',/,' --------'
     3  ,'-------------------------------------------------------')
 4110 FORMAT('  DRAWDOWN CONSTRAINTS READ FROM FILE: ',A120)
 4200 FORMAT(1X,/1X,'HEAD-DIFFERENCE CONSTRAINTS:',/,T24,'FIRST CELL',
     1  T38,'SECOND CELL',T53,'RIGHT-HAND',T65,'STRESS',/,T2,'NUMBER',
     2  T10,'NAME',T22,'LAY',T27,'ROW',T32,'COL',T37,'LAY',T42,'ROW',
     3  T47,'COL',T56,'SIDE',T65,'PERIOD',/,' --------------',
     3  '-------------------------------------------------------')
 4210 FORMAT(1X,I5,T10,A10,T21,I4,T26,I4,T31,I4,T36,I4,T41,I4,T46,I4,
     1  T52,ES11.4,T64,I5)
 4220 FORMAT('  HEAD-DIFFERENCE CONSTRAINTS READ FROM FILE: ',A120)
 4300 FORMAT(1X,/1X,'GRADIENT CONSTRAINTS:',/,T24,'FIRST CELL',
     1  T38,'SECOND CELL',T78,'STRESS',/,T2,'NUMBER',
     2  T10,'NAME',T22,'LAY',T27,'ROW',T32,'COL',T37,'LAY',T42,'ROW',
     3  T47,'COL',T56,'LEN',T68,'GRAD',T78,'PERIOD',/,' --------',
     4  '-------------------------------------------------------',
     5  '-------------------')
 4310 FORMAT(1X,I5,T10,A10,T21,I4,T26,I4,T31,I4,T36,I4,T41,I4,T46,I4,
     1  T52,ES11.4,T65,ES11.4,T77,I5)
 4320 FORMAT('  GRADIENT CONSTRAINTS READ FROM FILE: ',A120)
 5000 FORMAT(/,1X,I8,' BYTES OF MEMORY ALLOCATED TO STORE DATA FOR',
     1               ' HEAD CONSTRAINTS')
 6000 FORMAT(/,1X,'CLOSING HEAD CONSTRAINTS FILE',/)
C
      RETURN
C
C-----FOR ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),IOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1HDC3AR)')
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
     &2X,'-- STOP EXECUTION (GWM1HDC3AR)')
      CALL GSTOP(' ')
C
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1HDC3AR)')
      CALL GSTOP(' ')
C
  993 CONTINUE
C-----ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
      WRITE(IOUT,9930)
 9930 FORMAT(/,1X,'*** ERROR DEALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1HDC3AR)')
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
     &2X,'-- STOP EXECUTION (GWM1HDC3AR)')
      CALL GSTOP(' ')
C
      END SUBROUTINE GWM1HDC3AR
C
C
C***********************************************************************
      SUBROUTINE GWM1HDC3OS(IGRID,KPER,IPERT,HDRY,NDEP,DEPVALS)
C***********************************************************************
C     VERSION: 21JAN2010
C     PURPOSE: ASSIGN COMPUTED HEAD STATE TO STATE ARRAY
C-----------------------------------------------------------------------
      USE GWM1RMS3, ONLY :  DEWATER
      USE GLOBAL,      ONLY: NCOL,NROW,NLAY,HNEW
      USE GWFMNW2MODULE, ONLY:MNW2
      INTEGER(I4B),INTENT(IN)::IGRID,KPER,IPERT,NDEP
      REAL(DP),INTENT(IN)::HDRY
      DOUBLE PRECISION, DIMENSION(NDEP), INTENT(IN), OPTIONAL :: DEPVALS
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,K
      REAL(DP)::STATE1,STATE2
      REAL(DP) :: EPS = 1.0D-5
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(NDEP.EQ.0)THEN                          ! CALL IS NOT FROM GWM-VI
C---- FOR HEAD BOUND CONSTRAINTS ASSIGN HEAD VALUE TO HDCSTATE
      DO 100 I=1,NHB+NDD   
        IF(IGRID.EQ.GRDLOCHDC(I))THEN                       ! CONSTRAINT ON GRID
          IF(KPER.EQ.HDCSP(I)) THEN                         ! ACTIVE THIS SP
            IF(HDCKLOC(I,1).EQ.0)THEN ! THIS IS A FLAG; ITS AN MNW WELL HEAD
              STATE1 = MNW2(17,HDCJLOC(I,1)) ! VALUE OF HEAD AT MNW WELL
            ELSE
              STATE1 = HNEW(HDCJLOC(I,1),HDCILOC(I,1),HDCKLOC(I,1))
            ENDIF
            IF(FUZZY_EQUALS(STATE1,HDRY,EPS))DEWATER(IPERT)=.TRUE. ! CELL IS DEWATERED
            HDCSTATE(I) = STATE1
          ENDIF
        ENDIF
  100 ENDDO
C
C---- FOR HEAD DIFFERENCE CONSTRAINTS ASSIGN HEAD DIFFERENCE TO HDCSTATE
      DO 200 I=NHB+NDD+1,NHB+NDD + NDF+NGD
        IF(IGRID.EQ.GRDLOCHDC(I))THEN                       ! CONSTRAINT ON GRID
          IF(KPER.EQ.HDCSP(I)) THEN                         ! ACTIVE THIS SP
            STATE1 = HNEW(HDCJLOC(I,1),HDCILOC(I,1),HDCKLOC(I,1))
            STATE2 = HNEW(HDCJLOC(I,2),HDCILOC(I,2),HDCKLOC(I,2))
            IF(FUZZY_EQUALS(STATE1,HDRY,EPS))DEWATER(IPERT)=.TRUE. ! CELL IS DEWATERED
            IF(FUZZY_EQUALS(STATE2,HDRY,EPS))DEWATER(IPERT)=.TRUE. ! CELL IS DEWATERED
            HDCSTATE(I) = STATE1 - STATE2
          ENDIF
        ENDIF
  200 ENDDO
C
      ELSEIF(NDEP.GT.0 .AND. PRESENT(DEPVALS))THEN ! CALL IS FROM GWM-VI
C
C---- FOR HEAD BOUND CONSTRAINTS ASSIGN HEAD VALUE TO HDCSTATE
        DO 300 I=1,NHB+NDD   
          STATE1 = DEPVALS(I)
          IF(FUZZY_EQUALS(STATE1,HDRY,EPS))DEWATER(IPERT)=.TRUE.  ! CELL HAS DEWATERED
          HDCSTATE(I) = STATE1
  300   ENDDO
C
C---- FOR HEAD DIFFERENCE CONSTRAINTS ASSIGN HEAD DIFFERENCE TO HDCSTATE
        K = NHB+NDD
        DO 400 I=NHB+NDD+1,NHB+NDD + NDF+NGD
          K=K+1
          STATE1 = DEPVALS(K)
          K=K+1
          STATE2 = DEPVALS(K)
          IF(FUZZY_EQUALS(STATE1,HDRY,EPS))DEWATER(IPERT)=.TRUE.  ! CELL HAS DEWATERED
          IF(FUZZY_EQUALS(STATE2,HDRY,EPS))DEWATER(IPERT)=.TRUE.  ! CELL HAS DEWATERED
          HDCSTATE(I) = STATE1 - STATE2
  400   ENDDO
C
      ENDIF
C
      RETURN
      END SUBROUTINE GWM1HDC3OS
C
C
C***********************************************************************
      SUBROUTINE GWM1HDC3FP(RSTRT,IPERT)
C***********************************************************************
C     VERSION: 16JULY2009
C     PURPOSE - USE SIMULATION RESULTS TO COMPUTE RESPONSE MATRIX AND 
C            AUGMENTED RIGHT HAND SIDE
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : RMFILE
      USE GWM1DCV3, ONLY : FVBASE
      USE GWM1OBJ3, ONLY : SOLNTYP
      USE GWM1RMS3, ONLY : SLPITCNT,IBASE,IREF,IRM,DELINC,
     &                     RHSIN,RHSINF,RANGENAME,RANGENAMEF,CONTYP,NDV
C-----AMAT HAS LOCAL NAME RESMAT
      USE GWM1RMS3, ONLY : RHS,RESMAT => AMAT
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::IPERT
C-----LOCAL VARIABLES
      REAL(DP)::STATEN
      INTEGER(I4B)::ROW,IHDC,I
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IPERT.EQ.-1)THEN
C-------THIS IS A REFERENCE SIMULATION
        CALL SGWM1HDC3FP                         ! MODIFY CONSTRAINTS
      ELSEIF(IPERT.EQ.0)THEN
C-------THIS IS A BASE SIMULATION
        IF(SLPITCNT.GT.0)THEN
C---------THIS IS A BASE SIMULATION FOR A LATER SLP ITERATION
        ELSEIF(.NOT.(IREF.EQ.1 .AND. IBASE.EQ.1))THEN 
C---------THIS IS EITHER THE FIRST SLP OR BASE LP SIMULATION    
          CALL SGWM1HDC3FP                       ! MODIFY CONSTRAINTS                   
        ENDIF
        IHDC = 0
        DO 100 ROW=RSTRT,RSTRT+HDCNUM-1
          IHDC = IHDC+1
C
C---------RETRIEVE THE INITIAL STATE FROM MODFLOW OUTPUT
          HDCSTATE0(IHDC) = HDCSTATE(IHDC)
C
C---------IF SAVING RESPONSE MATRIX WRITE THE INITIAL STATE TO FILE
          IF(IRM.EQ.1 .OR. IRM.EQ.4)WRITE(RMFILE)HDCSTATE0(IHDC)
C
C---------COMPUTE THE DIFFERENCE BETWEEN INITIAL STATE AND RHS
          RHS(ROW) = HDCRHS(IHDC) - HDCSTATE0(IHDC)
100     ENDDO
C
      ELSEIF(IPERT.GT.0)THEN
C-----  THIS IS A PERTURBATION SIMULATION
C
C-------COMPUTE THE RESPONSE MATRIX VALUE FOR EACH CONSTRAINT
        IHDC = 0
        DO 200 ROW=RSTRT,RSTRT+HDCNUM-1
          IHDC = IHDC+1
C
C---------COMPUTE STATE DIFFERENCE (PERTURBED STATE MINUS INITIAL STATE)
          STATEN = HDCSTATE(IHDC) - HDCSTATE0(IHDC)
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
      RSTRT = RSTRT+HDCNUM
C
      RETURN
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1HDC3FP
C***********************************************************************
C
C  PURPOSE: USE RESULTS OF THE FIRST SIMULATION TO MODIFY CONSTRAINTS
C-----------------------------------------------------------------------
C
C-----STORE THE INITIAL RHS FOR LATER OUTPUT
      IHDC = 0
      DO 100 ROW=RSTRT,RSTRT+HDCNUM-1
        IHDC = IHDC+1
        RHSIN(ROW)= HDCRHS(IHDC)                  ! STORE THE INPUT RHS
        RANGENAME(ROW+NDV-1)=HDCNAME(IHDC)        ! LOAD FOR RANGE ANALYSIS OUTPUT
        RHSINF(ROW)= HDCRHS(IHDC)                 ! STORE THE INPUT RHS
        RANGENAMEF(ROW+NDV-1)=HDCNAME(IHDC)       ! LOAD FOR RANGE ANALYSIS OUTPUT
        IF(IHDC.GT.NHB.AND.IHDC.LE.NHB+NDD)THEN   ! DRAWDOWN CONSTRAINT
          CONTYP(ROW)=2                           ! THIS IS TRANSFORMED CONSTRAINT
        ELSE                                      ! NON-DRAWDOWN CONSTRAINT
          CONTYP(ROW)=1                           ! SET FOR RANGE OUTPUT
        ENDIF
  100 ENDDO
C
C-----FOR HEAD CONSTRAINTS, CONVERT DRAWDOWN TO HEAD BOUND
      DO 200 I=NHB+1,NHB+NDD
        IF(IRM.EQ.1 .OR. IRM.EQ.4)WRITE(RMFILE)HDCSTATE(I)! WRITE REFERENCE STATE
        HDCRHS(I) = HDCSTATE(I) - HDCRHS(I)       ! HDCRHS IS NOW MINIMUM HEAD
C------SWITCH DIRECTION OF INEQUALITY: UPPER BOUND ON DRAWDOWN
C      BECOMES LOWER BOUND ON HEAD
        IF(HDCDIR(I).EQ.1)THEN
          HDCDIR(I) = 2
        ELSEIF(HDCDIR(I).EQ.2)THEN
          HDCDIR(I) = 1
        ENDIF
  200 ENDDO
C
      RETURN
      END SUBROUTINE SGWM1HDC3FP
C
      END SUBROUTINE GWM1HDC3FP
C
C
C***********************************************************************
      SUBROUTINE GWM1HDC3FPR(RSTRT,IREAD,COL)
C***********************************************************************
C     VERSION: 19MAR2008
C     PURPOSE - READ RESPONSE MATRIX AND AUGMENTED RIGHT HAND SIDE
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : RMFILE
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,FVBASE
      USE GWM1RMS3, ONLY : RHSIN,RHSINF,RANGENAME,RANGENAMEF,CONTYP,NDV
C-----AMAT HAS LOCAL NAME RESMAT 
      USE GWM1RMS3, ONLY : RESMAT => AMAT
      USE GWM1RMS3, ONLY : RHS
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::IREAD
      INTEGER(I4B),INTENT(IN)::COL
C-----LOCAL VARIABLES
      INTEGER(I4B)::IHDC,ROW,I
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IREAD.EQ.0)THEN          
        CALL SGWM1HDC3FPR                        ! READ THE REFERENCE STATE             
        IHDC = 0
        DO 100 ROW=RSTRT,RSTRT+HDCNUM-1 
          IHDC = IHDC+1
          READ(RMFILE)HDCSTATE0(IHDC)            ! READ THE BASE STATE FROM A FILE
          RHS(ROW) = HDCRHS(IHDC) - HDCSTATE0(IHDC)
          HDCSTATE(IHDC) = HDCSTATE0(IHDC)       ! ASSIGN BASE STATE FOR OUTPUT
  100   ENDDO
      ELSEIF(IREAD.EQ.1)THEN                     ! READ RESPONSE MATRIX
        DO 200 ROW=RSTRT,RSTRT+HDCNUM-1
          READ(RMFILE)RESMAT(ROW,COL)            ! READ EACH RESPONSE COEFFICIENT
          RHS(ROW) = RHS(ROW) + RESMAT(ROW,COL)*FVBASE(COL)
  200   ENDDO
      ENDIF
      RSTRT = RSTRT+HDCNUM                       ! SET NEXT STARTING LOCATION
C
      RETURN
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1HDC3FPR
C***********************************************************************
C
C  PURPOSE: READ RESULTS OF THE FIRST SIMULATION TO MODIFY CONSTRAINTS
C-----------------------------------------------------------------------
C
C-----STORE THE INITIAL RHS FOR LATER OUTPUT
      IHDC = 0
      DO 100 ROW=RSTRT,RSTRT+HDCNUM-1
        IHDC = IHDC+1
        RHSIN(ROW)= HDCRHS(IHDC)                  ! STORE THE INPUT RHS
        RANGENAME(ROW+NDV-1)=HDCNAME(IHDC)        ! LOAD FOR RANGE ANALYSIS OUTPUT
        RHSINF(ROW)= HDCRHS(IHDC)                 ! STORE THE INPUT RHS
        RANGENAMEF(ROW+NDV-1)=HDCNAME(IHDC)       ! LOAD FOR RANGE ANALYSIS OUTPUT
        IF(IHDC.GT.NHB.AND.IHDC.LE.NHB+NDD)THEN   ! DRAWDOWN CONSTRAINT
          CONTYP(ROW)=2                           ! THIS IS TRANSFORMED CONSTRAINT
        ELSE                                      ! NON-DRAWDOWN CONSTRAINT
          CONTYP(ROW)=1                           ! SET FOR RANGE OUTPUT
        ENDIF
  100 ENDDO
C
      DO 200 I=NHB+1,NHB+NDD
        READ(RMFILE)HDCSTATE(I)                   ! READ THE REFERENCE STATE FROM A FILE
        HDCRHS(I) = HDCSTATE(I) - HDCRHS(I)       ! HDCRHS IS NOW MINIMUM HEAD
        IF(HDCDIR(I).EQ.1)THEN                    ! SWITCH DIRECTION OF INEQUALITY
          HDCDIR(I) = 2
        ELSEIF(HDCDIR(I).EQ.2)THEN
          HDCDIR(I) = 1
        ENDIF
 200  ENDDO
C
      RETURN
      END SUBROUTINE SGWM1HDC3FPR
C
      END SUBROUTINE GWM1HDC3FPR
C
C
C***********************************************************************
      SUBROUTINE GWM1HDC3FM(RSTRT,NSLK)
C***********************************************************************
C     VERSION: 19MAR2008
C     PURPOSE: PLACE SLACK COEFFICIENTS IN LP MATRIX STARTING IN ROW RSTRT
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE
      USE GWM1RMS3, ONLY : AMAT
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::NSLK
C-----LOCAL VARIABLES
      INTEGER(I4B)::ROW,IHDC
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----SET HEAD CONSTRAINT SLACKS/SURPLUSES
      IHDC = 0
      DO 100 ROW=RSTRT,RSTRT+NHB+NDD-1
        IHDC = IHDC+1
        IF(HDCDIR(IHDC).EQ.1) THEN
          AMAT(ROW,NSLK) =  ONE
        ELSE
          AMAT(ROW,NSLK) = -ONE
        ENDIF
  100 ENDDO
C
C-----SET HEAD DIFFERENCE CONSTRAINT SLACKS
      DO 200 ROW=RSTRT+NHB+NDD,RSTRT+HDCNUM-1
        AMAT(ROW,NSLK) = -ONE
  200 ENDDO
C
C-----SET NEXT STARTING LOCATION
      RSTRT = RSTRT+HDCNUM
C
      RETURN
      END SUBROUTINE GWM1HDC3FM
C
C
C***********************************************************************
      SUBROUTINE GWM1HDC3OT(RSTRT,IFLG)
C***********************************************************************
C     VERSION: 25JAN2011
C     PURPOSE - WRITE STATUS OF CONSTRAINT
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,ONE
      USE GWM1RMS3, ONLY : RHS,RANGENAME,NDV,HCLOSEG
      USE GWM1BAS3, ONLY : GWM1BAS3CS
      INTEGER(I4B),INTENT(INOUT)::RSTRT
      INTEGER(I4B),INTENT(IN)::IFLG
C-----LOCAL VARIABLES
      CHARACTER(LEN=25)::CTYPE
      INTEGER(I4B)::ROW,IHDC,DIRR
      REAL(DP)::DIFF
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----WRITE STATUS FOR A FLOW PROCESS SIMULATION
      IF(IFLG.EQ.1.OR.IFLG.EQ.4)THEN            
        DO 100 IHDC=1,HDCNUM
          CALL SGWM1HDC3OT
          DIFF = HDCSTATE(IHDC) - HDCRHS(IHDC)
          DIRR = HDCDIR(IHDC)             
          IF(IFLG.EQ.1)THEN
            CALL GWM1BAS3CS(CTYPE,HDCNAME(IHDC),0.0,0.0,DIFF,DIRR,1)
          ELSEIF(IFLG.EQ.4)THEN
            CALL GWM1BAS3CS(CTYPE,HDCNAME(IHDC),
     &                      HDCSTATE(IHDC),HDCRHS(IHDC),DIFF,DIRR,2)
          ENDIF
  100   ENDDO
C
C-----WRITE STATUS FOR THE LINEAR PROGRAM OUTPUT
      ELSEIF(IFLG.EQ.2)THEN
        IHDC = 0
        DO 200 ROW=RSTRT,RSTRT+HDCNUM-1 
          IHDC = IHDC+1
          IF(ABS(RHS(ROW)).GT.ZERO)THEN          ! DUAL VARIABLE IS NON-ZERO 
                                                 ! CONSTRAINT IS BINDING         
            CALL SGWM1HDC3OT
            CALL GWM1BAS3CS(CTYPE,HDCNAME(IHDC),0.0,0.0,RHS(ROW),0,0) 
          ENDIF
  200   ENDDO
        RSTRT = RSTRT+HDCNUM                     ! SET NEXT STARTING LOCATION
      ENDIF
C
      RETURN
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1HDC3OT
C***********************************************************************
C
C  PURPOSE - ASSIGN CONSTRAINT TYPE
C---------------------------------------------------------------------------
C
      IF(IHDC.LE.NHB)THEN
        CTYPE = 'Head Bound'
      ELSEIF(IHDC.LE.NHB+NDD)THEN
        CTYPE = 'Head Drawdown'
      ELSEIF(IHDC.LE.NHB+NDD+NDF)THEN
        CTYPE = 'Head Difference'
      ELSEIF(IHDC.LE.NHB+NDD+NDF+NGD)THEN
        CTYPE = 'Head Gradient'
      ENDIF
C
      RETURN
      END SUBROUTINE SGWM1HDC3OT
C 
      END SUBROUTINE GWM1HDC3OT
C
C
      END MODULE GWM1HDC3