      MODULE GWM1STA3
C     VERSION: 06AUG2015
      USE GWM_SUBS, ONLY: FUZZY_EQUALS, IGETUNIT
      USE GWM_STOP, ONLY:   GSTOP
      USE MF2005_UTLS, ONLY: U2DINT, URDCOM, URWORD
      IMPLICIT NONE
      PRIVATE
      PUBLIC::STANUM,SVNAME,SVBASE,SVILOC,SVJLOC,
     &        SVKLOC,SVSP,FLOWTYPE,STASTATE,STASTATE0,STATERES,STARHS,
     &        NHVAR,NRVAR,NSVAR,NDVAR,NSTADEP,GWMSTADAT
      PUBLIC::GWM1STA3AR,GWM1STA3OS,GWM1STA3FP,GWM1STA3FPR,GWM1STA3OT
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
C
C-----VARIABLES FOR STATE VARIABLES
      INTEGER(I4B),SAVE::STANUM,NHVAR,NRVAR,NSVAR,NDVAR,NSTADEP
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::SVNAME(:)
      REAL(DP),SAVE,ALLOCATABLE::SVBASE(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::SVILOC(:),SVJLOC(:),SVKLOC(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::SVSP(:,:)
      lOGICAL     ,SAVE,ALLOCATABLE::FLOWTYPE(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::GRDLOCSTA(:),NSSVGZL(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::SSVGLOC(:,:),SSVZLOC(:,:)
      REAL(DP),SAVE,ALLOCATABLE::STASTATE(:),STASTATE0(:),STATERES(:,:),
     &        STARHS(:)
      INTEGER(I4B),SAVE,POINTER::SVZONE(:,:,:,:)
      INTEGER(I4B),SAVE,POINTER::MDCELLOC(:,:)
C
      TYPE MDCEL_IJK
        INTEGER(I4B),POINTER,DIMENSION(:,:)::MDCELLOC_TARGET
      END TYPE
      TYPE, PUBLIC :: GWMSTATYPE
        INTEGER(I4B),POINTER,DIMENSION(:,:,:,:)::SVZONE_TARGET
        TYPE(MDCEL_IJK),POINTER,DIMENSION(:)::MDCELVAR
      END TYPE
      TYPE(GWMSTATYPE), SAVE  ::GWMSTADAT(10) 
C
C      STANUM   -total number of state variables     
C      NHVAR    -number of head state variables
C      NRVAR    -number of streamflow state variables
C      NSVAR    -number of storage state variables
C      NDVAR    -number of drain state variables
C      SVNAME   -name of state variable (10 digits)
C      SVBASE   -current value of state variable during solution process 
C      SVILOC(I), SVJLOC(I), SVKLOC(I)
C               -i,j,k location in MODFLOW grid of the ith head state variable
C                  or segment (SVILOC) and reach (SVJLOC) for ith streamflow SV
C                  or number of cells in drain state variable (SVILOC).
C      SVZONE   - indicator array for cells to be included in calculation of storage change
c                  for storage state variables
C      SVSP     -2-d array recording information about the stress period(s) in which the
C                  state variable is to be evaluated; for head and streamflow state 
C                  variables, the second entry is ignored; for storage state variables the 
C                  two entries record start and end times for state variable evaluation;
C                  for drain state variables, if the second entry is negative, the first 
C                  entry is the stress period to evaluate flow; if the second entry is
C                  positive the two entries record start and end times for volume evaluation.
C      FLOWTYPE -logical array indicating if constraint is flow-type (TRUE) or leakage type.
C      MDCELLOC -2-d array storing the layer, row, column indices for the set of drain cells 
C                  assoicated with the ith Drain state variable listed in MDCELVAR
C      GRDLOCSTA-grid number on which state variable is located
C      NSSVGZL  -vector containing the number of grids associated with a storage state variable
C      SSVGLOC  -vector holding storage state variable grid number information 
C      SSVZLOC  -vector holding storage state variable zone number information 
C      STASTATE -vector containing state values resulting from most 
C                 recent flow process simulation
C      STASTATE0-vector used to store base state values during pertubation
C                 for response matrix generation
C      STATERES -matrix response coefficient for the state variable to pumping
C      STARHS   -stores non-variable information from Taylor series eqn for state variables
C
C-----ERROR HANDLING
      INTEGER(I2B)::ISTAT  
      CHARACTER(LEN=200)::FLNM
      CHARACTER(LEN=20)::FILACT,FMTARG,ACCARG 
      INTEGER(I4B)::NDUM
      REAL(SP)::RDUM
C
      CONTAINS
C*************************************************************************
      SUBROUTINE GWM1STA3AR(FNAMEN,IOUT,NPER,NGRIDS)
C*************************************************************************
C     VERSION: 29DEC2011
C     PURPOSE: READ INPUT FROM THE STATE VARIABLE FILE (STAVAR)
C---------------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWM1BAS3PS,ZERO,CUTCOM
      USE GLOBAL,   ONLY : NCOL,NROW,NLAY
      USE GWM1DCV3, ONLY : NFVAR,FVNCELL,FVMNW,GRDLOCDCV
      INTEGER(I4B),INTENT(IN)::IOUT,NPER,NGRIDS
      CHARACTER(LEN=200),INTENT(IN),DIMENSION(NGRIDS)::FNAMEN
C-----LOCAL VARIABLES  
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::SVNM(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::NG(:),NLMX(:)
      INTEGER(I4B),DIMENSION(NGRIDS)::NHVARG,NRVARG,NSVARG,NDVARG
      INTEGER(I4B),DIMENSION(NGRIDS)::NUNOPN,NSVARGU
      INTEGER(I4B)::LOCAT,ISTART,ISTOP,IPRN,IPRNG,LLOC,INMS,INMF,BYTES
      INTEGER(I4B)::JSVROW,NSEG,NRCH,SVSPL,NMDCEL,DRNTYP
      INTEGER(I4B)::ICZS,ICZF,SPSTRT,SPEND,NSVZL,NGMAX,SEQNUM
      INTEGER(I4B)::G,I,J,K,L,II,JJ,KK,IR,IC,IL,JE,JS,JZ,JESV,JSR,IA
      INTEGER(I4B)::MNWCNT,FLORLK,NSLOC
      CHARACTER(LEN=200)::LINE
      CHARACTER(LEN=10)::TCTNM
      CHARACTER(LEN=20)::TCMNW
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----LOOP OVER ALL GRIDS TO OPEN FILES AND COUNT STATE VARIABLES
      NHVARG = 0
      NRVARG = 0
      NSVARG = 0
      NDVARG = 0
      NSVARGU= 0
      IPRN = -1
      DO 100 G=1,NGRIDS
        IF(FNAMEN(G).NE.' ')THEN
C---------OPEN FILE
          NUNOPN(G)=IGETUNIT(10,200)             
          LOCAT=NUNOPN(G)
          WRITE(IOUT,1000,ERR=990)LOCAT,FNAMEN(G)
          OPEN(UNIT=LOCAT,FILE=FNAMEN(G),ACTION='READ',ERR=999)
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
C---------READ NHVAR, NRVAR, NSVAR AND NDVAR
          READ(LOCAT,'(A)',ERR=991)LINE
          CALL CUTCOM(LINE,200)                  ! REMOVE COMMENTS FROM LINE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NHVARG(G),RDUM,IOUT,
     &                LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRVARG(G),RDUM,IOUT,
     &                LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSVARG(G),RDUM,IOUT,
     &                LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NDVARG(G),RDUM,IOUT,
     &                LOCAT)
          IF(NHVARG(G)+NRVARG(G)+NSVARG(G)+NDVARG(G).LE.0)THEN
            CALL GWM1BAS3PS(' PROGRAM STOPPED: NO STATE
     &                    VARIABLES ON GRID',G)
            CALL GSTOP(' ')
          ENDIF
        ENDIF
  100 ENDDO
C
C-----SUM UP THE STATE VARIABLES
      NHVAR=SUM(NHVARG)                          ! TOTAL NUMBER OF HEAD SV
      NRVAR=SUM(NRVARG)                          ! TOTAL NUMBER OF STRM SV
      NSVAR=SUM(NSVARG)     ! HERE, NSVAR=# OF SV READ FROM FILES
      IF(NSVAR.GT.0)THEN
        ALLOCATE(SVNM(NSVAR),NG(NSVAR)) ! TEMPORARY STORAGE
        ! ASSEMBLE LIST OF STATE VARIABLES
        CALL SVREAD(SVNM,NG)
        ! SEARCH OVER STORAGE STATE VARIABLE NAMES; FIND NUMBER UNIQUE  
        CALL SVSRCH(SVNM,NG,NSVAR,NSVARGU,NGMAX)
        DEALLOCATE (SVNM,NG)                     ! ERASE TEMP STORAGE
        NSVAR=SUM(NSVARGU)  ! NOW, NSVAR=NUMBER OF UNIQUE STATE VARIABLES
      ENDIF
      NDVAR=SUM(NDVARG)                          ! TOTAL NUMBER OF DRAIN SV
      STANUM=NHVAR+NRVAR+NSVAR+NDVAR             ! TOTAL NUMBER OF SV
C
C-----WRITE HEADER AND ALLOCATE SPACE FOR STATE VARIABLE INFORMATION
      IF(STANUM.GT.0)THEN
        WRITE(IOUT,3000,ERR=990)NHVAR,NRVAR,NSVAR,NDVAR
        IF(NSVAR.NE.SUM(NSVARG))WRITE(IOUT,3010,ERR=990)SUM(NSVARG)
        ALLOCATE (SVNAME(STANUM),SVBASE(STANUM),
     &      SVSP(STANUM,2),FLOWTYPE(STANUM),GRDLOCSTA(STANUM), 
     &      SVILOC(STANUM),SVJLOC(STANUM),SVKLOC(STANUM),
     &      STASTATE(STANUM),STASTATE0(STANUM),
     &      STATERES(STANUM,NFVAR),STARHS(STANUM),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992 
        BYTES = 10*STANUM + 8*STANUM +  3*4*STANUM + 4*3*STANUM +
     &          2*8*STANUM + 8*STANUM*NFVAR + 8*STANUM
        ALLOCATE(NSSVGZL(NSVAR),SSVGLOC(NSVAR,NGMAX),
     &           SSVZLOC(NSVAR,NGMAX),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992 
        BYTES = BYTES + 4*NSVAR + 2*4*NSVAR*NGMAX
        DO 120 G=1,NGRIDS                 ! ALLOCATE ON EACH GRID
          IF(NSVARG(G).GT.0)THEN
            CALL SGWF2BAS7PNT(G)                 ! POINT TO GRID DIMENSIONS
         ALLOCATE (GWMSTADAT(G)%SVZONE_TARGET(NCOL,NROW,NLAY,NSVARG(G)),
     &                STAT=ISTAT)
            IF(ISTAT.NE.0)GOTO 992 
            BYTES = BYTES + 4*NCOL*NROW*NLAY*NSVARG(G)
          ENDIF
 120    ENDDO
      ELSE
        RETURN                                   ! NO STATE VARIABLES
      ENDIF
C
C-----READ HEAD STATE VARIABLE INFORMATION
      JSVROW=0
      DO 200 G=1,NGRIDS
        IF(FNAMEN(G).NE.' '.AND.NHVARG(G).GT.0)THEN
C---------LOOP OVER HEAD STATE VARIABLES IN EACH ACTIVE GRID FILE
          LOCAT=NUNOPN(G)  
          DO 180 J=1+JSVROW,NHVARG(G)+JSVROW
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
C              TRACK START AND END COL FOR SVNAME
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            TCTNM=LINE(INMS:INMF)   
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL,RDUM,IOUT,LOCAT)
            IF(IL.GT.0)THEN ! THIS IS A CELL CONSTRAINT; READ ROW,COL
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC,RDUM,IOUT,LOCAT)
            ELSE             ! THIS IS AN MNW WELL CONSTRAINT;READ WELLID
              CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            ENDIF
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SVSPL,RDUM,IOUT,LOCAT)
C-----------PROCESS NAME OF STATE VARIABLE; CHECK IF NAME UNIQUE
            CALL CHKNAME(TCTNM,G,J,-1)
            SVNAME(J)=TCTNM                      ! SAVE THE NAME
C-----------PROCESS ROW, COLUMN, LAYER LOCATIONS FOR THIS VARIABLE
            CALL LOADSV(G)                       ! TEST INPUT AND STORE LOCATION
C-----------STORE STRESS PERIOD AND GRID FOR EVALUATION OF THIS STATE VARIABLE
            SVSP(J,1)=SVSPL
            GRDLOCSTA(J)=G
 180      ENDDO
          JSVROW=JSVROW+NHVARG(G)
        ENDIF
 200  ENDDO
C
C-----READ STREAMFLOW STATE VARIABLE INFORMATION
      DO 300 G=1,NGRIDS
        IF(FNAMEN(G).NE.' '.AND.NRVARG(G).GT.0)THEN
C---------LOOP OVER STREAMFLOW STATE VARIABLES IN EACH ACTIVE GRID FILE
          LOCAT=NUNOPN(G)  
          DO 280 J=1+JSVROW,NRVARG(G)+JSVROW
            READ(LOCAT,'(A)',ERR=991)LINE
            CALL CUTCOM(LINE,200)                ! REMOVE COMMENTS FROM LINE
            LLOC=1
C              TRACK START AND END COL FOR SVNAME
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSEG,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRCH,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SVSPL,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,FLORLK,RDUM,IOUT,LOCAT)
C-----------PROCESS NAME OF STATE VARIABLE; CHECK IF NAME UNIQUE
            CALL CHKNAME(LINE(INMS:INMF),G,J,-1)
            SVNAME(J)=LINE(INMS:INMF)            ! SAVE THE NAME
C-----------PROCESS SEGMENT AND REACH NUMBER 
            IF(NSEG.LT.0)CALL GSTOP('PROGRAM STOPPED: NSEG IS NEGATIVE')
            SVILOC(J)=NSEG                       ! STORE SEGMENT NUMBER 
            IF(NRCH.LT.0)CALL GSTOP('PROGRAM STOPPED: NRCH IS NEGATIVE')
            SVJLOC(J)=NRCH                       ! STORE REACH NUMBER 
C-----------STORE STRESS PERIOD AND GRID FOR EVALUATION OF THIS STATE VARIABLE
            SVSP(J,1)=SVSPL
            GRDLOCSTA(J)=G
C-----------DETERMINE IF THIS IS A FLOW OR LEAK TYPE CONSTRAINT
            IF(FLORLK.EQ.0)THEN
              FLOWTYPE(J) = .TRUE.
            ELSE
              FLOWTYPE(J) = .FALSE.
            ENDIF
 280      ENDDO
          JSVROW=JSVROW+NRVARG(G)
        ENDIF
 300  ENDDO
C
C-----WRITE INFORMATION TO OUTPUT FILE
      IF(IPRN.EQ.1)THEN
C-------WRITE HEAD SV INFO
        IF(NHVAR.GT.0)THEN
          WRITE(IOUT,5000,ERR=990)
          DO 500 J=1,NHVAR
            IF(SVKLOC(J).EQ.0)THEN ! THIS IS A FLAG; ITS AN MNW WELL HEAD
              MNWCNT = 0         ! COUNT THE NUMBER OF MNW-TYPE VARIABLES
              LOOP_FVW : DO JJ=1,NFVAR     ! LOOP OVER GWM FLOW VARIABLES
                IF(FVNCELL(JJ).LT.0)THEN   ! THIS IS MNW-TYPE FLOW VARIABLE
                  DO K=1,ABS(FVNCELL(JJ))  ! LOOP OVER CELLS FOR THIS VARIABLE
                    MNWCNT = MNWCNT + 1
                    IF(MNWCNT.EQ.SVILOC(J))THEN ! MATCH FOUND 
                      TCMNW=FVMNW(JJ)%FVWELLID(K)
                      EXIT LOOP_FVW
                    ENDIF               
                  ENDDO
                ENDIF
              ENDDO LOOP_FVW
              WRITE(IOUT,5110,ERR=990)J,SVNAME(J),TCMNW
            ELSE
              WRITE(IOUT,5100,ERR=990)J,SVNAME(J),
     &             SVKLOC(J),SVILOC(J),SVJLOC(J)
            ENDIF
            WRITE(IOUT,5200,ERR=990)SVSP(J,1)
  500     ENDDO
        ENDIF
C-------WRITE STREAMFLOW SV INFO
        IF(NRVAR.GT.0)THEN
          WRITE(IOUT,6000,ERR=990)
           DO 600 J=NHVAR+1,NHVAR+NRVAR
             IF(FLOWTYPE(J))THEN
               WRITE(IOUT,6100,ERR=990)J,SVNAME(J),
     &                SVILOC(J),SVJLOC(J)
             ELSE
               WRITE(IOUT,6102,ERR=990)J,SVNAME(J),
     &                SVILOC(J),SVJLOC(J)
             ENDIF
             WRITE(IOUT,6200,ERR=990)SVSP(J,1)
  600      ENDDO
        ENDIF
C-------WRITE STORAGE SV HEADER
        IF(NSVAR.GT.0)WRITE(IOUT,6500,ERR=990)
      ENDIF
C
C-----READ STORAGE STATE VARIABLE INFORMATION
      JSR = JSVROW  ! SET INDEX TO END OF HEAD AND STRM SV
      JS  = JSR     ! INITIALIZE COUNTER FOR NUMBER OF UNIQUE STATE VARIABLES
      DO 700 G=1,NGRIDS
        IF(FNAMEN(G).NE.' '.AND.NSVARG(G).GT.0)THEN
C---------LOOP OVER STORAGE STATE VARIABLES IN EACH ACTIVE GRID FILE
          LOCAT=NUNOPN(G)  
          SVZONE=>GWMSTADAT(G)%SVZONE_TARGET   ! POINT TO CORRECT SVZONE ARRAY 
          CALL SGWF2BAS7PNT(G)                 ! POINT TO GRID DIMENSIONS
          DO 680 JZ=1+JSVROW,NSVARG(G)+JSVROW
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
C-----------INTERPRET LINE
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SPSTRT,RDUM,IOUT,LOCAT) 
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SPEND,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ICZS,ICZF,1,NDUM,RDUM,IOUT,LOCAT)
C-----------PROCESS NAME OF STATE VARIABLE; CHECK THAT NAME IS UNIQUE
            JE = 0  ! TELL CHKNAME THAT THIS IS A STORAGE STATE VARIABLE
            CALL CHKNAME(LINE(INMS:INMF),G,JS,JE)  
            IF(JE.EQ.0)THEN                      ! NAME IS UNIQUE
              JS = JS + 1                        ! INCREMENT # UNIQUE SV
              SVNAME(JS)=LINE(INMS:INMF)         ! SAVE THE NAME
              SVSP(JS,1)=SPSTRT                  ! STORE START STRESS PERIODS
              SVSP(JS,2)=SPEND                   ! STORE END STRESS PERIOS
              GRDLOCSTA(JS)=G
              JESV = JS-JSR                      ! SHIFT TO SSV INDEXING
              NSSVGZL(JESV)= 1                   ! INITIALIZE NUMBER OF GRIDS
              SSVGLOC(JESV,NSSVGZL(JESV)) = G    ! STORE GRID NUMBER
              SSVZLOC(JESV,NSSVGZL(JESV)) = JZ   ! STORE ZONE NUMBER
              IF(IPRN.EQ.1)                      ! WRITE INFORMATION         
     &          WRITE(IOUT,6600,ERR=990)
     &          JZ,SVNAME(JS),SVSP(JS,1),SVSP(JS,2)
            ELSEIF(JE.GT.0)THEN                  ! NAME EXISTS ON ANOTHER GRID
C             JE IS THE SEQUENCE NUMBER FOR THE EXISTING STORAGE STATE VARIABLE
              JESV = JE-JSR                      ! SHIFT TO SSV INDEXING
              NSSVGZL(JESV)= NSSVGZL(JESV)+1     ! INCREMENT NUMBER OF GRIDS
              SSVGLOC(JESV,NSSVGZL(JESV)) = G    ! STORE GRID NUMBER
              SSVZLOC(JESV,NSSVGZL(JESV)) = JZ   ! STORE ZONE NUMBER
C             CHECK THAT SPSTRT AND SPEND READ ON THIS GRID ARE VALID
              IF(SVSP(JE,1).NE.SPSTRT.OR.SVSP(JE,2).NE.SPEND)THEN
                CALL GWM1BAS3PS(' PROGRAM STOPPED: SPSTRT/SPEND DO NOT M
     &ATCH FOR STORAGE STATE VARIABLE ON GRID',G)
                CALL GSTOP(' ')
              ENDIF     
              IF(IPRN.EQ.1)                      ! WRITE INFORMATION         
     &        WRITE(IOUT,6605,ERR=990)JZ,SVNAME(JE),GRDLOCSTA(JE)
            ENDIF
C
C-----------PROCESS THE ZONE FOR STORAGE VARIABLE 
            IF(LINE(ICZS:ICZF).EQ.'ALL')THEN     ! ALL DOMAIN IN ZONE
              DO 620 KK=1,NLAY
              DO 620 II=1,NROW
              DO 620 JJ=1,NCOL
                SVZONE(JJ,II,KK,JZ)=1
  620         ENDDO 
              IF(IPRN.EQ.1)WRITE(IOUT,6610,ERR=990)
            ELSEIF(LINE(ICZS:ICZF).EQ.'ZONE')THEN ! READ ZONE
              SVZONE(:,:,:,JZ) = 0               ! EMPTY ZONE ARRAY
              READ(LOCAT,*,ERR=991)NSVZL         ! READ NUMBER OF LAYERS
              IF(IPRN.EQ.1.AND.NSVZL.EQ.1)WRITE(IOUT,6620,ERR=990)
              IF(IPRN.EQ.1.AND.NSVZL.GT.1)WRITE(IOUT,6622,ERR=990)NSVZL
              DO 640 K=1,NSVZL
                READ(LOCAT,*,ERR=991)KK          ! READ LAYER NUMBER (LNUM)
                CALL U2DINT(SVZONE(:,:,KK,JZ),
     &          '     STATE VARIABLE ZONE',NROW,NCOL,KK,LOCAT,IOUT)
  640         ENDDO 
            ELSE                                 ! INVALID ENTRY
              CALL GWM1BAS3PS(' PROGRAM STOPPED: INVALID ENTRY CZONE FOR
     & STORAGE STATE VARIABLE ON GRID',G)
              CALL GSTOP(' ')
            ENDIF
 680      ENDDO
          JSVROW=JSVROW+NSVARG(G) ! INCREMENT JSVROW FOR SV READ
        ENDIF
 700  ENDDO
      JSVROW = JSR + NSVAR  ! NOW, JSVROW IS NUMBER OF UNIQUE SV
C
C-----READ DRAIN STATE VARIABLE INFORMATION
C-----WRITE DRAIN SV HEADER
      IF(NDVAR.GT.0.AND.IPRN.EQ.1)WRITE(IOUT,8500,ERR=990)
      DO 800 G=1,NGRIDS
        IF(FNAMEN(G).NE.' '.AND.NDVARG(G).GT.0)THEN
C---------LOOP OVER STREAMFLOW STATE VARIABLES IN EACH ACTIVE GRID FILE
          LOCAT=NUNOPN(G) 
C---------ALLOCATE FOR NUMBER OF DRAINS SV ON THIS GRID 
          ALLOCATE (GWMSTADAT(G)%MDCELVAR(NDVARG(G))) 
          JE = 0
          DO 780 J=1+JSVROW,NDVARG(G)+JSVROW
            JE = JE+1
            READ(LOCAT,'(A)',ERR=991)LINE
            LLOC=1
C              TRACK START AND END COL FOR SVNAME
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NMDCEL,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,DRNTYP,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SVSPL,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SPSTRT,RDUM,IOUT,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SPEND,RDUM,IOUT,LOCAT)
C-----------PROCESS NAME OF STATE VARIABLE; CHECK IF NAME UNIQUE
            CALL CHKNAME(LINE(INMS:INMF),G,J,-1)
            SVNAME(J)=LINE(INMS:INMF)            ! SAVE THE NAME
C-----------PROCESS DRAIN STATE VARIABLE TYPE
            IF(DRNTYP.EQ.1)THEN
              IF(SVSPL.LE.0)CALL GSTOP('PROGRAM STOPPED: SVSP IS NOT SPE
     &CIFIED, BUT DRAIN TYPE IS 1')
              SPSTRT=SVSPL                      ! STORE STRESS PERIOD HERE
              SPEND=-1                          ! USE -1 AS A FLAG FOR TYPE 1
              IF(IPRN.EQ.1)                      ! WRITE INFORMATION         
     &          WRITE(IOUT,8600,ERR=990)
     &          J,SVNAME(J),SPSTRT
            ELSEIF(DRNTYP.EQ.2)THEN
              IF(SPSTRT.LE.0)CALL GSTOP('PROGRAM STOPPED: SPSTRT IS NOT 
     &SPECIFIED, BUT DRAIN TYPE IS 2')
              IF(SPEND.LE.0)CALL GSTOP('PROGRAM STOPPED: SPEND IS NOT SP
     &ECIFIED, BUT DRAIN TYPE IS 2')
              IF(IPRN.EQ.1)                      ! WRITE INFORMATION         
     &          WRITE(IOUT,8610,ERR=990)J,SVNAME(J),SPSTRT,SPEND
            ELSE
              CALL GSTOP('PROGRAM STOPPED: INVALID VALUE OF DRAIN TYPE')
            ENDIF
            SVSP(J,1)=SPSTRT                  ! STORE START STRESS PERIODS
            SVSP(J,2)=SPEND                   ! STORE END STRESS PERIOS
C-----------PROCESS NUMBER OF DRAIN CELLS IN THIS STATE VARIABLE
            SVILOC(J)=NMDCEL                  ! STORE NUMBER OF CELLS
C-----------ALLOCATE FOR NUMBER OF CELLS ON THIS DRAIN SV
          ALLOCATE (GWMSTADAT(G)%MDCELVAR(JE)%MDCELLOC_TARGET(NMDCEL,4))     ! ALLOCATE TEMP SPACE
            BYTES = BYTES + 4*3*NMDCEL
            DO 760 II=1,NMDCEL
              READ(LOCAT,'(A)',ERR=991)LINE
              LLOC=1
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC,RDUM,IOUT,LOCAT)
              GWMSTADAT(G)%MDCELVAR(JE)%MDCELLOC_TARGET(II,1)=IL ! PLACE L,R,C INTO
              GWMSTADAT(G)%MDCELVAR(JE)%MDCELLOC_TARGET(II,2)=IR !  STORAGE 
              GWMSTADAT(G)%MDCELVAR(JE)%MDCELLOC_TARGET(II,3)=IC !
              ! READ OPTIONAL AUXILIARY VARIABLE
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IA,RDUM,IOUT,LOCAT)
              IF(IA.EQ.0)THEN
                GWMSTADAT(G)%MDCELVAR(JE)%MDCELLOC_TARGET(II,4)=-1
                IF(IPRN.EQ.1)WRITE(IOUT,8620,ERR=990)IL,IR,IC ! WRITE INFO    
              ELSE
                GWMSTADAT(G)%MDCELVAR(JE)%MDCELLOC_TARGET(II,4)=IA 
                IF(IPRN.EQ.1)WRITE(IOUT,8622,ERR=990)IL,IR,IC,IA ! WRITE INFO    
              ENDIF
 760        ENDDO           
C-----------STORE STRESS PERIOD AND GRID FOR EVALUATION OF THIS STATE VARIABLE
            GRDLOCSTA(J)=G
 780      ENDDO
          JSVROW=JSVROW+NDVARG(G)
        ENDIF
 800  ENDDO
C
C-----COUNT THE NUMBER OF DEPENDENTS IMPLIED BY STATE VARIABLES
      NSTADEP = NHVAR               ! Add the head state variables
      NSTADEP = NSTADEP + NRVAR          ! Just a single simulated value
      NSTADEP = NSTADEP + NSVAR        ! Add number of storage state variables
      NSLOC = NHVAR+NRVAR+NSVAR
      DO I=NSLOC+1,NSLOC+NDVAR  ! Loop over drain state variables
        NSTADEP = NSTADEP + SVILOC(I)   ! Add # drains included in state variable
      ENDDO
C
C-----CLOSE FILE
      DO I=1,NGRIDS
        CLOSE(NUNOPN(I))
      ENDDO
      WRITE(IOUT,7000,ERR=990)BYTES
      WRITE(IOUT,7100,ERR=990)
C
      RETURN
C
 1000 FORMAT(1X,/1X,'OPENING STATE VARIABLE FILE ON UNIT ',I4,':',
     1  /1X,A200)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. IPRN MUST BE EQUAL TO 0 OR 1: ',
     1  I4)
 2100 FORMAT(1X,/1X,'OPTIMAL STATE VARIABLE VALUES WILL BE WRITTEN ',
     1  'TO UNIT NUMBER:',I4)
 3000 FORMAT(1X,/1X,'NO. OF HEAD STATE VARIABLES (NHVAR):',
     1  T45,I8,/1X,'NO. OF STREAMFLOW STATE VARIABLES (NRVAR):',
     2  T45,I8,/1X,'NO. OF STORAGE STATE VARIABLES (NSVAR):',
     3  T45,I8,/1X,'NO. OF DRAIN STATE VARIABLES (NDVAR):',
     4  T45,I8)
 3010 FORMAT(1X,'NO. OF STORAGE VARIABLE RECORDS IN INPUT:',T45,I8)
 5000 FORMAT(/,T2,'HEAD TYPE STATE VARIABLES:',/,
     1  T30,'MNW WELLID or',/,
     1  T3,'NUMBER',T14,'NAME',T30,'LAY',T35,'ROW',
     2  T40,'COL',/,1X,'----------------------',
     3  '------------------------------------')
 5100 FORMAT(T3,I5,T14,A10,T26,3I5)
 5110 FORMAT(T3,I5,T14,A10,T26,A)
 5200 FORMAT('   AVAILABLE IN STRESS PERIOD: ',I5,/)
 6000 FORMAT(/,T2,'STREAMFLOW TYPE STATE VARIABLES:',/,
     1  T3,'NUMBER',T14,'NAME',T30,'SEG',T40,'REACH',T50,'TYPE',/,
     2  1X,'----------------------',
     3  '------------------------------------')
 6100 FORMAT(T3,I5,T14,A10,T26,I5,T40,I5,T48,'FLOW-TYPE')
 6102 FORMAT(T3,I5,T14,A10,T26,I5,T40,I5,T48,'LEAK-TYPE')
 6200 FORMAT('   AVAILABLE IN STRESS PERIOD: ',I5,/)
 6500 FORMAT(/,T2,'STORAGE TYPE STATE VARIABLES:')
 6600 FORMAT(/,' NUMBER:',I5,' NAME: ',A10,' SPSTRT:',I5,' SPEND:',I5)
 6605 FORMAT(/,' NUMBER:',I5,' NAME: ',A10,'MATCHES EXISTING STATE ',
     &       'VARIABLE ON GRID',I5)
 6610 FORMAT(T2,'STORAGE CHANGE CALCULATED FOR ',
     1  'ALL ACTIVE CELLS IN GRID')
 6620 FORMAT(T2,'STORAGE CHANGE CALCULATED FOR ',
     1  'ONE LAYER IN GRID')
 6622 FORMAT(T2,'STORAGE CHANGE CALCULATED FOR',I5,
     1  ' LAYERS IN GRID')
 7000 FORMAT(/,1X,I8,' BYTES OF MEMORY ALLOCATED TO STORE DATA FOR',
     1               ' STATE VARIABLES')
 7100 FORMAT(/,1X,'CLOSING STATE VARIABLE FILE',/)
 8500 FORMAT(/,T2,'DRAIN TYPE STATE VARIABLES:')
 8600 FORMAT(/,' NUMBER:',I5,', NAME: ',A,
     &  ', FLOW AT END OF STRESS PERIOD:',I5,/
     &  ,T2,'CALCULATED AT CELL(S):',
     &  /,T30,'LAY',T35,'ROW',T40,'COL',T45,'AUX'
     &  /,T30,'---',T35,'---',T40,'---',T45,'---')
 8610 FORMAT(/,' NUMBER:',I5,', NAME: ',A,', VOLUME FROM SPSTRT:',I5,
     &          ' TO SPEND:',I5,
     &  /,T2,'CALCULATED AT CELL(S):'
     &  /,T30,'LAY',T35,'ROW',T40,'COL',T45,'AUX'
     &  /,T30,'---',T35,'---',T40,'---',T45,'---')
 8620 FORMAT(T28,3I5,' NONE')
 8622 FORMAT(T28,4I5)
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
     &2X,'-- STOP EXECUTION (GWM1STA3AR)')
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
     &2X,'-- STOP EXECUTION (GWM1STA3AR)')
      CALL GSTOP(' ')
C
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1STA3AR)')
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
     &2X,'-- STOP EXECUTION (GWM1STA3AR)')
      CALL GSTOP(' ')
C
	CONTAINS
C
C***********************************************************************
      SUBROUTINE SVREAD(SVNM,NG)
C***********************************************************************
C     ASSEMBLE A TEMPORARY LIST OF STATE VARIABLE NAMES
C       SVNM(NSVAR)- THE READ-ORDER LIST OF SVNAMES FROM ALL GRID FILES
C       NG(NSVAR)  - THE GRID FROM WHICH THE ITH SVNM WAS READ
C
      INTEGER(I4B),INTENT(INOUT)::NG(NSVAR)
      CHARACTER(LEN=10),INTENT(INOUT),DIMENSION(NSVAR)::SVNM
      INTEGER(I4B)::IOTMP
      INTEGER(I4B),ALLOCATABLE::DUMMY(:,:)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IOTMP=IGETUNIT(10,200)                      ! TEMP OUTPUT UNIT  
      OPEN(UNIT=IOTMP,STATUS='SCRATCH')          ! CREATE SCRATCH FILE
      JJ=1
      DO 200 G=1,NGRIDS
        IF(FNAMEN(G).NE.' ')THEN
          CALL SGWF2BAS7PNT(G)                   ! POINT TO GRID DIMENSIONS
          LOCAT=NUNOPN(G)                        ! OBTAIN FILE #
          JSVROW=0                               ! START ROW COUNTER
          DO 120 J=1+JSVROW,NHVARG(G)+JSVROW
            READ(LOCAT,'(A)',ERR=991)LINE        ! READ PAST HEAD SV
 120      ENDDO
          JSVROW=JSVROW+NHVARG(G)
C
          DO 140 J=1+JSVROW,NRVARG(G)+JSVROW
            READ(LOCAT,'(A)',ERR=991)LINE        ! READ PAST STRFLW SV
 140      ENDDO
          JSVROW=JSVROW+NRVARG(G)
C
          DO 180 JZ=1+JSVROW,NSVARG(G)+JSVROW
            READ(LOCAT,'(A)',ERR=991)LINE        ! READ A STORAGE SV
            LLOC=1
            CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOTMP,LOCAT)
           CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SPSTRT,RDUM,IOTMP,LOCAT)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,SPEND,RDUM,IOTMP,LOCAT)
            CALL URWORD(LINE,LLOC,ICZS,ICZF,1,NDUM,RDUM,IOTMP,LOCAT)
            SVNM(JJ)= LINE(INMS:INMF)            ! SAVE THE NAME
            NG(JJ)  = G                          ! SAVE THE GRID IT IS ON
            JJ=JJ+1
C
C-----------READ PAST THE ZONE INFORMATION; COLLECT LAYER MAX
            IF(LINE(ICZS:ICZF).EQ.'ALL')THEN     ! ALL DOMAIN IN ZONE
            ELSEIF(LINE(ICZS:ICZF).EQ.'ZONE')THEN ! READ ZONE
              READ(LOCAT,*,ERR=991)NSVZL         ! READ NUMBER OF LAYERS
              DO 160 K=1,NSVZL
                READ(LOCAT,*,ERR=991)KK          ! READ LAYER NUMBER
                ALLOCATE (DUMMY(NROW,NCOL))      ! CREATE ARRAY FOR READ
                CALL U2DINT(DUMMY(:,:),
     &        '              DUMMY READ',NROW,NCOL,0,LOCAT,IOTMP)
                DEALLOCATE (DUMMY)               ! ERASE ARRAY
  160         ENDDO 
            ELSE                                 ! INVALID ENTRY
            CALL GWM1BAS3PS(' PROGRAM STOPPED: INVALID ENTRY CZONE FOR
     & STORAGE STATE VARIABLE ON GRID',G)
            CALL GSTOP(' ')
            ENDIF
  180     ENDDO          
C
          REWIND(LOCAT)                          ! REWIND THE FILE
          CALL URDCOM(LOCAT,IOTMP,LINE)          ! READ FIRST LINES 
          READ(LOCAT,'(A)',ERR=991)LINE          ! POSITION TO READ
        ENDIF
  200 ENDDO
C
      CLOSE(IOTMP)                               ! CLOSE/ERASE TEMP FILE
      RETURN
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
     &2X,'-- STOP EXECUTION (GWM1STA3AR)')
      CALL GSTOP(' ')
      END SUBROUTINE SVREAD
C
C***********************************************************************
      SUBROUTINE SVSRCH(SVNM,NG,NSVAR,NSVARGU,NGMAX)
C***********************************************************************
C      DETERMINE THE NUMBER OF UNIQUE STATE VARIABLE NAMES ON EACH GRID
C      IF A NAME LISTED IN ONE GRID APPEARED IN A PRIOR GRID LIST 
C      IT IS NOT UNIQUE AND INSTEAD IS COUNTED AS AN ADDITIONAL ENTRY 
C      ON THE FIRST GRID LIST.
C      INPUT:  SVNM(NSVAR)- THE READ-ORDER LIST OF SVNAMES FROM ALL GRID FILES
C              NG(NSVAR)  - THE GRID FROM WHICH THE ITH SVNM WAS READ
C              NSVAR      - THE TOTAL NUMBER OF SVNAMES READ
C
C      OUTPUT: NSVARGU(G) - THE NUMBER OF UNIQUE VARIABLE NAMES ON GRID G
C              NGMAX - THE MAXIMUM NUMBER OF UNIQUE VARIABLE NAMES ON ANY GRID
C
      INTEGER(I4B),INTENT(IN)::NSVAR,NG(NSVAR)
      CHARACTER(LEN=10),INTENT(IN),DIMENSION(NSVAR)::SVNM
      INTEGER(I4B),INTENT(INOUT),DIMENSION(NGRIDS)::NSVARGU
      INTEGER(I4B),INTENT(OUT)::NGMAX
      INTEGER(I4B)::IFLG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        NSVARGU = 0         ! ZERO OUT ENTIRE VECTOR
        NSVARGU(NG(1)) = 1  ! INITIALIZE FIRST UNIQUE STATE VARIABLES
        DO 200 I=2,NSVAR
          IFLG = 0
          DO 120 JJ=1,I-1
            IF(SVNM(I).EQ.SVNM(JJ))THEN       ! NAME MATCHES EARLIER NAME
              IF(NG(I).EQ.NG(JJ))THEN         ! THIS IS AN INVALID MATCH
                WRITE(IOUT,4000,ERR=990)SVNM(I),NG(I)
                CALL GSTOP(' ')
              ENDIF
              IFLG = 1
              ! INCREMENT COUNTER ON GRID WHERE NAME FIRST APPEARS
              NSVARGU(NG(JJ))=NSVARGU(NG(JJ))+1  
              GO TO 180
            ENDIF  
 120      ENDDO
          IF(IFLG.EQ.0)THEN                      ! SVNM(I) IS NEW
            NSVARGU(NG(I))=NSVARGU(NG(I))+1      ! INCREMENT GRID COUNTER
          ENDIF
 180      CONTINUE
 200    ENDDO 
C
        NGMAX = 0                                ! FIND NGMAX
        DO 300 G=1,NGRIDS
          NGMAX = MAX(NGMAX,NSVARGU(G))       
 300    ENDDO
C
      RETURN   
C      
 4000 FORMAT(1X,/1X,'PROGRAM STOPPED. STATE VARIABLE NAME ',A,
     1  ' ON GRID',I4,' HAS ALREADY BEEN USED.') 
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),IOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1STA3AR)')
      CALL GSTOP(' ')
C    
      END SUBROUTINE SVSRCH
C      
C***********************************************************************
      SUBROUTINE LOADSV(G)
C***********************************************************************
C     LOAD HEAD-TYPE STATE VARIABLES 
      INTEGER(I4B),INTENT(IN)::G
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CALL SGWF2BAS7PNT(G)                       ! POINT TO GRID DIMENSIONS
      IF(IL.EQ.0)THEN  ! THIS IS AN MNW WELL CONSTRAINT
        TCMNW=LINE(INMS:INMF) ! THIS IS THE WELLID NAME
        CALL UPCASE(TCMNW)  !     convert to uppercase
        MNWCNT = 0            ! COUNT THE NUMBER OF MNW-TYPE VARIABLES
        IC=-1
        LOOP_FV : DO J=1,NFVAR      ! LOOP OVER GWM FLOW VARIABLES
          IF(FVNCELL(J).LT.0)THEN   ! THIS IS MNW-TYPE FLOW VARIABLE
            DO K=1,ABS(FVNCELL(J))  ! LOOP OVER CELLS FOR THIS VARIABLE
              MNWCNT = MNWCNT + 1
              IF(FVMNW(J)%FVWELLID(K).EQ.TCMNW)THEN ! MATCH FOUND 
                IF(G.NE.GRDLOCDCV(J))THEN ! FLOW VARIABLE NOT ON GRID
                  WRITE(IOUT,3010,ERR=990)TRIM(TCMNW),G  
                  CALL GSTOP(' ')
                ENDIF
                IR = MNWCNT ! STORE THE INDEX IN IR1
                IC = IR     ! THIS INDEX GETS CHANGED IN GWM1DCV3ARW
                EXIT LOOP_FV
              ENDIF               
            ENDDO
          ENDIF
        ENDDO LOOP_FV
        IF(IC.EQ.-1)WRITE(IOUT,3020,ERR=990)TRIM(TCMNW)  
        IF(IC.EQ.-1)CALL GSTOP(' ')
      ELSE
C-------CHECK THAT CELL IS ON GRID      
        IF(IR.LT.1 .OR. IR.GT.NROW)THEN
          WRITE(IOUT,1000,ERR=990)IR
          CALL GSTOP(' ')
        ENDIF
        IF(IC.LT.1 .OR. IC.GT.NCOL)THEN
          WRITE(IOUT,2000,ERR=990)IC
          CALL GSTOP(' ')
        ENDIF
        IF(IL.LT.1 .OR. IL.GT.NLAY)THEN
          WRITE(IOUT,3000,ERR=990)IL
          CALL GSTOP(' ')
        ENDIF
      ENDIF
C
C-----LOAD VARIABLES 
      SVILOC(J)=IR                   
      SVJLOC(J)=IC
      SVKLOC(J)=IL
C
      RETURN
 1000 FORMAT(1X,/1X,'PROGRAM STOPPED. ROW NUMBER FOR STATE VARIABLE',
     1  ' IS OUT OF BOUNDS: ',I5)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. COLUMN NUMBER FOR STATE VARIABLE',
     1  ' IS OUT OF BOUNDS: ',I5)
 3000 FORMAT(1X,/1X,'PROGRAM STOPPED. LAYER NUMBER FOR STATE VARIABLE',
     1  ' IS OUT OF BOUNDS: ',I5)
 3010 FORMAT(/1X,'PROGRAM STOPPED READING STAVAR. MNW WELLID ',A,
     1  ' NOT FOUND ON GRID:',I5)
 3020 FORMAT(/1X,'PROGRAM STOPPED READING STAVAR. MNW WELLID ',A,
     1  ' NOT FOUND')
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
     &2X,'-- STOP EXECUTION (LOADSV)')
      CALL GSTOP(' ')
C
      END SUBROUTINE LOADSV
C
C***********************************************************************
      SUBROUTINE CHKNAME(SVNMT,G,J,JE)
C***********************************************************************
C     CHECK THAT STATE VARIABLE NAME HAS NOT BEEN USED
C     EXCEPT FOR STORAGE STATE VARIABLES
      USE GWM1DCV3, ONLY : FVNAME,EVNAME,BVNAME
      USE GWM1DCV3, ONLY : NFVAR,NBVAR,NEVAR
      CHARACTER(LEN=*),INTENT(IN)::SVNMT
      INTEGER(I4B),INTENT(IN)::J,G
      INTEGER(I4B)::JE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       DO 100 JJ=1,J-1
         IF(SVNMT.EQ.SVNAME(JJ))THEN
           IF(JE.EQ.0 .AND. GRDLOCSTA(JJ).NE.G)THEN
             ! A STORAGE STATE VARIABLE THAT HAS REDUNDANT NAME ON ANOTHER GRID
             JE = JJ                        ! RETURN THE SV NUMBER
             GOTO 110                       ! THERE ARE NO MORE REDUNDANCIES
           ELSE                             ! AN IMPROPER REDUNDANCY
             WRITE(IOUT,4000,ERR=990)SVNMT
             CALL GSTOP(' ')
           ENDIF
         ENDIF
  100  ENDDO
  110  CONTINUE
       DO 200 II=1,NFVAR
         IF(SVNMT.EQ.FVNAME(II))THEN        ! A FLOW VARIABLE NAME
           WRITE(IOUT,4000,ERR=990)SVNMT
           CALL GSTOP(' ')
         ENDIF
  200  ENDDO
       DO 300 II=1,NEVAR
         IF(SVNMT.EQ.EVNAME(II))THEN        ! AN EXTERNAL VARIABLE NAME
           WRITE(IOUT,4000,ERR=990)SVNMT
           CALL GSTOP(' ')
         ENDIF
  300  ENDDO
       DO 400 II=1,NBVAR                   
         IF(SVNMT.EQ.BVNAME(II))THEN        ! A BINARY VARIABLE NAME
           WRITE(IOUT,4000,ERR=990)SVNMT
           CALL GSTOP(' ')
         ENDIF
  400  ENDDO

      RETURN
 4000 FORMAT(1X,/1X,'PROGRAM STOPPED. STATE VARIABLE NAME ',A,
     1  ' HAS ALREADY BEEN USED.') 
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
     &2X,'-- STOP EXECUTION (CHKNAME)')
      CALL GSTOP(' ')
C
      END SUBROUTINE CHKNAME
C
	END SUBROUTINE GWM1STA3AR
C
C***********************************************************************
      SUBROUTINE GWM1STA3OS(IGRID,KPER,IPERT,KSTP,HDRY,ITIME,STRON,
     &                      SFRON,BCFON,LPFON,HUFON,IC1,NDEP,DEPVALS)
C***********************************************************************
C     VERSION: 24JAN2012
C     PURPOSE: ASSIGN COMPUTED STATE VARIABLE VALUES TO STATE ARRAY
C-----------------------------------------------------------------------
      USE GLOBAL,       ONLY: NCOL,NROW,NLAY,ISSFLG,IBOUND,HNEW,HOLD,
     1                        BOTM,LBOTM
      USE GWM1RMS3,     ONLY: DEWATER
      USE GWM1BAS3,     ONLY: GWM1BAS3PS,ZERO,ONE
      USE GWFMNW2MODULE,ONLY: MNW2
      USE GWFBASMODULE, ONLY: DELT
      USE GLOBAL,       ONLY: NPER,NSTP ! Only accessed for calls from GWM-VI
      ! use gwm1bas3, only: gwmout  ! add for debugging
      INTEGER(I4B),INTENT(IN)::IGRID,KPER,IPERT,KSTP,ITIME
      INTEGER(I4B),INTENT(IN)::STRON,SFRON,BCFON,LPFON,HUFON
      REAL(DP),INTENT(IN)::HDRY
C-----VARIABLES USED WHEN ROUTINE IS CALLED FROM GWM-VI
      INTEGER(I4B),INTENT(IN)::IC1,NDEP
      DOUBLE PRECISION, DIMENSION(NDEP), INTENT(IN), OPTIONAL :: DEPVALS
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K,T,SP,ISEG,IREC,LOC,NSLOC,SGRID,SZONE,SVN
      INTEGER(I4B)::DRPAUX
      REAL(DP)::STATE1,QTOT,DELTA,STATOT 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(STANUM.EQ.0)RETURN                      ! NO STATE VARIABLES
      IF(NDEP.EQ.0)THEN                          ! CALL IS NOT FROM GWM-VI
!
C-----INITIALIZE STATE VARIABLE STORAGE ARRAY AT FIRST CALL
      IF(KPER.EQ.1.AND.KSTP.EQ.1)THEN            ! THIS IS FIRST TIME THROUGH
        DO I=1,STANUM
          STASTATE(I) = ZERO                
        ENDDO
      ENDIF
C
C-----HEADS AND STREAMFLOWS EVALUATED AT THE END OF A STRESS PERIOD      
      IF(ITIME.EQ.1)THEN                         ! AT END OF STRESS PERIOD
        IF(NHVAR.GT.0)CALL GWM1STA3OSH           ! STORE HEAD VALUES
        IF(NRVAR.GT.0)CALL GWM1STA3OSR           ! STORE STREAMFLOW VALUES
C
C-----STORAGE CHANGE EVALUATED AT END OF EACH TIME STEP
      ELSEIF(ITIME.EQ.0)THEN                     ! AT END OF TIME STEP
        IF(NSVAR.GT.0)CALL GWM1STA3OSS           ! STORE STORAGE CHANGE
      ENDIF
C
C-----DRAIN STATE VARIABLES MAY BE EVALUATED AS FLOWS OR VOLUMES
      IF(NDVAR.GT.0)CALL GWM1STA3OSD             ! STORE DRAIN VALUES
!
      ELSEIF(NDEP.GT.0 .AND. PRESENT(DEPVALS))THEN ! CALL IS FROM GWM-VI
        J = IC1                                  ! POINTER TO NEXT VALUE IN DEPVAL
        DO I=1,NHVAR+NRVAR+NSVAR                 ! LOAD HEAD, STREAM AND STORAGE
          STASTATE(I) = DEPVALS(J)               ! MOVE INTO GWM ARRAY
          J=J+1                                  ! INCREMENT POINTER
        ENDDO
        NSLOC = NHVAR+NRVAR+NSVAR
        DO I=NSLOC+1,NSLOC+NDVAR                 ! LOOP OVER DRAIN STATE VARIABLES
          STATOT=0.0D0
          DO K=1,SVILOC(I)                       ! LOOP OVER NUMBER OF CELLS
            ! By convention, drain flow is negative. 
            ! Change sign to conform to GWM usage.
            STATOT=STATOT-DEPVALS(J)             ! ADD NEXT LEAKAGE
            J=J+1                                ! INCREMENT POINTER
          ENDDO
          STASTATE(I) = STATOT                   ! MOVE INTO GWM ARRAY
        ENDDO
      ENDIF
!     add for debugging
!      if(present(depvals) .or. (kper==nper .and. kstp==nstp(kper)))then
!        write(gwmout,7770)
!        write(gwmout,7771)(i,stastate(i),i=1,stanum)
! 7770   format('in GWM1STA3OS, stastate contains:')
! 7771   format(i3,2x,g25.16)
!      endif
      RETURN
C
	CONTAINS
C***********************************************************************
      SUBROUTINE GWM1STA3OSH
C***********************************************************************
C  PURPOSE: ASSIGN MODFLOW RESULTS FOR EACH HEAD TYPE STATE VARIABLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        REAL(DP) :: EPS = 1.0D-5
        DO I=1,NHVAR
          IF(IGRID.EQ.GRDLOCSTA(I))THEN          ! STATE VARIABLE ON GRID
            IF(SVSP(I,1).EQ.KPER) THEN           ! STATE VAR ACTIVE IN THIS SP
              IF(SVKLOC(I).EQ.0)THEN ! THIS IS A FLAG; ITS AN MNW WELL HEAD
                STATE1 = MNW2(17,SVJLOC(I))      ! VALUE OF HEAD AT MNW WELL
              ELSE
                STATE1 = HNEW(SVJLOC(I),SVILOC(I),SVKLOC(I))
              ENDIF
              IF(FUZZY_EQUALS(STATE1,HDRY,EPS))DEWATER(IPERT)=.TRUE. ! CELL IS DEWATERED
              STASTATE(I) = STATE1
            ENDIF
          ENDIF
        ENDDO
      RETURN
      END SUBROUTINE GWM1STA3OSH
C***********************************************************************
      SUBROUTINE GWM1STA3OSR
C***********************************************************************
C  PURPOSE: ASSIGN MODFLOW RESULTS FOR EACH STREAMFLOW TYPE STATE VARIABLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO I=NHVAR+1,NHVAR+NRVAR
          IF(IGRID.EQ.GRDLOCSTA(I))THEN          ! STATE VARIABLE ON GRID
            IF(SVSP(I,1).EQ.KPER) THEN           ! STATE VAR ACTIVE IN THIS SP
              IF(STRON.GT.0)THEN                 ! STR PACKAGE
                CALL GWM1STA3OSSTR(IGRID,I)
              ELSEIF(SFRON.GT.0)THEN             ! SFR PACKAGE
                CALL GWM1STA3OSSFR(IGRID,I)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      RETURN
      END SUBROUTINE GWM1STA3OSR
C***********************************************************************
      SUBROUTINE GWM1STA3OSS
C***********************************************************************
C  PURPOSE: ASSIGN MODFLOW RESULTS FOR STORAGE CHANGES
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        NSLOC = NHVAR+NRVAR
        DO 300 SVN=1,NSVAR
          J=SVN+NSLOC
          IF(KPER.GE.SVSP(J,1).AND.KPER.LE.SVSP(J,2))THEN
C-----------STORAGE CHANGE IS NEEDED IN THIS STRESS PERIOD 
            DO 280 I=1,NSSVGZL(SVN)              ! LOOP OVER GRIDS 
              SGRID = SSVGLOC(SVN,I)             ! OBTAIN GRID NUMBER
              SZONE = SSVZLOC(SVN,I)             ! OBTAIN ZONE NUMBER 
              SVZONE=>GWMSTADAT(SGRID)%SVZONE_TARGET ! POINT TO CORRECT ZONE
              IF(BCFON.GT.0)THEN                 ! BCF IS ACTIVE
                CALL GWM1STA3OSB(NSLOC,SGRID,SZONE,SVN)
              ELSEIF(LPFON.GT.0)THEN             ! LPF IS ACTIVE
                CALL GWM1STA3OSL(NSLOC,SGRID,SZONE,SVN)
              ELSEIF(HUFON.GT.0)THEN             ! HUF IS ACTIVE
                CALL GWM1BAS3PS('Program Stopped: Storage state variables can
     & not be used with HUF Package',0)          ! HUF NOT IMPLEMENTED
                CALL GSTOP(' ')                  ! WITH STORAGE SV
              ENDIF
 280        ENDDO
          ENDIF                              
 300    ENDDO
      RETURN
      END SUBROUTINE GWM1STA3OSS
C***********************************************************************
      SUBROUTINE GWM1STA3OSD
C***********************************************************************
C  PURPOSE: ASSIGN MODFLOW RESULTS FOR EACH DRAIN-FLOW TYPE STATE VARIABLE
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      USE GWFDRNMODULE,ONLY:DRNAUX
      NSLOC = NHVAR+NRVAR+NSVAR
      DRPAUX = -1                     ! SET INDEX FOR AUX LOCATION  
      DO I=1,20
        IF(DRNAUX(I).EQ.'GWM-DR')THEN
          DRPAUX = I+5                ! AUX LOCATION FOUND
        ENDIF
      ENDDO
      DO I=NSLOC+1,NSLOC+NDVAR
        IF(IGRID.EQ.GRDLOCSTA(I))THEN            ! STATE VARIABLE ON GRID
          IF(SVSP(I,2).LT.0.AND.ITIME.EQ.1)THEN  ! FLOW-TYPE DRAIN SV
            IF(SVSP(I,1).EQ.KPER)THEN            ! STATE VAR ACTIVE IN THIS SP
              CALL SGWM1STA3OSDQ(I-NSLOC,QTOT)
              STASTATE(I) = QTOT                 ! STORE THE FLOW VALUE
            ENDIF
          ELSEIF(SVSP(I,2).GT.0.AND.ITIME.EQ.0)THEN ! VOLUME-TYPE DRAIN SV
            IF(KPER.GE.SVSP(I,1).AND.KPER.LE.SVSP(I,2))THEN
              CALL SGWM1STA3OSDQ(I-NSLOC,QTOT)
              STASTATE(I) = STASTATE(I) + QTOT*DELT ! ACCUMULATE THE VOLUME
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE GWM1STA3OSD
C
C***********************************************************************
      SUBROUTINE GWM1STA3OSSTR(IGRID,ISVAR)
C***********************************************************************
C  PURPOSE: ASSIGN COMPUTED STREAMFLOW VALUES FROM STR PACKAGE
C-----------------------------------------------------------------------
      USE GWFSTRMODULE, ONLY: NSTREM, STRM, ISTRM
      INTEGER(I4B),INTENT(IN)::IGRID,ISVAR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CALL SGWF2STR7PNT(IGRID)                   ! POINT TO CORRECT GRID
      ISEG = SVILOC(ISVAR)                       ! STREAM SEGMENT FOR SV
      IREC = SVJLOC(ISVAR)                       ! STREAM REACH FOR SV
      DO 100 J = 1,NSTREM
        IF(ISTRM(4,J).EQ.ISEG .AND. ISTRM(5,J).EQ.IREC)THEN
          LOC = J                                ! THIS SV AT LOC
          GOTO 200
        ENDIF
  100 ENDDO
      RETURN                                     ! STATE VARIABLE NOT FOUND
  200 CONTINUE
      IF(FLOWTYPE(ISVAR))THEN
        STASTATE(ISVAR) = STRM(9,LOC)            ! LOAD STREAMFLOW VALUE
      ELSE
        STASTATE(ISVAR) = STRM(11,LOC)           ! LOAD LEAKAGE VALUE
      ENDIF
      RETURN
      END SUBROUTINE GWM1STA3OSSTR
C
C***********************************************************************
      SUBROUTINE GWM1STA3OSSFR(IGRID,ISVAR)
C***********************************************************************
C  PURPOSE: ASSIGN COMPUTED STREAMFLOW VALUES FROM SFR PACKAGE
C-----------------------------------------------------------------------
      USE GWFSFRMODULE, ONLY: NSTRM, STRM, ISTRM
      INTEGER(I4B),INTENT(IN)::IGRID,ISVAR
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CALL SGWF2SFR7PNT(IGRID)                   ! POINT TO CORRECT GRID
      ISEG = SVILOC(ISVAR)                       ! STREAM SEGMENT FOR SV
      IREC = SVJLOC(ISVAR)                       ! STREAM REACH FOR SV
      DO 100 J = 1,NSTRM
        IF(ISTRM(4,J).EQ.ISEG .AND. ISTRM(5,J).EQ.IREC)THEN
          LOC = J                                ! THIS SV AT LOC
          GOTO 200
        ENDIF
  100 ENDDO
      RETURN                                     ! STATE VARIABLE NOT FOUND
  200 CONTINUE
      IF(FLOWTYPE(ISVAR))THEN
        STASTATE(ISVAR) = STRM(9,LOC)            ! LOAD STREAMFLOW VALUE
      ELSE
        STASTATE(ISVAR) = STRM(11,LOC)           ! LOAD LEAKAGE VALUE
      ENDIF
      RETURN
      END SUBROUTINE GWM1STA3OSSFR
C
C***********************************************************************
      SUBROUTINE GWM1STA3OSB(SVARCT,SGRID,SZONE,SVN)
C***********************************************************************
C     PURPOSE: ASSEMBLE CHANGE-IN-STORAGE STATE VARIABLE VALUE
C      BUILT FROM SUBROUTINE GWF2BCF7BDS
C-----------------------------------------------------------------------
      USE GWFBCFMODULE,ONLY:LAYCON,SC1,SC2
      INTEGER(I4B),INTENT(IN)::SVARCT,SGRID,SZONE,SVN
      LOGICAL, ALLOCATABLE :: LAYFLG(:)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      CALL SGWF2BCF7PNT(SGRID)                   ! POINT TO CORRECT GRID
C2------IF STEADY STATE, STORAGE TERM IS ZERO
      IF(ISSFLG(KPER).NE.0) RETURN
C
      ALLOCATE(LAYFLG(NLAY))
      DO 100 K=1,NLAY
        IF(LAYCON(K).EQ.3 .OR. LAYCON(K).EQ.2)THEN
          LAYFLG(K)= .TRUE.
        ELSE
          LAYFLG(K)= .FALSE.
        ENDIF
 100  ENDDO
C        
      CALL SGWM1STA3OSS(NCOL,NROW,NLAY,SVARCT,LAYFLG,SC1,SC2,SZONE,SVN)

      RETURN
      END SUBROUTINE GWM1STA3OSB
C
C***********************************************************************
      SUBROUTINE GWM1STA3OSL(SVARCT,SGRID,SZONE,SVN)
C***********************************************************************
C     PURPOSE: GRAB RELEVANT LPF INFO
C      BUILT FROM SUBROUTINE GWF2LPF7BDS
C-----------------------------------------------------------------------
      USE GWFLPFMODULE,ONLY:LAYTYP,SC1,SC2
      INTEGER(I4B),INTENT(IN)::SVARCT,SGRID,SZONE,SVN
      LOGICAL, ALLOCATABLE :: LAYFLG(:)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      
      CALL SGWF2LPF7PNT(SGRID)                   ! POINT TO CORRECT GRID
C2------IF STEADY STATE, STORAGE TERM IS ZERO
      IF(ISSFLG(KPER).NE.0) RETURN
C
      ALLOCATE(LAYFLG(NLAY))
      DO 100 K=1,NLAY
        IF(LAYTYP(K).NE.0)THEN
          LAYFLG(K)= .TRUE.
        ELSE
          LAYFLG(K)= .FALSE.
        ENDIF
 100  ENDDO
C      
      CALL SGWM1STA3OSS(NCOL,NROW,NLAY,SVARCT,LAYFLG,SC1,SC2,SZONE,SVN)
C
      RETURN
      END SUBROUTINE GWM1STA3OSL
C
C***********************************************************************
      SUBROUTINE SGWM1STA3OSS(NCOL,NROW,NLAY,SVARCT,LAYFLG,SC1,SC2,
     &                       SZONE,SVN)
C***********************************************************************
C     PURPOSE: ASSEMBLE CHANGE-IN-STORAGE STATE VARIABLE VALUE
C      BUILT FROM SUBROUTINE GWF2LPF7BDS/GWF2BCF7BDS
C-----------------------------------------------------------------------
      INTEGER(I4B),INTENT(IN)::NCOL,NROW,NLAY,SVARCT,SZONE,SVN
      LOGICAL,     INTENT(IN)::LAYFLG(NLAY)
      REAL(DP),    INTENT(IN)::SC1(NCOL,NROW,NLAY),SC2(NCOL,NROW,NLAY)
C-----LOCAL VARIABLES
      INTEGER(I4B)::I2,I3,KT
      REAL(DP)::STRG,RHO,RHO1,RHO2,SOLD,SNEW,TP,HSING,TLED
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C1------INITIALIZE 1/DELT.
      TLED=ONE/DELT
C
C5----LOOP THROUGH EVERY CELL IN THE GRID.
      KT=0
      DO 300 K=1,NLAY
      IF(LAYFLG(K)) KT=KT+1
      DO 300 I=1,NROW
      DO 300 J=1,NCOL
      IF(SVZONE(J,I,K,SZONE).LE.0)GO TO 300   ! CELL NOT IN ZONE
C
C6------SKIP NO-FLOW AND CONSTANT-HEAD CELLS.
      IF(IBOUND(J,I,K).LE.0) GO TO 300
      HSING=HNEW(J,I,K)
C
C7-----CHECK LAYER TYPE TO SEE IF ONE STORAGE CAPACITY OR TWO.
      IF(.NOT.LAYFLG(K)) GO TO 285
C
C7A----TWO STORAGE CAPACITIES.
      TP=BOTM(J,I,LBOTM(K)-1)
      RHO2=SC2(J,I,KT)*TLED
      RHO1=SC1(J,I,K)*TLED
      SOLD=RHO2
      IF(HOLD(J,I,K).GT.TP) SOLD=RHO1
      SNEW=RHO2
      IF(HSING.GT.TP) SNEW=RHO1
      STRG=SOLD*(HOLD(J,I,K)-TP) + SNEW*TP - SNEW*HSING
      GO TO 288
C
C7B----ONE STORAGE CAPACITY.
  285 RHO=SC1(J,I,K)*TLED
      STRG=RHO*HOLD(J,I,K) - RHO*HSING
C
C-----ACCUMULATE CELL-BY-CELL STORAGE CHANGE IN STORAGE
  288 CONTINUE
      STRG = STRG/TLED                           ! COMPUTE VOLUME
      STASTATE(SVN) = STASTATE(SVN) + STRG       ! ACCUMULATE STORAGE
C
  300 CONTINUE
C
      RETURN
      END SUBROUTINE SGWM1STA3OSS
C
C***********************************************************************
      SUBROUTINE SGWM1STA3OSDQ(INDEX,QTOT)
C***********************************************************************
C     PURPOSE: ASSEMBLE DRAIN STATE VARIABLE VALUE (FROM ROUTINE GWF2DRN7BD)
      USE GWFDRNMODULE,ONLY:NDRAIN,DRAI
      INTEGER(I4B),INTENT(IN)::INDEX
      REAL(DP),INTENT(OUT)::QTOT
C-----LOCAL VARIABLES
      INTEGER(I4B)::IL,IR,IC,L,II,SEQNUM,DRAUX
      REAL(DP)::HHNEW,EL,C,Q 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SGWF2DRN7PNT(IGRID)
      MDCELLOC=>GWMSTADAT(IGRID)%MDCELVAR(INDEX)%MDCELLOC_TARGET ! POINT TO STORAGE
      QTOT=ZERO
C-----LOOP OVER THE NUMBER OF CELLS FOR DRAIN STATE VARIABLE I   
      DO 200 II=1,SVILOC(I)
        IL= MDCELLOC(II,1)  ! RETRIEVE THE L,R,C
        IR= MDCELLOC(II,2) 
        IC= MDCELLOC(II,3) 
        DRAUX = MDCELLOC(II,4)
C-------FIND THE SEQUENCE NUMBER IN THE DRAIN PARAMETER LIST ASSOCIATED WITH THIS CELL
        SEQNUM=-1
        DO 150 L=1,NDRAIN
          IF(IL.EQ.DRAI(1,L))THEN
            IF(IR.EQ.DRAI(2,L))THEN
              IF(IC.EQ.DRAI(3,L))THEN
                IF(DRAUX.LT.0)THEN          ! THIS IS ONLY DRAIN IN CELL
                  SEQNUM=L
                  GOTO 180
                ELSEIF(DRPAUX.NE.-1)THEN    ! AUX VARIABLES AVAILABLE
                  IF(DRAUX.EQ.INT(DRAI(DRPAUX,L)))THEN 
                    SEQNUM=L         ! MATCH TO DRAIN IN CELL FOUND
                    GOTO 180
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
  150   ENDDO
C
        IF(SEQNUM.LT.0)CALL GSTOP('PROGRAM STOPPED: DRAIN STATE VA
     &RIABLE CELL NOT IN THE LIST OF CELLS IN DRAIN PACKAGE')
C5C-----GET DRAIN PARAMETERS FROM DRAIN LIST.
  180   EL=DRAI(4,SEQNUM)
        C=DRAI(5,SEQNUM)
        HHNEW=HNEW(IC,IR,IL)
C5D-----IF HEAD HIGHER THAN DRAIN, CALCULATE Q=C*(HHNEW-EL).
        IF(HHNEW.GT.EL) THEN
          Q=C*(HHNEW-EL)
        ELSE
          Q=ZERO
        END IF
        QTOT=QTOT+Q
 200  ENDDO
C
      RETURN
      END SUBROUTINE SGWM1STA3OSDQ
C      
      END SUBROUTINE GWM1STA3OS
C
C***********************************************************************
      SUBROUTINE GWM1STA3FP(IPERT)
C***********************************************************************
C     VERSION: 11SEPT2009
C     PURPOSE - USE SIMULATION RESULTS TO COMPUTE RESPONSE MATRIX
C               AND RIGHT HAND SIDE VECTOR
C-----------------------------------------------------------------------
      USE GWM1DCV3, ONLY : FVBASE
      USE GWM1BAS3, ONLY : RMFILE
      USE GWM1RMS3, ONLY : IRM,DELINC
      INTEGER(I4B),INTENT(IN)::IPERT
C-----LOCAL VARIABLES
      REAL(DP)::STATEN
      INTEGER(I4B)::ROW
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IPERT.EQ.0)THEN
C-----THIS IS A BASE SIMULATION
        DO 100 ROW=1,STANUM       
C         
C---------RETRIEVE THE INITIAL STATE FROM MODFLOW OUTPUT
          STASTATE0(ROW) = STASTATE(ROW)
C---------IF SAVING RESPONSE MATRIX WRITE THE INITIAL STATE TO FILE
          IF(IRM.EQ.1 .OR. IRM.EQ.4)WRITE(RMFILE)STASTATE0(ROW)
C
C---------SET STARHS TO INITIAL STATE
          STARHS(ROW) = -STASTATE0(ROW)
100     END DO
C
      ELSEIF(IPERT.GT.0)THEN
C-----THIS IS A PERTURBATION SIMULATION
C-------COMPUTE THE RESPONSE MATRIX VALUE FOR EACH STATE VARIABLE
        DO 200 ROW=1,STANUM
C
C---------COMPUTE STATE DIFFERENCE (PERTURBED STATE MINUS INITIAL STATE)
          STATEN = STASTATE(ROW) - STASTATE0(ROW)
C
C---------DIVIDE STATE DIFFERENCE BY THE PERTURBATION SIZE 
          STATERES(ROW,IPERT) = STATEN/DELINC(IPERT)
C
C---------AUGMENT RIGHT HAND SIDE WITH KNOWN PERTURBATION RESPONSE TERM
          STARHS(ROW)=STARHS(ROW)+STATERES(ROW,IPERT)*FVBASE(IPERT)
C
C---------WRITE RESPONSE MATRIX ELEMENT TO DISK
          IF(IRM.EQ.1 .OR. IRM.EQ.4) WRITE(RMFILE)STATERES(ROW,IPERT)
C
  200   ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE GWM1STA3FP

C***********************************************************************
      SUBROUTINE GWM1STA3FPR(IREAD,IFV)
C***********************************************************************
C  VERSION: 11SEPT2009
C  PURPOSE - READ THE STATE VARIABLE RESPONSE MATRIX (STATERES)
C-----------------------------------------------------------------------
      USE GWM1DCV3, ONLY : FVBASE
      USE GWM1BAS3, ONLY : RMFILE
      INTEGER(I4B),INTENT(IN)::IREAD,IFV
C-----LOCAL VARIABLES
      INTEGER(I4B)::ROW
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----STATE VARIABLES DO NOT HAVE REFERENCE SIMULATION
      IF(IREAD.EQ.0)THEN                         ! READ THE BASE STATE
        DO 100 ROW=1,STANUM
          READ(RMFILE)STASTATE0(ROW)             ! READ THE INITIAL STATE FROM A FILE
          STARHS(ROW) =  -STASTATE0(ROW)
          STASTATE(ROW) = STASTATE0(ROW)         ! ASSIGN BASE STATE FOR OUTPUT
  100   ENDDO
      ELSEIF(IREAD.EQ.1)THEN                     ! READ RESPONSE MATRIX
        DO 200 ROW=1,STANUM 
          READ(RMFILE)STATERES(ROW,IFV)          ! READ EACH RESPONSE COEFFICIENT
          STARHS(ROW) = STARHS(ROW) + STATERES(ROW,IFV)*FVBASE(IFV)
  200   ENDDO
      ENDIF
C
      RETURN
      END SUBROUTINE GWM1STA3FPR
C
C
C***********************************************************************
      SUBROUTINE GWM1STA3OT
C***********************************************************************
C     VERSION: 14JAN2010
C     PURPOSE - WRITE STATUS OF CONSTRAINT
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWM1BAS3CS
C-----LOCAL VARIABLES
      CHARACTER(LEN=25)::CTYPE
      INTEGER(I4B)::ROW
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----WRITE STATUS FOR A FLOW PROCESS SIMULATION
        DO 100 ROW=1,STANUM
          IF(ROW.LE.NHVAR)THEN
            CTYPE = 'Head '
          ELSEIF(ROW.LE.NHVAR+NRVAR)THEN
            IF(FLOWTYPE(ROW))THEN
              CTYPE = 'Streamflow'
            ELSE
              CTYPE = 'Streamflow-Leak'
            ENDIF
          ELSEIF(ROW.LE.NHVAR+NRVAR+NSVAR)THEN
            CTYPE = 'Change in Storage'
          ELSEIF(ROW.LE.NHVAR+NRVAR+NSVAR+NDVAR)THEN
            CTYPE = 'Drain'
          ENDIF
          CALL GWM1BAS3CS(CTYPE,SVNAME(ROW),0.0,0.0,STASTATE(ROW),0,-1)
  100   ENDDO
C
      RETURN
C 
      END SUBROUTINE GWM1STA3OT
C
	END MODULE GWM1STA3