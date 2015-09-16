!    Structure of file GWM1DCV3.for
!
!    MODULE GWM1DCV3
!      define storage for all decision variable data
!      CONTAINS
!      SUBROUTINE GWM1DCV3AR
!        CONTAINS
!          routines related to reading GWM input for WEL-type and MNW-type decision variables
!      END SUBROUTINE GWM1DCV3AR
!      FM and BD routines 
!    END MODULE GWM1DCV3

      MODULE GWM1DCV3
C     VERSION: 27AUG2013
      USE GWM_SUBS, ONLY: IGETUNIT
      USE GWM_STOP, ONLY:   GSTOP
      USE MF2005_UTLS, ONLY: URDCOM, URWORD

      IMPLICIT NONE
      PRIVATE  
      PUBLIC::NFVAR,FVNAME,FVMIN,FVMAX,FVBASE,FVINI,FVRATIO,FVILOC,EVSP,
     &        FVJLOC,FVKLOC,FVNCELL,FVON,FVDIR,FVSP,EVNAME,EVMIN,EVMAX,
     &        EVBASE,EVDIR,NEVAR,BVNAME,BVBASE,BVNLIST,BVLIST,NBVAR,
     &        GRDLOCDCV,FVCURRENT,FVWELL,FVMNW,MNW2HERE
      PUBLIC::GWM1DCV3AR,GWM1DCV3FM,GWF2DCV3FM,GWF2DCV3BD,
     &        GWF2DCV3RP,GWM1DCV3FVCPNT
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
C
C-----STORAGE FOR FLOW DECISION VARIABLE DATA
      INTEGER(I4B),SAVE            ::NFVAR
      LOGICAL(LGT),SAVE            ::MNW2HERE
      INTEGER(I4B),SAVE,ALLOCATABLE::FVNCELL(:),FVON(:),FVDIR(:)
      REAL(DP),    SAVE,ALLOCATABLE::FVINI(:),FVCURRENT(:)
      REAL(DP),    SAVE,ALLOCATABLE::FVMIN(:),FVMAX(:),FVBASE(:)
      LOGICAL(LGT),SAVE,ALLOCATABLE::FVSP(:,:)
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::FVNAME(:)
      
      INTEGER(I4B),    SAVE, DIMENSION(:),POINTER::FVILOC,FVJLOC,FVKLOC
      DOUBLE PRECISION,SAVE, DIMENSION(:),POINTER::FVRATIO
      TYPE NCMAX
        INTEGER(I4B),    POINTER,DIMENSION(:)::FVILOC,FVJLOC,FVKLOC
        DOUBLE PRECISION,POINTER,DIMENSION(:)::FVRATIO
      END TYPE
      TYPE(NCMAX),ALLOCATABLE,SAVE  ::FVCELL(:) 
C
C      NFVAR    -number of flow-rate variables
C      FVNAME   -name of flow-rate variable (10 digits)
C      FVMIN    -lower bound for flow-rate variable
C      FVMAX    -upper bound for flow-rate variable
C      FVBASE   -current value of flow-rate variable during solution process 
C      FVCURRENT-flow-rate variable storage used by GWM-VI
C      FVINI    -initial value of flow-rate variables
C      FVCELL   -derived type array that holds cell information for each decision variable
C      FVRATIO  -fraction of total flow applied at each cell associated with
C                  the flow-rate variable in FVCELL(I)
C      FVILOC(J), FVJLOC(J), FVKLOC(J)
C               -i,j,k location in MODFLOW grid of the jth cell associated  
C                  with the ith flow-rate variable in FVCELL(I)
C      FVNCELL  -number of cells associated with the flow-rate variable
C                  if FVNCELL>0, this is WEL-type; if <0, this is MNW-type
C      FVON     -flag for availability of flow-rate variable (>0-yes, =<0-no) 
C      FVDIR    -flag for direction of flow for flow-rate variable 
C                 1 => inject, 2 => withdrawal
C      FVSP     -logical array (NFVAR x # stress periods) which records if a
C                  flow-rate variable is a candidate during a stress period
C
C-----ADDITIONAL STORAGE FOR MNW-TYPE FLOW DECISION VARIABLES 
C       DATA STRUCTURE MIMICS GWFMNW2MODULE
      TYPE MNWWELL
        DOUBLE PRECISION, POINTER,DIMENSION(:,:)::FVMNWINT
        DOUBLE PRECISION, POINTER,DIMENSION(:,:)::FVMNWNOD
        INTEGER                                 ::AUXVALUE
      END TYPE
      TYPE FVWELL
        TYPE(MNWWELL),    POINTER,DIMENSION(:)  ::FVCELLS
        DOUBLE PRECISION, POINTER,DIMENSION(:,:)::FVMNW2
        CHARACTER(LEN=20),POINTER,DIMENSION(:)  ::FVWELLID
        INTEGER          ,POINTER,DIMENSION(:)  ::FVWTOMNW
      END TYPE
      TYPE(FVWELL),ALLOCATABLE,SAVE  ::FVMNW(:)
C
C     Flow Variables may be WEL or MNW-type (this is indicated by the sign of 
C       FVNCELL).   For either type the flow variable data structure above is used.  
C       If MNW-type, then there may be multiple MNW wells in a single Flow Variable.  
C       For each MNW well in a Flow Variable, multiple nodes can be defined.  
C       The additional structure to store information for MNW-type flow variables 
C       is held in the allocatable, derived-type array FVMNW, which holds MNW 
C       information for each of the NFVAR flow variables.  If the Ith flow variable
C       is not MNW-type then the MNW information arrays are not allocated.
C
C	For the Ith flow variable that is MNW-type, two storage arrays are defined;
C		FVMNW(I)%FVWTOMNW(NWELLS)
C		FVMNW(I)%FVWELLID(NWELLS)
C		FVMNW(I)%FVMNW2(31,NWELLS)
C			Where NWELLS is the number of MNW wells associated with the ith flow 
C           decision variable, FVWTOMNW maps from the flow variable and well to the
C           index number of the MNW well, FVWELLID holds the name of each well and FVMNW2 
C           holds data about the well. The value of NWELLS is stored in FVNCELL(I).

C	And a derived-type variable, FVCELLS, is defined for each well.
C		FVMNW(I)%FVCELLS(NWELLS)
C
C		For the Kth well in the Ith flow variable, two arrays and a scalar are defined:
C         FVMNW(I)%FVCELLS(K)%FVMNWINT(NNODES)
C	      FVMNW(I)%FVCELLS(K)%FVMNWNOD(NNODES)
C         FVMNW(I)%FVCELLS(K)%AUXVALUE 
C				Where NNODES is the number of nodes associate with the well and 
C               FVMNWINT and FVMNWNOD hold data about the node. The value of NNODES 
C               for well K is stored in FVMNW2(2,K). Integer variable AUXVALUE holds a 
C               unique integer associated with the MNW2 well and is used by GWM-VI to 
C               uniquely identify MNW wells in binary data written by MNW2 
C               to enable MMProc to compare Qnet to Qdes.
      ! KMNWAUXVAL is a counter for AUXVALUE unique identifiers
      INTEGER :: KMNWAUXVAL = 0 
C
C-----VARIABLES FOR EXTERNAL DECISION VARIABLES
      INTEGER(I4B),SAVE::NEVAR
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::EVNAME(:)
      REAL(DP),SAVE,ALLOCATABLE::EVMIN(:),EVMAX(:)
      REAL(DP),SAVE,ALLOCATABLE::EVBASE(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::EVDIR(:)
      LOGICAL(LGT),SAVE,ALLOCATABLE::EVSP(:,:)
C
C      NEVAR    -number of external variables
C      EVNAME   -name of external variable (10 digits)
C      EVMIN    -lower bound for external variable
C      EVMAX    -upper bound for external variable
C      EVBASE   -current value of external variable during solution process 
C      EVDIR    -code for the type of external variable
C      EVSP     -logical array (NEVAR x # stress periods) which records if a
C                  external variable is a candidate during a stress period
C
C-----VARIABLES FOR BINARY DECISION VARIABLES
      INTEGER(I4B),SAVE::NBVAR
      CHARACTER(LEN=10),SAVE,ALLOCATABLE::BVNAME(:)
      REAL(DP),SAVE,ALLOCATABLE::BVBASE(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::BVNLIST(:),BVLIST(:,:)
C
C      NBVAR    -number of binary variables
C      BVNAME   -name of binary variable (10 digits)
C      BVBASE   -current value of binary variable during solution process 
C      BVNLIST  -number of variables associated with the binary variable
C      BVLIST(I,J)-index of the jth variable associated with the ith constraint
C
C-----VARIALBLE FOR MULTIPLE GRIDS
      INTEGER(I4B),SAVE,ALLOCATABLE::GRDLOCDCV(:)
C
C     GRDLOCDCV -grid on which the ith decision variable resides
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
      SUBROUTINE GWM1DCV3AR(FNAMEN,IOUT,NPER,NFVAR,NEVAR,
     &                      NBVAR,NDV,NGRIDS)
C*************************************************************************
C     VERSION: 17MAY2011
C     PURPOSE: READ INPUT FROM THE DECISION-VARIABLE FILE 
C---------------------------------------------------------------------------
      USE GWM1BAS3, ONLY: ONE,ZERO,GWM1BAS3PS,CUTCOM,GWMWFILE
      USE GLOBAL,   ONLY: NCOL,NROW,NLAY
      INTEGER(I4B),INTENT(IN)::IOUT,NPER,NGRIDS
      INTEGER(I4B),INTENT(OUT)::NFVAR,NEVAR,NBVAR,NDV
      CHARACTER(LEN=200),INTENT(IN),DIMENSION(NGRIDS)::FNAMEN
C-----LOCAL VARIABLES
      CHARACTER(LEN=10)::TFVNAME
      REAL(SP)::TRAT
      INTEGER(I4B)::I,II,J,JJ,IR,IC,IL,NPVMAX,BYTES,NCABS,ncmax
      INTEGER(I4B)::LLOC,INMS,INMF,INMVS,INMVF,N,TNPV
      INTEGER(I4B)::LOCAT,ISTART,ISTOP,IPRN,IPRNG
      INTEGER(I4B),DIMENSION(NGRIDS)::NUNOPN
      INTEGER(I4B)::NC
      INTEGER(I4B)::IPSS,IPSF,IPTS,IPTF,ITYPES,ITYPEF,LINEST,LINEEN
      CHARACTER(LEN=1)::FTYPE,FSTAT
      CHARACTER(LEN=2)::ETYPE
      LOGICAL(LGT)::NFOUND
      CHARACTER(LEN=200)::FNAME,LINE
      INTEGER(I4B),DIMENSION(NGRIDS)::NFVARG,NEVARG,NBVARG
      INTEGER(I4B)::JFVROW,JEVROW,JBVROW,G
      CHARACTER(LEN=10)::ETYPED(7)
C-----ALLOCATE TEMPORARY STORAGE UNTIL SIZE CAN BE DETERMINED
      INTEGER(I4B),ALLOCATABLE::TBVLIST(:,:)
      CHARACTER(LEN=10000),ALLOCATABLE::WSP(:),ESP(:)
      INTEGER(I4B),ALLOCATABLE::GRDLOCFV(:),GRDLOCEV(:),
     &             GRDLOCBV(:)
      DATA ETYPED /'  Import  ','  Export  ','  Head    ',
     &    '  Strmflow','  Storage ','  General ','  Drain   '/
C-----LOCAL VARIABLES RELATED TO MNW
      INTEGER(I4B)::Qlimit,PUMPLOC,PPFLAG,PUMPCAP,LOSSTYPEINDEX
      INTEGER(I4B)::NNA,GMNWMAX,NNODES,K,NCNT,MNWID,INODE,IINT,JFV
      REAL(DP)::RATIO,Rw,Rskin,Kskin,B,C,P,CWC,RwNode,RskinNode
      REAL(DP)::KskinNode,BNode,CNode,PNode,CWCNode,Ztop,Zbotm,PP
      CHARACTER(LEN=20)::LOSSTYPE,WELLNM
      CHARACTER(LEN=20)::WELLNAME,FVDIRNAME
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----LOOP OVER ALL GRIDS TO OPEN FILES AND COUNT DECISION VARIABLES
C
      MNW2HERE = .FALSE.
      NFVARG = 0
      NEVARG = 0
      NBVARG = 0
      IPRN = -1
      ALLOCATE (GWMWFILE(NGRIDS),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
      BYTES = 4*(NGRIDS)
      DO 100 G=1,NGRIDS
        GWMWFILE(G)=0
        IF(FNAMEN(G).NE.' ')THEN
C---------OPEN FILE
          NUNOPN(G)=IGETUNIT(10,200)
          LOCAT=NUNOPN(G)
          WRITE(IOUT,1000,ERR=990)LOCAT,FNAMEN(G)
          FLNM=FNAMEN(G)
          FILACT='READ                '
          ACCARG='SEQUENTIAL          '
          FMTARG='FORMATTED           '
          OPEN(UNIT=LOCAT,FILE=FLNM,ACTION=FILACT,ACCESS=ACCARG,
     &         FORM=FMTARG,ERR=999)
C
C---------CHECK FOR COMMENT LINES
          CALL URDCOM(LOCAT,IOUT,LINE)
C
C---------READ IPRN
          CALL CUTCOM(LINE,200)                  ! REMOVE COMMENTS FROM LINE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPRNG,RDUM,IOUT,LOCAT)
          IF(IPRNG.EQ.0 .OR. IPRNG.EQ.1)THEN
            IPRN = MAX(IPRN,IPRNG)               ! USE MOST DETAILED ECHO
          ELSE
            WRITE(IOUT,2000,ERR=990)IPRNG        ! INVALID VALUE OF IPRN
            CALL GSTOP(' ')
          ENDIF
C
C---------READ GWMWFILE IF PRESENT
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,GWMWFILE(G),RDUM,IOUT,
     &                LOCAT)
          IF(GWMWFILE(G).GT.0)WRITE(IOUT,2100,ERR=990)GWMWFILE(G)
C
C---------READ NFVAR, NEVAR, AND NBVAR
          READ(LOCAT,'(A)',ERR=991)LINE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NFVARG(G),RDUM,IOUT,
     &                LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NEVARG(G),RDUM,IOUT,
     &                LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NBVARG(G),RDUM,IOUT,
     &                LOCAT)
        ENDIF
  100 ENDDO
C-----ADD VARIABLES FROM ALL GRIDS
      NFVAR=SUM(NFVARG)
      NEVAR=SUM(NEVARG)
      NBVAR=SUM(NBVARG)
      IF(NFVAR.LE.0)RETURN                       ! NO FLOW VARIABLES DEFINED
      WRITE(IOUT,3000,ERR=990)NFVAR,NEVAR
      IF(NBVAR.GE.1)THEN
        WRITE(IOUT,3010,ERR=990)NBVAR
      ELSE
        WRITE(IOUT,3020,ERR=990)
      ENDIF
C
C-------READ FLOW-RATE VARIABLE INFORMATION 
C 
      IF(NFVAR.GT.0)THEN
        NDV = NFVAR+NEVAR+NBVAR+1                ! NUMBER OF DECISION VARIABLES
        ALLOCATE (GRDLOCDCV(NDV-1),STAT=ISTAT)   ! GRID LOCATION STORAGE
        BYTES = BYTES + 4*(NDV-1)
        IF(ISTAT.NE.0)GOTO 992 
C
C-------ALLOCATE SPACE FOR FLOW-RATE VARIABLE INFORMATION
        ALLOCATE (FVNAME(NFVAR),FVMIN(NFVAR),FVMAX(NFVAR),FVBASE(NFVAR),
     &            FVINI(NFVAR),FVNCELL(NFVAR),FVON(NFVAR),FVDIR(NFVAR),
     &            FVSP(NFVAR,NPER),FVCURRENT(NFVAR),STAT=ISTAT) 
        IF(ISTAT.NE.0)GOTO 992 
        ALLOCATE (FVMNW(NFVAR),FVCELL(NFVAR),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992 
        BYTES = BYTES + 10*NFVAR + 8*5*NFVAR + 4*3*NFVAR +(NFVAR*NPER)/8
        ALLOCATE (WSP(NFVAR),GRDLOCFV(NFVAR),STAT=ISTAT) ! TEMPORARY ARRAYS
        IF(ISTAT.NE.0)GOTO 992 
        WSP = ' '
        JFVROW=0
        DO 240 G=1,NGRIDS
          IF(FNAMEN(G).NE.' ')THEN
C-----------OPEN FILE
            LOCAT=NUNOPN(G)  
            DO 230 J=1+JFVROW,NFVARG(G)+JFVROW
              READ(LOCAT,'(A)',ERR=991)LINE
              LLOC=1
              CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NC,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IL,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IR,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IC,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,IPTS,IPTF,0,NDUM,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,IPSS,IPSF,0,NDUM,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ITYPES,ITYPEF,0,NDUM,RDUM,
     &                    IOUT,LOCAT)
C
C------------CHECK THAT DECISION-VARIABLE NAME HAS NOT BEEN USED
              DO 200 JJ=1,J-1
                IF(LINE(INMS:INMF).EQ.FVNAME(JJ))THEN
                  WRITE(IOUT,4000,ERR=990)FVNAME(JJ)
                  CALL GSTOP(' ')
                ENDIF
  200         ENDDO
              FVNAME(J)=LINE(INMS:INMF)
              GRDLOCFV(J)=G  
C-------------PROCESS NC, THE NUMBER OF CELLS FOR THIS VARIABLE
              IF(NC.GT.-1 .AND. NC.LT.1)THEN
                WRITE(IOUT,4010,ERR=990)NC       ! VALUE OF NC IS INVALID
                CALL GSTOP(' ')
              ENDIF
              FVNCELL(J)=NC                      ! STORE THE NUMBER OF CELLS
              NCABS=ABS(NC)                      ! ALLOCATE SPACE
              ALLOCATE (FVCELL(J)%FVILOC(NCABS),FVCELL(J)%FVJLOC(NCABS),
     &                  FVCELL(J)%FVKLOC(NCABS),FVCELL(J)%FVRATIO(NCABS)
     &                  ,STAT=ISTAT)
              IF(ISTAT.NE.0)GOTO 992 
              BYTES = BYTES + 4*3*NCABS + 8*NCABS
              CALL GWM1DCV3FVCPNT(J)             ! POINT TO ARRAYS
C
C-------------PROCESS ROW, COLUMN, LAYER LOCATIONS FOR THIS VARIABLE
              IF(NC.LE.-1)THEN                   ! THIS IS AN MNW-TYPE FV
                MNW2HERE = .TRUE.
                JFV = J
                CALL GWM1DCV3ARM(IOUT)           ! READ ADDITIONAL INFO
              ELSEIF(NC.EQ.1)THEN
                CALL LOADDV(1)                   ! TEST INPUT AND STORE LOCATION
                FVRATIO(1)=1.0                   ! ONLY ONE CELL
              ELSE                               ! MULTIPLE WELL-TYPE CELLS
                TRAT=0.0
                DO 210 JJ=1,NC                   ! LOOP OVER THE MULTIPLE CELLS
                  READ(LOCAT,*,ERR=991,END=991)FVRATIO(JJ),IL,IR,IC
                  CALL LOADDV(JJ)                ! TEST INPUT AND STORE LOCATION
                  TRAT=FVRATIO(JJ)+TRAT
  210           ENDDO
                IF(REAL(TRAT,DP).NE.ONE)THEN     ! FRACTIONS DO NOT ADD TO ONE
                  DO 220 JJ=1,NC
                    FVRATIO(JJ)=FVRATIO(JJ)/TRAT ! NORMALIZE FRACTIONS
  220             ENDDO
                ENDIF
              ENDIF
              CALL GWM1DCV3FVCSV(J)              ! SAVE CELL INFO IN ARRAYS
C
C-------------PROCESS FTYPE, THE FLOW-RATE VARIABLE DIRECTION 
              FTYPE=LINE(IPTS:IPTF)
              IF(FTYPE.EQ.'i'.OR.FTYPE.EQ.'I')THEN
                FVDIR(J)=1
              ELSEIF(FTYPE.EQ.'w'.OR.FTYPE.EQ.'W')THEN
                FVDIR(J)=2
              ELSE
                WRITE(IOUT,4020,ERR=990)FTYPE    ! VALUE OF FTYPE IS INVALID 
                CALL GSTOP(' ')
              ENDIF
C
C-------------PROCESS FSTAT, THE FLAG TO INDICATE FLOW-RATE VARIABLE IS ACTIVE
              FSTAT=LINE(IPSS:IPSF)
              IF(FSTAT.EQ.'y'.OR.FSTAT.EQ.'Y')THEN
                FVON(J) = 1
              ELSEIF(FSTAT.EQ.'n'.OR.FSTAT.EQ.'N')THEN
                FVON(J) = -1
              ELSE
                WRITE(IOUT,4030,ERR=990)FSTAT    ! VALUE OF FSTAT IS INVALID
                CALL GSTOP(' ')
              ENDIF
C
C-------------ASSIGN WELLSP TO IDENTIFY STRESS PERIODS VARIABLE IS ACTIVE
              LINEST = 1
              LINEEN = LINEST + ITYPEF-ITYPES
              WSP(J)(LINEST:LINEEN)=LINE(ITYPES:ITYPEF)
              DO WHILE (LINE(ITYPEF:ITYPEF).EQ.'&')  ! Continuation line present
                READ(LOCAT,'(A)',ERR=991)LINE        ! Read the next line
                LLOC=1
                CALL URWORD(LINE,LLOC,ITYPES,ITYPEF,0,NDUM,RDUM,
     &                    IOUT,LOCAT)                ! Find location of start/end
                LINEST = LINEEN                      ! Start of array; Overwrite the &
                LINEEN = LINEST + ITYPEF-ITYPES      ! End of array
                WSP(J)(LINEST:LINEEN)=LINE(ITYPES:ITYPEF)
              ENDDO
              CALL SPARRAY(WSP(J),FVSP,NFVAR,NPER,J)
  230       ENDDO
C
            JFVROW=JFVROW+NFVARG(G)
          ENDIF
  240   ENDDO
        GRDLOCDCV(1:NFVAR)=GRDLOCFV              ! STORE GRID LOCATIONS
        DEALLOCATE (GRDLOCFV,STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 993
      ENDIF
C
C-----READ EXTERNAL VARIABLES INFORMATION
C
      IF(NEVAR.GE.1)THEN
C
C-------ALLOCATE SPACE FOR EXTERNAL VARIABLES
        ALLOCATE (EVNAME(NEVAR),EVMAX(NEVAR),EVMIN(NEVAR),EVBASE(NEVAR),
     &            EVDIR(NEVAR),EVSP(NEVAR,NPER),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992 
        BYTES = BYTES + 10*NEVAR + 8*3*NEVAR + 4*NEVAR +(NEVAR*NPER)/8
        ALLOCATE (ESP(NEVAR),GRDLOCEV(NEVAR),STAT=ISTAT)  ! TEMPORARY ARRAYS
        IF(ISTAT.NE.0)GOTO 992 
        ESP = ' '
C
        JEVROW=0
        DO 320 G=1,NGRIDS
          IF(NEVARG(G).GE.1)THEN
            LOCAT=NUNOPN(G)                      ! SET DCV FILE UNIT
            DO 310 J=1+JEVROW,NEVARG(G)+JEVROW
              READ(LOCAT,'(A)',ERR=991)LINE
              LLOC=1
              CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,IPTS,IPTF,0,NDUM,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ITYPES,ITYPEF,0,NDUM,RDUM,
     &                    IOUT,LOCAT)
C
C-------------CHECK THAT EXTERNAL-VARIABLE NAME HAS NOT BEEN USED
              DO 300 JJ=1,J-1
                IF(LINE(INMS:INMF).EQ.EVNAME(JJ))THEN
                  WRITE(IOUT,5000,ERR=990)EVNAME(JJ)
                  CALL GSTOP(' ')
                ENDIF
  300         ENDDO
              EVNAME(J)=LINE(INMS:INMF)
              GRDLOCEV(J)=G  
C-------------PROCESS ETYPE, THE EXTERNAL VARIABLE DIRECTION 
              ETYPE=LINE(IPTS:IPTF)
              IF(ETYPE.EQ.'IM'.OR.ETYPE.EQ.'im'.OR.ETYPE.EQ.'Im')THEN
                EVDIR(J)=1
             ELSEIF(ETYPE.EQ.'EX'.OR.ETYPE.EQ.'ex'.OR.ETYPE.EQ.'Ex')THEN
                EVDIR(J)=2
             ELSEIF(ETYPE.EQ.'HD'.OR.ETYPE.EQ.'hd'.OR.ETYPE.EQ.'Hd')THEN
                EVDIR(J)=3
             ELSEIF(ETYPE.EQ.'SF'.OR.ETYPE.EQ.'sf'.OR.ETYPE.EQ.'Sf')THEN
                EVDIR(J)=4
             ELSEIF(ETYPE.EQ.'ST'.OR.ETYPE.EQ.'st'.OR.ETYPE.EQ.'St')THEN
                EVDIR(J)=5
             ELSEIF(ETYPE.EQ.'GN'.OR.ETYPE.EQ.'gn'.OR.ETYPE.EQ.'Gn')THEN
                EVDIR(J)=6
             ELSEIF(ETYPE.EQ.'DR'.OR.ETYPE.EQ.'dr'.OR.ETYPE.EQ.'Dr')THEN
                EVDIR(J)=7
              ELSE
                WRITE(IOUT,4020,ERR=990)ETYPE    ! VALUE OF ETYPE IS INVALID 
                CALL GSTOP(' ')
              ENDIF
C
C-------------ASSIGN EVSP TO IDENTIFY STRESS PERIODS VARIABLE IS ACTIVE
              LINEST = 1
              LINEEN = LINEST + ITYPEF-ITYPES
              ESP(J)(LINEST:LINEEN)=LINE(ITYPES:ITYPEF)
              DO WHILE (LINE(ITYPEF:ITYPEF).EQ.'&')  ! Continuation line present
                READ(LOCAT,'(A)',ERR=991)LINE        ! Read the next line
                LLOC=1
                CALL URWORD(LINE,LLOC,ITYPES,ITYPEF,0,NDUM,RDUM,
     &                    IOUT,LOCAT)                ! Find location of start/end
                LINEST = LINEEN                      ! Start of array; Overwrite the &
                LINEEN = LINEST + ITYPEF-ITYPES      ! End of array
                ESP(J)(LINEST:LINEEN)=LINE(ITYPES:ITYPEF)
              ENDDO
              CALL SPARRAY(ESP(J),EVSP,NEVAR,NPER,J)
  310       ENDDO
C
          JEVROW=JEVROW+NEVARG(G)
          ENDIF
  320   ENDDO
C
        GRDLOCDCV(NFVAR+1:NFVAR+NEVAR)=GRDLOCEV  ! STORE GRID LOCATIONS
        DEALLOCATE (GRDLOCEV,STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 993
      ENDIF
C
C----READ IN BINARY VARIABLE INFORMATION
C
      IF(NBVAR.GE.1)THEN
C
C-------ALLOCATE SPACE FOR BINARY VARIABLES
        ALLOCATE (BVNAME(NBVAR),BVBASE(NBVAR),BVNLIST(NBVAR),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992 
        BYTES = BYTES + 10*NBVAR + 8*NBVAR + 4*NBVAR 
        ALLOCATE (TBVLIST(NBVAR,NFVAR+NEVAR),GRDLOCBV(NBVAR),STAT=ISTAT) !TEMP
        IF(ISTAT.NE.0)GOTO 992 
        NPVMAX = 0
C
        JBVROW=0
        DO 450 G=1,NGRIDS
          IF(NBVARG(G).GE.1)THEN
            LOCAT=NUNOPN(G)
            DO 440 J=1+JBVROW,NBVARG(G)+JBVROW
              READ(LOCAT,'(A)',ERR=991)LINE
              LLOC=1
              CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,TNPV,RDUM,IOUT,LOCAT)
C
C-------------CHECK THAT BINARY-VARIABLE NAME HAS NOT BEEN USED
              DO 400 JJ=1,J-1
                IF(LINE(INMS:INMF).EQ.BVNAME(JJ))THEN
                  WRITE(IOUT,6000,ERR=990)BVNAME(JJ)
                  CALL GSTOP(' ')
                ENDIF
  400         ENDDO
              BVNAME(J)=LINE(INMS:INMF)          ! STORE THE VARIABLE NAME
              GRDLOCBV(J)=G  
C
C-------------PROCESS THE NUMBER OF VARIABLES ASSOCIATED WITH THE BINARY
              IF(TNPV.LT.1 .OR. TNPV.GT.NFVAR+NEVAR)THEN
                WRITE(IOUT,6010,ERR=990)BVNAME(J)
                CALL GSTOP(' ')
              ENDIF
              BVNLIST(J)=TNPV                    ! STORE THE NUMBER OF VARIABLES
              NPVMAX = MAX(NPVMAX,TNPV)          ! FIND THE MAXIMUM NUMBER 
C
C-------------READ DECISION VARIABLES ASSOCIATED WITH BINARY VARIABLE
              DO 430 JJ=1,TNPV
                CALL URWORD(LINE,LLOC,INMVS,INMVF,0,NDUM,RDUM,
     &                      IOUT,LOCAT)
                TFVNAME=LINE(INMVS:INMVF)        ! GRAB THE NAME
                IF (TFVNAME.EQ.'&') THEN         ! IF THE NAME IS &
                  READ(LOCAT,'(A)',ERR=991)LINE  ! BVLIST CONTINUES
                  LLOC=1                         ! TO A NEW LINE
                  CALL URWORD(LINE,LLOC,INMVS,INMVF,0,NDUM,RDUM,
     &                        IOUT,LOCAT)
	            TFVNAME=LINE(INMVS:INMVF)      ! GRAB NAME ON NEW LINE
                ENDIF	   
C
C---------------PROCESS THE ASSOCIATED VARIABLE
                NFOUND=.TRUE.
                DO 410 I=1,NFVAR
                  IF(FVNAME(I).EQ.TFVNAME)THEN
                    TBVLIST(J,JJ)=I              ! STORE THE VARIABLE NUMBER
                    NFOUND=.FALSE.
                    EXIT
                  ENDIF
  410           ENDDO
                DO 420 I=1,NEVAR
                  IF(EVNAME(I).EQ.TFVNAME)THEN
                    TBVLIST(J,JJ)=NFVAR+I        ! STORE THE VARIABLE NUMBER
                    NFOUND=.FALSE.
                    EXIT
                  ENDIF
  420           ENDDO
                IF(NFOUND)THEN
                  WRITE(IOUT,6020,ERR=990)TFVNAME
                  CALL GSTOP(' ')
                ENDIF
  430         ENDDO
  440       ENDDO
C
          JBVROW=JBVROW+NBVARG(G)
          ENDIF
  450   ENDDO
C   
C-------MOVE TEMPORARY BVLIST INTO PERMANENT STORAGE
        ALLOCATE (BVLIST(NBVAR,NPVMAX),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992
        BYTES = BYTES + 4*NBVAR*NPVMAX
        DO 510 J=1,NBVAR
          DO 500 JJ=1,BVNLIST(J)
            BVLIST(J,JJ) = TBVLIST(J,JJ)
  500     ENDDO
  510   ENDDO
        DEALLOCATE (TBVLIST,STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 993
C
        GRDLOCDCV(NFVAR+NEVAR+1:NFVAR+NEVAR+NBVAR)=GRDLOCBV ! STORE LOCATIONS
        DEALLOCATE (GRDLOCBV,STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 993
      ENDIF
C
C-----WRITE INFORMATION TO OUTPUT FILE
      IF(IPRN.EQ.1)THEN
        IF(NFVAR.GT.0)THEN                       ! WRITE PUMPING-VARIABLE INFO
        IF(MAXVAL(FVNCELL).GE.1)WRITE(IOUT,7000,ERR=990)
        DO 610 J=1,NFVAR
          IF(FVNCELL(J).GT.0)THEN                ! WRITE WEL-TYPE
            CALL GWM1DCV3FVCPNT(J)               ! POINT TO ARRAYS
            IF(NGRIDS.GT.1)THEN                  ! WRITE DECVAR FILE NAME
              IF(J.EQ.1)THEN                     ! WRITE FILE NAME FOR THIS GRID
                WRITE(IOUT,7005,ERR=990)FNAMEN(GRDLOCDCV(J))
              ELSEIF(GRDLOCDCV(J).NE.GRDLOCDCV(J-1))THEN
                WRITE(IOUT,7005,ERR=990)FNAMEN(GRDLOCDCV(J))
              ENDIF
            ENDIF
            DO 600 JJ=1,FVNCELL(J)
              IF(FVDIR(J).EQ.1)THEN
                IF(JJ.EQ.1)THEN
                  WRITE(IOUT,7010,ERR=990)J,FVNAME(J),'INJECTION',
     1                 FVKLOC(1),FVILOC(1),FVJLOC(1),FVRATIO(JJ)
                ELSE
                  WRITE(IOUT,7020,ERR=990)FVKLOC(JJ),FVILOC(JJ),
     1                            FVJLOC(JJ),FVRATIO(JJ)
                ENDIF
              ELSEIF(FVDIR(J).EQ.2)THEN
                IF(JJ.EQ.1)THEN
                  WRITE(IOUT,7010,ERR=990)J,FVNAME(J),'WITHDRAWAL',
     1                 FVKLOC(1),FVILOC(1),FVJLOC(1),FVRATIO(JJ)
                ELSE
                  WRITE(IOUT,7020,ERR=990)FVKLOC(JJ),FVILOC(JJ),
     1                          FVJLOC(JJ),FVRATIO(JJ)
                ENDIF
              ENDIF
  600       ENDDO
            IF(FVON(J).GT.0)THEN
              WRITE(IOUT,7030,ERR=990)TRIM(WSP(J))
            ELSE
              WRITE(IOUT,7040,ERR=990)
            ENDIF
          ENDIF
  610   ENDDO
 ! Begin output of MNW type wells
        IF(MINVAL(FVNCELL).LE.-1)WRITE(IOUT,7001,ERR=990)
        DO 630 J=1,NFVAR
          IF(FVNCELL(J).LT.0)THEN
            CALL GWM1DCV3FVCPNT(J)             ! POINT TO ARRAYS
            IF(NGRIDS.GT.1)THEN                ! WRITE DECVAR FILE NAME
              IF(J.EQ.1)THEN                   ! WRITE FILE NAME FOR THIS GRID
                WRITE(IOUT,7005,ERR=990)FNAMEN(GRDLOCDCV(J))
              ELSEIF(GRDLOCDCV(J).NE.GRDLOCDCV(J-1))THEN
                WRITE(IOUT,7005,ERR=990)FNAMEN(GRDLOCDCV(J))
              ENDIF
            ENDIF
            DO 620 JJ=1,-FVNCELL(J)
              WELLNAME = FVMNW(J)%FVWELLID(JJ) 
              NNODES =   FVMNW(J)%FVMNW2(2,JJ)
              IF(FVDIR(J).EQ.1)THEN
                FVDIRNAME ='INJECTION'
              ELSEIF(FVDIR(J).EQ.2)THEN
                FVDIRNAME ='WITHDRAWAL'
              ENDIF
              IF(JJ.EQ.1)THEN
                WRITE(IOUT,7012,ERR=990)J,FVNAME(J),FVDIRNAME,1,
     1                 WELLNAME,FVRATIO(JJ)
                CALL WRMNW_ECHO(J,JJ)
              ELSE
                WRITE(IOUT,7014,ERR=990)JJ,WELLNAME,FVRATIO(JJ)
                CALL WRMNW_ECHO(J,JJ)
              ENDIF
  620       ENDDO
            IF(FVON(J).GT.0)THEN
              WRITE(IOUT,7030,ERR=990)WSP(J)
            ELSE
              WRITE(IOUT,7040,ERR=990)
            ENDIF
          ENDIF
  630   ENDDO
        ENDIF
C
        IF(NEVAR.GE.1)THEN                       ! WRITE EXTERNAL-VARIABLE INFO
          WRITE(IOUT,7050,ERR=990)
          DO 700 J=1,NEVAR
            IF(NGRIDS.GT.1)THEN                  ! WRITE DECVAR FILE NAME
              JJ=NFVAR
              IF(J.EQ.1)THEN                     ! WRITE FILE NAME FOR THIS GRID
                WRITE(IOUT,7055,ERR=990)FNAMEN(GRDLOCDCV(J+JJ))
              ELSEIF(GRDLOCDCV(J+JJ).NE.GRDLOCDCV(J-1+JJ))THEN
                WRITE(IOUT,7055,ERR=990)FNAMEN(GRDLOCDCV(J+JJ))
              ENDIF
            ENDIF
            I=EVDIR(J)
            WRITE(IOUT,7060,ERR=990)J,EVNAME(J),ETYPED(I),TRIM(ESP(J))
  700     ENDDO
        ENDIF
C
        IF(NBVAR.GE.1)THEN                       ! WRITE BINARY-VARIABLE INFO
          WRITE(IOUT,7070,ERR=990)
          DO 800 I=1,NBVAR
            IF(NGRIDS.GT.1)THEN                  ! WRITE DECVAR FILE NAME
              JJ=NFVAR+NEVAR
              IF(I.EQ.1)THEN                     ! WRITE FILE NAME FOR THIS GRID
                WRITE(IOUT,7075,ERR=990)FNAMEN(GRDLOCDCV(I+JJ))
              ELSEIF(GRDLOCDCV(I+JJ).NE.GRDLOCDCV(I-1+JJ))THEN
                WRITE(IOUT,7075,ERR=990)FNAMEN(GRDLOCDCV(I+JJ))
              ENDIF
            ENDIF
            TNPV=BVNLIST(I)
            IF(BVLIST(I,1).LE.NFVAR)THEN
              TFVNAME=FVNAME(BVLIST(I,1))
            ELSE
              TFVNAME=EVNAME(BVLIST(I,1)-NFVAR)
            ENDIF
            WRITE(IOUT,7080,ERR=990)I,BVNAME(I),TNPV,TFVNAME
            DO 810 II=2,TNPV
              IF(BVLIST(I,II).LE.NFVAR)THEN
                TFVNAME=FVNAME(BVLIST(I,II))
              ELSE
                TFVNAME=EVNAME(BVLIST(I,II)-NFVAR)
              ENDIF
              WRITE(IOUT,7090,ERR=990)TFVNAME
  810       ENDDO
  800     ENDDO
        ENDIF
      ENDIF
C 
C----CLOSE FILE
      DO I=1,NGRIDS
        CLOSE(NUNOPN(I))
      ENDDO
      WRITE(IOUT,7100,ERR=990)BYTES
      WRITE(IOUT,7110,ERR=990)
C
C-----DEALLOCATE LOCAL ALLOCATABLE ARRAYS
      IF(ALLOCATED(WSP))THEN
        DEALLOCATE(WSP,STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 993
      ENDIF
      IF(ALLOCATED(ESP))THEN
        DEALLOCATE(ESP,STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 993
      ENDIF
C
 1000 FORMAT(1X,/1X,'OPENING DECISION-VARIABLE FILE ON UNIT ',I4,':',
     1  /1X,A200)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. IPRN MUST BE EQUAL TO 0 OR 1: ',
     1  I4)
 2100 FORMAT(1X,/1X,'OPTIMAL FLOW VARIABLE VALUES WILL BE WRITTEN ',
     1  'TO UNIT NUMBER:',I4)
 3000 FORMAT(1X,/1X,'NO. OF FLOW-RATE DECISION VARIABLES (NFVAR):',
     1  T45,I8,/1X,'NO. OF EXTERNAL DECISION VARIABLES (NEVAR):',
     2  T45,I8)
 3010 FORMAT(1X,'BINARY VARIABLES ARE ACTIVE. NBVAR:',T45,I8)
 3020 FORMAT(1X,'BINARY VARIABLES ARE NOT ACTIVE.')
 4000 FORMAT(1X,/1X,'PROGRAM STOPPED. DECISION-VARIABLE NAME ',A,
     1  ' HAS ALREADY BEEN USED.') 
 4010 FORMAT(1X,/1X,'PROGRAM STOPPED. FVNCELL MUST BE GREATER THAN OR',
     1  ' EQUAL TO 1, BUT IS: ',I5)
 4020 FORMAT(1X,/1X,'PROGRAM STOPPED. FTYPE MUST BE W OR I, BUT IS: ',
     1  A)
 4030 FORMAT(1X,/1X,'PROGRAM STOPPED. FSTAT MUST BE Y OR N, BUT IS: ',
     1  A)
 5000 FORMAT(1X,/1X,'PROGRAM STOPPED. EXTERNAL-VARIABLE NAME ',A,
     1  ' HAS ALREADY BEEN USED.')
 6000 FORMAT(1X,/1X,'PROGRAM STOPPED. BINARY-VARIABLE NAME ',A,
     1  ' HAS ALREADY BEEN USED.')
 6010 FORMAT(1X,/1X,'PROGRAM STOPPED. NUMBER OF DECISION VARIABLES',
     1  ' ASSOCIATED WITH BINARY VARIABLE ',A,' IS NOT VALID.')
 6020 FORMAT(1X,/1X,'PROGRAM STOPPED.',A,' WAS NOT DEFINED AS A', 
     1  ' VARIABLE NAME (FVNAME)')
 7000 FORMAT(/,T2,'FLOW-RATE VARIABLES: WEL-TYPE',/,T52,'FRACTION',/,
     1  T3,'NUMBER',T14,'NAME',T25,'TYPE',T35,'LAY',T41,'ROW',
     2  T47,'COL',T53,'OF FLOW',/,1X,'----------------------',
     3  '------------------------------------')
 7001 FORMAT(/,T2,'FLOW-RATE VARIABLES: MNW2-TYPE',
     1  /,T38,'MNW WELL',T59,'FRACTION',/,
     1  T3,'NUMBER',T14,'NAME',T25,'TYPE',T35,'NUMBER',T44,
     2  'WELLID',T60,'OF FLOW',/,1X,'----------------------',
     3  '-------------------------------------------')
C
 7005 FORMAT('  FLOW-RATE VARIABLES READ FROM FILE: ',A120,/)
 7010 FORMAT(I5,6X,A10,1X,A10,1X,3I5,4X,F6.4)
 7012 FORMAT(/,I5,6X,A10,1X,A10,1X,I5,1X,A20,1X,F6.4)
 7014 FORMAT(  33X,                I5,1X,A20,1X,F6.4)
 7020 FORMAT(T34,3I5,4X,F6.4)
 7030 FORMAT('   AVAILABLE IN STRESS PERIODS: ',A,/)
 7040 FORMAT('   UNAVAILABLE IN ALL STRESS PERIODS',/)
 7050 FORMAT(/,T2,'EXTERNAL VARIABLES:',/,/,T3,'NUMBER',T14,'NAME',
     1  T24,'TYPE',/,1X'------------------------------',/)
 7055 FORMAT('  EXTERNAL VARIABLES READ FROM FILE: ',A120,/)
 7060 FORMAT(I5,8X,A10,A10,/,'   AVAILABLE IN STRESS PERIODS: ',A,/)
 7070 FORMAT(/,T2,'BINARY VARIABLES:',/,T24,'NUMBER OF',T46,
     1  'NAME OF',/,T3,'NUMBER',T14,'NAME',T20,
     2  'ASSOCIATED VARIABLES',T42,'ASSOCIATED VARIABLES',/,
     3  1X'-------------------------------------------------------',/)
 7075 FORMAT('  BINARY VARIABLES READ FROM FILE: ',A120,/)
 7080 FORMAT(I5,8X,A10,T25,I5,T46,A)
 7090 FORMAT(T46,A)
 7100 FORMAT(/,1X,I8,' BYTES OF MEMORY ALLOCATED TO STORE DATA FOR',
     1               ' DECISION VARIABLES')
 7110 FORMAT(/,1X,'CLOSING DECISION-VARIABLE FILE',/)
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM)
      WRITE(*,9900)TRIM(FLNM),IOUT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
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
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')
C
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')
C
  993 CONTINUE
C-----ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
 9930 FORMAT(/,1X,'*** ERROR DEALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')
C
  999 CONTINUE
C-----FILE-OPENING ERROR
      WRITE(*,9990)TRIM(FLNM),LOCAT,'OLD',FMTARG,ACCARG,FILACT
      WRITE(IOUT,9990)TRIM(FLNM),LOCAT,'OLD',FMTARG,ACCARG,
     &                 FILACT
 9990 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE STATUS: ',A,/
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')
C
      CONTAINS
C***********************************************************************
      SUBROUTINE LOADDV(JJ)
C***********************************************************************
      INTEGER(I4B)::JJ
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----CHECK THAT CELL IS ON GRID      
      CALL SGWF2BAS7PNT(G)                       ! CHANGE POINTERS TO THIS GRID
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
C
C-----LOAD VARIABLES 
      FVILOC(JJ)=IR                   
      FVJLOC(JJ)=IC
      FVKLOC(JJ)=IL
C
 1000 FORMAT(1X,/1X,'PROGRAM STOPPED. ROW NUMBER FOR WELL IS OUT OF',
     1  ' BOUNDS: ',I5)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. COLUMN NUMBER FOR WELL IS OUT',
     1  ' OF BOUNDS: ',I5)
 3000 FORMAT(1X,/1X,'PROGRAM STOPPED. LAYER NUMBER FOR WELL IS OUT',
     1  ' OF BOUNDS: ',I5)
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM)
      WRITE(*,9900)TRIM(FLNM),IOUT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (LOADDV)')
      CALL GSTOP(' ')
C
      END SUBROUTINE LOADDV
C
C***********************************************************************
      SUBROUTINE SPARRAY(SPSTRNG,SARRAY,NV,NPER,J)
C***********************************************************************
C
C  PURPOSE - CONVERT CHARACTER STRING, SPSTRNG TO LOGICAL ARRAY, SARRAY
C-----------------------------------------------------------------------
      CHARACTER(LEN=10000),INTENT(IN)::SPSTRNG
      INTEGER(I4B),INTENT(IN)::NV,NPER,J
      LOGICAL(LGT),INTENT(OUT)::SARRAY(NV,NPER)
C-----LOCAL VARIABLES
      INTEGER(I4B)::IP,IFLG,IN,IT,NUM,NUMOLD,K
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO 50 K=1,NPER
        SARRAY(J,K) = .FALSE.
   50 ENDDO
      IP = 0
      IFLG = 0
C
C-----SET POINTERS TO ELEMENTS OF THE STRESS PERIOD CHARACTER STRING
  100 IP = IP + 1
      IN=IP
C
C-----FIND NEXT PUNCTUATION MARK IN THE STRESS PERIOD STRING
      CALL CHREAD(SPSTRNG,IP,IN,IT,NUM)
      IF(NUM.LT.1 .OR. NUM.GT.NPER)THEN
        CALL GSTOP('PROGRAM STOPPED: INVALID STRESS PERIOD NUMBER')
      ENDIF
C
C-----PROCESS THE PUNCTUATION MARK AND MOST RECENT NUMBER
      IF(IFLG.EQ.1) THEN                         ! PRIOR MARK WAS A DASH 
        DO 200 K=NUMOLD,NUM                      ! FOR ALL VARIABLES IN SEQUENCE
          SARRAY(J,K)=.TRUE.                     ! SET VARIABLE AS ACTIVE
  200   ENDDO
        IFLG = 0
      ENDIF
      IF(IT.EQ.1) THEN                           ! PUNCTUATION IS A DASH
        IFLG = 1                                 ! SET FLAG FOR DASH
        NUMOLD = NUM                             ! SAVE FIRST NUMBER IN SEQUENCE
      ELSEIF(IT.EQ.2) THEN                       ! PUNCTUATION IS A COLON
        SARRAY(J,NUM)=.TRUE.                     ! SET VARIABLE AS ACTIVE
      ELSEIF(IT.EQ.3) THEN                       ! PUNCTUATION IS A BLANK
        SARRAY(J,NUM)=.TRUE.                     ! SET VARIABLE AS ACTIVE
      ENDIF

      IF(IT.NE.3)GOTO 100                        ! CURRENT CHARACTER NOT A BLANK
C
      RETURN
      END SUBROUTINE SPARRAY
C
C***********************************************************************
      SUBROUTINE CHREAD(WSP,IP,IN,IT,NUM)
C***********************************************************************
C
C  PURPOSE - RETURN THE NEXT PUNCTUATION MARK AND THE PREVIOUS NUMERAL
C           VALUE IN A CHARACTER STRING OF STRESS PERIODS
C
C    INPUT:
C      WSP     CHARACTER STRING WITH LIST OF STRESS PERIODS
C      IP, IN  LOCATION OF FIRST DIGIT IN NEXT NUMERAL
C    OUTPUT:
C      IP      LOCATION OF NEXT PUNCTUATION MARK
C      IT      INTEGER WHICH INDICATES VALUE OF NEXT PUNCTUATION MARK
C      NUM     VALUE OF NEXT NUMERAL
C-----------------------------------------------------------------------
      CHARACTER(LEN=10000),INTENT(IN)::WSP
      INTEGER(I4B),INTENT(INOUT)::IP
      INTEGER(I4B),INTENT(IN)::IN
      INTEGER(I4B),INTENT(OUT)::IT,NUM
      CHARACTER(LEN=1)::CT
      INTEGER(I4B)::II
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----FIND NEXT PUNCTUATION MARK
  100 IF(WSP(IP:IP).EQ.'-') THEN
        IT = 1                                   ! MARK IS A DASH
      ELSEIF(WSP(IP:IP).EQ.':') THEN
        IT = 2                                   ! MARK IS A COLON
      ELSEIF(IP.EQ.INDEX(WSP,' ')) THEN
        IT = 3                                   ! MARK IS A SPACE
      ELSE
        IP = IP + 1
        GOTO 100
      ENDIF
C
C--- EVALUATE CURRENT NUMERAL
      NUM=0
      DO 200 II=IN,IP-1
        CT = WSP(II:II)
        IF(.NOT. (CT.EQ.'0' .OR. CT.EQ.'1' .OR. CT.EQ.'2' .OR.
     &     CT.EQ.'3' .OR. CT.EQ.'4' .OR. CT.EQ.'5' .OR. CT.EQ.'6'
     &     .OR. CT.EQ.'7' .OR.  CT.EQ.'8' .OR. CT.EQ.'9'))
     &  CALL GSTOP('STRESS PERIOD STRING CONTAINS NON-DIGIT')
        NUM = NUM + (ICHAR(CT)-48)*10**(IP-1-II)
  200 ENDDO
C
      RETURN
      END SUBROUTINE CHREAD
C
C***********************************************************************
      SUBROUTINE GWM1DCV3ARM(IOUT)
C***********************************************************************
C     VERSION: 12JAN2012
C     PURPOSE: READ MANAGED MNW INFO AND LOAD INTO SPECIAL STORAGE
C     THIS SUBROUTINE COPIES MANY PARTS OF GWF2MNW27RP
C---------------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,ONE
      IMPLICIT NONE
      INTEGER(I4B),INTENT(IN)::IOUT
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      GMNWMAX = ABS(FVNCELL(JFV))
      ALLOCATE (FVMNW(JFV)%FVMNW2(22,GMNWMAX),STAT=ISTAT)  
      ALLOCATE (FVMNW(JFV)%FVWTOMNW(GMNWMAX),STAT=ISTAT)
      ALLOCATE (FVMNW(JFV)%FVWELLID(GMNWMAX),STAT=ISTAT)
      ALLOCATE (FVMNW(JFV)%FVCELLS(GMNWMAX),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992 
      BYTES = BYTES + 4*GMNWMAX + 8*22*GMNWMAX + 20*GMNWMAX
      FVMNW(JFV)%FVWTOMNW = 0
      FVMNW(JFV)%FVWELLID =' '
      FVMNW(JFV)%FVMNW2   =ZERO
C-----LOOP OVER MNW WELLS ASSOCIATED WITH FLOW VARIABLE JFV
      DO 100 MNWID=1,GMNWMAX
C-------ITEM 3B (SIMILAR TO ITEM 2A IN MNW2)
        IF(GMNWMAX.EQ.1)THEN
          READ(LOCAT,*,ERR=991)WELLNM,NNODES
          RATIO = ONE
        ELSE
          READ(LOCAT,*,ERR=991)RATIO,WELLNM,NNODES
        ENDIF
        CALL ASSIGN3B
C-------ITEM 3C (SIMILAR TO ITEM 2B IN MNW2)
        READ(LOCAT,*,ERR=991)LOSSTYPE,PPFLAG
        CALL ASSIGN3C
C-------NOW THAT NUMBER OF NODES IS KNOWN FOR THIS WELL, ALLOCATE SPACE
        NNA = ABS(NNODES)
        ALLOCATE(FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT(11,NNA),STAT=ISTAT)
        ALLOCATE(FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD(19,NNA),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992 
        BYTES = BYTES + 8*11*NNA + 8*19*NNA
        FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT = ZERO
        FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD = ZERO
        FVMNW(JFV)%FVCELLS(MNWID)%AUXVALUE = GET_NEXT_MNW2_AUXVALUE() ! GWM-VI
C-------ITEM 3D (SIMILAR TO ITEM 2C IN MNW2)
        CALL READSET2C ! These are well-wide parameters; they are stored node-wise
C-------LOOP OVER NODES IN THIS MNW WELLS
        DO NCNT=1,NNA
C---------ITEM 3E (SIMILAR TO ITEM 2D IN MNW2)
          CALL READSET2D
          CALL ASSIGN3E
        ENDDO
c GWM DOES NOT SUPPORT MNW FUNCTIONS READ ON LINE 2E THROUGHT 2H        
c
  100 ENDDO
C GWM DOES NOT READ ANY MNW DATA FOR EACH STRESS PERIOD
      RETURN
!      
  991 CONTINUE ! FILE-READING ERROR
      INQUIRE(LOCAT,NAME=FLNM)
      WRITE(IOUT,9910)TRIM(FLNM),LOCAT
 9910 FORMAT(/,1X,'*** ERROR READING FILE "',A,'" ON UNIT ',I5,/,
     & 2X,'-- STOP EXECUTION (GWM1DCV3ARM)')
      CALL GSTOP(' ')
  992 CONTINUE ! ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1DCV3ARM)')
      CALL GSTOP(' ')

      END SUBROUTINE GWM1DCV3ARM
c
C***********************************************************************
      SUBROUTINE ASSIGN3B
      ! Store the ratio for MNW well of flow variable JFV
        FVRATIO(MNWID)=RATIO
        CALL UPCASE(WELLNM)  !     convert to uppercase
        FVMNW(JFV)%FVWELLID(MNWID)=WELLNM
        FVMNW(JFV)%FVMNW2(2,MNWID)=NNODES
      END SUBROUTINE ASSIGN3B
C***********************************************************************
      SUBROUTINE ASSIGN3C
        CALL UPCASE(LOSSTYPE)
C       STORE INTEGER CODE FOR LOSSTYPE
        if(LOSSTYPE.EQ.'NONE') then
c     for none, NNODES must be 1
          if(NNODES.NE.1) then
          write(iout,*) '***ERROR***  OPTION: NONE   REQUIRES NNODES=1'
            CALL GSTOP('MNW2 ERROR - OPTION: NONE  REQUIRES NNODES=1')
	    end if
          FVMNW(JFV)%FVMNW2(3,MNWID)=0
        elseif(LOSSTYPE.EQ.'THIEM') then
          FVMNW(JFV)%FVMNW2(3,MNWID)=1 
        elseif(LOSSTYPE.EQ.'SKIN') then
          FVMNW(JFV)%FVMNW2(3,MNWID)=2
        elseif(LOSSTYPE.EQ.'GENERAL') then
          FVMNW(JFV)%FVMNW2(3,MNWID)=3
        elseif(LOSSTYPE.EQ.'SPECIFYCWC') then
          FVMNW(JFV)%FVMNW2(3,MNWID)=4	
	  end if
C    GWM DOES NOT SUPPORT PUMPLOC,QLIMIT OR PUMPCAP
        PUMPLOC=0
        PUMPCAP=0
        Qlimit=0
        FVMNW(JFV)%FVMNW2(11,MNWID)=PUMPLOC
        FVMNW(JFV)%FVMNW2( 6,MNWID)=Qlimit
        FVMNW(JFV)%FVMNW2(22,MNWID)=PUMPCAP
        FVMNW(JFV)%FVMNW2(19,MNWID)=PPFLAG
      END SUBROUTINE ASSIGN3C
C***********************************************************************
      SUBROUTINE ASSIGN3E
        IF(NNODES.GT.0)THEN
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 1,NCNT)=IL             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 2,NCNT)=IR            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 3,NCNT)=IC           
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 4,NCNT)=0.0           
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 5,NCNT)=RwNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 6,NCNT)=RskinNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 7,NCNT)=KskinNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 8,NCNT)=BNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD( 9,NCNT)=CNode            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD(10,NCNT)=PNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD(11,NCNT)=CWCNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWNOD(19,NCNT)=PP  
        ELSE
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 1,NCNT)=Ztop            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 2,NCNT)=Zbotm            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 3,NCNT)=IR            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 4,NCNT)=IC            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 5,NCNT)=RwNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 6,NCNT)=RskinNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 7,NCNT)=KskinNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 8,NCNT)=BNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT( 9,NCNT)=CNode            
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT(10,NCNT)=PNode             
          FVMNW(JFV)%FVCELLS(MNWID)%FVMNWINT(11,NCNT)=CWCNode             
        ENDIF   
             
      END SUBROUTINE ASSIGN3E
C***********************************************************************
      SUBROUTINE READSET2C
c     read Data Set 2c, depending on LOSSTYPE (FVMNW2(3,MNWID)
        SELECT CASE (INT(FVMNW(JFV)%FVMNW2(3,MNWID)))
          CASE (1)
            READ(LOCAT,*,ERR=991) Rw
c     don't allow Rw = 0
            if(Rw.eq.0.0) then
              write(iout,*) '***ERROR*** Rw=0.0; Rw read=',Rw
              CALL GSTOP('MNW2 ERROR - Rw')
            endif
          CASE (2)
            READ(LOCAT,*,ERR=991) Rw,Rskin,Kskin
            if(Rw.eq.0.0) then
              write(iout,*) '***ERROR*** Rw=0.0; Rw read=',Rw
              CALL GSTOP('MNW2 ERROR - Rw')
            endif
            if(Rskin.eq.0.0) then
              write(iout,*) '***ERROR*** Rskin=0.0; Rskin read=',Rskin
              CALL GSTOP('MNW2 ERROR')
            endif
            if(Kskin.eq.0.0) then
              write(iout,*) '***ERROR*** Kskin=0.0; Kskin read=',Kskin
              CALL GSTOP('MNW2 ERROR - Kskin')
            endif
          CASE (3)
            READ(LOCAT,*,ERR=991) Rw,B,C,P
            if(Rw.eq.0.0) then
              write(iout,*) '***ERROR*** Rw=0.0; Rw read=',Rw
              CALL GSTOP('MNW2 ERROR - Rw')
            endif
            if(P.gt.0.0.and.(P.lt.1.0.or.P.gt.3.5)) then
              write(iout,*) '***ERROR*** P=',P,' exceeds 1 <= P <=3.5'
              CALL GSTOP('MNW2 ERROR - P')
            endif
          CASE (4)
            READ(LOCAT,*,ERR=991) CWC            
        END SELECT
c     end read Data Set 2c
        RETURN
!      
  991 CONTINUE ! FILE-READING ERROR
      INQUIRE(LOCAT,NAME=FLNM)
      WRITE(IOUT,9910)TRIM(FLNM),LOCAT
 9910 FORMAT(/,1X,'*** ERROR READING FILE "',A,'" ON UNIT ',I5,/,
     & 2X,'-- STOP EXECUTION (GWM1DCV3ARM)')
      CALL GSTOP(' ')
      END SUBROUTINE READSET2C
C***********************************************************************
      SUBROUTINE READSET2D
        IF(NNODES.GT.0) THEN  
          CALL READSET2D1
        ELSE   ! if nnodes<0, read in Ztop and Zbot which define intervals
          CALL READSET2D2
        END IF
      RETURN  
      END SUBROUTINE READSET2D
C***********************************************************************
      SUBROUTINE READSET2D1
c     read Data Set 2D-1 FOR NODE 
c     If PPFLAG=0, don't read PP variable
           PPFLAG=INT(FVMNW(JFV)%FVMNW2(19,mnwid))
           IF(PPFLAG.eq.0) then           
c     access the LOSSTYPE
            SELECT CASE (INT(FVMNW(JFV)%FVMNW2(3,MNWID)))
c     LOSSTYPE=NONE, read IL,IR,IC only
              CASE (0)
                READ(LOCAT,*,ERR=991) IL,IR,IC
c
c     LOSSTYPE=THIEM, read IL,IR,IC,{Rw}
              CASE (1)
                IF(Rw.GT.0.0) THEN            
                  READ(LOCAT,*,ERR=991) IL,IR,IC
c     Rw at each node is the same if Rw>0.0  
                  RwNode=Rw         
                ELSE
c     If Rw<0, read in separate Rw for each node
                  READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode
                END IF
c
c     LOSSTYPE=SKIN, read IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
                IF(Rw.GT.0.0) THEN            
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC
                      RwNode=Rw             
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991) IL,IR,IC,KskinNode
                      RwNode=Rw             
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RskinNode
                      RwNode=Rw             
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RskinNode,KskinNode
                      RwNode=Rw             
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,KskinNode
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,RskinNode
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991)
     &                  IL,IR,IC,RwNode,RskinNode,KskinNode
                    ENDIF
                  ENDIF
                END IF
c
c     LOSSTYPE=GENERAL, read IL,IR,IC,{Rw B C P}
              CASE (3)
                IF(Rw.GT.0.0) THEN            
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,PNode
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,CNode
                        RwNode=Rw             
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,CNode,PNode
                        RwNode=Rw             
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode
                        RwNode=Rw             
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode,PNode
                        RwNode=Rw             
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode,CNode
                        RwNode=Rw             
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode,CNode,PNode
                        RwNode=Rw             
                      END IF
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,Rwnode,PNode
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,Rwnode,CNode
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)IL,IR,IC,Rwnode,CNode,PNode
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,Rwnode,BNode
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)IL,IR,IC,Rwnode,BNode,PNode
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991)IL,IR,IC,Rwnode,BNode,CNode
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)
     &                    IL,IR,IC,Rwnode,BNode,CNode,PNode
                      END IF
                    ENDIF
                  ENDIF
                END IF
c
c     LOSSTYPE=SPECIFYcwc, read IL,IR,IC,{CWC}
              CASE (4)
                IF(CWC.GT.0.0) THEN            
                  READ(LOCAT,*,ERR=991) IL,IR,IC
                  CWCNode=CWC         
                ELSE
                  READ(LOCAT,*,ERR=991) IL,IR,IC,CWCNode
                END IF
            END SELECT
c    ELSE if PPFLAG NE 0, read PP flag
           ELSE
c     access the LOSSTYPE
            SELECT CASE (INT(FVMNW(JFV)%FVMNW2(3,MNWID)))
c     LOSSTYPE=NONE, read IL,IR,IC only
              CASE (0)
                READ(LOCAT,*,ERR=991) IL,IR,IC,PP
c
c     LOSSTYPE=THIEM, read IL,IR,IC,{Rw}
              CASE (1)
                IF(Rw.GT.0.0) THEN           
                  READ(LOCAT,*,ERR=991) IL,IR,IC,PP
c     Rw at each node is the same if Rw>0.0  
                  RwNode=Rw         
                ELSE
c     If Rw<0, read in separate Rw for each node
                  READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,PP
                END IF
c
c     LOSSTYPE=SKIN, read IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
                IF(Rw.GT.0.0) THEN            
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,PP
                      RwNode=Rw             
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991) IL,IR,IC,KskinNode,PP
                      RwNode=Rw             
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RskinNode,PP
                      RwNode=Rw             
                      KskinNode=Kskin
                    ELSE             
                    READ(LOCAT,*,ERR=991)IL,IR,IC,RskinNode,KskinNode,PP
                      RwNode=Rw             
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,PP
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,KskinNode,PP
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,RskinNode,PP
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991)
     &                 IL,IR,IC,RwNode,RskinNode,KskinNode,PP
                    ENDIF
                  ENDIF
                END IF
c
c     LOSSTYPE=GENERAL, read IL,IR,IC,{Rw B C P}
              CASE (3)
                IF(Rw.GT.0.0) THEN            
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,PP
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,PNode,PP
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,CNode,PP
                        RwNode=Rw             
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,CNode,PNode,PP
                        RwNode=Rw             
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode,PP
                        RwNode=Rw             
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode,PNode,PP
                        RwNode=Rw             
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,BNode,CNode,PP
                        RwNode=Rw             
                        PNode=P
                      ELSE
                      READ(LOCAT,*,ERR=991)IL,IR,IC,BNode,CNode,PNode,PP
                        RwNode=Rw             
                      END IF
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,RwNode,PP
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) IL,IR,IC,Rwnode,PNode,PP
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,Rwnode,CNode,PP
                        BNode=B            
                        PNode=P
                      ELSE
                     READ(LOCAT,*,ERR=991)IL,IR,IC,Rwnode,CNode,PNode,PP
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) IL,IR,IC,Rwnode,BNode,PP
                        CNode=C
                        PNode=P
                      ELSE
                     READ(LOCAT,*,ERR=991)IL,IR,IC,Rwnode,BNode,PNode,PP
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                     READ(LOCAT,*,ERR=991)IL,IR,IC,Rwnode,BNode,CNode,PP
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)
     &                    IL,IR,IC,Rwnode,BNode,CNode,PNode,PP
                      END IF
                    ENDIF
                  ENDIF
                END IF
c
c     LOSSTYPE=SPECIFYcwc, read IL,IR,IC,{CWC}
              CASE (4)
                IF(CWC.GT.0.0) THEN            
                  READ(LOCAT,*,ERR=991) IL,IR,IC,PP
                  CWCNode=CWC         
                ELSE
                  READ(LOCAT,*,ERR=991) IL,IR,IC,CWCNode,PP
                END IF
            END SELECT            
           END IF
         RETURN
!      
  991   CONTINUE ! FILE-READING ERROR
        INQUIRE(LOCAT,NAME=FLNM)
        WRITE(IOUT,9910)TRIM(FLNM),LOCAT
 9910   FORMAT(/,1X,'*** ERROR READING FILE "',A,'" ON UNIT ',I5,/,
     &   2X,'-- STOP EXECUTION (GWM1DCV3ARM)')
        CALL GSTOP(' ')
        END SUBROUTINE READSET2D1
C***********************************************************************
        SUBROUTINE READSET2D2
c     read Data Set 2D-2 FOR INTERVAL
c     access the LOSSTYPE
            SELECT CASE (INT(FVMNW(JFV)%FVMNW2(3,MNWID)))
c     LOSSTYPE=NONE, read Ztop,Zbotm,IR,IC only
              CASE (0)
                READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC
c
c     LOSSTYPE=THIEM, read Ztop,Zbotm,IR,IC,{Rw}
              CASE (1)
                IF(Rw.GT.0.0) THEN            
                  READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC
c     Rw at each node is the same if Rw>0.0  
                  RwNode=Rw         
                ELSE
c     If Rw<0, read in separate Rw for each node
                  READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,RwNode
                END IF
c
c     LOSSTYPE=SKIN, read Ztop,Zbotm,IR,IC,{Rw Rskin Kskin}
              CASE (2)
                IF(Rw.GT.0.0) THEN            
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC
                      RwNode=Rw             
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,KskinNode
                      RwNode=Rw             
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,RskinNode
                      RwNode=Rw             
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991)
     &                  Ztop,Zbotm,IR,IC,RskinNode,KskinNode
                      RwNode=Rw             
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,RwNode
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991)
     &                  Ztop,Zbotm,IR,IC,RwNode,KskinNode
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(LOCAT,*,ERR=991)
     &                  Ztop,Zbotm,IR,IC,RwNode,RskinNode
                        KskinNode=Kskin
                    ELSE             
                      READ(LOCAT,*,ERR=991)
     &                  Ztop,Zbotm,IR,IC,RwNode,RskinNode,KskinNode
                    ENDIF
                  ENDIF
                END IF
c
c     LOSSTYPE=GENERAL, read Ztop,Zbotm,IR,IC,{Rw B C P}
              CASE (3)
                IF(Rw.GT.0.0) THEN            
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,PNode
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,CNode
                        RwNode=Rw             
                        BNode=B            
                        PNode=P
                      ELSE
                       READ(LOCAT,*,ERR=991)Ztop,Zbotm,IR,IC,CNode,PNode
                        RwNode=Rw             
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,BNode
                        RwNode=Rw             
                        CNode=C
                        PNode=P
                      ELSE
                       READ(LOCAT,*,ERR=991)Ztop,Zbotm,IR,IC,BNode,PNode
                        RwNode=Rw             
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                       READ(LOCAT,*,ERR=991)Ztop,Zbotm,IR,IC,BNode,CNode
                        RwNode=Rw             
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)
     &                    Ztop,Zbotm,IR,IC,BNode,CNode,PNode
                        RwNode=Rw             
                      END IF
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,RwNode
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                      READ(LOCAT,*,ERR=991)Ztop,Zbotm,IR,IC,Rwnode,PNode
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                      READ(LOCAT,*,ERR=991)Ztop,Zbotm,IR,IC,Rwnode,CNode
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991) 
     &                   Ztop,Zbotm,IR,IC,Rwnode,CNode,PNode
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                      READ(LOCAT,*,ERR=991)Ztop,Zbotm,IR,IC,Rwnode,BNode
                        CNode=C
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)
     &                    Ztop,Zbotm,IR,IC,Rwnode,BNode,PNode
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(LOCAT,*,ERR=991)
     &                    Ztop,Zbotm,IR,IC,Rwnode,BNode,CNode
                        PNode=P
                      ELSE
                        READ(LOCAT,*,ERR=991)
     &                    Ztop,Zbotm,IR,IC,Rwnode,BNode,CNode,PNode
                      END IF
                    ENDIF
                  ENDIF
                END IF
c
c     LOSSTYPE=SPECIFYcwc, read Ztop,Zbotm,IR,IC,{CWC}
              CASE (4)
                IF(CWC.GT.0.0) THEN            
                  READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC
                  CWCNode=CWC         
                ELSE
                  READ(LOCAT,*,ERR=991) Ztop,Zbotm,IR,IC,CWCNode
                END IF
            END SELECT
c     Set vars for interval
         RETURN
!      
  991   CONTINUE ! FILE-READING ERROR
        INQUIRE(LOCAT,NAME=FLNM)
        WRITE(IOUT,9910)TRIM(FLNM),LOCAT
 9910   FORMAT(/,1X,'*** ERROR READING FILE "',A,'" ON UNIT ',I5,/,
     &   2X,'-- STOP EXECUTION (GWM1DCV3ARM)')
        CALL GSTOP(' ')
        END SUBROUTINE READSET2D2
C***********************************************************************
      SUBROUTINE FIND_LOSSTYPE2(INDEX,LTYPE)
      INTEGER(I4B),INTENT(IN)::INDEX
      CHARACTER(LEN=20),INTENT(OUT)::LTYPE
        IF    (INDEX.EQ.0)THEN
          LTYPE='NONE'
        ELSEIF(INDEX.EQ.1)THEN
          LTYPE='THIEM'
        ELSEIF(INDEX.EQ.2)THEN
          LTYPE='SKIN'
        ELSEIF(INDEX.EQ.3)THEN
          LTYPE='GENERAL'
        ELSEIF(INDEX.EQ.4)THEN
          LTYPE='SPECIFYCWC'
	  ENDIF
      END SUBROUTINE FIND_LOSSTYPE2
 !*****************************************************************************
      SUBROUTINE WRMNW_ECHO(J,K)
       INTEGER(I4B),INTENT(IN)::J,K
       LOSSTYPEINDEX=INT(FVMNW(J)%FVMNW2(3,K))
        CALL FIND_LOSSTYPE2(LOSSTYPEINDEX,LOSSTYPE)
        PPFLAG =FVMNW(J)%FVMNW2(19,K)
        WRITE(IOUT,9010,ERR=990) ABS(NNODES),LOSSTYPE,PPFLAG
        IF(NNODES.GT.0)THEN       ! ITEM 2D-1 
          CALL WRMNW_HEAD1
          DO INODE =1,NNODES            
          IL       =FVMNW(J)%FVCELLS(K)%FVMNWNOD(1,INODE)           
          IR       =FVMNW(J)%FVCELLS(K)%FVMNWNOD(2,INODE)           
          IC       =FVMNW(J)%FVCELLS(K)%FVMNWNOD(3,INODE)           
          RwNode   =FVMNW(J)%FVCELLS(K)%FVMNWNOD(5,INODE)             
          RskinNode=FVMNW(J)%FVCELLS(K)%FVMNWNOD(6,INODE)             
          KskinNode=FVMNW(J)%FVCELLS(K)%FVMNWNOD(7,INODE)             
          BNode    =FVMNW(J)%FVCELLS(K)%FVMNWNOD(8,INODE)            
          CNode    =FVMNW(J)%FVCELLS(K)%FVMNWNOD(9,INODE)          
          PNode    =FVMNW(J)%FVCELLS(K)%FVMNWNOD(10,INODE)           
          CWCNode  =FVMNW(J)%FVCELLS(K)%FVMNWNOD(11,INODE)             
          PP       =FVMNW(J)%FVCELLS(K)%FVMNWNOD(19,INODE) 
          CALL WRMNW_DATA1
          ENDDO
        ELSE                     ! ITEM 2D-2
          CALL WRMNW_HEAD2
          DO IINT =1,ABS(NNODES)       
          Ztop     =FVMNW(J)%FVCELLS(K)%FVMNWINT(1,IINT)           
          Zbotm    =FVMNW(J)%FVCELLS(K)%FVMNWINT(2,IINT)          
          IR       =FVMNW(J)%FVCELLS(K)%FVMNWINT(3,IINT)           
          IC       =FVMNW(J)%FVCELLS(K)%FVMNWINT(4,IINT)           
          RwNode   =FVMNW(J)%FVCELLS(K)%FVMNWINT(5,IINT)           
          RskinNode=FVMNW(J)%FVCELLS(K)%FVMNWINT(6,IINT)             
          KskinNode=FVMNW(J)%FVCELLS(K)%FVMNWINT(7,IINT)            
          BNode    =FVMNW(J)%FVCELLS(K)%FVMNWINT(8,IINT)            
          CNode    =FVMNW(J)%FVCELLS(K)%FVMNWINT(9,IINT)            
          PNode    =FVMNW(J)%FVCELLS(K)%FVMNWINT(10,IINT)            
          CWCNode  =FVMNW(J)%FVCELLS(K)%FVMNWINT(11,IINT)            
          CALL WRMNW_DATA2
          ENDDO
        ENDIF 
      RETURN
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM)
      WRITE(*,9900)TRIM(FLNM),IOUT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')

 9010 FORMAT(T12,'NNODES:',I5,
     &         T25,'LOSSTYPE: ',A20,T56,'PPFLAG:',I4)
            
      END SUBROUTINE WRMNW_ECHO

!*****************************************************************************
      SUBROUTINE WRMNW_HEAD1
c
      if(PPFLAG.EQ.0) then
c     write header depending on LOSSTYPE
        SELECT CASE (LOSSTYPEINDEX)
          CASE (0)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col'
          CASE (1)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col      Rw     '
          CASE (2)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col      Rw     ',
     &'Rskin     Kskin'
          CASE (3)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col      Rw     B
     &         C          P  '
          CASE (4)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col    spec.CWC'
        END SELECT
c     If PPFLAG>0 print PP input
      ELSE
c     write header depending on LOSSTYPE
        SELECT CASE (LOSSTYPEINDEX)
          CASE (0)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col        PP'
          CASE (1)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col      Rw        PP'
          CASE (2)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col      Rw     ',
     &'Rskin     Kskin        PP'
          CASE (3)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col      Rw     B
     &         C          P        PP'
          CASE (4)
      write(iout,9000,ERR=990) ' Node  Lay  Row  Col    spec.CWC  ',
     &'      PP'
        END SELECT
      ENDIF
c
      RETURN
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM)
      WRITE(*,9900)TRIM(FLNM),IOUT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')

 9000 FORMAT(T10,100A)
      END SUBROUTINE WRMNW_HEAD1
!*****************************************************************************
      SUBROUTINE WRMNW_HEAD2
c     write header depending on LOSSTYPE
            SELECT CASE (LOSSTYPEINDEX)
              CASE (0)
            write(iout,9000,ERR=990) ' Interval      Ztop        ',
     &'Zbotm     Row  Col'
              CASE (1)
            write(iout,9000,ERR=990) ' Interval      Ztop        ',
     &'Zbotm     Row  Col      Rw     '
              CASE (2)
            write(iout,9000,ERR=990) ' Interval      Ztop       ',
     &'Zbotm      Row  Col      Rw     Rskin    ',
     &' Kskin '
              CASE (3)
            write(iout,9000,ERR=990) ' Interval      Ztop       ',
     &'Zbotm      Row  Col      Rw     B         C         P  '
              CASE (4)
            write(iout,9000,ERR=990) ' Interval      Ztop       ',
     &'Zbotm      Row  Col      spec.CWC'
            END SELECT
      RETURN
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM)
      WRITE(*,9900)TRIM(FLNM),IOUT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')

 9000 FORMAT(T10,100A)
      END SUBROUTINE WRMNW_HEAD2
!*****************************************************************************
      SUBROUTINE WRMNW_DATA1
c
      if(PPFLAG.EQ.0) then
c     write data depending on LOSSTYPE
            SELECT CASE (LOSSTYPEINDEX)
              CASE (0)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4)')
     &          INODE,IL,IR,IC      
              CASE (1)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P1G10.4)')
     &          INODE,IL,IR,IC,RwNode       
              CASE (2)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P3G10.4)')
     &          INODE,IL,IR,IC,RwNode,RskinNode,KskinNode       
              CASE (3)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P4G10.4)')
     &          INODE,IL,IR,IC,RwNode,BNode,CNode,PNode     
              CASE (4)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P1G10.4)')
     &          INODE,IL,IR,IC,CWCNode     
            END SELECT
c     If PPFLAG>0 print PP input
      else
            SELECT CASE (LOSSTYPEINDEX)
              CASE (0)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,G10.3)')
     &          INODE,IL,IR,IC,PP       
              CASE (1)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P2G10.3)')
     &          INODE,IL,IR,IC,RwNode,PP       
              CASE (2)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P4G10.4)')
     &          INODE,IL,IR,IC,RwNode,RskinNode,KskinNode,PP     
              CASE (3)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P5G10.4)')
     &          INODE,IL,IR,IC,RwNode,BNode,CNode,PNode,PP      
              CASE (4)
          write(iout,'(T11,I4,1x,I4,1x,I4,1x,I4,2x,1P2G10.4)')
     &          INODE,IL,IR,IC,CWCNode,PP     
            END SELECT
      end if      
      RETURN
      END SUBROUTINE WRMNW_DATA1
!*****************************************************************************
      SUBROUTINE WRMNW_DATA2
c     write data depending on LOSSTYPE
        SELECT CASE (LOSSTYPEINDEX)
          CASE (0)
            write(iout,'(T11,I4,6x,1P2G12.5,1x,I4,1x,I4)')
     &            IINT,Ztop,Zbotm,IR,IC 
          CASE (1)
            write(iout,'(T11,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P1G10.4)')
     &            IINT,Ztop,Zbotm,IR,IC,RwNode
          CASE (2)
            write(iout,'(T11,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P3G10.4)')
     &            IINT,Ztop,Zbotm,IR,IC,RwNode,RskinNode,KskinNode
          CASE (3)
            write(iout,'(T11,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P4G10.4)')
     &            IINT,Ztop,Zbotm,IR,IC,RwNode,BNode,CNode,PNode
          CASE (4)
            write(iout,'(T11,I4,6x,1P2G12.5,1x,I4,1x,I4,4x,1P1G10.4)')
     &            IINT,Ztop,Zbotm,IR,IC,CWCNode
        END SELECT
      RETURN
      END SUBROUTINE WRMNW_DATA2

      END SUBROUTINE GWM1DCV3AR

C
C***********************************************************************
      SUBROUTINE GWM1DCV3FM
C***********************************************************************
C     VERSION: 21MAR2008
C     PURPOSE: FORMULATE THE BOUNDS ON DECISION VARIABLES
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE,SMALLEPS
      USE GWM1RMS3, ONLY : BNDS
C-----LOCAL VARIABLES
      INTEGER(I4B)::I
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO 100 I=1,NFVAR
        IF(FVON(I).GT.0)THEN
          BNDS(I) = FVMAX(I)                     ! LOAD FLOW VARIABLE BOUNDS 
        ELSEIF(FVON(I).LT.0)THEN
          BNDS(I) = SQRT(SMALLEPS)              ! MAKE VARIABLE UNAVAILABLE 
        ENDIF
  100 ENDDO
      DO 200 I=1,NEVAR
        BNDS(NFVAR+I) = EVMAX(I)                   ! LOAD EXTERNAL VARIABLE BOUNDS 
  200 ENDDO
      DO 300 I=NFVAR+NEVAR+1,NFVAR+NEVAR+NBVAR
        BNDS(I) = ONE                            ! LOAD BINARY VARIABLE BOUNDS
  300 ENDDO
C
      RETURN
      END SUBROUTINE GWM1DCV3FM
C
C***********************************************************************
      SUBROUTINE GWM1DCV3FVCSV(J)
      INTEGER(I4B),INTENT(IN)::J
        FVCELL(J)%FVILOC   => FVILOC
        FVCELL(J)%FVJLOC   => FVJLOC
        FVCELL(J)%FVKLOC   => FVKLOC
        FVCELL(J)%FVRATIO  => FVRATIO
      RETURN
      END SUBROUTINE GWM1DCV3FVCSV
C***********************************************************************
      SUBROUTINE GWM1DCV3FVCPNT(J)
      INTEGER(I4B),INTENT(IN)::J
         FVILOC => FVCELL(J)%FVILOC  
         FVJLOC => FVCELL(J)%FVJLOC  
         FVKLOC => FVCELL(J)%FVKLOC  
         FVRATIO=> FVCELL(J)%FVRATIO   
      RETURN
      END SUBROUTINE GWM1DCV3FVCPNT
C
C***********************************************************************
      SUBROUTINE GWF2DCV3FM(KPER,IGRID)
C***********************************************************************
C     VERSION: 21MAR2008
C     PURPOSE: FORMULATE THE MANAGED FLOW VARIABLES INTO THE SIMULATION RHS
C-----------------------------------------------------------------------
      USE GLOBAL,  ONLY: IBOUND,RHS
      INTEGER(I4B),INTENT(IN)::KPER,IGRID
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,K,IL,IR,IC
      REAL(SP)::Q 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      DO 100 I=1,NFVAR                      ! LOOP OVER GWM FLOW VARIABLES
        IF(IGRID.EQ.GRDLOCDCV(I))THEN       ! FLOW VARIABLE ON CURRENT GRID
        IF(FVSP(I,KPER)) THEN               ! FLOW VARIABLE ACTIVE IN STRESS PERIOD
        IF(FVNCELL(I).GT.0)THEN             ! THIS IS A WEL-TYPE FLOW VARIABLE
          CALL GWM1DCV3FVCPNT(I)            ! POINT TO CORRECT CELL INFO
          DO 110 K=1,FVNCELL(I)             ! LOOP OVER CELLS FOR THIS VARIABLE
            IL = FVKLOC(K)                  ! ASSIGN CELL LAYER
            IR = FVILOC(K)                  ! ASSIGN CELL ROW
            IC = FVJLOC(K)                  ! ASSIGN CELL COLUMN
            Q  = FVBASE(I)*FVRATIO(K)       ! ASSIGN FLOW RATE 
            IF(IBOUND(IC,IR,IL).GT.0)THEN   ! CELL IS ACTIVE 
              RHS(IC,IR,IL)=RHS(IC,IR,IL)-Q ! SUBTRACT FROM THE RHS ACCUMULATOR
            ENDIF  
  110     ENDDO
        ENDIF
        ENDIF
        ENDIF
  100 ENDDO
C
      RETURN
      END SUBROUTINE GWF2DCV3FM
C
C***********************************************************************
      SUBROUTINE GWF2DCV3BD(KSTP,KPER,IGRID)
C***********************************************************************
C     VERSION: 21MAR2008
C     PURPOSE: CALCULATE VOLUMETRIC BUDGET FOR THE MANAGED FLOW VARIABLES 
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY: ZERO 
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IBOUND,BUFF
      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,DELT,PERTIM,TOTIM,VBVL,VBNM
      USE GWFWELMODULE,ONLY:IWELCB
      INTEGER(I4B),INTENT(IN)::KSTP,KPER,IGRID
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,K,IL,IR,IC
      REAL(SP)::Q 
      REAL(DP)::RATIN,RATOUT,QQ,RIN,ROUT 
      CHARACTER(LEN=16)::TEXT='   MANAGED WELLS' 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
CGWM TO AVOID EXCESSIVE OUTPUT CELL-BY-CELL WRITING NOT IMPLEMENTED
C1------CLEAR RATIN AND RATOUT ACCUMULATORS, AND SET CELL-BY-CELL
C1------BUDGET FLAG.
      RATIN=ZERO
      RATOUT=ZERO
C      IBD=0
C      IF(IWELCB.LT.0 .AND. ICBCFL.NE.0) IBD=-1
C      IF(IWELCB.GT.0) IBD=ICBCFL
C      IBDLBL=0
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
C      IF(IBD.EQ.2) CALL UBDSV2(KSTP,KPER,TEXT,IWELCB,NCOL,NROW,NLAY,
C     1          NWELLS,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C
C3------CLEAR THE BUFFER.
   50 BUFF=ZERO
C
C4------IF THERE ARE NO WELLS, DO NOT ACCUMULATE FLOW.
      IF(NFVAR.EQ.0)GOTO 200
C
C5------LOOP THROUGH EACH WELL CALCULATING FLOW.
      DO 100 I=1,NFVAR                      ! LOOP OVER GWM FLOW VARIABLES
        IF(IGRID.EQ.GRDLOCDCV(I))THEN       ! FLOW VARIABLE IS ON CURRENT GRID
          IF(FVSP(I,KPER)) THEN             ! FLOW VARIABLE ACTIVE IN STRESS PERIOD
          IF(FVNCELL(I).GT.0)THEN           ! THIS IS A WEL-TYPE FLOW VARIABLE
            CALL GWM1DCV3FVCPNT(I)          ! POINT TO CORRECT CELL INFO
            DO 110 K=1,FVNCELL(I)           ! LOOP OVER CELLS FOR THIS VARIABLE
C
C5A-----GET LAYER, ROW & COLUMN OF CELL CONTAINING WELL.
              IL = FVKLOC(K)                ! ASSIGN CELL LAYER
              IR = FVILOC(K)                ! ASSIGN CELL ROW
              IC = FVJLOC(K)                ! ASSIGN CELL COLUMN
              Q=ZERO
C
C5B-----IF THE CELL IS NO-FLOW OR CONSTANT_HEAD, IGNORE IT.
              IF(IBOUND(IC,IR,IL).LE.0)GOTO 99
C
C5C-----GET FLOW RATE FROM WELL LIST.
              Q  = FVBASE(I)*FVRATIO(K)     ! ASSIGN FLOW RATE 
              QQ=Q
C
C5D-----PRINT FLOW RATE IF REQUESTED.
C      IF(IBD.LT.0) THEN
C         IF(IBDLBL.EQ.0) WRITE(IOUT,61) TEXT,KPER,KSTP
C   61    FORMAT(1X,/1X,A,'   PERIOD',I3,'   STEP',I3)
C         WRITE(IOUT,62) L,IL,IR,IC,Q
C   62    FORMAT(1X,'WELL',I4,'   LAYER',I3,'   ROW',I4,'   COL',I4,
C     1       '   RATE',1PG15.6)
C         IBDLBL=1
C      END IF
C
C5E-----ADD FLOW RATE TO BUFFER.
              BUFF(IC,IR,IL)=BUFF(IC,IR,IL)+Q
C
C5F-----SEE IF FLOW IS POSITIVE OR NEGATIVE.
              IF(Q) 90,99,80
C
C5G-----FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
   80         RATIN=RATIN+QQ
              GOTO 99
C
C5H-----FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
   90         RATOUT=RATOUT-QQ
C
C5I-----IF CELL-BY-CELL FLOWS ARE BEING SAVED AS A LIST, WRITE FLOW.
C5I-----OR IF RETURNING THE FLOW IN THE WELL ARRAY, COPY FLOW TO WELL.
   99         CONTINUE
CGWM   IF(IBD.EQ.2) CALL UBDSVA(IWELCB,NCOL,NROW,IC,IR,IL,Q,IBOUND,NLAY)
  110       ENDDO
          ENDIF
          ENDIF
        ENDIF
  100 ENDDO
C
C6------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A 3-D ARRAY,
C6------CALL UBUDSV TO SAVE THEM.
CGWM      IF(IBD.EQ.1) CALL UBUDSV(KSTP,KPER,TEXT,IWELCB,BUFF,NCOL,NROW,
CGWM     1                          NLAY,IOUT)
C
C7------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
  200 RIN=RATIN
      ROUT=RATOUT
      VBVL(3,MSUM)=RIN
      VBVL(4,MSUM)=ROUT
      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
      VBNM(MSUM)=TEXT
C
C8------INCREMENT BUDGET TERM COUNTER(MSUM).
      MSUM=MSUM+1
C
C9------RETURN
      RETURN
      END SUBROUTINE GWF2DCV3BD 

      subroutine gwf2dcv3rp(kper,igrid)
      ! Write managed-flow rates for current stress period.
      ! Logic is as in GWF2DCV3FM, except RHS is not modified.
      use global,  only: iout
      implicit none
      integer(i4b),intent(in)::kper,igrid
c-----local variables
      integer(i4b)::i,k,il,ir,ic,n
      real(sp)::q 
      double precision :: qtot
   10 format(/,'Flows at Managed Wells for Stress Period ',i4)
   20 format(/,1x,'Well no.  Layer   Row   Col   Stress Rate')
   25 format(1x,'----------------------------------------------------')
   30 format(1x,i8,1x,i6,1x,i5,1x,i5,1x,g24.16)
   40 format(22x,'Total:',1x,g24.16)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      write(iout,10)kper
      write(iout,20)
      write(iout,25)
      n = 0
      qtot = 0.0d0
      do 100 i=1,nfvar                      ! loop over gwm flow variables
        if(igrid.eq.grdlocdcv(i))then       ! flow variable on current grid
        if(fvsp(i,kper)) then               ! flow variable active in stress period
        if(fvncell(i).gt.0)then             ! this is a wel-type flow variable
          CALL GWM1DCV3FVCPNT(I)            ! POINT TO CORRECT CELL INFO
          do 110 k=1,fvncell(i)             ! loop over cells for this variable
            il = fvkloc(k)                ! assign cell layer
            ir = fviloc(K)                ! assign cell row
            ic = fvjloc(k)                ! assign cell column
            q  = fvbase(i)*fvratio(k)     ! assign flow rate 
            qtot = qtot + q
            n = n + 1
            write(iout,30)n,il,ir,ic,q
  110     enddo
        endif
        endif
        endif
  100 enddo
      write(iout,25)
      write(iout,40)qtot
c
      return
      end subroutine gwf2dcv3rp
C***********************************************************************      
      INTEGER FUNCTION GET_NEXT_MNW2_AUXVALUE()
        IMPLICIT NONE
        KMNWAUXVAL = KMNWAUXVAL + 1
        GET_NEXT_MNW2_AUXVALUE = KMNWAUXVAL
        RETURN
      END FUNCTION GET_NEXT_MNW2_AUXVALUE
C***********************************************************************
      END MODULE GWM1DCV3
