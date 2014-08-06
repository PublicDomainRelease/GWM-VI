!    Structure of file gwm1dcv3mnw.for
!
!    MODULE GWM1DCV3MNW
!      define storage for unmanaged MNW data 
!      CONTAINS
!      SUBROUTINE GWM1MNW27AR
!        CONTAINS
!          routines related to reading original MNW2 input file
!      END SUBROUTINE GWM1DCV3AR
!      SUBROUTINE GWM1DCV3ARW
!        CONTAINS
!          routines related to writing new MWN2 input file with managed and unmanaged flows
!      END SUBROUTINE GWM1DCV3ARW
!      Utility routines 
!    END MODULE GWM1DCV3MNW
!
      MODULE GWM1DCV3MNW
C     VERSION: 27AUG2013
      USE GWM1BAS3, ONLY: GWMOUT
      USE GWM1DCV3, ONLY: MNW2HERE
      USE GWM_SUBS, ONLY: IGETUNIT
      USE GWM_STOP, ONLY: GSTOP

      IMPLICIT NONE
      PRIVATE 
      PUBLIC  GWM1DCV3MNW2
      ! GWM-VI needs to be able to set NAUX=1  
      PUBLIC  NAUX_EXTERNAL
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
C
C-----VARIABLES FOR STORING UNMANAGED MNW WELL DATA - MIMICS GWFMNW2MODULE
C       This data read from the MNW file and stored for repeated rewriting of MNW file
C
C     Create storage for MNW2 well data at each stress period on each grid
c       This will store ITEMS 3 and 4 from MNW2 input file
      INTEGER,          SAVE,POINTER::ITMP
      INTEGER,          SAVE,POINTER::Qlimit,QCUT,WELLINDX
      DOUBLE PRECISION, SAVE,POINTER::Qdes,Qfrcmn,Qfrcmx
      DOUBLE PRECISION, SAVE,POINTER::CapMult,Cprime,Hlim
      ! CHARACTER(LEN=16),SAVE,DIMENSION(:),POINTER::GMNWAUX3
      CHARACTER(LEN=20),SAVE,POINTER::WELLID
      TYPE GTWELL
        INTEGER,POINTER           ::Qlimit,QCUT,WELLINDX
        DOUBLE PRECISION,POINTER  ::Qdes,Qfrcmn,Qfrcmx
        DOUBLE PRECISION,POINTER  ::CapMult,Cprime,Hlim
       ! CHARACTER(LEN=16),DIMENSION(:),POINTER ::GMNWAUX3
      END TYPE GTWELL
      TYPE GRID_TIME
        INTEGER,POINTER :: ITMP
        TYPE(GTWELL),DIMENSION(:), POINTER   ::GTW_TARGET
      END TYPE GRID_TIME
C
C     Create storage for each MNW2 well, mimic MNW2 storage structure       
c       This will store ITEMS 1 and 2 from MNW2 input file
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), POINTER     ::GMNW2
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), POINTER     ::GMNWNOD
      DOUBLE PRECISION, SAVE, DIMENSION(:,:), POINTER     ::GMNWINT
      DOUBLE PRECISION, SAVE, DIMENSION(:,:,:),POINTER    ::GCapTable
      CHARACTER(LEN=20),SAVE, DIMENSION(:),   POINTER     ::GWELLID
      CHARACTER(LEN=10) :: AUXVARNAME = 'MNWSEQNUM '
      INTEGER(I4B) :: NAUX_EXTERNAL = 0
      ! CHARACTER(LEN=16),SAVE, DIMENSION(:),   POINTER     ::GMNWAUX
      INTEGER,SAVE,POINTER  ::GMNWMAX,GNODTOT,GIWL2CB,GMNWPRNT
      TYPE GWMMNWTYPE
        INTEGER,POINTER  ::GMNWMAX,GNODTOT,GIWL2CB,GMNWPRNT
        DOUBLE PRECISION,  DIMENSION(:,:), POINTER    ::GMNW2_TARGET
        DOUBLE PRECISION,  DIMENSION(:,:), POINTER    ::GMNWNOD_TARGET
        DOUBLE PRECISION,  DIMENSION(:,:), POINTER    ::GMNWINT_TARGET
        DOUBLE PRECISION,  DIMENSION(:,:,:), POINTER  ::GCapTable_TARGET
        CHARACTER(LEN=20), DIMENSION(:),   POINTER    ::GWELLID_TARGET
        ! CHARACTER(LEN=16), DIMENSION(:), POINTER      ::GMNWAUX_TARGET
        TYPE(GRID_TIME),DIMENSION(:), POINTER   ::GT_TARGET
      END TYPE
      TYPE(GWMMNWTYPE), SAVE:: GWMMNWDAT(10)
C     The allocatable, derived-type array GWMMNWTYPE holds MNW information for each 
C       of ten possible LGR-grids.  
C
C	For the Ith grid, six storage arrays are defined that have the structure 
C       used in GWF2MNW2MODULE and one derived-type array that holds stress 
C       period information.
C		GWMMNWTYPE(I)%GT_TARGET(NPER)
C			Where NPER is the number of stress periods
C
C		For the Tth stress period on the Ith grid a derived-type array is defined.
C		    GWMMNWTYPE(I)%GT_TARGET(T)%GTW_TARGET(NWELLS)
C				Where NWELLs is the number of wells active in stress period T
C               on grid I.
C
C		    For the Kth well at the Tth stress period on the Ith grid a set of
C             scalar data is defined.  For example,
C		        GWMMNWTYPE(I)%GT_TARGET(T)%GTW_TARGET(K)%Qdes
C                    is the value Qdes for well K, stress period T, grid I.
C
      CONTAINS 
C      
C*************************************************************************
      SUBROUTINE GWM1DCV3MNW2(MNW2UNIT,FIRSTSIM,Iusip,
     1                       Iude4,Iupcg,Iugmg,VIFLG,MNWJTF)
C*************************************************************************
C     VERSION: 27AUG2013
C     PURPOSE: READ EXISTING MNW FILE, WRITE A NEW MNW FILE
C---------------------------------------------------------------------------
C
      INTEGER(I4B),INTENT(IN)::MNW2UNIT,Iusip,Iude4,Iupcg
      INTEGER(I4B),INTENT(IN)::Iugmg,VIFLG
      LOGICAL,INTENT(IN)::FIRSTSIM
      CHARACTER(LEN=2000),INTENT(IN),OPTIONAL::MNWJTF
C-----LOCAL VARIABLES
      CHARACTER(LEN=200)::FNAME
      INTEGER(I4B)      ::ISCRCH,LOCAT
      INTEGER(I4B),SAVE ::IGRID,VIFLGP
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----CHECK IF DECISION VARIABLES INCLUDE MNW-TYPE WELLS
      IF(.NOT.MNW2HERE)RETURN   
      IF(FIRSTSIM.AND.MNW2UNIT.NE.0)THEN
C-------READ THE MNW2 INPUT FILE INTO GWM MEMORY
        ISCRCH  = IGETUNIT(7,1000)       ! Open Scratch file for unneeded output
        LOCAT = ISCRCH
        FNAME = 'SCRATCH'
        OPEN(UNIT=ISCRCH,STATUS='SCRATCH',ACTION='WRITE',ERR=999)
        IF(VIFLG.EQ.0)THEN   ! THIS IS A CALL FROM GWM-VI
          IGRID = 1          ! GWM-VI DOES NOT SUPPORT MULTI-GRID (LGR)
          VIFLGP = 0
          FNAME=MNWJTF       ! Write a Jupiter Template File for GWM-VI
        ELSE
          IGRID  = VIFLG
          VIFLGP = 1
          FNAME='GWM_MNW_Combined_Input.txt'
        ENDIF
C-------READ ORIGINAL MNW2 INPUT FILE
        CALL GWM2MNW27AR(MNW2UNIT,Iusip,Iude4,Iupcg,Iugmg,IGRID,ISCRCH)
C-------USE THE DATA TO COMPLETE CONSTRUCTION OF DATA STORAGE FOR MANAGED MNW WELLS
        CALL GWM1DCV3ARS(IGRID) 
C-------CLOSE THE MNW2 INPUT FILE TO END LINK TO THE ORIGINAL MNW FILE
        CLOSE(MNW2UNIT,STATUS='KEEP')
C-------OPEN NEW FILE ON ORIGINAL MNW2 INPUT FILE NUMBER
        LOCAT = MNW2UNIT
        OPEN(UNIT=MNW2UNIT,FILE=FNAME,STATUS='REPLACE',ERR=999)
        IF(VIFLGP.EQ.0)WRITE(MNW2UNIT,9100)
      ELSEIF(FIRSTSIM.AND.MNW2UNIT.EQ.0)THEN
C-------ERROR - THERE MUST BE A MNW2 INPUT FILE IN THE NAM FILE
        WRITE(GWMOUT,9000)   
        CALL GSTOP(' ')
      ELSEIF(.NOT.FIRSTSIM.AND.MNW2UNIT.NE.0)THEN
C-------REWRITE THE MNW2 FILE AGAIN
        REWIND(UNIT=MNW2UNIT)
      ENDIF 
C
C-----WRITE THE NEW MNW2 FILE
      CALL GWM1DCV3ARW(IGRID,MNW2UNIT,VIFLGP)
      REWIND(MNW2UNIT)
      RETURN
 9000 FORMAT(/,'PROGRAM STOPPED. AN MNW2 INPUT FILE MUST BE INCLUDED',/,
     &         '  IN THE NAM FILE IF MNW2-TYPE DECISION VARIABLES',/,
     &         '  ARE PRESENT IN THE GWM INPUT')
 9100 FORMAT('jtf %')
 999  CONTINUE ! FILE-OPENING ERROR
      WRITE(*,9990)TRIM(FNAME),LOCAT
 9990 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,
     &2X,'-- STOP EXECUTION (GWM1DCV3MNW2)')
      CALL GSTOP(' ')
      END SUBROUTINE GWM1DCV3MNW2
C
C*************************************************************************
      SUBROUTINE GWM1DCV3ARS(IGRID)
C*************************************************************************
C     VERSION: 20FEB2012
C     PURPOSE: COMPLETE ASSIGNMENT OF INDICES FOR LOCATING MNW DATA
C---------------------------------------------------------------------------
      USE GWM1DCV3
      USE GWM1HDC3,     ONLY:HDCNUM,HDCKLOC,HDCILOC,HDCJLOC
      USE GWM1STA3,     ONLY:NHVAR,SVKLOC,SVILOC,SVJLOC
      INTEGER(I4B),INTENT(IN)::IGRID
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K,MNWID,MNWMAN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL GWM1DCV3MNW2PNT12(IGRID)
      MNWID = ABS(GMNWMAX)
      MNWMAN = 0
      DO 50 I=1,NFVAR                            ! LOOP OVER GWM FLOW VARIABLES
        IF(IGRID.EQ.GRDLOCDCV(I).AND.            ! FLOW VARIABLE ON CURRENT GRID
     &       FVNCELL(I).LT.0)THEN                ! THIS IS A MNW-TYPE FV
          DO K=1,ABS(FVNCELL(I))                 ! LOOP OVER CELLS FOR FV
            MNWID = MNWID + 1
            MNWMAN=MNWMAN + 1
            ! Save sequence number in MNW well list for flow variable/cell
            FVMNW(I)%FVWTOMNW(K)=MNWID           ! SAVE LOCATION IN MNW ARRAYS
            ! Save sequence number in MNW well list for head constraint at MNW well
            LOOP_HEDCON : DO J=1,HDCNUM
              IF(HDCKLOC(J,1).EQ.0)THEN ! THIS IS A CONSTRAINT ON HEAD IN MNW WELL
                IF(HDCILOC(J,1).EQ.MNWMAN)THEN ! MATCH IS FOUND
                  HDCJLOC(J,1) = MNWID ! STORE INDEX TO WELL LOCATION IN FULL LIST
                ENDIF
              ENDIF
            ENDDO LOOP_HEDCON
            ! Save sequence number in MNW well list for head state variable at MNW well
            LOOP_STAVAR : DO J=1,NHVAR
              IF(SVKLOC(J).EQ.0)THEN ! THIS IS A STATE VARIABLE OF HEAD IN MNW WELL
                IF(SVILOC(J).EQ.MNWMAN)THEN ! MATCH IS FOUND
                  SVJLOC(J) = MNWID ! STORE INDEX TO WELL LOCATION IN FULL LIST
                ENDIF
              ENDIF
            ENDDO LOOP_STAVAR
          ENDDO
        ENDIF  
  50  ENDDO
      RETURN
      END SUBROUTINE GWM1DCV3ARS
C
C*************************************************************************
      SUBROUTINE GWM2MNW27AR(IN,Iusip,Iude4,Iupcg,Iugmg,IGRID,ISCRCH)
C*************************************************************************
C     VERSION: 15AUG2013
C     PURPOSE: READ THE MNW2 FILE (USING STRUCTURE OF GWF2MNW27AR)
C              STORE INFORMATION IN GWM STORAGE THAT MIMICS MNW2 STORAGE
C---------------------------------------------------------------------------
      USE GLOBAL,       ONLY:NLAY,NPER,IOUT
      USE GWFMNW2MODULE, ONLY:NMNW2,MNWMAX,NMNWVL,IWL2CB,MNWPRNT,
     1                       NODTOT,INTTOT,MNW2,MNWNOD,MNWINT,
     2                       CapTable,SMALL,NTOTNOD,WELLID
      INTEGER(I4B),INTENT(IN)::IN,Iusip,Iude4,Iupcg
      INTEGER(I4B),INTENT(IN)::Iugmg,IGRID,ISCRCH
C-----LOCAL VARIABLES
      CHARACTER*200 LINE
      INTEGER(I4B)::LLOC,ISTART,ISTOP,KKPER,N,NAUX,ISTAT
      INTEGER(I4B)::MNWID,NODNUM1,NODNUM2,INTFIRST,INTLAST
      INTEGER(I4B)::TMP_NNODES,NEW_NNODES,ICNT,IOUT_ORIG,MNWMAX_S
      REAL  R
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C1------Allocate scalar variables, which makes it possible for multiple
C1------grids to be defined.
      ALLOCATE(NMNW2,MNWMAX,NTOTNOD,IWL2CB,MNWPRNT,NODTOT,INTTOT,SMALL,
     1 NMNWVL,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992 
      NMNW2=0
      NTOTNOD=0
      NODTOT=0
C
C3------READ MAXIMUM NUMBER OF MNW2 WELLS, UNIT OR FLAG FOR
C3------CELL-BY-CELL FLOW TERMS, AND PRINT FLAG
      CALL URDCOM(IN,ISCRCH,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MNWMAX,R,IOUT,IN)
      MNWMAX_S = MNWMAX    ! Store original MNWMAX here
      IF (MNWMAX.LT.0) THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NODTOT,R,IOUT,IN)
        MNWMAX=-MNWMAX
      ENDIF
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IWL2CB,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MNWPRNT,R,IOUT,IN)
C
C4------READ AUXILIARY VARIABLES 
!      ALLOCATE (GMNWAUX(20))
      NAUX=0
!     GWM does not read AUX variables that may be present
!   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
!      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
!     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
!         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
!         IF(NAUX.LT.5) THEN
!            NAUX=NAUX+1
!            GMNWAUX(NAUX)=LINE(ISTART:ISTOP)
!         END IF
!         GO TO 10
!      END IF
C
C5------ALLOCATE SPACE FOR MNW2 ARRAYS.
      NMNWVL=30+NAUX
      IF (NODTOT.EQ.0)NODTOT=(MNWMAX*NLAY)+(10*NLAY)+25
      ALLOCATE (MNW2(NMNWVL,MNWMAX),STAT=ISTAT)
      ALLOCATE (MNWNOD(34,NODTOT),STAT=ISTAT)
      ALLOCATE (MNWINT(11,NODTOT),STAT=ISTAT)
      ALLOCATE (CapTable(mnwmax,27,2),STAT=ISTAT)
      ALLOCATE (WELLID(mnwmax+1),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992 
C
C7------SAVE POINTERS TO DATA 
      CALL SGWF2MNW2PSV(IGRID)
C-----READ THE REST OF THE MNW2 FILE AND FILL UP THE ARRAYS 
      ! GWF2MNW27RP writes to IOUT using the value in module GLOBAL 
      IOUT_ORIG = IOUT                 ! Save the correct value of IOUT
      IOUT = ISCRCH                    ! Reassign IOUT to scratch file
      KKPER=1
      CALL GWF2MNW27RP(IN,KKPER,Iusip,Iude4,0,Iupcg,0,Iugmg,0,IGRID)
      DO MNWID=1,MNWMAX
        TMP_NNODES=MNW2(2,MNWID)
        IF(TMP_NNODES.LT.0)THEN        ! This is a node count for intervals
          NODNUM1 = MNW2(4,MNWID)      ! Index for first node in well
          IF(MNWID.LT.MNWMAX)THEN
            NODNUM2=MNW2(4,MNWID+1)-1  ! Index for last node in well
          ELSEIF(MNWID.EQ.MNWMAX)THEN
            NODNUM2=NTOTNOD            ! Index for last node in well
          ENDIF
          INTFIRST=MNWNOD(12,NODNUM1)  ! First interval in first node of well
          INTLAST =MNWNOD(13,NODNUM2)  ! Last interval in last node of well  
          NEW_NNODES=INTLAST-INTFIRST+1
          MNW2(2,MNWID)=-NEW_NNODES    ! Replace with interval count
        ENDIF
      ENDDO
C-----STORE WELL DATA FOR UNMANAGED WELLS
      MNWMAX = MNWMAX_S    ! Put value with original sign back for storage
      CALL GWM1DCV3MNW2PSV12(IGRID)
      MNWMAX = ABS(MNWMAX) ! Return to the postive value
C-----ALLOCATE SPACE FOR MNW DATA OVER STRESS PERIODS
      ALLOCATE(GWMMNWDAT(IGRID)%GT_TARGET(NPER),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992 
      CALL COPY_ITEM4(1)           ! SAVE ALL DATA AT THIS STRESS PERIOD 
C-----LOOP OVER REMAINING STRESS PERIODS
      DO KKPER=2,NPER
        CALL GWF2MNW27RP(IN,KKPER,Iusip,Iude4,0,Iupcg,0,Iugmg,0,IGRID)
        CALL COPY_ITEM4(KKPER)     ! SAVE ALL DATA AT THIS STRESS PERIOD 
      ENDDO
      CLOSE(UNIT=ISCRCH)               ! Erase unneeded output
      IOUT = IOUT_ORIG                 ! Restore correct value of IOUT
      MNWMAX = MNWMAX_S    ! Put value with original sign back for storage
C
      RETURN
  992 CONTINUE ! ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')
      CONTAINS
C***********************************************************************
      SUBROUTINE COPY_ITEM4(T)
      INTEGER(I4B),INTENT(IN)::T
C     Copy all the data for Item 4 in MNW2 input
      ICNT = 0
      DO MNWID=1,MNWMAX
        IF(MNW2(1,MNWID).NE.0.0)THEN ! This well is active
          ICNT = ICNT+1
        ENDIF
      ENDDO
C-----ALLOCATE SPACE FOR MNW DATA OVER ALL WELLS AT THIS STRESS PERIOD
      ALLOCATE(GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(ICNT),
     & STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992 
      ICNT = 0
      DO MNWID=1,MNWMAX
        IF(MNW2(1,MNWID).NE.0.0)THEN ! This well is active
          ICNT = ICNT+1
C---------ALLOCATE SPACE FOR MNW DATA AT EACH WELL DURING STRESS PERIOD
          ALLOCATE(ITMP,Qlimit,QCUT,Qdes,Qfrcmn,Qfrcmx,
     &          CapMult,Cprime,Hlim,WELLINDX,STAT=ISTAT)
          IF(ISTAT.NE.0)GOTO 992 
          Qdes    =MNW2(5,MNWID)          
          CapMult =MNW2(24,MNWID)
          Cprime  =MNW2(12,MNWID)
          Qlimit  =MNW2(6,MNWID)
          Hlim    =MNW2(7,MNWID)
          QCUT    =MNW2(8,MNWID)
          Qfrcmn  =MNW2(9,MNWID)
          Qfrcmx  =MNW2(10,MNWID)
          WELLINDX = MNWID
          CALL GWM1DCV3MNW2PSV34(IGRID,T,ICNT)
        ENDIF
      ENDDO

      ALLOCATE(ITMP,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992 
      ITMP = ICNT
      CALL GWM1DCV3MNW2PSV34(IGRID,T,0)  ! SAVE ITMP
!
      RETURN
  992 CONTINUE ! ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1DCV3AR)')
      CALL GSTOP(' ')
      END SUBROUTINE COPY_ITEM4

      END SUBROUTINE GWM2MNW27AR
C  
C***********************************************************************
      SUBROUTINE GWM1DCV3ARW(IGRID,IO,VIFLG)
C***********************************************************************
C     VERSION: 15AUG2013
C     PURPOSE: WRITE A NEW MNW FILE COMBINING THE ORIGINAL MNW FILE
C              AND THE MNW-TYPE FLOW VARIABLES FROM GWM
C-----------------------------------------------------------------------
      USE GWM1DCV3
      USE GLOBAL,       ONLY:NPER
      USE GWM1BAS3,     ONLY:ONE,ZERO
      IMPLICIT NONE
      INTEGER,INTENT(IN)::IGRID,IO,VIFLG
!    Local variables used to hold data from GMNW2 or GWELLID 
      INTEGER NNODES,LOSSTYPEINDEX,QlimitL,NODNUM,INTNUM,QCUTL
      INTEGER PUMPLAY,PUMPROW,PUMPCOL,PUMPLOC,PPFLAG,PUMPCAP
      DOUBLE PRECISION HliftL,LIFTq0L,LIFTqdesL,HWtolL,LiftnL,QnL 
      DOUBLE PRECISION HlimL,QfrcmnL,QfrcmxL
      CHARACTER*20 WELLNAME,LOSSTYPE
!    Local variables used to hold data from GMNWNOD or GMNWINT
      INTEGER IL,IR,IC
      DOUBLE PRECISION Ztop,Zbotm,RwNode,RskinNode,KskinNode,
     &                 BNode,CNode,PNode,CWCNode,PP,GQDES
!    Local variables used for indices, counters
      INTEGER I,J,T,K,INODE,IINT,index,MNWID,ITMPG,NFVMNW
      INTEGER NAUX,IAUX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL GWM1DCV3MNW2PNT12(IGRID)
C-----ITEM 1: MNWMAX,[NODTOT],IWL2CB,MNWPRNT,OPTION
      ! For GWM-2005, NAUX_EXTERNAL has been pre-set to zero
      ! For GWM-VI, NAUX_EXTERNAL will be set to one externally
      NAUX = NAUX_EXTERNAL                       
      NFVMNW=0                                   ! COUNT NUMBER OF MNW-TYPE FV
      DO I=1,NFVAR                               ! LOOP OVER GWM FLOW VARIABLES
        IF(IGRID.EQ.GRDLOCDCV(I).AND.            ! FLOW VARIABLE ON CURRENT GRID
     &       FVNCELL(I).LT.0)THEN                ! THIS IS A MNW-TYPE FV
          NFVMNW=NFVMNW+ABS(FVNCELL(I))   
        ENDIF  
      ENDDO
      CALL WRITE_ITEM1
C-----LOOP OVER MULTINODE WELLS
      DO  MNWID=1,ABS(GMNWMAX)
C-----ITEM 2: FOR UNMANAGED MNW WELLS
        CALL WRITE_ITEM2_UNM
      ENDDO
      DO 50 I=1,NFVAR                            ! LOOP OVER GWM FLOW VARIABLES
        IF(IGRID.EQ.GRDLOCDCV(I).AND.            ! FLOW VARIABLE ON CURRENT GRID
     &       FVNCELL(I).LT.0)THEN                ! THIS IS A MNW-TYPE FV
          DO K=1,ABS(FVNCELL(I))                 ! LOOP OVER CELLS FOR FV
C-----ITEM 2: FOR MANAGED MNW WELLS ASSOCIATED WITH AN MNW-TYPE FLOW VARIABLE
            CALL WRITE_ITEM2_MAN
          ENDDO
        ENDIF  
  50  ENDDO
C-----LOOP OVER STRESS PERIODS
      DO 200 T=1,NPER
        CALL GWM1DCV3MNW2PNT34(IGRID,T,0)   ! RETRIEVE ITMP VALUE FOR MNW DATA
C-------COUNT NUMBER OF ACTIVE MNW WELLS IN THIS STRESS PERIOD
        ITMPG=0                             ! DETERMINE ITMP FOR MNW-TYPE FV
        DO 150 I=1,NFVAR                    ! LOOP OVER GWM FLOW VARIABLES
          IF(IGRID.EQ.GRDLOCDCV(I).AND.     ! FLOW VARIABLE ON CURRENT GRID
     &       FVSP(I,T)  .AND.            ! FLOW VARIABLE ACTIVE IN STRESS PERIOD
     &       FVNCELL(I).LT.0)THEN                ! THIS IS A MNW-TYPE FV
            ITMPG = ITMPG+ABS(FVNCELL(I))
          ENDIF  
  150   ENDDO
C-------ITEM 3 AND 4: ITMP AND WELLID and QDES 
        CALL WRITE_ITEM4
  200 ENDDO

      RETURN
      CONTAINS
c
C***********************************************************************
      SUBROUTINE WRITE_ITEM1
      CHARACTER(LEN=15) :: AUXTEXT
      AUXTEXT = ' AUX ' // AUXVARNAME
C-----WRITE TEXT LINE
      WRITE(IO,9000)
C-----WRITE ITEM 1  
      ! Add NODTOT to the line with negative MNWMAX as implemented in MF-2005 V1.10
      ! GWM does not write the AUX variables in original MNW file when 
      ! called from GWM-2005, but it does write one AUX variable when
      ! called from GWM-VI.
      IF (NAUX_EXTERNAL == 0) THEN   ! Call from GWM-2005
        GMNWMAX = ABS(GMNWMAX) ! Regardless of original sign of MNWMAX write NODTOT
        WRITE(IO,9100)-GMNWMAX-NFVMNW,GNODTOT,GIWL2CB,GMNWPRNT  
      ELSE                           ! Call from GWM-VI
C        IF(GMNWMAX.GT.0)THEN       ! Don't use NODTOT      
        IF(GMNWMAX.GE.0)THEN       ! Don't use NODTOT - Bug fix 5/5/2014 ERB
          WRITE(IO,9110) GMNWMAX+NFVMNW,        GIWL2CB,GMNWPRNT,AUXTEXT
        ELSEIF(GMNWMAX.LT.0)THEN   ! Use NODTOT and force MNWMAX to negative    
          WRITE(IO,9115) GMNWMAX-NFVMNW,GNODTOT,GIWL2CB,GMNWPRNT,AUXTEXT
        ENDIF  
      ENDIF  
      RETURN
 9000 FORMAT('# MNW2 input file written by GWM',/,
     & '#   Combines unmanaged MNW2 wells read from original MNW file'
     & ,/,'#     with managed MNW2 wells read from DECVAR file'
     & ,/,'#     This file can be discarded')
 9100 FORMAT(4I10,    T70,'# MNWMAX,NODTOT,IWL2CB,MNWPRNT')
 9110 FORMAT(3I10,A15,T70,'# MNWMAX,IWL2CB,MNWPRNT,OPTION')
 9115 FORMAT(4I10,A15,T70,'# MNWMAX,NODTOT,IWL2CB,MNWPRNT,OPTION')
      END SUBROUTINE WRITE_ITEM1
C***********************************************************************
      SUBROUTINE WRITE_ITEM2_UNM
      ! ITEM 2A
        WELLNAME = GWELLID(MNWID)   
        NNODES =   GMNW2(2,MNWID)
        WRITE(IO,9200)WELLNAME,NNODES
      ! ITEM 2B
        LOSSTYPEINDEX=INT(GMNW2(3,MNWID))
        CALL FIND_LOSSTYPE(LOSSTYPEINDEX,LOSSTYPE)
        PUMPLOC=GMNW2(11,MNWID)
        PUMPLOC=ABS(PUMPLOC)  ! If negative, reset to postive since well indices already computed
        QlimitL =GMNW2( 6,MNWID)
        PPFLAG =GMNW2(19,MNWID)
        PUMPCAP=GMNW2(22,MNWID)
        WRITE(IO,9210) LOSSTYPE,PUMPLOC,QlimitL,PPFLAG,PUMPCAP
      ! ITEM 2C - ALL LOSSTYPE INFO READ NODE-WISE, SO SET 2C ITEMS TO -1
        IF(LOSSTYPEINDEX.NE.0)WRITE(IO,9220)-1.,-1.,-1.,-1.  
      ! ITEM 2D
        IF(NNODES.GT.0)THEN       ! ITEM 2D-1 
          NODNUM = GMNW2( 4,MNWID)
          DO INODE =1,NNODES            
          IL       =GMNWNOD(1,NODNUM+INODE-1)           
          IR       =GMNWNOD(2,NODNUM+INODE-1)           
          IC       =GMNWNOD(3,NODNUM+INODE-1)           
          RwNode   =GMNWNOD(5,NODNUM+INODE-1)             
          RskinNode=GMNWNOD(6,NODNUM+INODE-1)             
          KskinNode=GMNWNOD(7,NODNUM+INODE-1)             
          BNode    =GMNWNOD(8,NODNUM+INODE-1)            
          CNode    =GMNWNOD(9,NODNUM+INODE-1)          
          PNode    =GMNWNOD(10,NODNUM+INODE-1)           
          CWCNode  =GMNWNOD(11,NODNUM+INODE-1)             
          PP       =GMNWNOD(19,NODNUM+INODE-1) 
          CALL WRITE_2D1
          ENDDO
        ELSE                     ! ITEM 2D-2
          INTNUM = GMNW2(13,MNWID)
          DO IINT =1,ABS(NNODES)       
          Ztop     =GMNWINT(1,INTNUM+IINT-1)           
          Zbotm    =GMNWINT(2,INTNUM+IINT-1)          
          IR       =GMNWINT(3,INTNUM+IINT-1)           
          IC       =GMNWINT(4,INTNUM+IINT-1)           
          RwNode   =GMNWINT(5,INTNUM+IINT-1)           
          RskinNode=GMNWINT(6,INTNUM+IINT-1)             
          KskinNode=GMNWINT(7,INTNUM+IINT-1)            
          BNode    =GMNWINT(8,INTNUM+IINT-1)            
          CNode    =GMNWINT(9,INTNUM+IINT-1)            
          PNode    =GMNWINT(10,INTNUM+IINT-1)            
          CWCNode  =GMNWINT(11,INTNUM+IINT-1)            
          CALL WRITE_2D2
          ENDDO
        ENDIF 
      ! ITEMS 2E THROUGH 2H  
        CALL WRITE_THE_REST
      RETURN
 9200 FORMAT(A20,I10, T70,'# WELLID,NNDODE - UnManaged')
 9210 FORMAT(A20,4I10,T70,'# LOSSTYPE,PUMPLOC,ETC. - UnManaged')
 9220 FORMAT(4D10.2,  T70,'# LOSSTYPE ITEMS - UnManaged')
      END SUBROUTINE WRITE_ITEM2_UNM
C***********************************************************************
      SUBROUTINE WRITE_ITEM2_MAN
      ! ITEM 2A
        WELLNAME = FVMNW(I)%FVWELLID(K) 
        NNODES =   FVMNW(I)%FVMNW2(2,K)
        WRITE(IO,9200)WELLNAME,NNODES
      ! ITEM 2B
        LOSSTYPEINDEX=INT(FVMNW(I)%FVMNW2(3,K))
        CALL FIND_LOSSTYPE(LOSSTYPEINDEX,LOSSTYPE)
        PUMPLOC=FVMNW(I)%FVMNW2(11,K)
        PUMPLOC=ABS(PUMPLOC)  ! If negative, reset to postive since well indices already computed
        QlimitL=FVMNW(I)%FVMNW2( 6,K)
        PPFLAG =FVMNW(I)%FVMNW2(19,K)
        PUMPCAP=FVMNW(I)%FVMNW2(22,K)
        WRITE(IO,9210) LOSSTYPE,PUMPLOC,QlimitL,PPFLAG,PUMPCAP
      ! ITEM 2C - ALL LOSSTYPE INFO READ NODE-WISE, SO SET 2C ITEMS TO -1
        WRITE(IO,9220)-1.,-1.,-1.,-1.
      ! ITEM 2D
        IF(NNODES.GT.0)THEN       ! ITEM 2D-1 
          DO INODE =1,NNODES            
          IL       =FVMNW(I)%FVCELLS(K)%FVMNWNOD(1,INODE)           
          IR       =FVMNW(I)%FVCELLS(K)%FVMNWNOD(2,INODE)           
          IC       =FVMNW(I)%FVCELLS(K)%FVMNWNOD(3,INODE)           
          RwNode   =FVMNW(I)%FVCELLS(K)%FVMNWNOD(5,INODE)             
          RskinNode=FVMNW(I)%FVCELLS(K)%FVMNWNOD(6,INODE)             
          KskinNode=FVMNW(I)%FVCELLS(K)%FVMNWNOD(7,INODE)             
          BNode    =FVMNW(I)%FVCELLS(K)%FVMNWNOD(8,INODE)            
          CNode    =FVMNW(I)%FVCELLS(K)%FVMNWNOD(9,INODE)          
          PNode    =FVMNW(I)%FVCELLS(K)%FVMNWNOD(10,INODE)           
          CWCNode  =FVMNW(I)%FVCELLS(K)%FVMNWNOD(11,INODE)             
          PP       =FVMNW(I)%FVCELLS(K)%FVMNWNOD(19,INODE) 
          CALL WRITE_2D1
          ENDDO
        ELSE                     ! ITEM 2D-2
          DO IINT =1,ABS(NNODES)       
          Ztop     =FVMNW(I)%FVCELLS(K)%FVMNWINT(1,IINT)           
          Zbotm    =FVMNW(I)%FVCELLS(K)%FVMNWINT(2,IINT)          
          IR       =FVMNW(I)%FVCELLS(K)%FVMNWINT(3,IINT)           
          IC       =FVMNW(I)%FVCELLS(K)%FVMNWINT(4,IINT)           
          RwNode   =FVMNW(I)%FVCELLS(K)%FVMNWINT(5,IINT)           
          RskinNode=FVMNW(I)%FVCELLS(K)%FVMNWINT(6,IINT)             
          KskinNode=FVMNW(I)%FVCELLS(K)%FVMNWINT(7,IINT)            
          BNode    =FVMNW(I)%FVCELLS(K)%FVMNWINT(8,IINT)            
          CNode    =FVMNW(I)%FVCELLS(K)%FVMNWINT(9,IINT)            
          PNode    =FVMNW(I)%FVCELLS(K)%FVMNWINT(10,IINT)            
          CWCNode  =FVMNW(I)%FVCELLS(K)%FVMNWINT(11,IINT)            
          CALL WRITE_2D2
          ENDDO
        ENDIF 
      ! ITEMS 2E THROUGH 2H  
      !   GWM DOES NOT SUPPORT CAPABILITIES IN THESE ITEMS
      RETURN
 9200 FORMAT(A20,I10, T70,'# WELLID,NNDODE - Managed')
 9210 FORMAT(A20,4I10,T70,'# LOSSTYPE,PUMPLOC,ETC. - Managed')
 9220 FORMAT(4D10.2,  T70,'# LOSSTYPE ITEMS - Managed')     
      END SUBROUTINE WRITE_ITEM2_MAN
C***********************************************************************
      SUBROUTINE WRITE_THE_REST
    ! ITEM 2E
        IF(PUMPLOC.NE.0)THEN
          PUMPLAY = GMNW2(14,MNWID)
          PUMPROW = GMNW2(15,MNWID)
          PUMPCOL = GMNW2(16,MNWID)
          WRITE(IO,9230) PUMPLAY,PUMPROW,PUMPCOL 
        ENDIF 
    ! ITEM 2F  
        IF(Qlimit.GT.0.0) THEN
          HlimL = GMNW2(7,MNWID)
          QCUTL = GMNW2(8,MNWID)
          IF(QCUT.NE.0) THEN
            QfrcmnL=GMNW2(9,MNWID)
            QfrcmxL=GMNW2(10,MNWID)
            WRITE(IO,9240) HlimL,QCUTL,QfrcmnL,QfrcmxL
          ELSE
            WRITE(IO,9240) HlimL,QCUTL
          END IF
        ENDIF
    ! ITEM 2G
        IF(PUMPCAP.GT.0.0) THEN
          HliftL   = GMNW2(23,MNWID)
          LIFTq0L  = GCapTable(MNWID,1,1)
          LIFTqdesL= GCapTable(MNWID,PUMPCAP+2,1)
          HWtolL   = GMNW2(28,MNWID)
          WRITE(IO,9250) HliftL,LIFTq0L,LIFTqdesL,HWtolL   
    ! ITEM 2H
          DO index=2,PUMPCAP+1
            LIFTnL= GCapTable(MNWID,index,1)
            QnL   = GCapTable(MNWID,index,2)
            WRITE(IO,9250) LiftnL,QnL              
          ENDDO
	  END IF
      RETURN
 9230 FORMAT(4I10)
 9240 FORMAT(D20.12,I20,2D20.12)
 9250 FORMAT(4D20.12)
      END SUBROUTINE WRITE_THE_REST
C***********************************************************************
      SUBROUTINE WRITE_2D1
        IF(PPFLAG.GT.0)THEN
        SELECT CASE (LOSSTYPEINDEX)
          CASE (0)
            WRITE(IO,9240) IL,IR,IC
          CASE (1)
            WRITE(IO,9240) IL,IR,IC,Rwnode,PP
          CASE (2)
            WRITE(IO,9240) IL,IR,IC,Rwnode,RskinNode,KskinNode,PP
          CASE (3)
            WRITE(IO,9240) IL,IR,IC,Rwnode,BNode,CNode,PNode,PP
          CASE (4)
            WRITE(IO,9240) IL,IR,IC,CWCNode,PP
        END SELECT
        ELSE
        SELECT CASE (LOSSTYPEINDEX)
          CASE (0)
            WRITE(IO,9240) IL,IR,IC
          CASE (1)
            WRITE(IO,9240) IL,IR,IC,Rwnode
          CASE (2)
            WRITE(IO,9240) IL,IR,IC,Rwnode,RskinNode,KskinNode
          CASE (3)
            WRITE(IO,9240) IL,IR,IC,Rwnode,BNode,CNode,PNode
          CASE (4)
            WRITE(IO,9240) IL,IR,IC,CWCNode
        END SELECT
        ENDIF
      RETURN
 9240 FORMAT(3I5,4D20.12)
      END SUBROUTINE WRITE_2D1
C***********************************************************************
      SUBROUTINE WRITE_2D2
        SELECT CASE (LOSSTYPEINDEX)
          CASE (1)
            WRITE(IO,9240) Ztop,Zbotm,IR,IC,Rwnode
          CASE (2)
            WRITE(IO,9240) Ztop,Zbotm,IR,IC,Rwnode,RskinNode,KskinNode
          CASE (3)
            WRITE(IO,9240) Ztop,Zbotm,IR,IC,Rwnode,BNode,CNode,PNode
          CASE (4)
            WRITE(IO,9240) Ztop,Zbotm,IR,IC,CWCNode
        END SELECT
      RETURN
 9240 FORMAT(2D20.12,2I5,4D20.12)
      END SUBROUTINE WRITE_2D2
C***********************************************************************
      SUBROUTINE WRITE_ITEM4
      INTEGER :: AUXVALUE
C-----ITEM 3: ITMP
      WRITE(IO,9300)ITMP+ITMPG
C-----WRITE ITEM 4 FOR UNMANAGED MNW WELLS
      ! For GWM-VI calls (VIFLG=0), add an AUX variable value of zero
      DO 100 MNWID=1,ITMP
        CALL GWM1DCV3MNW2PNT34(IGRID,T,MNWID)    ! RETRIEVE MNW DATA
        PUMPCAP=GMNW2(22,MNWID)
        WELLNAME = GWELLID(WELLINDX)   
        IF(PUMPCAP.EQ.0)THEN ! Don't write CapMult
          IF(Cprime.LE.ZERO)THEN
            IF(VIFLG.EQ.1)WRITE(IO,9400) WELLNAME,Qdes
            IF(VIFLG.EQ.0)WRITE(IO,9402) WELLNAME,Qdes,0
          ELSE
            IF(VIFLG.EQ.1)WRITE(IO,9410) WELLNAME,Qdes,Cprime
            IF(VIFLG.EQ.0)WRITE(IO,9412) WELLNAME,Qdes,Cprime,0
          ENDIF  
        ELSE                 ! Write CapMult
          IF(Cprime.LE.ZERO)THEN
            IF(VIFLG.EQ.1)WRITE(IO,9420) WELLNAME,Qdes,CapMult
            IF(VIFLG.EQ.0)WRITE(IO,9422) WELLNAME,Qdes,CapMult,0
          ELSE
            IF(VIFLG.EQ.1)WRITE(IO,9430) WELLNAME,Qdes,CapMult,Cprime
            IF(VIFLG.EQ.0)WRITE(IO,9432) WELLNAME,Qdes,CapMult,Cprime,0
          ENDIF  
        ENDIF
        IF(Qlimit.LT.0.0) THEN
          IF(QCUT.NE.0) THEN
            WRITE(IO,9450) Hlim,QCUT,Qfrcmn,Qfrcmx
          ELSE
            WRITE(IO,9460) Hlim,QCUT
          END IF
        ENDIF
 100  ENDDO
C-----WRITE ITEM 4 FOR MNW-TYPE FLOW VARIABLES
      DO 220 I=1,NFVAR                           ! LOOP OVER GWM FV
        IF(IGRID.EQ.GRDLOCDCV(I))THEN            ! FV ON CURRENT GRID
        IF(FVSP(I,T)) THEN                       ! FV ACTIVE IN STRESS PERIOD
        IF(FVNCELL(I).LT.0)THEN                  ! THIS IS A MNW-TYPE FV
          CALL GWM1DCV3FVCPNT(I)                 ! POINT TO CORRECT CELL INFO
          DO 110 K=1,ABS(FVNCELL(I))             ! LOOP OVER CELLS FOR FV
            AUXVALUE = FVMNW(I)%FVCELLS(K)%AUXVALUE
            WELLNAME = FVMNW(I)%FVWELLID(K)      ! RETRIEVE NAME
            IF(VIFLG.EQ.1)THEN                   ! THIS IS A GWM-2005 CALL
              GQDES = FVBASE(I)*FVRATIO(K)! ASSIGN FLOW RATE 
              WRITE(IO,9400) WELLNAME,GQdes
            ELSEIF(VIFLG.EQ.0)THEN               ! THIS IS A GWM-VI CALL
              WRITE(IO,9404) WELLNAME,('%'//ADJUSTL(WELLNAME)//'%'),
     1                       AUXVALUE
            ENDIF
 110      ENDDO
        ENDIF
        ENDIF
        ENDIF
 220  ENDDO
     
      RETURN
 9300 FORMAT(I5,T70,'# ITMP')
 9400 FORMAT(A20,D20.12,          T70,'# WELLNAME,Qdes')
 9402 FORMAT(A20,D20.12,I10,      T70,'# WELLNAME,Qdes,Auxvalue')
 9404 FORMAT(A20,A     ,I10,      T70,'# WELLNAME,Qdes,Auxvalue')
 9410 FORMAT(A20,2D20.12,    T72,'# WELLNAME,Qdes,Cprime')
 9412 FORMAT(A20,2D20.12,F10.0,    T72,
     1   '# WELLNAME,Qdes,Cprime,Auxvalue')
 9420 FORMAT(A20,2D20.12,    T72,'# WELLNAME,Qdes,CapMult')
 9422 FORMAT(A20,2D20.12,F10.0,    T72,
     1   '# WELLNAME,Qdes,CapMult,Auxvalue')
 9430 FORMAT(A20,3D20.12,    T82,'# WELLNAME,Qdes,CapMult,Cprime')
 9432 FORMAT(A20,3D20.12,F10.0,    T82,
     1   '# WELLNAME,Qdes,CapMult,Cprime,Auxvalue')
 9450 FORMAT(D20.12,2X,I5,2D20.12,T70,'# Hlim,QCUT,Qfrcmn,Qfrcmx')
 9460 FORMAT(D20.12,2X,I5,        T70,'# Hlim,QCUT')
      END SUBROUTINE WRITE_ITEM4

      END SUBROUTINE GWM1DCV3ARW
      
C***********************************************************************
      SUBROUTINE FIND_LOSSTYPE(INDEX,LTYPE)
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
      END SUBROUTINE FIND_LOSSTYPE
C***********************************************************************
      SUBROUTINE GWM1DCV3MNW2PNT12(IGRID)
      INTEGER(I4B),INTENT(IN)::IGRID
C     Retrieve MNW data for a grid well 
      GMNW2    => GWMMNWDAT(IGRID)%GMNW2_TARGET
      GMNWNOD  => GWMMNWDAT(IGRID)%GMNWNOD_TARGET
      GMNWINT  => GWMMNWDAT(IGRID)%GMNWINT_TARGET
      GCapTable=> GWMMNWDAT(IGRID)%GCapTable_TARGET
      GWELLID  => GWMMNWDAT(IGRID)%GWELLID_TARGET
      GMNWMAX  => GWMMNWDAT(IGRID)%GMNWMAX
      GNODTOT  => GWMMNWDAT(IGRID)%GNODTOT
      GIWL2CB  => GWMMNWDAT(IGRID)%GIWL2CB
      GMNWPRNT => GWMMNWDAT(IGRID)%GMNWPRNT
      RETURN
      END SUBROUTINE GWM1DCV3MNW2PNT12
C***********************************************************************
      SUBROUTINE GWM1DCV3MNW2PSV12(IGRID)
      USE GWFMNW2MODULE, ONLY:NMNW2,MNWMAX,NMNWVL,IWL2CB,MNWPRNT,
     1                       NODTOT,INTTOT,MNW2,MNWNOD,MNWINT,
     2                       CapTable,SMALL,NTOTNOD,WELLID
      INTEGER(I4B),INTENT(IN)::IGRID
C     Save MNW data for a grid and well
      GWMMNWDAT(IGRID)%GMNW2_TARGET     =>MNW2    
      GWMMNWDAT(IGRID)%GMNWNOD_TARGET   =>MNWNOD   
      GWMMNWDAT(IGRID)%GMNWINT_TARGET   =>MNWINT   
      GWMMNWDAT(IGRID)%GCapTable_TARGET =>CapTable 
      GWMMNWDAT(IGRID)%GWELLID_TARGET   =>WELLID   
      ! GWMMNWDAT(IGRID)%GMNWAUX_TARGET   =>GMNWAUX
      GWMMNWDAT(IGRID)%GMNWMAX  =>MNWMAX
      GWMMNWDAT(IGRID)%GNODTOT  =>NODTOT
      GWMMNWDAT(IGRID)%GIWL2CB  =>IWL2CB
      GWMMNWDAT(IGRID)%GMNWPRNT =>MNWPRNT
      RETURN
      END SUBROUTINE GWM1DCV3MNW2PSV12
C***********************************************************************
      SUBROUTINE GWM1DCV3MNW2PNT34(IGRID,T,N)
      INTEGER(I4B),INTENT(IN)::IGRID,T,N
C     Retrieve MNW data for a grid, stress period and node 
      IF(N.EQ.0)THEN
        ITMP   =>GWMMNWDAT(IGRID)%GT_TARGET(T)%ITMP
      ELSE
        Qlimit =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qlimit
        QCUT   =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%QCUT
        Qdes   =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qdes
        Qfrcmn =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qfrcmn
        Qfrcmx =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qfrcmx
        CapMult=>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%CapMult
        Cprime =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Cprime
        Hlim   =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Hlim
        ! GMNWAUX3=>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%GMNWAUX3
        WELLINDX =>GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%WELLINDX
      ENDIF
      RETURN
      END SUBROUTINE GWM1DCV3MNW2PNT34
C***********************************************************************
      SUBROUTINE GWM1DCV3MNW2PSV34(IGRID,T,N)
      INTEGER(I4B),INTENT(IN)::IGRID,T,N
C     Save MNW data for a grid, stress period and node
      IF(N.EQ.0)THEN
        GWMMNWDAT(IGRID)%GT_TARGET(T)%ITMP                  =>ITMP  
      ELSE 
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qlimit =>Qlimit
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%QCUT   =>QCUT
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qdes   =>Qdes
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qfrcmn =>Qfrcmn
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Qfrcmx =>Qfrcmx
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%CapMult=>CapMult
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Cprime =>Cprime
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%Hlim   =>Hlim
        ! GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%GMNWAUX3=>GMNWAUX3
        GWMMNWDAT(IGRID)%GT_TARGET(T)%GTW_TARGET(N)%WELLINDX=>WELLINDX
      ENDIF
      RETURN
      END SUBROUTINE GWM1DCV3MNW2PSV34
     

      END MODULE GWM1DCV3MNW
