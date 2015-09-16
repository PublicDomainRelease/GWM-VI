MODULE MKMMPROC
  USE UTILITIES
  USE GWM1HDC3, ONLY: HDCNAME,HDCNUM,HDCSP,NHB,NDD,NDF,NGD,&
                    HDCILOC,HDCJLOC,HDCKLOC
  USE GWM1STC3, ONLY: STCNAME,STCNUM,STCSP,NSF,NSD,NLK,STCSLOC,STCRLOC
  USE GWM1STA3, ONLY: SVNAME,STANUM,NHVAR,NRVAR,NSVAR,NDVAR,NSTADEP,SVSP,&
                    FLOWTYPE,SVILOC,SVJLOC,SVKLOC,GWMSTADAT
  USE GWM1DCV3, ONLY: FVILOC,FVJLOC,FVKLOC,FVNAME,FVNCELL,FVSP,NFVAR,FVMNW
  USE MKMMPROC_RDNAM, ONLY: READ_NAME_FILE
  USE MKMMPROC_RDNAM, ONLY: LISTFILE,HCLOSE, &
                            IHEDUNFILE,IFLOWCBFILE,SFRFILE,DRNFILE,WELFILE, &
                            SFRB1FILE,DRNBFILE,SFRB2FILE,GWMWFILENAME
  USE GLOBAL,         ONLY: NCOL,NROW,NLAY,NPER,NSTP,PERLEN, &
                            TSMULT,BOTM,IFREFM
  USE GWFBASMODULE,   ONLY: HDRY,HNOFLO
  USE GWM_SUBS,       ONLY: EPSQNET
  IMPLICIT NONE 
  PRIVATE
  PUBLIC :: WRMMPROCIN,GWMWFILENAME
!      Local variables
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  CHARACTER(LEN=2000) :: FILENAME
  CHARACTER(LEN=20)   :: TEMPNAME,WELLNAME
  INTEGER:: IUMMP,IUJIF,LOCAT,ISTAT,ISTRT,NABF,NMQWEL,NMQMNW,NCELL,NNODES
  INTEGER:: NSV,I,K,II,JJ,KK,T,IL,IR,IC,INODE,IINT,DRN 
  INTEGER(I4B),SAVE,POINTER::SVZONE(:,:,:,:)
  INTEGER(I4B),SAVE,POINTER::MDCELLOC(:,:),WINC(:)
  DOUBLE PRECISION Ztop,Zbotm
CONTAINS
!
!***********************************************************************
      SUBROUTINE WRMMPROCIN(NAMFILE,SIMVALJIF_FILE,MMPROCJTF_FILE, &
                            MNW2JTF_FILE,GWMOUT)
!***********************************************************************
!
! Create three files 
!     1) The template file, MMProc.in.jtf, for the MODFLOW Management Pre- and Post-Processor
!     2) The JUPITER instruction file, SimulatedValues.jif
!     3) The JUPITER template file, mnw2_input.jtf
!
!-----------------------------------------------------------------------
      USE UTILITIES, ONLY: UTL_SAMENAME, UTL_STOP
      USE GWMMOD, ONLY: SIMCOMMAND
      USE GWM1DCV3, ONLY : GWM1DCV3FVCPNT
      USE GWM_SUBS, ONLY: TO_LONGNAME
      IMPLICIT NONE
      ! Arguments
      CHARACTER(LEN=2000),INTENT(IN)::NAMFILE,SIMVALJIF_FILE, &
                                MMPROCJTF_FILE,MNW2JTF_FILE
      INTEGER, INTENT(IN) :: GWMOUT
      ! Local variables
      INTEGER :: LEN
      CHARACTER(LEN=2000) :: SIMCOMMANDTEMP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Open file SimulatedValues.jif, which will be a JUPITER
        ! instruction file for extracting values from SimulatedValues.out
        IUJIF = UTL_GETUNIT(980,1000)
        FILENAME = SIMVALJIF_FILE
        LOCAT=IUJIF
        OPEN(LOCAT,FILE=FILENAME,STATUS='REPLACE',IOSTAT=ISTAT) 
        IF(ISTAT.NE.0)GOTO 999
        ! Open file MMProc.in.jtf, which will be a JUPITER template file
        ! for creating MMProc.in
        IUMMP = UTL_GETUNIT(980,1000)
        FILENAME = MMPROCJTF_FILE 
        LOCAT=IUMMP
        OPEN(LOCAT,FILE=FILENAME,STATUS='REPLACE',IOSTAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 999
        !
        WRITE(IUMMP,9100)
        ! Item 1: Command -- Text string to be used by MMProc to invoke a model run.
        IF (SIMCOMMAND=='') THEN
          WRITE(GWMOUT,9900)
          CALL UTL_STOP('Error: No SIMCOMMAND entry found in OPTIONS input block')
        ENDIF
        CALL UTL_ADDQUOTE(SIMCOMMAND,SIMCOMMANDTEMP)
        LEN = LEN_TRIM(SIMCOMMANDTEMP)
        IF (LEN<54) THEN
          WRITE(IUMMP,9120)TRIM(SIMCOMMANDTEMP)
        ELSE
          WRITE(IUMMP,9125)TRIM(SIMCOMMANDTEMP)
        ENDIF
        ! Items 2-6:  NAM file name and basic model information
        CALL WRMMPROCBD  
!      Items 7-10:  Binary file information
        CALL WRMMPROCBF
!      Items 11-16:  Managed well information
        CALL WRMMPROCMW
!      Items 17-20:  Simulated value information
        CALL WRMMPROCSIM
!      Close files     
        CLOSE(IUMMP)
        CLOSE(IUJIF)
        
        CALL MKMNWINPUT(MNW2JTF_FILE)

      RETURN
  999 CONTINUE
!-----FILE-OPENING ERROR
        WRITE(*,9990)TRIM(FILENAME),LOCAT
        STOP
9100 FORMAT('jtf %')
9120 FORMAT(A,T55,'# Item 1: Command')
9125 FORMAT(A,2X,'# Item 1: Command')
9900 FORMAT(/,1X,'Error: No command with PURPOSE=FORWARD found')
9990 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,&
      2X,'-- STOP EXECUTION (WRMMPROCIN)')!
   CONTAINS
!      
!***********************************************************************
      SUBROUTINE WRMMPROCBD
!***********************************************************************
!
        CALL READ_NAME_FILE(NAMFILE)
!
!  Item 2: Name-file -- Modflow/Seawat/CFP/NWT/FMP name file
        WRITE(IUMMP,9000)TRIM(NAMFILE), '     # ITEM 2: NAMEFILE NAME'
  
!  Item 2a: LISTFILE -- Name of Modflow LIST output file
        WRITE(IUMMP,9000)TRIM(LISTFILE),'     # ITEM 2A: LISTFILE NAME'
  
!  Item 2b: NLAY,NROW,NCOL,NPER -- From DIS file
        WRITE(IUMMP,*)NLAY,NROW,NCOL,NPER  
     
!  Item 2c: NSTP -- From DIS file
        DO I=1,NPER
          WRITE(IUMMP,*)NSTP(I)
        ENDDO  

!  Item 3: HNOFLO -- From BAS file
        WRITE(IUMMP,9010)HNOFLO,'     # ITEM 3: HNOFLO'
  
!  Item 4: IFREFM -- If BAS file has FREE option, IFREFM=1; otherwise, 0
        WRITE(IUMMP,9020)IFREFM,'     # ITEM 4: IFREFM'
  
!  Item 5: HDRY -- From flow-package file
        WRITE(IUMMP,9010)HDRY,  '     # ITEM 5:  HDRY'
  
!  Item 6: HCLOSE -- From solver-package file
        WRITE(IUMMP,9010)HCLOSE,'     # ITEM 6: HCLOSE'
! 
      RETURN    
 9000 FORMAT(2A)
 9010 FORMAT(G24.16,T50,A)
 9020 FORMAT(I20,T50,A)

      END SUBROUTINE WRMMPROCBD
 !
!***********************************************************************
      SUBROUTINE WRMMPROCBF
!***********************************************************************
 
 ! Item 7: Name of binary head file  
        IF(IHEDUNFILE.NE.'')THEN               
          WRITE(IUMMP,9000)TRIM(IHEDUNFILE),'     # ITEM 7: BINARY FILE FOR HEADS'
        ELSE
          WRITE(IUMMP,9100)' ITEM 7: BINARY FILE FOR HEADS'
        ENDIF
 ! Item 8: Name of binary budget file for flow package ("none" if not needed)
        IF(NSVAR.GT.0)THEN
          IF(IFLOWCBFILE.NE.'')THEN               
            WRITE(IUMMP,9000)TRIM(IFLOWCBFILE),'     # ITEM 8: BINARY FILE FOR FLOW'
          ELSE
            WRITE(IUMMP,9100)' ITEM 8: BINARY FILE FOR FLOW'
          ENDIF
        ELSE
          WRITE(IUMMP,9000)'none     # ITEM 8: BINARY FILE FOR FLOW'
        ENDIF
!  Item 9: NABF -- Number of additional binary budget files.  I
        NABF = 0
        IF (STCNUM.GT.0.OR.NRVAR.GT.0) THEN
          IF(SFRB1FILE.NE.'') NABF = NABF + 1
          IF(SFRB2FILE.NE.'') NABF = NABF + 1
        ENDIF
        IF(NDVAR.GT.0)      NABF = NABF + 1
        WRITE(IUMMP,9020)NABF,'     # ITEM 9: NABF - Number of other binary files'
!  Item 10: SFR-LEAK  SFR-budget-file-name for SFR-LEAK (may be same as any other budget file)
        IF(STCNUM.GT.0.OR.NRVAR.GT.0)THEN                
          IF(SFRB1FILE.NE.'')THEN               
            WRITE(IUMMP,9002)TRIM(SFRB1FILE),'     # ITEM 10: BINARY FILE FOR SFR-LEAK'
          ELSE
!            WRITE(IUMMP,9100)' ITEM 10: BINARY FILE FOR SFR-LEAK'
          ENDIF
        ENDIF
!  Item 10: SFR-FLOW  SFR-budget-file-name for SFR-FLOW (may be same as any other budget file)
        IF(STCNUM.GT.0.OR.NRVAR.GT.0)THEN                
          IF(SFRB2FILE.NE.'')THEN               
            WRITE(IUMMP,9004)TRIM(SFRB2FILE),'     # ITEM 10: BINARY FILE FOR SFR-FLOW'
          ELSE
!            WRITE(IUMMP,9100)' ITEM 10: BINARY FILE FOR SFR-FLOW'
          ENDIF
        ENDIF
!  Item 10: DRN  DRN-budget-file-name (may be same as any other budget file)
        IF(NDVAR.GT.0)THEN                
          IF(DRNBFILE.NE.'')THEN               
            WRITE(IUMMP,9006)TRIM(DRNBFILE),'     # ITEM 10: BINARY FILE FOR DRAINS'
          ELSE
            WRITE(IUMMP,9100)' ITEM 10: BINARY FILE FOR DRAINS'
          ENDIF
        ENDIF
! 
      RETURN    
 9000 FORMAT(2A)
 9002 FORMAT('SFR-LEAK',2X,2A)
 9004 FORMAT('SFR-FLOW',2X,2A)
 9006 FORMAT('DRN',2X,2A)
 9010 FORMAT(G24.16,T50,A)
 9020 FORMAT(I20,T50,A)
 9100 FORMAT('none   #no file found for',A)
    END SUBROUTINE WRMMPROCBF
!
!***********************************************************************
      SUBROUTINE WRMMPROCMW
!***********************************************************************
        USE GWM_SUBS, ONLY: TO_LONGNAME
        CHARACTER(LEN=22) :: VALFIELD  ! Need to provide for maximum precision
        INTEGER :: AUXVALUE, LENVALFIELD
        LENVALFIELD = 22
!
!  Item 11:  Number of cells associated with managed flows 
!      Calculate number of managed wells in this stress period 
        NMQWEL = 0
        DO I=1,NFVAR
          IF(FVNCELL(I).GT.0)THEN
            NMQWEL = NMQWEL+FVNCELL(I) 
          ENDIF  
        ENDDO
        WRITE(IUMMP,*)NMQWEL,'     # ITEM 11: NMQWEL -- Number of managed WEL cells'
!  Item 12: WEL-input-file-name, read if NMQWEL > 0
        IF (NMQWEL.GT.0) THEN ! ERB 2/19/2014 CORRECTION - Added IF statement
          IF(WELFILE.NE.'')THEN               
            WRITE(IUMMP,9000)TRIM(WELFILE), '     # ITEM 12: WELL FILE NAME'
          ELSE
            WRITE(IUMMP,9000)'none','     # ITEM 12: WELL FILE NAME'
          ENDIF
!  
!  Items 13: NMQWEL pairs of records of the form:
!       13a:  Flow-identifier  Layer  Row  Column  Q
!       13b:  IACTIVE - Integer array of length NPER: 0 = inactive, >0 = active.  
          DO 200 I=1,NFVAR
            TEMPNAME = FVNAME(I)
            IF(FVNCELL(I).GT.0)THEN              ! This is WEL-type
              NCELL = FVNCELL(I)   
              CALL GWM1DCV3FVCPNT(I)               ! POINT TO ARRAYS
              DO 180 IC=1,NCELL
                IF (NCELL>1) CALL TO_LONGNAME(FVNAME(I),TEMPNAME,IC)
                VALFIELD = ' '
                VALFIELD = ADJUSTL(TEMPNAME)
                WRITE(IUMMP,9040)TRIM(TEMPNAME),FVKLOC(IC),FVILOC(IC),FVJLOC(IC),&
                               ('%'//VALFIELD//'%')
                DO 150 T=1,NPER
                  IF(FVSP(I,T))THEN
                    WRITE(IUMMP,9030)1
                  ELSE
                    WRITE(IUMMP,9030)0
                  ENDIF
  150           ENDDO
  180         ENDDO
            ENDIF
  200     ENDDO
        ENDIF ! ERB 2/19/2014 CORRECTION - Added ENDIF
!         
!  Item 14:  Number of MNW wells associated with managed flows 
        NMQMNW = 0
        DO I=1,NFVAR
          IF(FVNCELL(I).LT.0)THEN
            NMQMNW = NMQMNW+ABS(FVNCELL(I))
          ENDIF  
        ENDDO
        WRITE(IUMMP,*)NMQMNW,'     # ITEM 14: NMQMNW -- Number of managed MNW wells'
!         
!  Item 15a and 15b:  List of all cells in which the MNW well is located
        ! Iterate through list of flow variables
        DO I=1,NFVAR
          IF(FVNCELL(I).LT.0)THEN
          ! Iterate through list of MNW2 wells included in current flow variable
          DO K=1,ABS(FVNCELL(I))                 ! LOOP OVER CELLS FOR FV
            WELLNAME = FVMNW(I)%FVWELLID(K) 
            NNODES =   FVMNW(I)%FVMNW2(2,K)
            AUXVALUE = FVMNW(I)%FVCELLS(K)%AUXVALUE
            VALFIELD = ' '
            VALFIELD = '%' // ADJUSTL(WELLNAME)
            VALFIELD(LENVALFIELD:LENVALFIELD) = '%'
            IF(NNODES.GT.0)THEN       ! ITEM 2D-1 
              WRITE(IUMMP,9050)WELLNAME,NNODES,VALFIELD,AUXVALUE
              ! Iterate through list of cells included in current MNW2 well
              DO INODE =1,NNODES            
                IL =FVMNW(I)%FVCELLS(K)%FVMNWNOD(1,INODE)           
                IR =FVMNW(I)%FVCELLS(K)%FVMNWNOD(2,INODE)           
                IC =FVMNW(I)%FVCELLS(K)%FVMNWNOD(3,INODE)           
                WRITE(IUMMP,*)IL,IR,IC,'     # ITEM 15B: LAY,ROW,COL'
              ENDDO
            ELSE                     ! ITEM 2D-2
              ALLOCATE (WINC(NLAY))
              WINC = 0
              ! Iterate through list of cells included in current MNW2 well
              DO IINT =1,ABS(NNODES)       
                Ztop =FVMNW(I)%FVCELLS(K)%FVMNWINT(1,IINT)           
                Zbotm=FVMNW(I)%FVCELLS(K)%FVMNWINT(2,IINT)          
                IR   =FVMNW(I)%FVCELLS(K)%FVMNWINT(3,IINT)           
                IC   =FVMNW(I)%FVCELLS(K)%FVMNWINT(4,IINT)  
                DO KK=1,NLAY
                  IF( Ztop .GT.BOTM(IC,IR,KK-1))THEN  ! Interval top above layer top
                    IF(Zbotm.LT.BOTM(IC,IR,KK-1))THEN ! Interval bottom below layer top
                      WINC(KK)=1     ! Interval and layer intersect
                    ENDIF
                  ELSE                               ! Interval top below layer top  
                    IF(Ztop.GT.BOTM(IC,IR,KK))THEN   ! Interval top above layer bottom
                      WINC(KK)=1     ! Interval and layer intersect
                    ENDIF
                  ENDIF  
                ENDDO       
              ENDDO
              NCELL= SUM(WINC)
              WRITE(IUMMP,9050)WELLNAME,NCELL,VALFIELD,AUXVALUE
              DO KK=1,NLAY
                IF(WINC(KK).EQ.1)WRITE(IUMMP,*)KK,IR,IC,'     # ITEM 15B: LAY,ROW,COL'
              ENDDO
              DEALLOCATE (WINC)
            ENDIF 
!         
!  Item 15c:  List active stress periods for this MNW well
            DO T=1,NPER
              IF(FVSP(I,T))THEN
                WRITE(IUMMP,9030)1
              ELSE
                WRITE(IUMMP,9030)0
              ENDIF
            ENDDO
          ENDDO
          ENDIF
        ENDDO  
      RETURN    
 9000 FORMAT(2A)
 9010 FORMAT(G24.16,T50,A)
 9020 FORMAT(I20,T50,A)
 9030 FORMAT(10I5)
 9040 FORMAT(A,3I10,2X,A,'  # Item 13A: Name Lay Row Col dummyQ')
 9050 FORMAT(A20,1X,I5,1X,A,1X,I5,'  # ITEM 15A: NAME, number of cells, Qdes, and SeqNum')
      END SUBROUTINE WRMMPROCMW
!***********************************************************************
      SUBROUTINE WRMNWTPL
!***********************************************************************


      END SUBROUTINE WRMNWTPL
!***********************************************************************
      SUBROUTINE WRMMPROCSIM
!***********************************************************************
!    Item 17 - NSV - number of model simulated records
        NSV = NHB+NDD+2*NDF+2*NGD       ! Number of head constraints
        NSV = NSV + STCNUM              ! Number of streamflow constraints
        NSV = NSV + NSTADEP             ! Number of dependents related to state variables
        WRITE(IUMMP,9020)NSV,'   # ITEM 17: NSV - number of simulated values'
!      Write JUPITER instruction file header lines
        WRITE(IUJIF,9100)
        WRITE(IUJIF,9110)NSV

!    Item 18 - Simulated Value Information
!    Heads
        ISTRT = 0
        CALL WRMMPROCSH
!    Streamflows
        ISTRT = ISTRT+NHVAR
        CALL WRMMPROCSR
!    Storage Change
        ISTRT = ISTRT+NRVAR
        CALL WRMMPROCSS
!    Drains
        ISTRT = ISTRT+NSVAR
        CALL WRMMPROCSD
        
!  Item 19: SFR input-file-name 
        IF(STCNUM.GT.0.OR.NRVAR.GT.0)THEN                
          IF(SFRFILE.NE.'')THEN               
            WRITE(IUMMP,9200)TRIM(SFRFILE),'     # ITEM 19: SFR INPUT FILE'
          ELSE
            WRITE(IUMMP,9210)' ITEM 19: BINARY FILE FOR STREAMFLOW'
          ENDIF
        ENDIF
!  Item 20: DRN aux-variable name 
        IF(NDVAR.GT.0)THEN                
          WRITE(IUMMP,9200)'GWM-DR','     # ITEM 20: DRN-Aux-variable-name'
        ENDIF
!  Item 21: EPSQNET
        WRITE(IUMMP,9220)EPSQNET 
!
      RETURN    
 9020 FORMAT(I20,T50,A)
 9100 FORMAT('jif %')
 9110 FORMAT('StandardFile 0  2',2X,I10)
 9200 FORMAT(2A)
 9210 FORMAT('none   #no file found for',A)
 9220 FORMAT(G11.4,13X,'# ITEM 21: EPSQNET')
     END SUBROUTINE WRMMPROCSIM
!
!***********************************************************************
      SUBROUTINE WRMMPROCSH
!***********************************************************************
!
!      Write head constraint information 
        IC = 0
        IF (NHB>0 .OR. NDD>0) THEN
!      First write for head bound and drawdown constraints
          DO 100 I=1,NHB+NDD
            IC = IC+1
            WRITE(IUMMP,9050)HDCNAME(IC)
            WRITE(IUMMP,9060)HDCKLOC(IC,1),HDCILOC(IC,1), &
                  HDCJLOC(IC,1),HDCSP(IC)
            WRITE(IUJIF,9070)HDCNAME(IC)
  100     ENDDO
        ENDIF
!      Second, write for head difference and head gradient constraints
        IF (NDF>0 .OR. NGD>0) THEN
          DO 200 I=1,NDF+NGD
            IC = IC+1
!          Need to define names for each head in the constraint
            CALL TO_LONGNAME(HDCNAME(IC),TEMPNAME,1)
            WRITE(IUMMP,9050)TEMPNAME
            WRITE(IUMMP,9060)HDCKLOC(IC,1),HDCILOC(IC,1), &
                  HDCJLOC(IC,1),HDCSP(IC)
            WRITE(IUJIF,9070)TEMPNAME
            CALL TO_LONGNAME(HDCNAME(IC),TEMPNAME,2)
            WRITE(IUMMP,9050)TEMPNAME
            WRITE(IUMMP,9060)HDCKLOC(IC,2),HDCILOC(IC,2), &
                  HDCJLOC(IC,2),HDCSP(IC)
            WRITE(IUJIF,9070)TEMPNAME
200      ENDDO
        ENDIF
!      Third, write for head state variables
        DO 300 I=ISTRT+1,ISTRT+NHVAR
          WRITE(IUMMP,9050)SVNAME(I)
          WRITE(IUMMP,9060)SVKLOC(I),SVILOC(I),SVJLOC(I),&
                           SVSP(I,1)
          WRITE(IUJIF,9070)SVNAME(I)
 300    ENDDO
!
      RETURN
 9050 FORMAT(A,'  HEAD')
 9060 FORMAT(4I10,T50,'# Item 18B: LAY,ROW,COL,SP')
 9070 FORMAT(A)
      END SUBROUTINE WRMMPROCSH
!
!***********************************************************************
      SUBROUTINE WRMMPROCSR
!***********************************************************************
!      Write the name, colum, row, layer and stress period for GWM constraint
!
!      First, write streamflow constraint information
        DO 100 I=1,NSF+NSD               ! Flow-type
          WRITE(IUMMP,9050)STCNAME(I)
          WRITE(IUMMP,9060)STCSLOC(I),STCRLOC(I),STCSP(I)
          WRITE(IUJIF,9070)STCNAME(I)
  100   ENDDO
        DO 110 I=NSF+NSD+1,NSF+NSD+NLK   ! Leak-type
          WRITE(IUMMP,9052)STCNAME(I)
          WRITE(IUMMP,9060)STCSLOC(I),STCRLOC(I),STCSP(I)
          WRITE(IUJIF,9070)STCNAME(I)
  110   ENDDO
!      Second, write for streamflow state variables
        DO 200 I=ISTRT+1,ISTRT+NRVAR
          IF(FLOWTYPE(I))THEN
            WRITE(IUMMP,9050)SVNAME(I)
          ELSE
            WRITE(IUMMP,9052)SVNAME(I)
          ENDIF
          WRITE(IUMMP,9060)SVILOC(I),SVJLOC(I),SVSP(I,1)
          WRITE(IUJIF,9070)SVNAME(I)
 200    ENDDO
! 
      RETURN    
 9050 FORMAT(A,'  SFR-FLOW')
 9052 FORMAT(A,'  SFR-LEAK')
 9060 FORMAT(4I10,5X,'# Item 18:BEGINNING AND ENDING STRESS PERIOD')
 9070 FORMAT(A)
      END SUBROUTINE WRMMPROCSR
!
!***********************************************************************
      SUBROUTINE WRMMPROCSS
!***********************************************************************
!
        SVZONE=>GWMSTADAT(1)%SVZONE_TARGET   ! POINT TO CORRECT SVZONE ARRAY 
        DO 200 I=ISTRT+1,ISTRT+NSVAR
          WRITE(IUMMP,9050)SVNAME(I)
          WRITE(IUMMP,9060)SVSP(I,1),SVSP(I,2)
          WRITE(IUJIF,9070)SVNAME(I)
          DO 120 KK=1,NLAY
          DO 120 II=1,NROW
          DO 120 JJ=1,NCOL
            WRITE(IUMMP,'(10(I5,1X))')SVZONE(JJ,II,KK,I)
  120     ENDDO 
  200   ENDDO 
!  
      RETURN
 9050 FORMAT(A,'  STORAGE-CHANGE')
 9060 FORMAT(2I10,2X,'ZONE',5X,'# Item 18:BEGINNING AND ENDING STRESS PERIOD')
 9070 FORMAT(A)
      END SUBROUTINE WRMMPROCSS
!
!
!***********************************************************************
      SUBROUTINE WRMMPROCSD
!***********************************************************************
!
       DO 200 I=ISTRT+1,ISTRT+NDVAR
          MDCELLOC=>GWMSTADAT(1)%MDCELVAR(I-ISTRT)%MDCELLOC_TARGET ! POINT TO STORAGE
          KK = 0
          DO 120 II=1,SVILOC(I)
            KK = KK + 1
            CALL TO_LONGNAME(SVNAME(I),TEMPNAME,KK)
            IF(SVSP(I,2).LT.0)THEN         ! This is flow-type
              WRITE(IUMMP,9050)TEMPNAME
              WRITE(IUMMP,9060)MDCELLOC(II,1),MDCELLOC(II,2),MDCELLOC(II,3),&
                               MDCELLOC(II,4),SVSP(I,1)
            ELSEIF(SVSP(I,2).GT.0)THEN     ! This is volume-type
              WRITE(IUMMP,9150)TEMPNAME
              WRITE(IUMMP,9160)MDCELLOC(II,1),MDCELLOC(II,2),MDCELLOC(II,3),&
                               MDCELLOC(II,4),SVSP(I,1),SVSP(I,2)
            ENDIF
            WRITE(IUJIF,9070)TEMPNAME
  120     ENDDO 
  200   ENDDO 
!  
      RETURN
 9050 FORMAT(A,'  DRN-FLOW'  ,T50,'# Item 18A: Simulated-value-identifier Type')
 9060 FORMAT(5I9,T50,'# Item 18B: LAY ROW COL aux-value stress-period')
 9150 FORMAT(A,'  DRN-VOLUME',T50,'# Item 18A: Simulated-value-identifier Type')
 9160 FORMAT(6I8,T50,'# Item 18B:LAY ROW COL Aux-value Begin-stress-period End-stress-period')
 9070 FORMAT(A)
      END SUBROUTINE WRMMPROCSD
!
      END SUBROUTINE WRMMPROCIN
  
END MODULE MKMMPROC

