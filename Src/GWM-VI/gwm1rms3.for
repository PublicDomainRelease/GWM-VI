
C
      MODULE GWM1RMS3LP
C     VERSION: 21MAR2008 - package references converted to version 3
      USE GWM_SUBS, ONLY: IGETUNIT
      USE GWM_STOP, ONLY:   GSTOP
      USE MF2005_UTLS, ONLY: URDCOM, URWORD
      IMPLICIT NONE
      PRIVATE
      PUBLIC::GWM1SIMPLEX1,FNDFES,FACTOR,SLVLP,PRICE,REDCST,
     1        RATEST,RATIOC,SOLUTN,CSTRNG,RHSRNG,MAXXI,MINXI,BSOLVE,
     2        FATHOM,BRANCH,BINWRT,GWM1PSIMPLEX1,SLVBAS
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
C
C
C-----VARIABLES USED BY SIMPLEX ALGORITHM IN GWM1SIMPLEX1
      INTEGER(I4B),SAVE,ALLOCATABLE::LBNB(:),IRTIC(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::IAAO(:),JAAO(:)
      REAL(DP),SAVE,ALLOCATABLE::AAO(:),CAUG(:),W1(:),W2(:)
C       LBNB     -indidator list for all variables
C                  0  = basic variable
C                  1  = non-basic variable at lower bound
C                 -1  = non-basic variable at upper bound
C       IRTIC    -pointer list of basic column to full matrix column
C       AAO      -vector of non-zero value in original A matrix
C       IAAO     -row index for each non-zero value in A matrix
C       JAAO     -pointer to elements of AA0 and IAA0 that are the 
C                  beginning of columns of AA0
C       CAUG     -augmented cost vector
C       W1,W2    -work arrays
C
C-----VARIABLES USED TO CREATE TEMPORARY LP IN GWM1PSIMPLEX1
      REAL(DP),SAVE,ALLOCATABLE::TAMAT(:,:),TAMAT2(:,:),TCST(:),
     1                           TBNDS(:),TRHS(:)
      INTEGER(I4B),SAVE,ALLOCATABLE::ROWCHK(:),IDROP(:)
C
C-----VARIABLES USED BY LAPACK ROUTINES
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
      INTEGER(I4B),SAVE,ALLOCATABLE::IPIV(:),IWORK(:)
      REAL(DP),SAVE,ALLOCATABLE::A(:,:),AF(:,:),X(:,:),BLA(:,:)
      REAL(DP),SAVE,ALLOCATABLE::R(:),C(:),FERR(:),BERR(:),WORK(:)
      INTEGER(I4B),SAVE::NRHS, LDA, LDB, LDAF, LDX
      REAL(DP),SAVE::RCOND
      CHARACTER(LEN=1),SAVE::FACT,TRANS,EQUED
C
C-----FOR ERROR HANDLING
      INTEGER(I2B)::ISTAT  
      CHARACTER(LEN=200)::FLNM
      CHARACTER(LEN=20)::FILACT,FMTARG,ACCARG,FILSTAT
      INTEGER(I4B)::LLOC,ISTART,ISTOP,NDUM,FUNIT
      REAL(SP)::RDUM
C
      CONTAINS
C***********************************************************************
      SUBROUTINE GWM1SIMPLEX1(M,NV,NDV,AMAT,COST,BNDS,RHS,OBJ,
     &                        IFLG,LPITMAX,IPRT)
C***********************************************************************
C    VERSION: 27AUG2013
C    PURPOSE: SOLVE LINEAR PROGRAMS OF THE FORM
C
C                  MINIMIZE CX
C                   
C               SUCH THAT AX = B
C
C                     0 <= X <= U
C
C             AND PERFORM RANGE ANALYSIS
C
C   FEASIBLE SOLUTION OBTAINED USING THE PHASE 1/PHASE 2 METHOD
C   BASIS FACTORIZATION PERFORMED USING LAPACK ROUTINES
C
C   SEE: "ADVANCED LINEAR PROGRAMMING", BY BRUCE MURTAGH, MCGRAW-HILL, 1981
C   FOR DETAILS ON METHODS FOR IMPLEMENTATION
C
C-----------------------------------------------------------------------
C
C       INPUT 
C             M      - NUMBER OF CONSTRAINTS
C             NV     - DIMENSION OF X - THIS INCLUDES ALL DECISION VARIABLES
C                          AND ALL SLACK AND SURPLUS VARIABLES
C             NDV    - NUMBER OF COLUMNS IN AMAT 
C                          THIS WILL BE EQUAL TO THE NUMBER OF DECISION 
C                          VARIABLE PLUS 1.  
C             AMAT(M,NDV) - DENSE COEFFICIENT MATRIX
C                          GWM1SIMPLEX1 ASSUMES THAT COEFFICIENT MATRIX IS DENSE
C                          BUT THAT THE LAST COLUMN OF AMAT CONTAINS THE 
C                          COEFFICIENT (EITHER -1,0 OR 1) ON THE SLACK/SURPLUS
C                          VARIABLE FOR THE CORRESPONDING ROW
C             COST(NV) - COST COEFFICIENT VECTOR
C             BNDS(NV) - UPPER BOUND VECTOR
C             RHS(M)   - RIGHT HAND SIDE VECTOR
C             LPITMAX  - MAXIMUM NUMBER OF ITERATIONS ALLOWED
C             BIGINF   - VALUE OF MACHINE INFINITY; USED BY RATEST IN UNBOUNDED
C                          TO DETERMINE IF AN UPPER BOUND VALUE IS SET TO INFINITY
C             IPRT     - PRINT FLAG (1=PRINT PROGRESS, 0=NO PRINT)
C
C        OUTPUT
C             COST(NV) - VALUE OF X SOLUTION VECTOR (VALID ONLY IF IFLG=0)
C             OBJ      - VALUE OF OBJECTIVE FUNCTION (VALID ONLY IF IFLG=0)
C             RHS(M)   - VALUE OF DUAL VARIABLES (VALID ONLY IF IFLG=0)
C             BNDS(NV) - REDUCED COSTS OF ALL VARIABLES
C             IFLG     - OUTPUT STATUS FLAG
C                        0 = OPTIMAL SOLUTION FOUND
C                        1 = PROBLEM IS INFEASIBLE
C                        2 = PROBLEM IS UNBOUNDED
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE,ZERO,BIGINF,GWMOUT
      USE GWM1BAS3, ONLY : GWM1BAS3PS
      INTEGER(I4B),INTENT(IN)::M,NV,NDV,LPITMAX,IPRT
      INTEGER(I4B),INTENT(OUT)::IFLG
      REAL(DP),INTENT(IN)::AMAT(M,NDV)
      REAL(DP),INTENT(OUT)::OBJ
      REAL(DP),INTENT(INOUT)::BNDS(NV),RHS(M),COST(NV)
C-----LOCAL VARIABLES
      INTEGER(I4B)::NVA,NAC,I,N
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!      for debugging, write contents of AMAT,COST,BNDS,RHS
!      integer :: j,k
!      write(gwmout,10)
!   10 format('at top of GWM1SIMPLEX1, amat contains:')
!   15 format(10(g23.16,1x))
!      do j=1,m
!        write(gwmout,15)(amat(j,k),k=1,ndv)
!      enddo
!      write(gwmout,20)
!   20 format('at top of GWM1SIMPLEX1, cost contains:')
!      write(gwmout,15)(cost(j),j=1,nv)
!      write(gwmout,25)
!   25 format('at top of GWM1SIMPLEX1, bnds contains:')
!      write(gwmout,15)(bnds(j),j=1,nv)
!      write(gwmout,30)
!   30 format('at top of GWM1SIMPLEX1, rhs contains:')
!      write(gwmout,15)(rhs(j),j=1,m)
C
C-----CHECK THAT NO BOUNDS EXCEED BIGINF, THE VALUE OF MACHINE INFINITY 
      DO 100 I=1,NV
        IF(BNDS(I).GT.BIGINF)THEN
          CALL GSTOP('INPUT ERROR, BOUNDS EXCEED BIGINF')
        ENDIF
  100 ENDDO
C
C-----DEFINE SIZE AND ALLOCATE WORK SPACE FOR LINEAR PROGRAM
      NVA = NV+M
      NAC = NDV*M+M
      ALLOCATE (LBNB(NVA), IRTIC(M),
     &          IAAO(NAC),JAAO(NVA+1),
     &          AAO(NAC),CAUG(NVA),W1(M),W2(M),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
      N = M
      NRHS = 1
      LDA  = N
      LDB  = N
      LDAF = N
      LDX  = N
      ALLOCATE (IPIV(N),IWORK(N),
     &          A(LDA,N),AF(LDAF,N),X(LDX,NRHS),BLA(LDB,NRHS),
     &          R(N),C(N),FERR(NRHS),BERR(NRHS),WORK(4*N),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
C
C-----PERFORM PHASE 1 - SOLVE FROM AN ALL ARTIFICIAL START
	CALL FNDFES(AMAT,COST,BNDS,RHS,M,NV,NDV,LPITMAX,IFLG,IPRT)
      IF(IFLG.EQ.0)THEN
C-------SOLUTION FOUND TO PHASE 1.
        IF(IPRT.EQ.1)CALL GWM1BAS3PS('    Feasible Solution Found',0)   
C
C-------PERFORM PHASE 2 - SOLVE FROM THE CURRENT FEASIBLE POINT
        CALL SLVLP (RHS,BNDS,M,NV,LPITMAX,IFLG,IPRT)
C
        IF(IFLG.EQ.0)THEN
C-------SOLUTION FOUND TO PHASE 2 - COMPUTE SOLUTION VECTOR AND OBJECTIVE  
          CALL SOLUTN(RHS,COST,BNDS,W2,M,NV,OBJ)
        ENDIF
      ENDIF
C
C-----DEALLOCATE SPACE USED BY LINEAR PROGRAM SOLVER
      DEALLOCATE(LBNB,IRTIC,IAAO,JAAO,AAO,CAUG,W1,W2,
     &           IPIV,IWORK, A,AF,X,BLA,
     &           R,C,FERR,BERR,WORK,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 993
C
      RETURN
C
C-----ERROR HANDLING
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (LPSUB)')
      CALL GSTOP(' ')
C
  993 CONTINUE
C-----ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
 9930 FORMAT(/,1X,'*** ERROR DEALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (LPSUB)')
      CALL GSTOP(' ')
C
      END SUBROUTINE GWM1SIMPLEX1
C
C
C***********************************************************************
      SUBROUTINE FNDFES(AMAT,COST,BNDS,RHS,M,NV,NDV,LPITMAX,IFLG,IPRT)
C***********************************************************************
C     VERSION: 20FEB2005
C     PURPOSE - PHASE 1 OF PHASE 1/PHASE 2 METHOD:
C               FIND A FEASIBLE SOLUTION USING ARTIFICIAL VARIABLES
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE,ZERO
      INTEGER(I4B),INTENT(IN)::M,NV,NDV,LPITMAX,IPRT
      INTEGER(I4B),INTENT(OUT)::IFLG
      REAL(DP),INTENT(IN)::AMAT(M,NDV),COST(NV),BNDS(NV),RHS(M)
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K,NDVAR,NSLACK
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----COPY FULL A MATRIX INTO SPARSE MATRIX STORAGE - FIRST THE VARIABLES
      K = 0
      JAAO(1) = 1
      NDVAR = NDV-1
      DO 180 J=1,NDVAR
        DO 160 I=1,M  
          IF(AMAT(I,J).NE.ZERO)THEN
            K = K + 1
            AAO(K) = AMAT(I,J)
            IAAO(K) = I
          ENDIF
  160   ENDDO
        JAAO(J+1) = K + 1
  180 ENDDO
C
C-----NEXT COPY SLACKS AND SURPLUS COEFFICIENTS INTO SPARSE MATRIX STORAGE
      NSLACK=0
      DO 200 I=1,M  
        IF(AMAT(I,NDV).NE.ZERO)THEN
          NSLACK = NSLACK + 1
          K = K + 1
          AAO(K) = AMAT(I,NDV)
          IAAO(K) = I
          JAAO(NSLACK+NDVAR+1) = K + 1
        ENDIF
 200  ENDDO
C
C-----TEST INPUT ON NUMBER OF VARIABLES
      IF((NDVAR + NSLACK).NE. NV) THEN
        CALL GSTOP('INPUT ERROR, NV DOES NOT MATCH NDVAR+NSLACK')
      ENDIF
C
C-----ADD ARTIFICIALS TO SPARSE MATRIX STORAGE - IGNORED IN PHASE 2 
      DO 300 I=1,M  
        K = K + 1
        IF(RHS(I).LT.ZERO)THEN
C-------NEGATIVE COEFFICIENTS NEEDED FOR NEGATIVE RHS
          AAO(K) = -ONE
        ELSE
          AAO(K) = ONE
        ENDIF
        IAAO(K) = I
        JAAO(I+NV+1) = K + 1
 300  ENDDO
C
C-----CREATE COST VECTOR AND THE BASIS INDICATOR AND POINTER LISTS 
C
C-----ALL ACTUAL VARIABLES HAVE ZERO COST AND ARE INITIALLY NON-BASIC
      DO 400 I=1,NV
        CAUG(I) = ZERO
        LBNB(I) = 1
 400  ENDDO
C
C-----ALL ARTIFICIAL VARIABLES HAVE COST OF ONE AND ARE INITIALLY BASIC
      K = 0
      DO 500 I=NV+1,NV+M
        CAUG(I) = ONE
        LBNB(I) = 0
        K = K + 1
        IRTIC(K) = I
 500  ENDDO
C
C-----SOLVE THE ALL-ARTIFICIAL PROBLEM TO GET THE INITIAL FEASIBLE SOLUTION
      CALL SLVLP (RHS,BNDS,M,NV,LPITMAX,IFLG,IPRT)
C-----IF SOLUTION FAILED THEN RETURN
      IF(IFLG.NE.0)RETURN
C
C-----CHECK THAT ALL ARTIFICIALS ARE OUT OF THE BASIS           
      DO 600 I=1,M       
        IF(IRTIC(I).GT.NV)THEN
C---------A BASIC IS ARTIFICIAL - PHASE 1 SOLUTION IS NOT FEASIBLE
          IFLG = 1
          RETURN
        ENDIF
 600  ENDDO
C
      IF(IFLG.EQ.0)THEN
C-------PHASE 1 SOLUTION IS FEASIBLE, RESTORE ACTUAL COST VECTOR
        DO 700 I=1,NV
          CAUG(I) = COST(I)
 700    ENDDO
C-------MATRIX STORAGE RETAINS ARTIFICIALS, BUT THESE WILL NEVER BE ACCESSED
      ENDIF
C
      RETURN
      END SUBROUTINE FNDFES
C
C
C***********************************************************************
      SUBROUTINE FACTOR(M)
C***********************************************************************
C     VERSION: 20FEB2005
C     PURPOSE - REASSEMBLE AND FACTORIZE THE BASIS MATRIX
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO
      INTEGER M
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----ZERO THE MATRIX STORAGE
      A=ZERO
C
C---- RECREATE THE BASIS MATRIX FROM SCRATCH 
      DO 200 K=1,M
C-------KTH MATRIX COLUMN IS ASSOCIATED WITH THE JTH VARIABLE
        J = IRTIC(K)
        DO 150 I=JAAO(J),JAAO(J+1)-1
C---------FILL KTH COLUMN OF BASIS WITH JTH COLUMN OF ORIGINAL MATRIX
          A(IAAO(I),K)=AAO(I)
  150   ENDDO
C-------FOR FACTORIZATION W1 IS ARBITRARY, SO ZERO IT
        W1(K) = ZERO
  200 ENDDO
C
C---- PERFORM BASIS FACTORIZATION
      FACT  = 'N'  
      TRANS = 'N'
      CALL SLVBAS(M,W1)
C
      RETURN
      END SUBROUTINE FACTOR
C
C
C***********************************************************************
      SUBROUTINE SLVLP (RHS,BNDS,M,NV,LPITMAX,IFLG,IPRT)
C***********************************************************************
C     VERSION: 25JUL2006
C     PURPOSE - SOLVE A LINEAR PROGRAM FROM THE CURRENT FEASIBLE BASIS
C               USING THE REVISED SIMPLEX METHOD WITH UPPER BOUNDING
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWM1BAS3PS
      INTEGER(I4B),INTENT(IN)::M,NV,LPITMAX,IPRT
      INTEGER(I4B),INTENT(OUT)::IFLG
      REAL(DP),INTENT(IN)::RHS(M),BNDS(NV)
C-----LOCAL VARIABLES
      LOGICAL(LGT)::OPT,UNBND
      INTEGER(I4B)::IT,IENTER,ILEAVE,ILCOL
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----BEGIN ITERATION LOOP
      IT = 0 
  100 CONTINUE
      IT = IT + 1
C
C-----TEST FOR ITERATION LIMIT
      IF(IT.GT.LPITMAX)THEN
        CALL GSTOP('PROGRAM STOPPED: LPITMAX REACHED')
      ENDIF
	IF(MOD(IT,100).EQ.0)THEN  ! PRINT MESSAGE TO SCREEN 
	  IF(IPRT.EQ.1)CALL GWM1BAS3PS('      Begin LP Iteration',IT)   
      ENDIF
C
C-----FACTORIZE BASIS MATRIX 
      CALL FACTOR(M)
C
C-----PRICE OUT VARIABLES TO FIND ENTERING VARIABLE
      CALL PRICE (M,NV,IENTER,OPT)
C
C-----TEST FOR OPTIMALITY
      IF(OPT)THEN
        IFLG = 0
        RETURN
      ENDIF
C 
C---- PERFORM RATIO TEST TO FIND LEAVING VARIABLE
      CALL RATEST(RHS,BNDS,M,NV,IENTER,ILEAVE,UNBND)
C
C-----TEST FOR UNBOUNDEDNESS
      IF(UNBND)THEN
        IFLG = 2
        RETURN
      ENDIF
C
C-----UPDATE BASIS POINTER LISTS
      IF(ILEAVE.EQ.0)THEN                     
C-----THE ENTERING VARIABLE GOES TO ITS BOUND FIRST, SWITCH INDICATOR
        LBNB(ABS(IENTER)) = -LBNB(ABS(IENTER))
      ELSEIF(ILEAVE.GT.0)THEN
C-----LEAVING VARIABLE GOES TO ITS LOWER BOUND
C-------FIND LEAVING BASIC VARIABLE COLUMN AND REPLACE WITH ENTERING VARIABLE
        ILCOL = IRTIC(ILEAVE)
        IRTIC(ILEAVE) = ABS(IENTER)
C-------CHANGE INDICATORS FOR LEAVING AND ENTERING VARIABLES
        LBNB(ILCOL) = 1
        LBNB(ABS(IENTER)) = 0
      ELSEIF(ILEAVE.LT.0)THEN
C-----LEAVING VARIABLE GOES TO ITS UPPER BOUND
C-------FIND LEAVING BASIC VARIABLE COLUMN AND REPLACE WITH ENTERING VARIABLE
        ILCOL = IRTIC(ABS(ILEAVE))
        IRTIC(ABS(ILEAVE)) = ABS(IENTER)
C-------CHANGE INDICATORS FOR LEAVING AND ENTERING VARIABLES
        LBNB(ILCOL) = -1
        LBNB(ABS(IENTER)) = 0
      ENDIF
C
C-----ITERATE
  200 GO TO 100
C
      RETURN
      END SUBROUTINE SLVLP
C
C
C***********************************************************************
      SUBROUTINE PRICE (M,NV,IENTER,OPT)
C***********************************************************************
C   VERSION: 20FEB2005
C   PURPOSE - PERFORMING A PRICING OPERATION TO FIND THE ENTERING COLUMN
C     OUTPUT:  IENTER - COLUMN NUMBER OF VARIABLES THAT WILL ENTER BASIS
C              OPT    - STATUS OF PROBLEM OPTIMALITY
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO
      INTEGER(I4B),INTENT(IN)::M,NV
      INTEGER(I4B),INTENT(OUT)::IENTER
      LOGICAL(LGT),INTENT(OUT)::OPT
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,IUP,IDN
      REAL(DP)::PMAX,PMIN,PRICEV
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----ASSEMBLE VECTOR OF COST COEFFICIENTS ON BASIC VARIABLES
      DO 100 I=1,M
        W2(I)= CAUG(IRTIC(I))
 100  ENDDO
C
C-----COMPUTE DUAL VARIABLES AS BASIC COSTS TIMES BASIS INVERSE
      TRANS = 'T'
      FACT  = 'F'  
      CALL SLVBAS(M,W2)
C
C-----LOOP OVER ALL NON-BASIC VARIABLES LOOKING FOR BEST REDUCED COST
      PMAX = ZERO
      PMIN = ZERO
C-----PRICE ALL ACTUAL VARIABLES - ARTIFICIAL VARIABLES NOT PRICED
      DO 200 I=1,NV
        IF(LBNB(I).EQ.1)THEN
C-------NON-BASIC VARIABLE IS AT ITS LOWER BOUND
          CALL REDCST(W2,M,I,PRICEV)
          IF(PRICEV.LT.PMIN)THEN
C---------INCREASING THIS VARIABLE WILL CAUSE LARGEST REDUCTION IN OBJECTIVE
            PMIN = PRICEV
            IUP = I
          ENDIF
        ELSEIF(LBNB(I).EQ.-1)THEN
C-------NON-BASIC VARIABLE IS AT ITS UPPER BOUND
          CALL REDCST(W2,M,I,PRICEV)
          IF(PRICEV.GT.PMAX)THEN
C---------DECREASING THIS VARIABLE WILL CAUSE LARGEST REDUCTION IN OBJECTIVE
            PMAX = PRICEV
            IDN = I
          ENDIF
        ENDIF
 200  ENDDO
C
      IF(PMAX.EQ.ZERO .AND. PMIN.EQ.ZERO)THEN
C-----NO ENTERING VARIABLE FOUND - CURRENT BASIS IS OPTIMAL
        OPT = .TRUE.
        IENTER = 0
      ELSE
C-----CURRENT BASIS IS NOT OPTIMAL - PICK BEST ENTERING VARIABLE
        OPT = .FALSE.
        IF(PMAX.GE.-PMIN)THEN
          IENTER = -IDN
        ELSEIF(PMAX.LT.-PMIN)THEN
          IENTER = IUP 
        ENDIF
      ENDIF
C
      RETURN
      END SUBROUTINE PRICE
C
C
C***********************************************************************
      SUBROUTINE REDCST(CB,M,I,PRICEV)
C***********************************************************************
C   VERSION: 14AUG2006
C   PURPOSE - COMPUTE THE REDUCED COST
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,SMALLEPS
      INTEGER(I4B),INTENT(IN)::M,I
      REAL(DP),INTENT(IN)::CB(M)
      REAL(DP),INTENT(OUT)::PRICEV
      REAL(DP)::CBN
      INTEGER(I4B)::II
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----COMPUTE PRODUCT OF CB AND ITH COLUMN OF NON-BASIC MATRIX N
      CBN = ZERO
      DO 100 II=JAAO(I),JAAO(I+1)-1
        CBN = CBN + AAO(II)*CB(IAAO(II))
 100  CONTINUE
C-----SUBTRACT FROM COST FOR THIS NON-BASIC VARIABLE
 	PRICEV = CAUG(I) - CBN
C-----SET PRICEV TO ZERO IF IT APPEARS TO BE ROUND-OFF JUNK
      IF(CAUG(I).NE.ZERO)THEN       !COMPARE PRICEV WITH THE ROUNDOFF IN CBN
C-------IF DIFFERENCE LESS THAN ROUND-OFF PORTION OF A TERM, ASSUME IT IS ZERO
        IF(ABS(PRICEV) .LT. ABS(SMALLEPS*CBN))PRICEV = ZERO
      ELSEIF(CAUG(I).EQ.ZERO)THEN   ! COMPARE PRICEV WITH ROUNDOFF OF ZERO
C-------IF DIFFERENCE LESS THAN ROUND-OFF, ASSUME IT IS ZERO
        IF(ABS(PRICEV).LT.SMALLEPS)PRICEV=ZERO
      ENDIF
C
      RETURN
      END SUBROUTINE REDCST
C
C
C*************************************************************************
      SUBROUTINE RATEST(RHS,BNDS,M,NV,IENTER,ILEAVE,UNBND)
C*************************************************************************
C   VERSION: 20SEPT2006
C   PURPOSE - PERFORM THE RATIO TEST TO DETERMINE THE EXITING COLUMN
C     OUTPUT:  ILEAVE - THE BASIC VARIABLE THAT WILL LEAVE BASIS
C              UNBND  - STATUS OF UNBOUNDEDNESS
C--------------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,BIGINF
      INTEGER(I4B),INTENT(IN)::M,NV,IENTER
      INTEGER(I4B),INTENT(OUT)::ILEAVE
      REAL(DP)::RHS(M),BNDS(NV)
      LOGICAL(LGT),INTENT(OUT)::UNBND
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,ILB
      REAL(DP)::XBMIN,XBNEW,BIGINF2,XBTST
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----LOAD CURRENT SOLUTION INTO ARRAY W1 AND ENTERING COLUMN INTO ARRAY W2
      CALL RATIOC(RHS,BNDS,M,NV,IENTER)
C
      BIGINF2 = BIGINF*2.0D0
C-----MAXIMUM NON-BASIC CHANGE ASSUMING NON-BASIC HITS ITS OWN BOUND
      XBMIN  = BNDS(ABS(IENTER))
      ILEAVE = 0
C
      IF(IENTER.GT.0)THEN
C-----NON-BASIC WILL INCREASE
C-------LOOP OVER ALL BASIC VARIABLES TO FIND THE LEAVING COLUMN        
        DO 100 I=1,M
          IF(W2(I).EQ.ZERO)GOTO 100  
          IF(W2(I).GT.ZERO) THEN
C---------MAXIMUM NON-BASIC CHANGE ASSUMING BASIC GOES TO LOWER BOUND
            XBNEW = W1(I)/W2(I)
            ILB   = 1
          ELSEIF(W2(I).LT.ZERO) THEN
C---------MAXIMUM NON-BASIC CHANGE ASSUMING BASIC GOES TO UPPER BOUND
            IF(IRTIC(I).LE.NV) THEN 
C-----------FOR ARTIFICIAL VARIABLES BNDS WILL BE HUGE SO DON'T BOTHER            
              XBNEW = -(BNDS(IRTIC(I)) - W1(I))/W2(I)
              ILB   = -1
C-----------FOR ARTIFICIAL VARIABLES BNDS NOT DEFINED BUT ARE INFINITY            
            ELSE
              XBNEW = BIGINF2
            ENDIF
          ENDIF
C---------COMPARE WITH CURRENT SMALLEST NON-BASIC CHANGE
          IF(XBNEW.LT.XBMIN)THEN
            XBMIN = XBNEW
            ILEAVE = I*ILB
          ENDIF
 100    ENDDO
C
      ELSEIF(IENTER.LT.0)THEN
C-----NON-BASIC WILL DECREASE
C-------LOOP OVER ALL BASIC VARIABLES TO FIND THE LEAVING COLUMN
        DO 200 I=1,M
          IF(W2(I).EQ.ZERO)GOTO 200
          IF(W2(I).LT.ZERO) THEN
C---------MAXIMUM NON-BASIC CHANGE IF BASIC GOES TO LOWER BOUND
            XBNEW = - W1(I)/W2(I)
            ILB   = 1
          ELSEIF(W2(I).GT.ZERO) THEN
C---------MAXIMUM NON-BASIC CHANGE IF BASIC GOES TO UPPER BOUND
            IF(IRTIC(I).LE.NV) THEN   
              XBNEW = (BNDS(IRTIC(I))-W1(I))/W2(I)
              ILB   = -1
C-----------FOR ARTIFICIAL VARIABLES BNDS NOT DEFINED BUT ARE INFINITY            
            ELSE
              XBNEW = BIGINF2
            ENDIF
          ENDIF
C---------COMPARE WITH CURRENT SMALLEST NON-BASIC CHANGE
          IF(XBNEW.LT.XBMIN)THEN
            XBMIN = XBNEW
            ILEAVE = I*ILB
          ENDIF
  200   ENDDO
      ENDIF
C
C-----CHECK FOR AN UNBOUNDED PROBLEM
      UNBND = .FALSE.
      IF(ILEAVE.NE.0)THEN
        IF(IRTIC(ABS(ILEAVE)).GT.NV)THEN
C-------THIS MUST BE AN ARTIFICIAL PROBLEM: IT CAN'T BE UNBOUNDED
          RETURN
        ENDIF
      ENDIF 
      IF(ILEAVE.LT.0)THEN
        IF(BNDS(IRTIC(ABS(ILEAVE))).EQ.BIGINF)THEN  
C-------A BASIC HAS GONE TO ITS UPPER BOUND WHICH IS MACHINE INFINITY
          UNBND = .TRUE.
        ENDIF
      ELSEIF(ILEAVE.EQ.0 .AND. LBNB(ABS(IENTER)).EQ.1 
     1                   .AND. BNDS(ABS(IENTER)).EQ.BIGINF)THEN
C-------THE NON-BASIC WILL GO TO ITS UPPER BOUND WHICH IS MACHINE INFINITY
        UNBND = .TRUE.
      ENDIF
C
      RETURN
      END SUBROUTINE RATEST
C
C
C***********************************************************************
      SUBROUTINE RATIOC(RHS,BNDS,M,NV,IENTER)
C***********************************************************************
C   VERSION: 20FEB2005
C   PURPOSE - COMPUTE CURRENT SOLUTION AND ENTERING COLUMN FOR RATIO TEST
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO
      INTEGER(I4B),INTENT(IN)::M,NV,IENTER
      REAL(DP),INTENT(IN)::RHS(M),BNDS(NV)
C-----LOCAL VARIABLES
      INTEGER I,J,IRHS
      REAL(DP)::BD,RHSTEST
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----FILL W1 WITH CURRENT BASIC SOLUTION
      DO 100 I=1,M  
        W1(I) = ZERO
  100 ENDDO
C
C-----COMPUTE NONBASIC PART - LOAD W1 WITH NXN
      DO 200 I=1,NV
        IF(LBNB(I).EQ.-1)THEN
          BD = BNDS(I)
          DO 150 J=JAAO(I), JAAO(I+1)-1
            W1(IAAO(J)) = W1(IAAO(J)) + AAO(J)*BD
  150     ENDDO
        ENDIF
  200 ENDDO
C           
C-----SUBTRACT FROM RHS - LOAD W1 WITH (B - NXN)
      DO 300 I=1,M 
        W1(I) = RHS(I) - W1(I)
  300 ENDDO
C
C-----MULTIPLY BY BASIS INVERSE - LOAD W1 WITH B-1(B - NXN)
      TRANS = 'N'
      FACT = 'F'  
      CALL SLVBAS(M,W1)
C
C-----FILL W2 WITH B-1(N) FOR ENTERING COLUMN
      DO 400 I=1,M  
        W2(I) = ZERO
  400 ENDDO
C
C-----LOAD W2 WITH ENTERING COLUMN OF N
      DO 500 I=JAAO(ABS(IENTER)), JAAO(ABS(IENTER)+1)-1
        W2(IAAO(I)) = AAO(I)
  500 ENDDO
C
C-----MULTIPLY BY BASIS INVERSE - LOAD W2 WITH B-1(N)
      TRANS = 'N'
      FACT = 'F'  
      CALL SLVBAS(M,W2)
C
      RETURN
      END SUBROUTINE RATIOC
C
C
C***********************************************************************
      SUBROUTINE SOLUTN(RHS,COST,BNDS,DUALS,M,NV,OBJ)
C***********************************************************************
C   VERSION: 20FEB2005
C   PURPOSE - COMPUTE THE SOLUTION VECTOR AND OBJECTIVE FUNCTION
C      INPUT   RHS - THE RIGHT HAND SIDE
C              COST - COST COEFFICIENTS FOR REAL VARIABLES
C      OUTPUT  COST - SOLUTION VECTOR
C              OBJ  - OBJECTIVE VALUE
C              RHS  - VALUE OF DUAL VARIABLES
C              BNDS - REDUCED COST FOR EACH VARIABLE
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO
      USE GWM1DCV3, ONLY : NBVAR
      USE GWM1RMS3, ONLY : CSTROR,NONLIN,RANGEFLG,LASTLP
      INTEGER(I4B),INTENT(IN)::M,NV
      REAL(DP),INTENT(OUT)::OBJ,DUALS(M)
      REAL(DP),INTENT(INOUT)::RHS(M),COST(NV),BNDS(NV)
      INTEGER(I4B)::I,IBC,J
      REAL(DP)::CI
      REAL(DP)::BD
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(RANGEFLG.EQ.1 .AND. LASTLP)THEN         ! PERFORM RANGE ANALYSIS
        CALL CSTRNG(RHS,COST,BNDS,DUALS,M,NV)
        CALL RHSRNG(RHS,COST,BNDS,DUALS,M,NV)
      ELSE                                       ! NO RANGE ANALYSIS
        DO 10 I=1,M
          W1(I) = ZERO                           ! COMPUTE BASIC VARIABLES
   10   ENDDO
        DO 40 I=1,NV
          IF(LBNB(I).EQ.-1)THEN                  ! PUT NON-BASIC PART IN W1
            BD = BNDS(I)
            DO 30 J=JAAO(I),JAAO(I+1)-1          ! LOOP OVER NON-BASIC I
              W1(IAAO(J)) = W1(IAAO(J)) + AAO(J)*BD
   30       ENDDO
          ENDIF
   40   ENDDO
        DO 50 I=1,M
          RHS(I) = RHS(I) - W1(I)                ! SUBTRACT NON-BASIC FROM RHS
   50   ENDDO
        TRANS = 'N'                              ! MULTIPLY BY BASIS INVERSE
        FACT = 'F'
        CALL SLVBAS(M,RHS)                       ! RHS STORES BASIC SOLUTION
      ENDIF
C
C-----COMPUTE THE OBJECTIVE FUNCTION
      OBJ = ZERO
C
C-----LOAD THE NON-BASIC INFO INTO OUTPUT VARIABLES
      DO 100 I=1,NV
        CI = COST(I)
        IF(LBNB(I).EQ.-1)THEN
          OBJ = OBJ + COST(I)*BNDS(I)
          COST(I) = BNDS(I)
        ELSEIF(LBNB(I).EQ.1)THEN
          COST(I) = ZERO
        ENDIF
        BNDS(I) = CSTROR(I) 
        CSTROR(I) = CI
  100 ENDDO
C
C---- ADD ON THE BASIC COLUMNS AND MOVE DUALS INTO THIS VECTOR
      DO 200 I=1,M
        IBC = IRTIC(I)
        IF(IBC.LE.NV)THEN  
          OBJ = OBJ + COST(IBC)*RHS(I)
          COST(IBC) = RHS(I)
        ENDIF
        RHS(I) = DUALS(I)
  200 ENDDO
C
      RETURN
      END SUBROUTINE SOLUTN
C
C
C***********************************************************************
      SUBROUTINE CSTRNG(RHS,COST,BNDS,DUALS,M,NV)
C***********************************************************************
C   VERSION: 20FEB2005
C   PURPOSE - COMPUTE COST COEFFICIENT RANGE ANALYSIS
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE,ZERO,SMALLEPS,BIGINF 
      USE GWM1RMS3, ONLY : CSTREL,CSTRLL,CSTREU,CSTRLU,CSTRLB,CSTRUB,
     &                     CSTROR
      INTEGER(I4B),INTENT(IN)::M,NV
      REAL(DP),INTENT(IN)::RHS(M),COST(NV),BNDS(NV)
      REAL(DP),INTENT(INOUT)::DUALS(M)
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K,KI,KJ,KH,IENL,IENU,IENTER,ILEAVE
      REAL(DP)::CKL,CKU,SUM,RJ,AKJ
      LOGICAL(LGT)::UNBND
c-----LOCAL ARRAY TO TEMPORARILY STORE THE DUALS
      REAL(DP),ALLOCATABLE::TEMPA(:)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----ALLOCATE SPACE FOR TEMPORARY ARRAY AND COPY DUALS
      ALLOCATE(TEMPA(M),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
      DO 10 I=1,M
        TEMPA(I) = DUALS(I)
   10 ENDDO
C
C-----GET A FRESH FACTORIZATION OF THE BASIS
      CALL FACTOR(M)
C
C-----COMPUTE COST COEFFICIENT RANGES FOR NON-BASIC VARIABLES
      DO 100 K=1,NV
        IF(LBNB(K).NE.0)THEN
C-------NON BASIC VARIABLES HAVE EASY RANGE TEST
          CKU = BIGINF
          CKL = -CKU
          SUM = ZERO
          DO 50 J=JAAO(K), JAAO(K+1)-1
            SUM  = SUM + AAO(J)*DUALS(IAAO(J))
   50     ENDDO
          IF(LBNB(K).EQ.1) THEN
            CKL = SUM - COST(K)
            CSTROR(K) = -CKL
            IENL = K
            IENU = 0
          ELSE
            CKU = SUM - COST(K)
            CSTROR(K) = -CKU
            IENL = 0 
            IENU = -K 
          ENDIF
          IF(CKL.EQ.-BIGINF)THEN
            CSTRLB(K) = CKL
          ELSE
            CSTRLB(K) = COST(K) + CKL
          ENDIF
          IF(CKU.EQ. BIGINF)THEN
            CSTRUB(K) = CKU
          ELSE
            CSTRUB(K) = COST(K) + CKU
          ENDIF
          CSTREL(K) = IENL
          CSTREU(K) = IENU
        ENDIF
  100 ENDDO
C
C-----COMPUTE COST COEFFICIENTS FOR BASIC VARIABLES
      KI = 0
      DO 200 K=1,NV
        IF(LBNB(K).EQ.0)THEN
          CSTROR(K) = ZERO
C
C-----LOAD W1 WITH THE KITH ROW OF B INVERSE
          DO 120 I=1,M
            W1(I) = ZERO
  120     ENDDO
          KI = KI + 1
C
C-----FIND ROW OF B-1 CORRESPONDING TO BASIC VAR. K
          DO 140 KJ=1,M
            IF(K .EQ. IRTIC(KJ)) THEN
              KH=KJ
            ENDIF
  140     ENDDO
C
          W1(KH) = ONE
          TRANS = 'T'
          FACT = 'F'  
          CALL SLVBAS(M,W1)
C
C-----SEARCH OVER ALL NON-BASICS FOR EXTREME COST CHANGES
          CKL = -BIGINF
          CKU = -CKL
          IENL = 0
          IENU = 0
          DO 180 J=1,NV
            IF(LBNB(J).EQ.0)GOTO 180
            RJ = CSTROR(J)
            AKJ = ZERO
            DO 160 I=JAAO(J), JAAO(J+1)-1
              AKJ  = AKJ + AAO(I)*W1(IAAO(I))
  160       ENDDO
            IF((AKJ.GT.SMALLEPS.AND.LBNB(J).EQ.1) .OR. (AKJ.LT.-SMALLEPS
     &              .AND. LBNB(J).EQ.-1))THEN
              CALL MINXI(CKU,RJ/AKJ,IENU,J)
              IF(IENU.EQ.J .AND. LBNB(J).EQ.-1)IENU = -IENU
            ELSEIF(AKJ.LT.-SMALLEPS.AND.LBNB(J).EQ.1 .OR.AKJ.GT.SMALLEPS
     &              .AND. LBNB(J).EQ.-1)THEN
              CALL MAXXI(CKL,RJ/AKJ,IENL,J)
              IF(IENL.EQ.J .AND. LBNB(J).EQ.-1)IENL = -IENL
            ENDIF
  180     ENDDO
          IF(CKL.EQ.-BIGINF)THEN
            CSTRLB(K) = CKL
          ELSE
            CSTRLB(K) = COST(K) + CKL
          ENDIF
          IF(CKU.EQ. BIGINF)THEN
            CSTRUB(K) = CKU
          ELSE
            CSTRUB(K) = COST(K) + CKU
          ENDIF
          CSTREL(K) = IENL
          CSTREU(K) = IENU
        ENDIF
  200 ENDDO
C
C-----DETERMINE LEAVING VARIABLES DUE TO CHANGE IN COST COEFFICIENT
      DO 300 K=1,NV
        IENTER = CSTREL(K)
        IF(IENTER.NE.0)THEN
          CALL RATEST(RHS,BNDS,M,NV,IENTER,ILEAVE,UNBND)
          IF(ILEAVE.NE.0)THEN
C---------BASIC VARIABLE HITS ITS BOUND
            CSTRLL(K) = IRTIC(ABS(ILEAVE))
          ELSEIF(.NOT. UNBND)THEN
C---------NON-BASIC ENTERING VARIABLE HITS ITS OWN BOUND
            CSTRLL(K) = ABS(IENTER)
          ENDIF
        ELSE
          CSTRLL(K) = 0           ! ALLOW FOR UNLIMITED COST COEFFICIENT CHANGES
        ENDIF
C
        IENTER = CSTREU(K)
        IF(IENTER.NE.0)THEN
          CALL RATEST(RHS,BNDS,M,NV,IENTER,ILEAVE,UNBND)
          IF(ILEAVE.NE.0)THEN
C---------BASIC VARIABLE HITS ITS BOUND
            CSTRLU(K) = IRTIC(ABS(ILEAVE))
          ELSEIF(.NOT. UNBND)THEN
C---------NON-BASIC ENTERING VARIABLE HITS ITS OWN BOUND
            CSTRLU(K) = ABS(IENTER)
          ENDIF
	  ELSE
          CSTRLU(K) = 0           ! ALLOW FOR UNLIMITED COST COEFFICIENT CHANGES
        ENDIF
        CSTREL(K) = ABS(CSTREL(K))
        CSTREU(K) = ABS(CSTREU(K))
  300 ENDDO
C
C-----COPY DUALS BACK TO ORIGINAL ARRAY AND DEALLOCATE SPACE
      DO 310 I=1,M
        DUALS(I) = TEMPA(I)
  310 ENDDO
      DEALLOCATE(TEMPA,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 993
C
      RETURN
C
C-----ERROR HANDLING
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (CSTRNG)')
      CALL GSTOP(' ')
C
  993 CONTINUE
C-----ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
 9930 FORMAT(/,1X,'*** ERROR DEALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (CSTRNG)')
      CALL GSTOP(' ')
C
      END SUBROUTINE CSTRNG
C
C
C***********************************************************************
      SUBROUTINE RHSRNG(RHS,COST,BNDS,DUALS,M,NV)
C***********************************************************************
C   VERSION: 14AUG2006
C   PURPOSE - COMPUTE THE RIGHT HAND SIDE RANGES
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ONE,ZERO,SMALLEPS,BIGINF
      USE GWM1RMS3, ONLY : CSTROR,RHSREL,RHSRLL,RHSREU,RHSRLU,RHSRLB,
     1                     RHSRUB,RHSROR   
      INTEGER(I4B),INTENT(IN)::M,NV
      REAL(DP),INTENT(IN)::COST(NV),BNDS(NV),DUALS(M)
      REAL(DP),INTENT(INOUT)::RHS(M)
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K,L,KI,KJ,KH,IENT,ILVL,ILVU,IUPRL,IUPRU,ILEAVE
      REAL(DP)::BD,BKL,BKU,RJ,ARJ,CKU
      LOGICAL(LGT)::UPRBND
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----COMPUTE CURRENT SOLUTION
      DO 100 I=1,M  
        W1(I) = ZERO
  100 ENDDO
C
C-----COMPUTE NONBASIC PART
      DO 200 I=1,NV
        IF(LBNB(I).EQ.-1)THEN
          BD = BNDS(I)
          DO 150 J=JAAO(I), JAAO(I+1)-1
            W1(IAAO(J)) = W1(IAAO(J)) + AAO(J)*BD
  150     CONTINUE
        ENDIF
  200 ENDDO
C           
C-----SUBTRACT FROM RHS
      DO 300 I=1,M  
        RHSROR(I) = RHS(I)
        RHS(I) = RHS(I) - W1(I)
  300 ENDDO
C
C---- MULTIPLY BY BASIS INVERSE - RHS NOW CONTAINS THE BASIC SOLUTION
C
      TRANS = 'N'
      FACT = 'F'  
      CALL SLVBAS(M,RHS)
C
C---- COMPUTE RIGHT HAND SIDE RANGES
C
      DO 500 K=1,M
C
C-----LOAD W1 WITH THE KTH COLUMN OF B INVERSE
        DO 350 I=1,M
          W1(I) = ZERO
  350   ENDDO
        W1(K) = ONE
        TRANS = 'N'
        FACT = 'F'  
        CALL SLVBAS(M,W1)
        BKL = -BIGINF
        BKU = -BKL
        ILVL = 0
        ILVU = 0
C
C-------FIND THE LEAVING VARIABLE
        DO 400 I=1,M
C---------FIRST DETERMINE LOWER LIMIT OF W1
          IF(W1(I) .GT. SMALLEPS)THEN
            CALL MAXXI(BKL,-RHS(I)/W1(I),ILVL,I)
            IF(ILVL .EQ. I)THEN
              IUPRL = 0
            ENDIF
          ELSEIF(W1(I) .LT. -SMALLEPS)THEN
            CALL MAXXI(BKL,(BNDS(IRTIC(I))-RHS(I))/W1(I),ILVL,I)
            IF(ILVL .EQ. I)THEN
              IUPRL = 1
            ENDIF
          ENDIF
C
C---------NOW DETERMINE UPPER LIMIT OF W1
          IF(W1(I) .LT. -SMALLEPS)THEN
            CALL MINXI(BKU,-RHS(I)/W1(I),ILVU,I)
            IF(ILVU .EQ. I)THEN
              IUPRU = 0
            ENDIF
          ELSEIF(W1(I) .GT. SMALLEPS)THEN
            CALL MINXI(BKU,(BNDS(IRTIC(I))-RHS(I))/W1(I),ILVU,I)
            IF(ILVU .EQ. I)THEN
              IUPRU = 1
            ENDIF
          ENDIF
  400   ENDDO
C 
        IF(BKL.EQ.-BIGINF)THEN
          RHSRLB(K) = BKL
        ELSE
          RHSRLB(K) = RHSROR(K) + BKL
        ENDIF
        IF(BKU.EQ. BIGINF)THEN
          RHSRUB(K) = BKU
        ELSE
          RHSRUB(K) = RHSROR(K) + BKU
        ENDIF
        RHSRLL(K) = IRTIC(ILVL)
        RHSRLU(K) = IRTIC(ILVU)
C
C-------COMPUTE THE ENTERING VARIABLE
        DO 480 KI=1,2
          UPRBND = .FALSE.
C
C     FIRST DETERMINE ENTERING VARIABLE WHEN LOWER BOUND IS EXCEEDED
C      THEN FIND ENTERING VARIABLE WHEN UPPER BOUND IS EXCEEDED
C
          IF(KI.EQ.1)THEN
            ILEAVE = RHSRLL(K)
            IF(IUPRL.EQ.1)UPRBND=.TRUE.
          ELSE
            ILEAVE = RHSRLU(K)
            IF(IUPRU.EQ.1)UPRBND=.TRUE.
          ENDIF
C
C---------EXTRACT ROW CONTAINING LEAVING VARIABLE FROM B INVERSE
          DO 410 J=1,M
            W1(J)=ZERO
 410      ENDDO
          DO 420 KJ = 1,M
            IF(ILEAVE .EQ. IRTIC(KJ)) THEN
              KH = KJ
            ENDIF
 420      ENDDO
          W1(KH)=ONE
          TRANS = 'T'
          FACT = 'F'  
          CALL SLVBAS(M,W1)
C
C---------PERFORM RATIO TEST ON ALL NONBASICS
          CKU = BIGINF
          IENT=0
          DO 460 J=1,NV
            IF(LBNB(J).EQ.0) GOTO 460
            RJ = CSTROR(J)
            ARJ = ZERO
            DO 440 L=JAAO(J),JAAO(J+1)-1
              ARJ = ARJ + AAO(L)*W1(IAAO(L))
 440        ENDDO
C-----------TEST IS DIFFERENT FOR BASIC LEAVING AT UPPER VS. LOWER BOUND
            IF(.NOT. UPRBND) THEN
C---------- BASIC LEAVES AT LOWER BOUND
              IF(ARJ .LT. -SMALLEPS .AND. LBNB(J) .EQ. 1) THEN
                CALL MINXI(CKU,-RJ/ARJ, IENT, J)
              ELSEIF(ARJ .GT. SMALLEPS .AND. LBNB(J) .EQ. -1) THEN
                CALL MINXI(CKU, -RJ/ARJ, IENT, J)
              ENDIF
            ELSEIF(UPRBND) THEN
C-----------BASIC LEAVES AT UPPER BOUND
              IF(ARJ .GT. SMALLEPS .AND. LBNB(J) .EQ. 1) THEN
                CALL MINXI(CKU, RJ/ARJ, IENT, J)
              ELSEIF(ARJ .LT. -SMALLEPS .AND. LBNB(J) .EQ. -1) THEN
                CALL MINXI(CKU, RJ/ARJ, IENT, J)
              ENDIF
            ENDIF
 460      ENDDO
          IF(KI .EQ. 1)THEN
            RHSREL(K) = IENT
          ELSE
            RHSREU(K) = IENT
          ENDIF
 480    ENDDO
 500  ENDDO
C
      RETURN
      END SUBROUTINE RHSRNG
C
C
C***********************************************************************
      SUBROUTINE MAXXI(XMAX,XNEW,INDEX,I)
C***********************************************************************
      USE GWM1BAS3, ONLY : SMALLEPS
      INTEGER(I4B),INTENT(IN)::I
      INTEGER(I4B),INTENT(OUT)::INDEX
      REAL(DP),INTENT(IN)::XNEW
      REAL(DP),INTENT(INOUT)::XMAX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(XNEW-XNEW*SMALLEPS.GE.XMAX)THEN
C       GIVE A BOOST TO XNEW TO COMPENSATE FOR POSSIBLE ROUND-OFF
        XMAX = XNEW
        INDEX = I
      ENDIF
      RETURN
      END SUBROUTINE MAXXI
C
C
C***********************************************************************
      SUBROUTINE MINXI (XMIN,XNEW,INDEX,I)
C***********************************************************************
      USE GWM1BAS3, ONLY : SMALLEPS
      INTEGER(I4B),INTENT(IN)::I
      INTEGER(I4B),INTENT(OUT)::INDEX
      REAL(DP),INTENT(IN)::XNEW
      REAL(DP),INTENT(INOUT)::XMIN
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(XNEW-XNEW*SMALLEPS.LE.XMIN)THEN
        XMIN = XNEW
        INDEX = I
      ENDIF
      RETURN
      END SUBROUTINE MINXI
C
C
C***********************************************************************
      SUBROUTINE BSOLVE(M,NV,NDV,AMAT,CST,BNDS,RHS,OBJ,IFLG)
C***********************************************************************
C  VERSION: 13JULY2005
C  PURPOSE - CALL AN LP SOLVER REPEATEDLY TO SOLVE MIXED BINARY PROBLEM
C
C    INPUT:
C      LP PROBLEM AND STORAGE -AMAT,CST,BNDS,RHS,XM1,LXM1
C      WELL CHARACTERISTICS -FVON,FVNAME,LWELLSP
C      WORK SPACE FOR MIXED-BINARY LP -DCST,DBND,DRHS,IRLX,IBIN,OBJS
C      BBITMAX     -DIMENSION OF WORK ARRAYS AND MAXIMUM NUMBER OF SUBPROBLEMS
C
C    OUTPUT:
C      IFLG    -SOLUTION STATUS FLAG
C      OBJ     -OPTIMAL SOLUTION TO THE MIXED-BINARY LINEAR PROGRAM
C
C    WORK ARRAYS AND INTERNAL VARIABLESFOR MIXED-BINARY SOLUTION
C      IBIN    -ARRAY WHICH STORES BINARY SETTINGS FOR EACH SOLUTION
C                 EACH COLUMN OF THE ARRAY STORES THE SETTING FOR A SOLUTION
C                 EACH ENTRY IN THE COLUMN GIVES THE WELL SETTING AS FOLLOWS:
C                  1  SET TO ONE 
C                  0  SET TO ZERO
C                 -1  NOT SET AND NOT USED AS BRANCHING VARIABLE
C                 -2  NOT SET BUT HAS BEEN USED AS A BRANCHING VARIABLE
C      IRLX    -STATUS OF EACH SUBPROBLEM
C                  1  UNFATHOMED
C                  0  NOT YET TESTED
C                 -1  FATHOMED BECAUSE IT IS INFEASIBLE 
C                 -2  FATHOMED BECAUSE IT IS RELAXED AND BOUND > INCUMBENT OBJ  
C                 -3  FATHOMED BECAUSE IT IS FEASIBLE INTEGER SOLN
C                 -4  FATHOMED BECAUSE ALL NON-BINARY VARIABLES HAVE BEEN USED
C                       AS A BRANCHING VARIABLE
C      OBJS    -ARRAY OF OBJECTIVE FUNCTION VALUES FOR SUBPROBLEMS
C      DCST,DBND,DRHS - TEMPORARY STORAGE ARRAYS FOR THE LP CST, BND AND RHS
C      INC     -INDEX OF THE BEST CURRENT SOLUTION (THE INCUMBENT SOLUTION)
C      NVAR    -THE WELL TO BE USED FOR BRANCHING IN NEXT PAIR OF SUBPROBLEMS
C
C-----------------------------------------------------------------------
      USE GWM1DCV3, ONLY : NFVAR,NBVAR,FVON,BVLIST,BVNLIST
      USE GWM1BAS3, ONLY : GWMOUT,BIGINF,ZERO
      USE GWM1RMS3, ONLY : BBITMAX,LPITMAX,BBITPRT,LASTLP
      INTEGER(I4B),INTENT(IN)::M,NV,NDV
      REAL(DP),INTENT(INOUT)::BNDS(NV),CST(NV),RHS(M)
      REAL(DP),INTENT(IN)::AMAT(M,NDV)
      REAL(DP),INTENT(OUT)::OBJ
      INTEGER(I4B),INTENT(OUT)::IFLG
C-----LOCAL VARIABLES
      REAL(DP), ALLOCATABLE :: OBJS(:),DCST(:),DBND(:)
      INTEGER(I4B),ALLOCATABLE::IRLX(:),IBIN(:,:)
      INTEGER(I4B)::I,J,II,ION,INC,NVAR,KK,K
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----ALLOCATE MEMORY FOR USE BY BRANCH AND BOUND ALGORITHM
      ALLOCATE (OBJS(BBITMAX),IRLX(BBITMAX),IBIN(NBVAR,BBITMAX),
     &          DCST(NV),DBND(NV),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
C
C-----WRITE HEADER FOR ITERATION FILE
      IF(BBITPRT.EQ.1)THEN
        WRITE(GWMOUT,*,ERR=990)'REPORT BRANCH AND BOUND ITERATIONS'
        WRITE(GWMOUT,*,ERR=990)
      ENDIF
C
C-----CLEAR AND INITIALIZE WORK ARRAYS
      DO 100 I=1,BBITMAX
        IRLX(I) = 0
        OBJS(I) = ZERO
  100 ENDDO
C
C-----SAVE ORIGINAL VALUES OF LP INPUT PARAMETERS
      DO 120 KK=1,NV
        DCST(KK) = CST(KK)
        DBND(KK) = BNDS(KK)
  120 ENDDO
C
C-----INITIALIZE THE INCUMBENT SUBPROBLEM INDEX AND ASSOCIATED OBJECTIVE VALUE
      INC = BBITMAX
      OBJS(INC) = BIGINF
C
C-----INITIALIZE THE COUNTER ON SUBPROBLEM NUMBER
      J = 0
C
C-----INITIALIZE WELL SETTINGS FOR FIRST TWO SUBPROBLEMS/PICK BRANCHING VARIABLE
      DO 160 I=NBVAR,1,-1
        ION = 0
        DO 140 K=1,BVNLIST(I)                    ! LOOP OVER ASSOCIATED VARIABLES
          IF(BVLIST(I,K).GT.NFVAR)THEN           ! AN EXTERNAL VARIABLE IS ASSOCIATED
            ION = 1
          ELSEIF(FVON(BVLIST(I,K)).GT.0)THEN     ! AN ASSOCIATED FLOW VARIABLE IS ON
            ION = 1
          ENDIF
  140   ENDDO
        IF(ION.EQ.0) THEN
C-------SET TO ZERO THOSE WELLS THAT ARE NOT ACTIVE CANDIDATES
          IBIN(I,1) = 0   
          IBIN(I,2) = 0
        ELSE
C-------RELAX ALL ACTIVE WELLS FOR FIRST TWO SUBPROBLEMS
          IBIN(I,1) = -1
          IBIN(I,2) = -1
C---------PICK THE FIRST ACTIVE CANDIDATE WELL AS THE FIRST BRANCHING VARIABLE
          NVAR = I
        ENDIF
  160 ENDDO
C
C-----BEGIN BOUNDING EVALUATION OF TWO SUBPROBLEMS
  200 CONTINUE
C
C-----LOOP OVER THE BRANCHING VARIABLE SETTING IT OFF AND ON 
      DO 300 II=0,1
C
C-------INCREMENT THE SUBPROBLEM NUMBER AND SET THE BRANCHING VARIABLE
        J = J + 1
        IBIN(NVAR,J) = II
C
C-------SOLVE THE LINEAR PROGRAM FOR THIS SUBPROBLEM
        LASTLP = .FALSE.
        CALL GWM1PSIMPLEX1(IBIN,J,M,NV,NDV,AMAT,CST,BNDS,RHS,
     1              OBJ,IFLG,LPITMAX,BBITPRT)
C
C-------PERFORM FATHOMING TEST ON THIS SUBPROBLEM
        CALL FATHOM(J,IFLG,OBJ,INC,OBJS,IRLX,IBIN,CST,NV)
C
C-------PRINT STATUS OF THIS SUBPROBLEM IF SOLUTION EXISTS
        IF(BBITPRT.NE.0 .AND. IFLG.EQ.0) CALL BINWRT(NV,CST,OBJS(J))
C
C-------RETURN LP PARAMETERS TO THEIR ORIGINAL VALUES
        DO 260 KK=1,NV
          BNDS(KK) = DBND(KK)
          CST(KK) = DCST(KK)
  260   ENDDO
C
  300 ENDDO
C    
C-----SELECT NEW BRANCH FOR BRANCH AND BOUND ALGORITHM
      CALL BRANCH(J,IRLX,IBIN,NVAR)
C
C-----NVAR IS THE NEW BRANCHING VARIABLE
      IF (NVAR.NE.0) THEN
C-----A NEW SUBPROBLEM IS POSSIBLE; RETURN TO BOUNDING EVALUATION LOOP
        GOTO 200
      ELSEIF(NVAR.EQ.0)THEN
C-----NO NEW SUBPROBLEMS HAVE BEEN FOUND AVAILABLE
        IF(INC.EQ.BBITMAX) THEN
C-------NO FEASIBLE SOLUTION WAS FOUND; SET IFLG AND EXIT
          IFLG = 1
          GOTO 400
        ELSEIF(INC.NE.BBITMAX)THEN
C-------THE INCUMBENT IS FEASIBLE AND BEST; RERUN OPTIMAL SOLUTION
          WRITE(GWMOUT,*,ERR=990)
          WRITE(GWMOUT,1000,ERR=990)J
          J = INC
C
C---------SOLVE THE LINEAR PROGRAM ONE LAST TIME
          LASTLP = .TRUE.
          CALL GWM1PSIMPLEX1(IBIN,J,M,NV,NDV,AMAT,CST,BNDS,RHS,
     1              OBJ,IFLG,LPITMAX,BBITPRT)
        ENDIF
      ENDIF
C
  400 CONTINUE
C
C-----DEALLOCATE MEMORY USED BY BRANCH AND BOUND ALGORITHM
      DEALLOCATE (OBJS,DCST,DBND,IRLX,IBIN,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 993
C
 1000 FORMAT(I12,' RELAXED PROBLEMS WERE SOLVED.')
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (BSOLVE)')
      CALL GSTOP(' ')
C
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (BSOLVE)')
      CALL GSTOP(' ')
C
  993 CONTINUE
C-----ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
 9930 FORMAT(/,1X,'*** ERROR DEALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (BSOLVE)')
      CALL GSTOP(' ')
C
C
      END SUBROUTINE BSOLVE
C
C
C***********************************************************************
      SUBROUTINE FATHOM(J,IFLG,OBJ,INC,OBJS,IRLX,IBIN,CST,NV)
C***********************************************************************
C  VERSION: 20FEB2005
C  PURPOSE - PERFORM FATHOMING TESTS FOR BRANCH AND BOUND ALGORITHM
C
C    INPUT:
C      J       -NUMBER OF SUBPROBLEM AND INDEX ON SUBPROBLEM STORAGE ARRAYS
C      IFLG    -STATUS FLAG FROM SIMPLEX SOLUTION
C      OBJ     -OBJECTIVE FUNCTION VALUE FROM SOLUTION OF SUBPROBLEM J
C      CST     -OPTIMAL PUMPING RATES FROM LP SOLUTION
C      INC     -THE NUMBER OF THE SUBPROBLEM THAT IS CURRENT BEST 
C
C    OUTPUT:
C      IRLX    -THE FATHOMING STATUS FOR SUBPROBLEM J
C      OBJS    -UPDATED STORAGE OF SUBPROBLEM OBJECTIVE VALUES
C      IBIN    -UPDATED STORAGE OF SUBPROBLEM BINARY SETTINGS
C      INC     -UPDATED NUMBER OF THE SUBPROBLEM THAT IS CURRENT BEST 
C
C   NOTE: AN ADDITIONAL FATHOMING TEST IS PERFORMED IN SUBROUTINE BRANCH 
C-----------------------------------------------------------------------
      USE GWM1DCV3, ONLY : NFVAR,NBVAR,BVLIST,BVNLIST
      USE GWM1BAS3, ONLY : ZERO,BIGINF,GWMOUT
      USE GWM1RMS3, ONLY : BBITMAX,BBITPRT
      INTEGER(I4B),INTENT(IN)::J,IFLG,NV
      INTEGER(I4B),INTENT(INOUT)::INC
      REAL(DP),INTENT(IN)::CST(NV)
      REAL(DP),INTENT(INOUT)::OBJ,OBJS(BBITMAX)
      INTEGER(I4B),INTENT(INOUT)::IBIN(NBVAR,BBITMAX)
      INTEGER(I4B),INTENT(OUT)::IRLX(BBITMAX)
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,K,ION,ISINT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----RESET THE BINARY INDICATORS FOR SUBPROBLEM J BASED ON ACTUAL SOLUTION
C      IF ALL RELAXED VARIABLES ARE ZERO 
      K = 0
C
C-----FIND IF ANY ASSOCIATED VARIABLES THEN SET ION TO 1
      ION = 0
      DO 120 I=1,NBVAR
        IF(IBIN(I,J).LT.0)THEN
C-------BINARY VARIABLE IS RELAXED
          DO 100 K=1,BVNLIST(I)
            IF(CST(BVLIST(I,K)).NE.ZERO)THEN
C-----------AN ASSOCIATED VARIABLE IS ACTIVE
              ION = 1
            ENDIF
  100     ENDDO
        ENDIF
  120 ENDDO
C
  140 IF(ION.EQ.0)THEN
C-----ALL BINARY VARIABLES HAVE ZERO ACTIVITY IN ASSOCIATED VARIABLES SO RESET
        DO 160 I=1,NBVAR
          IF(IBIN(I,J).LT.0) IBIN(I,J)=0
  160   ENDDO
      ENDIF
C
C-----DETERMINE IF THE SUBPROBLEM HAS A BINARY SOLUTION
C      ISINT WILL BE 0 IF SOLUTION IS FEASIBLE AND ALL BINARY
      ISINT = 0
      DO 180 I=1,NBVAR
        IF(IBIN(I,J).LT.0) ISINT = 1
  180 ENDDO
C
C-----FATHOM TEST 1: SUBPROBLEM IS FATHOMED IF IT IS INFEASIBLE
      IF(IFLG.NE.0) THEN
C-------SAVE OBJECTIVE AS INFINITY, RECORD AS FATHOMED AND EXIT
        OBJ = BIGINF
        OBJS(J) = BIGINF
        IRLX(J) = -1
        GOTO 300
      ENDIF
C
C-----FATHOM TEST 2: SUBPROBLEM IS FATHOMED IF IT IS RELAXED AND
C                  ITS BOUND IS GREATER THAN INCUMBENT
      IF(ISINT.EQ.1 .AND. OBJ.GT.OBJS(INC)) THEN
C-------SAVE OBJECTIVE AND RECORD AS FATHOMED
        OBJS(J) = OBJ
        IRLX(J) = -2
        GOTO 300
      ENDIF 
C
C-----FATHOM TEST 3: SUBPROBLEM IS FATHOMED IF IT IS ALL INTEGER
C
      IF(ISINT.EQ.0) THEN
C-----THIS IS AN INTEGER SOLUTION; SAVE OBJECTIVE AND RECORD AS FATHOMED
        OBJS(J) = OBJ
        IRLX(J) = -3
        IF(OBJ.LT.OBJS(INC)) THEN 
C-------THIS IS BEST SOLUTION SO FAR; UPDATE INDEX ON INCUMBENT SOLUTION
          INC = J
C---------IF INCUMBENT SOLUTION CHANGES, REDO FATHOM TEST 2 FOR ALL 
C         PREVIOUS UNFATHOMED SUBPROBLEMS
          DO 200 I=1,J-1
            IF(IRLX(I).EQ.1 .AND. OBJS(I).GE.OBJS(INC))THEN
              IRLX(I) = -2
              WRITE(GWMOUT,1000)I
            ENDIF
  200     ENDDO
        ENDIF
        GOTO 300
      ENDIF
C
C-----ALL FATHOMING TESTS HAVE FAILED:  THIS SUBPROBLEM IS UNFATHOMED
      OBJS(J) = OBJ
      IRLX(J) = 1
C
C-----PRINT STATUS OF SUBPROBLEM
  300 IF(BBITPRT.EQ.1)THEN
        IF(IRLX(J).EQ.1) WRITE(GWMOUT,1001)J
        IF(IRLX(J).EQ.-1)WRITE(GWMOUT,1002)J
        IF(IRLX(J).EQ.-2)WRITE(GWMOUT,1003)J
        IF(IRLX(J).EQ.-3 .AND. INC.EQ.J)WRITE(GWMOUT,1004)J
        IF(IRLX(J).EQ.-3 .AND. INC.NE.J)WRITE(GWMOUT,1005)J
      ENDIF
C
 1000 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - ',
     &                'IT IS RELAXED AND INFERIOR')
 1001 FORMAT(' SUBPROBLEM',I12,' IS UNFATHOMED - IT IS RELAXED')
 1002 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - IT IS INFEASIBLE')
 1003 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - IT IS RELAXED',
     &                ' AND INFERIOR')
 1004 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - ', 
     &                'IT IS INTEGER AND THE NEW BEST SOLUTION')
 1005 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - ',
     &                'IT IS INTEGER AND INFERIOR')
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (FATHOM)')
      CALL GSTOP(' ')
C
C
      END SUBROUTINE FATHOM
C
C
C***********************************************************************
      SUBROUTINE BRANCH(J,IRLX,IBIN,NVAR)
C***********************************************************************
C  VERSION: 20FEB2005
C  PURPOSE - SELECT NEW BRANCH FOR BRANCH AND BOUND ALGORITHM AND SET AS 
C              FATHOMED SUBPROBLEMS WHICH HAVE USED ALL BRANCHING VARIABLES 
C
C    INPUT:
C      J       -NUMBER OF SUBPROBLEM AND INDEX ON SUBPROBLEM STORAGE ARRAYS
C
C    OUTPUT:
C      IBIN    -TWO NEW SUBPROBLEMS COPIED FROM THE BRANCHING SUBPROBLEM
C      NVAR    -THE BINARY VARIABLE TO BE FIXED IN THE NEW SUBPROBLEMS
C                 IF NVAR = 0 THEN NO NEW SUBPROBLEM WAS FOUND 
C      IRLX    -UPDATED FATHOMING STATUS FOR SUBPROBLEMS THAT ARE USED UP
C---------------------------------------------------------------------------
      USE GWM1DCV3, ONLY : NBVAR
      USE GWM1BAS3, ONLY : GWMOUT,ONE
      USE GWM1RMS3, ONLY : BBITMAX
      INTEGER(I4B),INTENT(INOUT)::IRLX(BBITMAX),IBIN(NBVAR,BBITMAX)
      INTEGER(I4B),INTENT(IN)::J
      INTEGER(I4B),INTENT(OUT)::NVAR
C-----LOCAL VARIABLES
      INTEGER(I4B)::ISET,NSUB,I
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----INITIALIZE SUBPROBLEM COUNTER FOR INCREMENTING BACKWARDS
      NSUB = J + 1
C    
C-----FIND THE MOST RECENTLY GENERATED UNFATHOMED SUBPROBLEM
  100 NSUB = NSUB - 1
C
      IF(IRLX(NSUB).EQ.1) THEN
C-----SUBPROBLEM IS UNFATHOMED
C-------FIND THE FIRST UNBRANCHED NON-BINARY VALUED WELL AND JUMP OUT OF LOOP
        DO 120 NVAR=1,NBVAR
          IF(IBIN(NVAR,NSUB).EQ.-2) GO TO 140
          IF(IBIN(NVAR,NSUB).EQ.-1) GO TO 180
  120   ENDDO
C
C-------NO UNBRANCHED NON-BINARY VALUED WELLS AVAILABLE, 
C        SO THIS SUBPROBLEM IS ACTUALLY FATHOMED - RESET IRLX AND CONTINUE
        IRLX(NSUB) = -4
        WRITE(GWMOUT,1001,ERR=990)NSUB
        GOTO 160
  140   CONTINUE
C-------HAVE ALREADY BRANCHED FROM THIS SUBPROBLEM
        IRLX(NSUB) = -4
        WRITE(GWMOUT,1002,ERR=990)NSUB
      ENDIF

  160 IF(NSUB.GT.1)GOTO 100
C
C-----NO UNFATHOMED SUBPROBLEMS HAVE BEEN FOUND
      NVAR = 0
      GOTO 300
C
C-----NEW UNFATHOMED SUBPROBLEM HAS BEEN FOUND
  180 CONTINUE
C
C-----TEST THAT STORAGE SPACE IS AVAILABLE FOR NEXT SUBPROBLEM
      IF(J+2.GT.BBITMAX)
     &          CALL GSTOP('STOP: BRANCH AND BOUND ITERATIONS EXCEEDED')
C
C-----INDICATE THAT THIS NON-BINARY VALUED WELL HAS BEEN BRANCHED FROM
      IBIN(NVAR,NSUB) = -2
C
C-----COPY THE SETTINGS OF UNFATHOMED SUBPROBLEM TO THE TWO NEW SUBPROBLEMS
      DO 200 I=1,NBVAR
        ISET = IBIN(I,NSUB)
        IBIN(I,J+1) = ISET
        IBIN(I,J+2) = ISET
  200 ENDDO
      WRITE(GWMOUT,*,ERR=990)'NEXT: BRANCH FROM SUBPROBLEM',NSUB
      WRITE(GWMOUT,*,ERR=990)'                AND VARIABLE',NVAR
C
 1001 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - ',
     &                         'ALL BRANCHING VARIABLES USED')
 1002 FORMAT(' SUBPROBLEM',I12,' IS FATHOMED - ',
     &                         'ALREADY BRANCHED FROM IT ONCE')
C
  300 RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (BRANCH)')
      CALL GSTOP(' ')
C
C
      END SUBROUTINE BRANCH
C
C
C***********************************************************************
      SUBROUTINE BINWRT(NV,CST,OBJ)
C***********************************************************************
C   VERSION: 20FEB2005
C   PURPOSE - WRITE OUTPUT FOR THE CURRENT SUBPROBLEM
C-----------------------------------------------------------------------
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,FVNAME,EVNAME,BVNAME
      USE GWM1BAS3, ONLY : ZERO,GWMOUT
      USE GWM1BAS3, ONLY : GWM1BAS3PF
      INTEGER(I4B),INTENT(IN)::NV
      REAL(DP),INTENT(IN)::CST(NV),OBJ
C-----LOCAL VARIABLES
      INTEGER(I4B)::K,KK
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(GWMOUT,1000,ERR=990)
      WRITE(GWMOUT,2000,ERR=990)'Flow Variables'
      WRITE(GWMOUT,3000,ERR=990)(FVNAME(K),CST(K),K=1,NFVAR)
      IF(NEVAR.GT.0)THEN
        KK=NFVAR
        WRITE(GWMOUT,2000,ERR=990)'External Variables'
        WRITE(GWMOUT,3000,ERR=990)(EVNAME(K),CST(K+KK),K=1,NEVAR)
      ENDIF
      IF(NBVAR.GT.0)THEN
        KK=NFVAR+NEVAR
        WRITE(GWMOUT,2000,ERR=990)'Binary Variables'
        WRITE(GWMOUT,3000,ERR=990)(BVNAME(K),CST(K+KK),K=1,NBVAR)
      ENDIF
C
      CALL GWM1BAS3PF('    Objective Function =',1,OBJ)
      CALL GWM1BAS3PF(' ',0,ZERO)
C
 1000 FORMAT(3(4X'Name',10X,'Value',3X))
 2000 FORMAT(T1,A)
 3000 FORMAT(3(4X,A10,ES12.4))
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (BINWRT)')
      CALL GSTOP(' ')
C
C
      END SUBROUTINE BINWRT
C
C
C***********************************************************************
      SUBROUTINE GWM1PSIMPLEX1(IBIN,J,M,NV,NDV,AMAT,CST,BNDS,RHS,
     1                         OBJ,IFLG,LPITMAX,IPRT)
C***********************************************************************
C  VERSION: 13JULY2005
C  PURPOSE - BUILD A TEMPORARY LP PROBLEM THAT ELIMINATES BINARY 
C            VARIABLES THAT ARE SET IN BRANCH AND BOUND OR IMPLICITLY SET.   
C            IF A BINARY VARIABLE IS SET TO ZERO OR ONE IT IS REMOVED
C            ALONG WITH ITS DEFINITION CONSTRAINT OF FORM -Q + MX >=0
C---------------------------------------------------------------------------
      USE GWM1DCC3, ONLY : NLBR
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,BVNLIST,BVLIST
      USE GWM1BAS3, ONLY : ZERO,ONE,BIGINF
      USE GWM1RMS3, ONLY : BBITMAX,NCONF,NVF,
     &                     RHSIN,RHSINF,RANGENAME,RANGENAMEF,LASTLP
      INTEGER(I4B),INTENT(IN)::IBIN(NBVAR,BBITMAX)
      INTEGER(I4B),INTENT(INOUT)::J
      INTEGER(I4B),INTENT(IN)::M,NV,NDV,LPITMAX,IPRT
      INTEGER(I4B),INTENT(OUT)::IFLG
      REAL(DP),INTENT(IN)::AMAT(M,NDV)
      REAL(DP),INTENT(OUT)::OBJ
      REAL(DP),INTENT(INOUT)::CST(NV),BNDS(NV),RHS(M)
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,II,III,INUM,K,KK,NDROP,NDVTMP,MTMP,NVTMP,IT
      INTEGER(I4B)::NCVAR,ROWEMPTY,IINF,REMSLACK
      REAL(DP)    ::LSUM,CMIN,CMAX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----DETERMINE THE BINARY VARIABLES THAT NEED TO BE DROPPED
      ALLOCATE (IDROP(NBVAR),STAT=ISTAT)         ! CREATE WORK SPACE
      IF(ISTAT.NE.0)GOTO 992
C-----THE LP ENTERS WITH SOME BINARY VARIALBLES SET.  ALL-BINARY CONSTRAINTS 
C       MAY IMPLICITY SET SOME RELAXED BINARY VARIABLES.   THIS CODE LOOKS FOR THIS 
C       CONDITION AND SETS THE RELAXED BINARY VARIABLES TO ZERO.  
C       FOR CONSTRAINTS OF FORM SUM(a*X) <= N, IF THE SET VARIABLES CAUSE THE 
C         CONSTRAINT TO BIND AND THE a COEFFICIENTS ARE POSITIVE, RELAXED BV MUST BE ZERO
C       FOR CONSTRAINTS OF FORM SUM(a*X) >= N, IF THE SET VARIABLES CAUSE THE 
C         CONSTRAINT TO BIND AND THE a COEFFICIENTS ARE NEGATIVE, RELAXED BV MUST BE ZERO
C       FOR CONSTRAINTS OF FORM SUM(a*X) = N, IF THE SET VARIABLES SATISFY THE 
C         CONSTRAINT THE RELAXED BINARY VARIABLES CAN ONLY TAKE A VALUE OF ZERO
C 
      DO 100 I=1,NBVAR
        IDROP(I) = IBIN(I,J)                     ! COPY BINARY STATUS
  100 ENDDO
C-----LOOK AT ALL BUT LAST NBVAR CONSTRAINTS FOR ROWS WITH ONLY BINARY VARIABLES
      NCVAR  = NFVAR+NEVAR                       ! NUMBER NON-BINARY VARIABLES
      DO 150 I=1,M-NBVAR
        DO 110 K=1,NCVAR
          IF(AMAT(I,K).NE.ZERO)GOTO 150          ! THIS ROW HAS NON-BINARY
  110   ENDDO
C-------IF YOU GOT HERE THEN THE CONSTRAINT HAS ONLY BINARY VARIABLES
        LSUM = ZERO                              ! ACCUMULATE FIXED LHS
        CMIN = BIGINF                            ! FIND SMALLEST COEFFICIENT
        CMAX = -CMIN                             ! FIND LARGEST COEFFICIENT
        DO 120 K=1,NBVAR                         ! LOOP OVER BINARY COLUMNS
          KK = K + NCVAR
          IF(IDROP(K).EQ.1)LSUM = LSUM + AMAT(I,KK) ! ADD UP COEFFICIENTS 
          IF(IDROP(K).LT.0)THEN
            CMIN = MIN(CMIN,AMAT(I,KK))          ! FOR RELAXED BINARIES
            CMAX = MAX(CMAX,AMAT(I,KK))          ! FOR RELAXED BINARIES
          ENDIF
  120   ENDDO
C-------TEST FOR CONDITIONS THAT CAUSE RELAXED BV TO BE IMPLICITY SET
        III = 0
        IF(AMAT(I,NDV).EQ.ONE)THEN               ! CONSTRAINT IS <=
C---------IF CONSTRAINT IS MET AS = AND NO NEGATIVE COEFF. ON RELAXED BINARIES
          IF(LSUM.EQ.RHS(I).AND.CMIN.GE.ZERO)III=1  ! SET FLAG         
        ELSEIF(AMAT(I,NDV).EQ.-ONE)THEN          ! CONSTRAINT IS >=
C---------IF CONSTRAINT IS MET AS = AND NO POSITIVE COEFF. ON RELAXED BINARIES
          IF(LSUM.EQ.RHS(I).AND.CMAX.LE.ZERO)III=1  ! SET FLAG         
        ELSEIF(AMAT(I,NDV).EQ.ZERO)THEN          ! CONSTRAINT IS =
C---------IF CONSTRAINT IS MET AS = 
          IF(LSUM.EQ.RHS(I))III=1                ! SET FLAG
        ENDIF
        IF(III.EQ.1)THEN                       
C---------CONDITION IS MET: ALL RELAXED BINARIES IN THE CONSTRAINT MUST BE ZERO 
          DO 130 K=1,NBVAR 
            IF(AMAT(I,K+NCVAR).GT.ZERO .AND. IDROP(K).LT.0)THEN
              IDROP(K) = 0                     ! SET RELAXED BINARY TO ZERO
            ENDIF
  130     ENDDO
        ENDIF
  150 ENDDO
C
C-----COUNT THE NUMBER OF BINARY VARIABLES THAT ARE TO BE DROPPED
      NDROP = 0
      DO 180 I=1,NBVAR
        IF(IDROP(I).GE.0)NDROP = NDROP + 1
  180 ENDDO
C
C-----DEFINE DIMENSIONS OF TEMPORARY PROBLEM AND ALLOCATE SPACE
      NDVTMP = NDV - NDROP                       ! SHIFT LAST COLUMN
      MTMP   = M   - NDROP                       ! DROP DEFINITION CONSTRAINT
      NVTMP  = NV  - NDROP*2                     ! DROP BV & DEF CONST. SLACK
      ALLOCATE (TAMAT(MTMP,NDVTMP),TCST(NVTMP),TBNDS(NVTMP),TRHS(MTMP),
     &          ROWCHK(MTMP),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
      TAMAT = ZERO
      TCST  = ZERO
      TBNDS = ZERO
      TRHS  = ZERO
C
C-----LOAD TEMPORARY CST AND BNDS ARRAYS FROM FULL PROBLEM 
C
C-----LOOP OVER NON-BINARY VARIABLES:THESE ARE UNCHANGED
      DO 200 I=1,NCVAR
        TCST(I)=CST(I)
        TBNDS(I)=BNDS(I)
  200 ENDDO
C
C-----LOOP OVER BINARY VARIABLES SKIPPING THOSE THAT ARE SET
      K = NCVAR + 1
      DO 210 I=NCVAR+1,NCVAR+NBVAR
        IF(IDROP(I-NCVAR).LT.0)THEN        ! BINARY VARIABLE STAYS IN PROBLEM
          TCST(K)=CST(I)
          TBNDS(K)=BNDS(I)
          K = K + 1
        ENDIF
  210 ENDDO
C
C-----LOOP OVER ALL SLACK AND SURPLUS VARIABLES EXCEPT -Q + MX >=0 
      DO 220 I=NDV,NV-NBVAR
        TCST(K)=CST(I)
        TBNDS(K)=BNDS(I)
        K = K + 1
  220 ENDDO
C
C-----LOOP OVER SLACKS FOR -Q + MX >=0 CONSTRAINTS, SKIP SET CONSTRAINTS
      DO 225 I=NV-NBVAR + 1,NV
        IF(IDROP(I-(NV-NBVAR)).LT.0)THEN
C-------BINARY VARIABLE STAYS IN PROBLEM, KEEP ITS CONSTRAINT
          TCST(K)=CST(I)
          TBNDS(K)=BNDS(I)
          K = K + 1
        ENDIF
  225 ENDDO
C
C-----LOAD TEMPORARY MATRIX AND RIGHT HAND SIDE FROM FULL PROBLEM
C
C-----FIRST, FILL ALL BUT LAST NBVAR ROWS
      DO 320 I=1,M-NBVAR
        TRHS(I) = RHS(I)
C-------LOAD NON-BINARY VARIABLE COLUMNS
        DO 300 K=1,NCVAR
          TAMAT(I,K) = AMAT(I,K)
  300   ENDDO
C-------CONSIDER COLUMNS WITH BINARY VARIABLES
        KK = NCVAR + 1
        DO 310 K=NCVAR+1,NCVAR+NBVAR
          IF(IDROP(K-NCVAR).LT.0)THEN
C---------BINARY IS ACTIVE IN TEMPORARY PROBLEM
            TAMAT(I,KK) = AMAT(I,K)
            KK = KK + 1
          ELSEIF(IDROP(K-NCVAR).EQ.0)THEN
C---------BINARY IS TURNED OFF, DON'T DO ANYTHING
          ELSEIF(IDROP(K-NCVAR).EQ.1)THEN
C---------BINARY IS TURNED ON, PLACE ITS COEFFICIENT IN RHS
            TRHS(I) = TRHS(I) - AMAT(I,K)
          ENDIF
  310   ENDDO
C-------PLACE THE SLACK COEFFICIENT IN LAST COLUMN OF TEMPORARY PROBLEM
        TAMAT(I,KK) = AMAT(I,NDV)
  320 ENDDO
C
C-----SECOND, CONSIDER LAST NBVAR ROWS WHICH CONTAIN -Q + MX <=0 CONSTRAINTS
      II = M - NBVAR + 1
      DO 360 I=M-NBVAR + 1,M
        INUM = I - (M-NBVAR)
        IF(IDROP(INUM).LT.0)THEN           ! BINARY VARIABLE STAYS IN PROBLEM
          TRHS(II)=RHS(I)                   ! LOAD ITS CONSTRAINT
C---------LOAD NON-BINARY VARIABLE COLUMNS
          DO 330 K=1,NCVAR
            TAMAT(II,K) = AMAT(I,K)
  330     ENDDO
C---------CONSIDER COLUMNS WITH BINARY VARIABLES
          KK = NCVAR + 1
          DO 340 K=NCVAR+1,NCVAR+NBVAR
            IF(IDROP(K-NCVAR).LT.0)THEN
C-----------BINARY IS ACTIVE IN TEMPORARY PROBLEM
              TAMAT(II,KK) = AMAT(I,K)
              KK = KK + 1
            ELSEIF(IDROP(K-NCVAR).EQ.0)THEN
C-----------BINARY IS TURNED OFF, DON'T DO ANYTHING
            ELSEIF(IDROP(K-NCVAR).EQ.1)THEN
C-----------BINARY IS TURNED ON, PLACE ITS COEFFICIENT IN RHS
              TRHS(II) = TRHS(II) - AMAT(I,K)
            ENDIF
  340     ENDDO
C---------PLACE THE SLACK COEFFICIENT IN LAST COLUMN
          TAMAT(II,KK) = AMAT(I,NDV)
          II = II + 1
        ELSEIF(IDROP(INUM).EQ.0)THEN
C-------BINARY IS TURNED OFF, REPLACE CONSTRAINT WITH AN UPPER BOUND 
C         ON THIS VARIABLE TO TURN OFF WELL
          DO 350 III=1,BVNLIST(INUM)
            TBNDS(BVLIST(INUM,III)) = ZERO
  350    ENDDO
        ENDIF
  360 ENDDO
C
C-----ELIMINATE ROWS THAT ARE EMPTY
      ROWEMPTY = 0                               ! NUMBER OF EMPTY ROWS
      REMSLACK = 0                               ! NUMBER OF SLACKS
      IINF = 0
      DO 375 I=1,MTMP
        ROWCHK(I)=1
        DO 370 K=1,NDVTMP-1                      ! SEARCH NON-SLACK COLUMNS
          IF(TAMAT(I,K).NE.ZERO)GOTO 375
  370   ENDDO
C-----IF YOU GOT HERE THEN THE ROW IS EMPTY
        ROWEMPTY = ROWEMPTY+1
        ROWCHK(I)=0
C-------COUNT THE SLACK IF NOT AN EQUALITY CONSTRAINT
        IF(TAMAT(I,NDVTMP).NE.ZERO)REMSLACK = REMSLACK + 1
C-------CHECK IF EMPTY CONSTRAINT ROW IS FEASIBLE
        IF(TAMAT(I,NDVTMP).EQ.-ONE .AND. TRHS(I).GT.ZERO .OR.  ! 0 /> RHS
     &     TAMAT(I,NDVTMP).EQ.ZERO .AND. TRHS(I).NE.ZERO .OR.  ! 0 /= RHS
     &     TAMAT(I,NDVTMP).EQ. ONE .AND. TRHS(I).LT.ZERO)THEN  ! 0 /< RHS
          IINF = 1                               ! ROW IS NOT FEASIBLE
        ENDIF
  375 ENDDO
C
      IF(ROWEMPTY.GT.0.AND.IINF.EQ.0)THEN
C-------EMPTY CONSTRAINTS ROWS ARE PRESENT AND MET; REARRANGE MATRIX
        ALLOCATE(TAMAT2(MTMP-ROWEMPTY,NDVTMP),STAT=ISTAT)
        IF(ISTAT.NE.0)GOTO 992
        IT = 0 
        DO 390 I=1,MTMP
          IF(ROWCHK(I).EQ.1)THEN                 ! THIS ROW STAYS
	      IT = IT+1                            ! INCREMENT NEW ROW COUNTER
            DO 380 K=1,NDVTMP
              TAMAT2(IT,K)=TAMAT(I,K)            ! LOAD ROW COEFFICIENTS       
  380       ENDDO
            TRHS(IT) = TRHS(I)                   ! MOVE RHS INTO ROW
            RHSINF(IT) = RHSIN(I)                ! LOAD RHSIN FOR OUTPUT
            RANGENAMEF(IT+NDV-1)=RANGENAME(I+NDV-1)! LOAD RANGENAME 
          ENDIF
  390   ENDDO
        MTMP  = MTMP-ROWEMPTY                    ! DECREMENT ROW COUNTER
        NVTMP = NVTMP-REMSLACK                   ! REMOVE EMPTY ROW SLACKS
C-------SOLVE THE LINEAR PROGRAM 
        CALL GWM1SIMPLEX1(MTMP,NVTMP,NDVTMP,TAMAT2,TCST,
     &                TBNDS,TRHS,OBJ,IFLG,LPITMAX,IPRT)
C
      ELSEIF(ROWEMPTY.GT.0.AND.IINF.EQ.1)THEN
C-------ONE OR MORE EMPTY CONSTRAINT ROWS ARE INFEASIBLE; DON'T BOTHER SOLVING
        IFLG = 1                                 ! SET SIMPLEX FLAG 
C
      ELSEIF(ROWEMPTY.EQ.0)THEN
C-------SOLVE THE LINEAR PROGRAM 
        CALL GWM1SIMPLEX1(MTMP,NVTMP,NDVTMP,TAMAT,TCST,
     &                    TBNDS,TRHS,OBJ,IFLG,LPITMAX,IPRT)
      ENDIF
C
      IF(IFLG.EQ.0)THEN                          ! SOLUTION IS AVAILABLE                     
C-------LOAD SOLUTION OF TEMPORARY PROBLEM INTO FULL PROBLEM LOCATIONS
        DO 400 I=1,NCVAR                         ! LOAD NON-BINARY VARIABLES
          CST(I) = TCST(I)
          BNDS(I)= TBNDS(I)
  400   ENDDO
C
        K = NCVAR + 1
        DO 410 I=NCVAR+1,NCVAR+NBVAR             ! LOAD BINARY VARIABLES 
          IF(IDROP(I-NCVAR).LT.0)THEN           ! BINARY WAS IN TEMP PROBLEM
            CST(I) = TCST(K)
            K = K + 1
          ELSEIF(IDROP(I-NCVAR).EQ.1)THEN       ! BINARY VARIABLE WAS SET TO 1:   
            OBJ = OBJ + CST(I)                   ! ADD ITS COEFFICIENT TO OBJ
            CST(I) = ONE
          ELSEIF(IDROP(I-NCVAR).EQ.0)THEN       ! BINARY VARIABLE WAS SET TO 0
            CST(I) = ZERO
          ENDIF
  410   ENDDO
      ENDIF
C
C-----BNDS AND RHS ARE IGNORED DURING BRANCH AND BOUND ITERATIONS
C       BUT HAVE MEANING ON THE LAST CALL TO GWM1PSIMPLEX1
      IF(LASTLP)THEN 
        NCONF = MTMP-NLBR                        ! SET FINAL NUMBER OF CONSTRAINTS   
        NVF   = NDVTMP-1                         ! SET FINAL NUMBER OF VARIABLES   
        DO 420 I=NCVAR+1,NCVAR+NBVAR
          BNDS(I)=ZERO                           ! REDUCED COSTS ARE NOT CALCULATED
  420   ENDDO
        DO 430 I=NCVAR+NBVAR+1,NVTMP
          CST(I) = TCST(I-NBVAR)                 ! LOAD SLACK VALUES
  430   ENDDO
        I = 0
        DO 440 J=1,NCONF+ROWEMPTY                ! LOAD RHS FOR ALL NON-BINARY ROWS
          IF(ROWCHK(J).EQ.1)THEN                 ! ROW WAS IN THE LP
            I = I+1
            RHS(J)=TRHS(I)                       ! LOAD COMPUTED SHADOW PRICE
          ELSEIF(ROWCHK(J).EQ.0)THEN             ! ROW WAS NOT IN LP
            RHS(J)=BIGINF*2.0D0                  ! LOAD A FLAG TO INDICATE THIS 
          ENDIF
  440   ENDDO
      ENDIF
C 
      DEALLOCATE(TAMAT,TCST,TBNDS,TRHS,IDROP,ROWCHK,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 993
      IF(ROWEMPTY.GT.0.AND.IINF.EQ.0) DEALLOCATE(TAMAT2,STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 993
C
      RETURN
C
C-----ERROR HANDLING
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (LPSUB0)')
      CALL GSTOP(' ')
C
  993 CONTINUE
C-----ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
 9930 FORMAT(/,1X,'*** ERROR DEALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (LPSUB0)')
      CALL GSTOP(' ')
C
C
      END SUBROUTINE GWM1PSIMPLEX1
C
C
C***********************************************************************
      SUBROUTINE SLVBAS(N,BL)
C***********************************************************************
C   VERSION: 14AUG2006
C   PURPOSE - FACTOR OR SOLVE THE SYSTEM OF EQUATIONS DEFINED BY BASIS
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWMOUT,ZERO
      INTEGER(I4B),INTENT(IN)::N
      REAL(DP),INTENT(INOUT)::BL(N)
      INTERFACE 
        SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
     $                   WORK, IWORK, INFO )
        CHARACTER          EQUED, FACT, TRANS
        INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
        DOUBLE PRECISION   RCOND
        INTEGER            IPIV( * ), IWORK( * )
        DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), C( * ), FERR( * ), R( * ),
     $                   WORK( * ), X( LDX, * )
        END
C
      END INTERFACE
C-----LOCAL VARIABLE
      INTEGER(I4B)::I,INFO
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----IF MATRIX HAS NO ROWS RETURN
      IF(N.EQ.0)RETURN
      X = ZERO                                   ! BLANK WORK ARRAYS
      WORK = ZERO
      IWORK = 0
C-----LOAD RIGHT HAND SIDE INTO LAPACK STORAGE
      DO 100 I=1,N
        BLA(I,1) = BL(I)
  100 ENDDO
C
      IF(FACT.EQ.'N' .OR. FACT.EQ.'E')THEN
C-----FACTORIZE THE MATRIX USING LU DECOMPOSITION
        AF = ZERO                                ! BLANK WORK ARRAYS
        IPIV = 0
        CALL DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
     1             EQUED, R, C, BLA, LDB, X, LDX, RCOND, FERR, BERR,
     2             WORK, IWORK, INFO )
C-------FACTORIZATION RETURNS A SOLUTION BUT IGNORE IT
C
      ELSEIF(FACT.EQ.'F')THEN
c-----SOLVE THE FACTORIZED MATRIX FOR THE STORED RIGHT HAND SIDE
        CALL DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
     1             EQUED, R, C, BLA, LDB, X, LDX, RCOND, FERR, BERR,
     2             WORK, IWORK, INFO )
C
C-------MOVE THE SOLUTION FROM LAPACK STORAGE
        DO 200 I=1,N
          BL(I)=X(I,1)  
  200   ENDDO     
      ENDIF
C
C-----CHECK LAPACK SOLUTION STATUS FLAG
      IF(INFO.EQ.0) RETURN
C
      IF(INFO.GT.0 .AND. INFO.LE.N)THEN
        WRITE(GWMOUT,1000,ERR=990)INFO
      ELSEIF(INFO.GT.N)THEN
        WRITE(GWMOUT,2000,ERR=990)
      ELSEIF(INFO.LT.0)THEN
        WRITE(GWMOUT,3000,ERR=990)INFO
      ENDIF
C
      CALL GSTOP(' ')
C
 1000 FORMAT('LAPACK FAILURE: UPPER DIAGONAL IS EXACTLY ZERO IN ROW',I6
     &    ,/,'  LP BASIS MATRIX IS APPARENTLY SINGULAR.'
     &    ,/,'  CHECK THAT CONSTRAINTS ARE NOT REDUNDANT.'
     &    ,/,'  TRY REMOVING SELECTED CONSTRAINTS.')
 2000 FORMAT('LAPACK FAILURE: NONSINGULAR TO MACHINE PRECISION')
 3000 FORMAT('LAPACK FAILURE: ILLEGAL VALUE IN DGESVX ARGUMENT',I6)
C
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (SLVBAS)')
      CALL GSTOP(' ')
C
      END SUBROUTINE SLVBAS
C
C
      END MODULE GWM1RMS3LP

C
C
C
C
      MODULE GWM1RMS3SUBS
      USE GWM_SUBS, ONLY: IGETUNIT, URWORD2
      USE GWM_STOP, ONLY:   GSTOP
C     VERSION: 21MAR2008 - package references converted to version 3
      IMPLICIT NONE
      PRIVATE
      PUBLIC::GWM1RMS3AR,GWM1RMS3PL,GWM1RMS3PP,GWM1RMS3FP,GWM1RMS3FM,
     1        GWM1RMS3AP,GWM1RMS3MPS,GWM1RMS3SLP,GWM1RMS3LP_CHKSOL,
     2        GWM1RMS3OT,GWM1RMS3OS
C
      INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
      INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
      INTEGER, PARAMETER :: SP = KIND(1.0)
      INTEGER, PARAMETER :: DP = KIND(1.0D0)
      INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
C-----FOR ERROR HANDLING
      INTEGER(I2B)::ISTAT  
      CHARACTER(LEN=200)::FLNM
      CHARACTER(LEN=20)::FILACT,FMTARG,ACCARG,FILSTAT
      INTEGER(I4B)::NDUM,FUNIT
      REAL(SP)::RDUM
C
      CONTAINS
C***********************************************************************
      SUBROUTINE GWM1RMS3AR(FNAME,IOUT,NFVAR,NEVAR,NBVAR,NDV,NV,NCON)
C***********************************************************************
C   VERSION: 30MAY2011
C   PURPOSE: READ INPUT FROM THE SOLUTION AND OUTPUT-CONTROL FILE
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,RMFILE,RMFILEF,MPSFILE,CUTCOM
      USE GWM1DCV3, ONLY : FVINI,FVDIR,FVNAME,FVBASE,EVNAME,EVBASE,
     &                     BVNAME,BVBASE
      USE GWM1OBJ3, ONLY : SOLNTYP    
      USE GWM1DCC3, ONLY : NLBR
      USE GWM1RMS3, ONLY : SLPVCRIT,SLPZCRIT,BBITPRT,LPITMAX,NONLIN,
     1                     SLPITPRT,SLPITMAX,BBITMAX,NSIGDIG,NPGNMX,
     2                     AFACT,PGFACT,NINFMX,DELINC,FVOLD,CRITMFC,
     3                     MFCNVRG,AMAT,CST,BNDS,RHS,IBASE,
     4                     DINIT,DMIN,DSC,DELTA,IRM,RHSRLB,RHSRUB,
     5                     RHSREL,RHSROR,RHSRLL,RHSREU,RHSRLU,CSTROR,
     6                     CSTRLB,CSTRUB,CSTREL,CSTRLL,CSTREU,CSTRLU,
     7                     RANGEFLG,RHSIN,RHSINF,RANGENAME,RANGENAMEF,
     8                     CONTYP,IPGNA,NPGNA,DEWATER,DEWATERQ,VARBASE
      INTEGER(I4B),INTENT(IN)::NFVAR,NEVAR,NBVAR,NDV,NV,NCON,IOUT
      CHARACTER(LEN=200),INTENT(IN)::FNAME
C-----LOCAL VARIABLES
      REAL(SP)::TPBASE
      CHARACTER(LEN=10)::TFVNAME
      INTEGER(I4B)::LLOC,INMS,INMF,IKEYS,IKEYF,ISTART,ISTOP,BYTES
      INTEGER(I4B)::NUNOPN,LOCAT,I,J
      CHARACTER(LEN=200)::LINE,RMNAME,RMNAMEF,MPSNAME
      LOGICAL(LGT)::NFOUND
      INCLUDE 'openspec.inc' ! USE FORM TO DEFINE OPENING OF NON-FORMATTED FILES
C-----------------------------------------------------------------------
C
C1----OPEN FILE
      NUNOPN=IGETUNIT(10,200)
      LOCAT=NUNOPN
      WRITE(IOUT,1000,ERR=990)LOCAT,FNAME
      FUNIT=LOCAT; FILSTAT='OLD' 
      OPEN(UNIT=FUNIT,FILE=FNAME,STATUS=FILSTAT,ERR=999)
C
C-----CHECK FOR COMMENT LINES
      CALL URDCOM(LOCAT,IOUT,LINE)
C
C-----READ SOLUTION TYPE
      LLOC=1
      CALL URWORD(LINE,LLOC,IKEYS,IKEYF,1,NDUM,RDUM,IOUT,LOCAT)
      IF(LINE(IKEYS:IKEYF).NE.'NS'  .AND.
     1   LINE(IKEYS:IKEYF).NE.'MPS' .AND.
     2   LINE(IKEYS:IKEYF).NE.'LP'  .AND.
     3   LINE(IKEYS:IKEYF).NE.'SLP' .AND.
     4   LINE(IKEYS:IKEYF).NE.'FR') THEN
        WRITE(IOUT,2000,ERR=990)
        CALL GSTOP(' ')
      ENDIF
      SOLNTYP=LINE(IKEYS:IKEYF)
C
C-----IF SOLNTYP IS NS (NO SOLUTION) READ THE NAME
C      OF THE FILE TO WHICH THE RESPONSE MATRIX WILL BE WRITTEN
      IF(SOLNTYP.EQ.'NS')THEN
        WRITE(IOUT,3000,ERR=990)
        READ(LOCAT,*,ERR=991)DELTA
        READ(LOCAT,'(A)',ERR=991)LINE
        CALL CUTCOM(LINE,200)                   ! REMOVE COMMENTS FROM LINE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSIGDIG,RDUM,IOUT,LOCAT)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPGNMX, RDUM,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,PGFACT ,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,CRITMFC,IOUT,LOCAT)
        READ(LOCAT,'(A)',ERR=991)LINE
        CALL CUTCOM(LINE,200)                    ! REMOVE COMMENTS FROM LINE
        LLOC=1                                   ! READ FILENAME AND IRM
        CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
        RMNAME=LINE(INMS:INMF)
        RMFILE=IGETUNIT(980,1000)
        IRM=0                                    ! SET DEFAULT VALUE
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IRM,RDUM,IOUT,LOCAT)
        IF(IRM.EQ.1 .OR. IRM.EQ.0)THEN           ! SET TO UNFORMATTED OUTPUT
          IRM=1                                  ! IF ZERO, RESET TO ONE
          WRITE(IOUT,3010,ERR=990)RMFILE,RMNAME
          FMTARG='UNFORMATTED' 
          RMFILE=IGETUNIT(980,1000)
          FUNIT=RMFILE 
        ELSE                                     ! SET TO FORMATTED OUTPUT
          IRM=3                                  ! RESET TO INTERNAL FLAG
          WRITE(IOUT,3011,ERR=990)RMFILE,RMNAME
          FMTARG='FORMATTED' 
          RMFILEF=IGETUNIT(980,1000)
          FUNIT=RMFILE 
        ENDIF
        OPEN(FUNIT,FILE=RMNAME,FORM=FMTARG,ERR=999)
C
C-------PROCESS DELTA, THE PERTURBATION SIZE
        IF(ABS(DELTA).GT.50. .OR. ABS(DELTA).LT.5.0D-5)
     &                                          WRITE(IOUT,3020,ERR=990)
        WRITE(IOUT,3030,ERR=990)DELTA
      ENDIF
C
C-----IF SOLNTYP IS MPS (MPS FILE WILL BE WRITTEN) 
C      READ THE NAME OF THE MPS FILE TO WHICH
C      THE MANAGEMENT PROBLEM WILL BE WRITTEN IN MPS FORMAT 
      IF(SOLNTYP.EQ.'MPS')THEN
        WRITE(IOUT,3040,ERR=990)
        NONLIN=.FALSE.                           ! SOLVE AS A LINEAR PROBLEM
        IRM=2
        READ(LOCAT,*,ERR=991)DELTA
        READ(LOCAT,'(A)',ERR=991)LINE
        CALL CUTCOM(LINE,200)                   ! REMOVE COMMENTS FROM LINE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSIGDIG,RDUM,IOUT,LOCAT)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPGNMX, RDUM,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,PGFACT ,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,CRITMFC,IOUT,LOCAT)
        READ(LOCAT,'(A)',ERR=991)LINE
        LLOC=1
        CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
        MPSNAME=LINE(INMS:INMF)
        MPSFILE=IGETUNIT(980,1000)
        WRITE(IOUT,3010,ERR=990)MPSFILE,MPSNAME
        FUNIT=MPSFILE 
        OPEN(FUNIT,FILE=MPSNAME,ERR=999)
C
C-------PROCESS DELTA, THE PERTURBATION SIZE
        IF(ABS(DELTA).GT.50. .OR. ABS(DELTA).LT.5.0D-5)
     &                                          WRITE(IOUT,3020,ERR=990)
        WRITE(IOUT,3030)DELTA
      ENDIF
C
C-----IF SOLNTYP IS LP (LINEAR), READ IRM, LPITMAX, AND PERTI:
      IF(SOLNTYP.EQ.'LP')THEN
        WRITE(IOUT,3050,ERR=990)
        NONLIN=.FALSE.                           ! SOLVE AS A LINEAR PROBLEM
        READ(LOCAT,*,ERR=991)IRM
        READ(LOCAT,*,ERR=991)LPITMAX,BBITMAX
        READ(LOCAT,*,ERR=991)DELTA
        READ(LOCAT,'(A)',ERR=991)LINE
        CALL CUTCOM(LINE,200)                   ! REMOVE COMMENTS FROM LINE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSIGDIG,RDUM,IOUT,LOCAT)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPGNMX, RDUM,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,PGFACT ,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,CRITMFC,IOUT,LOCAT)
        READ(LOCAT,*,ERR=991)BBITPRT,RANGEFLG
C
C-------PROCESS IRM, THE RESPONSE MATRIX CALCULATION FLAG
        IF(IRM.EQ.0)THEN               ! READ EXISTING RM FILE
          WRITE(IOUT,3060,ERR=990)
          READ(LOCAT,'(A)',ERR=991)LINE  
		  LLOC=1        
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAME=LINE(INMS:INMF)
          RMFILE=IGETUNIT(980,1000)
          OPEN(RMFILE,FILE=RMNAME,STATUS='OLD',FORM=FORM,ERR=997)
          WRITE(IOUT,3012,ERR=990) RMFILE,RMNAME
        ELSEIF(IRM.EQ.1)THEN           ! COMPUTE RM AND SAVE TO FILE  
          WRITE(IOUT,3070,ERR=990)
          READ(LOCAT,'(A)',ERR=991)LINE   
		  LLOC=1        
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAME=LINE(INMS:INMF)
          RMFILE=IGETUNIT(980,1000)
          OPEN(RMFILE,FILE=RMNAME,STATUS='REPLACE',FORM=FORM,ERR=997)
          WRITE(IOUT,3010,ERR=990) RMFILE,RMNAME
        ELSEIF(IRM.EQ.2) THEN          ! COMPUTE RM FILE - DON'T RECORD IT    
          WRITE(IOUT,3080,ERR=990)
        ELSEIF(IRM.EQ.3)THEN           ! COMPUTE RM AND PRINT TO FILE   
          WRITE(IOUT,3082,ERR=990)
          READ(LOCAT,'(A)',ERR=991)LINE   
		  LLOC=1        
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAMEF=LINE(INMS:INMF)
          RMFILEF=IGETUNIT(980,1000)
          OPEN(RMFILEF,FILE=RMNAMEF,STATUS='REPLACE',
     &         FORM='FORMATTED',ERR=997)
          WRITE(IOUT,3014,ERR=990) RMFILEF,RMNAMEF
        ELSEIF(IRM.EQ.4)THEN           ! COMPUTE RM, SAVE AND PRINT IT
          WRITE(IOUT,3084,ERR=990)
          READ(LOCAT,'(A)',ERR=991)LINE    
		  LLOC=1        
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAME=LINE(INMS:INMF)
          RMFILE=IGETUNIT(980,1000)
          OPEN(RMFILE,FILE=RMNAME,STATUS='REPLACE',FORM=FORM,ERR=997)
          WRITE(IOUT,3010,ERR=990) RMFILE,RMNAME
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAMEF=LINE(INMS:INMF)
          RMFILEF=IGETUNIT(980,1000)
          OPEN(RMFILEF,FILE=RMNAMEF,STATUS='REPLACE',
     &         FORM='FORMATTED',ERR=997)
          WRITE(IOUT,3014,ERR=990) RMFILEF,RMNAMEF
        ELSEIF(IRM.EQ.5)THEN           ! READ EXISTING RM FILE AND PRINT IT
          WRITE(IOUT,3086,ERR=990)
          READ(LOCAT,'(A)',ERR=991)LINE      
		  LLOC=1        
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAME=LINE(INMS:INMF)
          RMFILE=IGETUNIT(980,1000)
          OPEN(RMFILE,FILE=RMNAME,STATUS='OLD',FORM=FORM,ERR=997)
          WRITE(IOUT,3012,ERR=990) RMFILE,RMNAME
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          RMNAMEF=LINE(INMS:INMF)
          RMFILEF=IGETUNIT(980,1000)
          OPEN(RMFILEF,FILE=RMNAMEF,STATUS='REPLACE',
     &         FORM='FORMATTED',ERR=997)
          WRITE(IOUT,3014,ERR=990) RMFILEF,RMNAMEF
        ELSE
          WRITE(IOUT,3090,ERR=990)IRM            ! INVALID VALUE OF IRM
          CALL GSTOP(' ')
        ENDIF
C
C-------PROCESS LPITMAX AND BBITMAX, MAXIMUM LP AND BRANCH AND BOUND ITERATIONS
        WRITE(IOUT,3100,ERR=990)LPITMAX,BBITMAX
        IF(NBVAR.GT.0.AND.BBITMAX.LT.2)CALL GSTOP('BBITMAX MUST BE >1')
C
C-------PROCESS DELTA, THE PERTURBATION SIZE
        IF(ABS(DELTA).GT.50. .OR. ABS(DELTA).LT.5.0D-5)
     &                                          WRITE(IOUT,3020,ERR=990)
        WRITE(IOUT,3030,ERR=990)DELTA
C
C-------PROCESS NPGNMX,AFACT,PGFACT, THE PERTURBATION ATTEMPT PARAMETERS
        WRITE(IOUT,3110,ERR=990)NPGNMX,PGFACT
C
C-------PROCESS BBITPRT
        IF(NBVAR.EQ.0)BBITPRT=0                  ! SET TO ZERO IF NO BINARIES
        IF(BBITPRT.EQ.0)THEN
          WRITE(IOUT,3120,ERR=990)
        ELSEIF(BBITPRT.EQ.1)THEN
          WRITE(IOUT,3130,ERR=990)
        ELSE
          WRITE(IOUT,3140,ERR=990)BBITPRT        ! INVALID ENTRY
          CALL GSTOP(' ')
        ENDIF
      ENDIF
C
C----IF SOLNTYP IS SLP 
      IF(SOLNTYP.EQ.'SLP')THEN
        IRM = 2                                  ! RESPONSE MATRIX MUST BE CALCULATED
        WRITE(IOUT,3150,ERR=990)
        NONLIN=.TRUE.                            ! SOLVE AS NONLINEAR PROBLEM
        READ(LOCAT,*,ERR=991)SLPITMAX,LPITMAX,BBITMAX
        READ(LOCAT,*,ERR=991)SLPVCRIT,SLPZCRIT,DINIT,DMIN,DSC
        READ(LOCAT,'(A)',ERR=991)LINE
        CALL CUTCOM(LINE,200)                   ! REMOVE COMMENTS FROM LINE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSIGDIG,RDUM,IOUT,LOCAT)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPGNMX, RDUM,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,PGFACT ,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,AFACT  ,IOUT,LOCAT)
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NINFMX, RDUM,IOUT,LOCAT)
        CALL URWORD2(LINE,LLOC,ISTART,ISTOP,3,NDUM,CRITMFC,IOUT,LOCAT)
        READ(LOCAT,*,ERR=991)SLPITPRT,BBITPRT,RANGEFLG
C
C-------PROCESS LPITMAX AND BBITMAX, MAXIMUM LP AND BRANCH AND BOUND ITERATIONS
        WRITE(IOUT,3100,ERR=990)LPITMAX,BBITMAX
C
C-------PROCESS SLPITMAX MAXIMUM SLP ITERATIONS
        WRITE(IOUT,3160,ERR=990)SLPITMAX
C
C-------PROCESS SLPVCRIT,SLPZCRIT
        WRITE(IOUT,3170,ERR=990)SLPVCRIT,SLPZCRIT
C
C-------PROCESS THE PERTURBATION PARAMETERS
        IF(ABS(DINIT).GT.50. .OR. ABS(DINIT).LT.5.0D-5)
     &                                          WRITE(IOUT,3020,ERR=990)
        IF(DSC.LE.ZERO)THEN
          WRITE(IOUT,3180,ERR=990)
          CALL GSTOP(' ')
        ENDIF
        IF(DINIT*DMIN.LT.ZERO)THEN            ! DINIT AND DMIN ARE NOT SAME SIGN
          WRITE(IOUT,3185,ERR=990)
          CALL GSTOP(' ')
        ENDIF
        WRITE(IOUT,3190,ERR=990)DINIT,DMIN,DSC
        DELTA = DINIT                              ! INITIALIZE WELL PERTURBATION
C
C-------PROCESS NPGNMX,AFACT,PGFACT, THE PERTURBATION ATTEMPT PARAMETERS
        WRITE(IOUT,3200,ERR=990)NPGNMX,AFACT,PGFACT
C
C-------PROCESS SLPITPRT
        IF(SLPITPRT.EQ.0)WRITE(IOUT,3210,ERR=990)
        IF(SLPITPRT.GE.1)WRITE(IOUT,3220,ERR=990)
        IF(SLPITPRT.EQ.2)WRITE(IOUT,3225,ERR=990)
C
C-------PROCESS BBITPRT
        IF(NBVAR.EQ.0)BBITPRT=0                  ! SET TO ZERO IF NO BINARIES
        IF(BBITPRT.EQ.0)THEN
          WRITE(IOUT,3120,ERR=990)
        ELSEIF(BBITPRT.EQ.1)THEN
          WRITE(IOUT,3130,ERR=990)
        ELSE
          WRITE(IOUT,3140,ERR=990)BBITPRT                  ! INVALID ENTRY
          CALL GSTOP(' ')
        ENDIF
C
      ENDIF
C
C-----IF SOLNTYP IS FR (FORWARD RUN) CALCULATE A FORWARD RUN, 
C       EVALUATE THE CONSTRAINTS AND OBJECTIVE FUNCTION AND THEN STOP.
      IF(SOLNTYP.EQ.'FR')THEN
        WRITE(IOUT,3300,ERR=990)
        CRITMFC = 0.0
        IRM = 2
      ENDIF
C
C-----ECHO INPUT PARAMETERS THAT ARE COMMON TO ALL SOLNTYP
C
C-------PROCESS CRITMFC, THE FLOW PROCESS ACCEPTABILITY CRITERIA
      IF(CRITMFC.LT.0.0)WRITE(IOUT,3230,ERR=990)CRITMFC
      IF(CRITMFC.EQ.0.0)WRITE(IOUT,3231,ERR=990)CRITMFC
      IF(CRITMFC.GT.0.0)WRITE(IOUT,3232,ERR=990)CRITMFC
C
C-----ASSIGN BASE PUMPING RATES
      READ(LOCAT,*,ERR=991)IBASE
      IF(IBASE.EQ.0)THEN                         ! BASE PUMPING ALREADY INPUT
        WRITE(IOUT,4000,ERR=990)                           ! WRITE HEADER 
      ELSEIF(IBASE.EQ.1)THEN
        DO 100 I=1,NFVAR
          READ(LOCAT,'(A)',ERR=991,END=991)LINE
          LLOC=1
          CALL URWORD(LINE,LLOC,INMS,INMF,0,NDUM,RDUM,IOUT,LOCAT)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,NDUM,TPBASE,IOUT,LOCAT)
C
C---------PROCESS WELL NAME
          TFVNAME=LINE(INMS:INMF)
          NFOUND=.TRUE.
          DO 110 J=1,NFVAR
            IF(FVNAME(J).EQ.TFVNAME)THEN         ! NAME OF VARIABLE IS VALID
              NFOUND=.FALSE.
              EXIT
            ENDIF
  110     ENDDO
          IF(NFOUND)THEN
            WRITE(IOUT,4010,ERR=990)TFVNAME      ! INVALID WELL NAME
            CALL GSTOP(' ')
          ELSE
            FVINI(J)=REAL(TPBASE,DP)             ! STORE PUMPING RATE
            IF(FVDIR(J).EQ.2)FVINI(J)=-FVINI(J)  ! SWITCH SIGN FOR WITHDRAWALS
          ENDIF
  100   ENDDO
        WRITE(IOUT,4020,ERR=990)                 ! WRITE HEADER 
        DO 120 I=1,NFVAR
          WRITE(IOUT,4030,ERR=990)I,FVNAME(I),ABS(FVINI(I)) ! WRITE BASE PUMPING RATES
  120   ENDDO
      ELSE
        WRITE(IOUT,4040,ERR=990)IBASE            ! INVALID VALUE FOR IBASE
        CALL GSTOP(' ')
      ENDIF
C
C-----INITIALIZE EXTERNAL AND BINARY BASE VALUES FOR CONVERGENCE TESTING
      EVBASE=ZERO
      BVBASE=ZERO

C-----ALLOCATE WORK ARRAYS USED IN PERTURBATION CALCULATIONS
      ALLOCATE (DELINC(NFVAR),IPGNA(1:NFVAR+1),NPGNA(NFVAR),
     &          MFCNVRG(-1:NFVAR+1),DEWATER(-1:NFVAR+1),
     &          DEWATERQ(-1:NFVAR+1),VARBASE(1:NFVAR+1),
C
C-----ALLOCATE SPACE FOR THE OLD PUMPING RATE AND SAVE IT
     &          FVOLD(NFVAR),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
      BYTES = 8*2*NFVAR + (2*4 + 8 + 3)*(NFVAR+2)
      DO 130 I=1,NFVAR
        IPGNA(I) = 0
        NPGNA(I) = 0
  130 ENDDO
      IF(IBASE.EQ.0)THEN                         ! BASE PUMPING IN FVBASE
        DO 300 I=1,NFVAR
          FVOLD(I) = FVBASE(I)
  300   ENDDO
      ELSEIF(IBASE.EQ.1)THEN                     ! BASE PUMPING IN FVINI
        DO 310 I=1,NFVAR
          FVOLD(I) = FVINI(I)
  310   ENDDO
      ENDIF
C
C-----ALLOCATE MEMORY FOR ARRAYS TO STORE THE LINEAR PROGRAM
      ALLOCATE (AMAT(NCON,NDV),CST(NV),BNDS(NV),RHS(NCON),
     1          RHSRLB(NCON),RHSRUB(NCON),RHSROR(NCON),
     2          RHSREL(NCON),RHSRLL(NCON),RHSREU(NCON),RHSRLU(NCON),
     3          CSTROR(NV),CSTRLB(NV),CSTRUB(NV),
     4          CSTREL(NV),CSTRLL(NV),CSTREU(NV),CSTRLU(NV),
     5          RHSIN(NCON),RHSINF(NCON),CONTYP(NCON),
     6          RANGENAME(NCON+NDV),RANGENAMEF(NCON+NDV),
     7          STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
      BYTES = BYTES + 8*(NCON*NDV + 2*NV + NCON) +
     1         8*3*NCON + 4*4*NCON + 8*3*NV + 4*4*NV +
     2         8*2*NCON + 10*2*(NCON+NDV)
C-----CLEAR MEMORY
      RHS  = ZERO
      AMAT = ZERO
      CST  = ZERO
      BNDS = ZERO
C
C-----ADD THE INTEGER AND DOUBLE PRECISION SPACE NEEDED BY SIMPLEX
      BYTES = BYTES + 4*(3*NCON + NV+NCON)
      BYTES = BYTES 
     1            + 8*(10*NCON + 2*NCON*NCON + NV+NCON + NCON*(NDV+1))
        WRITE(IOUT,*,ERR=990)
        WRITE(IOUT,*,ERR=990)'    PROBLEM SIZE'
        WRITE(IOUT,*,ERR=990)
        WRITE(IOUT,*,ERR=990)' NUMBER OF VARIABLES (INCLUDING SLACKS)',
     &                       NV
        WRITE(IOUT,*,ERR=990)' NUMBER OF CONSTRAINT EQUATIONS        ',
     &                       NCON
        WRITE(IOUT,*,ERR=990)
C
C-----CLOSE FILE
      CLOSE(UNIT=LOCAT)
      WRITE(IOUT,6000,ERR=990)BYTES
      WRITE(IOUT,7000,ERR=990)
C
 1000 FORMAT(1X,/1X,'OPENING SOLUTION FILE ON UNIT ',I4,':',
     1  /1X,A200)
 2000 FORMAT(1X,/1X,'PROGRAM STOPPED. USER MUST SPECIFY',
     1  ' SOLNTYP EITHER AS NS, MPS, LP, OR SLP.')
 3000 FORMAT(1X,/1X,'NO SOLUTION TO THE FORMULATION WILL BE FOUND.',
     1  ' GWM WILL',/,' CALCULATE THE RESPONSE MATRIX AND WRITE IT',
     2  ' TO THE RMNAME',/,' FILE AND THEN WILL STOP.') 
 3010 FORMAT(1X,/1X,'OPENING UNFORMATTED RESPONSE-MATRIX FILE FOR',
     1  ' SAVING ON UNIT ',I4,':',/1X,A200)
 3011 FORMAT(1X,/1X,'OPENING FORMATTED RESPONSE-MATRIX FILE FOR',
     1  ' SAVING ON UNIT ',I4,':',/1X,A200)
 3012 FORMAT(1X,/1X,'OPENING RESPONSE-MATRIX FILE FOR READING ON UNIT '
     1  ,I4,':',/1X,A200)
 3014 FORMAT(1X,/1X,'OPENING RESPONSE-MATRIX FILE FOR PRINTING ON UNIT '
     1  ,I4,':',/1X,A200)
 3020 FORMAT (1X,/1X,' WARNING:',/,
     1     'PERTURBATION PARAMETERS MAY BE TOO LARGE OR TOO SMALL',/,
     2     'COMPARED TO MAXIMUM PUMPING RATE.')
 3030 FORMAT(1X,/1X,'PERTURBATION VALUE: ',T45,D10.2)
 3040 FORMAT(1X,/1X,'NO SOLUTION TO THE FORMULATION WILL BE FOUND.',
     1  ' GWM WILL',/,' WRITE THE FORMULATION IN',
     2  ' MPS FORMAT TO FILE MPSNAME AND THEN STOP.')    
 3050 FORMAT(1X,/1X,'SOLNTYP IS LP: GWM WILL COMPLETE A SINGLE',
     1  ' ITERATION OF THE LINEAR PROBLEM.')
 3060 FORMAT(1X,/1X,'IRM EQUALS 0: RESPONSE MATRIX WILL BE READ',
     1  ' FROM FILE')
 3070 FORMAT(1X,/1X,'IRM EQUALS 1: RESPONSE MATRIX WILL BE',
     1  ' CALCULATED BY GWM',/,' AND SAVED TO FILE') 
 3080 FORMAT(1X,/1X,'IRM EQUALS 2: RESPONSE MATRIX WILL BE',
     1  ' CALCULATED BY GWM',/,' BUT NOT WRITTEN TO FILE')
 3082 FORMAT(1X,/1X,'IRM EQUALS 3: RESPONSE MATRIX WILL BE',
     1  ' CALCULATED BY GWM',/,' AND PRINTED TO FILE ')
 3084 FORMAT(1X,/1X,'IRM EQUALS 4: RESPONSE MATRIX WILL BE',
     1  ' CALCULATED BY GWM',/,' AND PRINTED AND SAVED')
 3086 FORMAT(1X,/1X,'IRM EQUALS 5: RESPONSE MATRIX WILL BE',
     1  ' READ FROM FILE',/,' AND PRINTED TO ANOTHER FILE')
 3090 FORMAT(1X,/1X,'PROGRAM STOPPED. IRM MUST BE BETWEEN 0 AND 5,',
     1  ' BUT WAS SET AT: ',I5)
 3100 FORMAT(1X,/1X,'MAXIMUM NUMBER OF LP ITERATIONS: ',T45,I8,/,
     &   1X,'MAXIMUM NUMBER OF BRANCH AND BOUND ITER: ',T45,I8)
 3110 FORMAT(/1X,'MAXIMUM NUMBER OF PERTURBATION ATTEMPTS: ',T45,I8,/,
     &   1X,'PERTURBATION ADJUSTMENT FACTOR (PGFACT): ',T45,F8.5)
 3120 FORMAT(1X,/1X,'OUTPUT FROM BRANCH-AND-BOUND ALGORITHM WILL',
     1  ' NOT BE PRINTED.')
 3130 FORMAT(1X,/1X,'OUTPUT FROM BRANCH AND BOUND ALGORITHM WILL',
     1  ' BE PRINTED.')
 3140 FORMAT(1X,/1X,'PROGRAM STOPPED. BBITPRT MUST BE EQUAL TO 0 OR 1: ',
     1  I4)
 3150 FORMAT(1X,/1X,'SOLNTYP IS SLP: GWM WILL USE SEQUENTIAL',
     1  ' ITERATION OF THE',/,' NONLINEAR PROBLEM UNTIL A',
     2  ' SOLUTION IS FOUND OR THE PROBLEM',/,' DOES NOT CONVERGE.')
 3160 FORMAT(1X,/1X,'MAXIMUM NUMBER OF SLP ITERATIONS: ',T45,I8)
 3170 FORMAT(1X,/1X,'SLP VARIABLE CONVERGENCE CRITERION',
     1  ' (SLPVCRIT):  ',ES13.5,
     2   /1X,'SLP OBJECTIVE CONVERGENCE CRITERION',
     3  ' (SLPZCRIT): ',ES13.5)
 3180 FORMAT(1X,/1X,'PROGRAM STOPPED. PERTURBATION SCALING',
     1  ' PARAMETER (PERTS)',/,' MUST BE GREATER THAN ZERO.')
 3185 FORMAT(1X,/1X,'PROGRAM STOPPED. PERTURBATION PARAMETERS ',
     1    'DINIT AND DMIN MUST HAVE SAME SIGN')
 3190 FORMAT(1X,/1X,'PERTURBATION VALUES (DINIT, DMIN, AND',
     1  ' DSC): ',3D10.2)
 3200 FORMAT(/1X,'MAXIMUM NUMBER OF PERTURBATION ATTEMPTS: ',T45,I8,/,
     &   1X,'BASE FLOW RATE RELAXATION PARAMETER (AFACT): ',T45,F8.5,/,
     &   1X,'PERTURBATION ADJUSTMENT FACTOR (PGFACT): ',T45,F8.5)
 3210 FORMAT(1X,/1X,'SLPITPRT=0: DO NOT PRINT SLP ITERATIONS.')
 3220 FORMAT(1X,/1X,'SLPITPRT>=1: PRINT SLP ITERATIONS.')
 3225 FORMAT(1X,/1X,'SLPITPRT=2: PRINT INTERMEDIATE LP SOLUTIONS.')
 3230 FORMAT(/1X,'CRITMFC SET TO ',T45,D11.3,/,
     &   1X,'GWM WILL ALWAYS ACCEPT FLOW PROCESS RESULTS')
 3231 FORMAT(/1X,'CRITMFC SET TO ',T45,D11.3,/,
     &   1X,'GWM WILL ACCEPT FLOW PROCESS RESULTS THAT MEET GWF',/,
     &      ' CONVERGENCE CRITERIA')
 3232 FORMAT(/1X,'CRITMFC SET TO ',T45,D11.3,/,
     &   1X,'GWM WILL ACCEPT FLOW PROCESS RESULTS IF PERCENT ',/,
     &      'DISCREPANCY OF WATER BALANCE FLOW RATE IS BELOW CRITMFC')
 3300 FORMAT(1X,/1X,'NO RESPONSE MATRIX OR SOLUTION WILL BE COMPUTED.',
     1  ' GWM WILL',/,' CALCULATE A FORWARD RUN, EVALUATE THE CONS',
     2  'TRAINTS AND',/,' OBJECTIVE FUNCTION AND THEN STOP.') 
 4000 FORMAT(/,T2,'BASE PUMPING RATES TAKEN FROM FVREF SPECIFIED IN',
     1  ' VARCON INPUT FILE')
 4010 FORMAT(1X,/1X,'PROGRAM STOPPED. ',A10,' WAS NOT DEFINED AS A', 
     1  ' VARIABLE NAME (FVNAME)',/,' IN THE DECISION-VARIABLE FILE.')
 4020 FORMAT(/,T2,'BASE PUMPING RATES:',/,T2,'VARIABLE',T11,
     1  'VARIABLE',T26,'BASE',/,T3,'NUMBER',T13,'NAME',T23,
     2  'PUMPING RATE',/,' ---------------------------------')
 4030 FORMAT(I5,6X,A10,T23,D11.3)
 4040 FORMAT(1X,/1X,'PROGRAM STOPPED. IBASE MUST EQUAL 0 OR 1,',
     1  ' BUT WAS SET AT: ',I5)
 6000 FORMAT(/,1X,I8,' BYTES OF MEMORY ALLOCATED FOR',
     1               ' RESPONSE MATRIX ALGORITHM')
 7000 FORMAT(/,1X,'CLOSING SOLUTION AND OUTPUT FILE',/)
C
      RETURN
C
C-----ERROR HANDLING. 
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(IOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),IOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1RMS3AR)')
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
     &2X,'-- STOP EXECUTION (GWM1RMS3AR)')
      CALL GSTOP(' ')
C
  992 CONTINUE
C-----ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
      WRITE(IOUT,9920)
 9920 FORMAT(/,1X,'*** ERROR ALLOCATING ARRAY(S)',
     &2X,'-- STOP EXECUTION (GWM1RMS3AR)')
      CALL GSTOP(' ')
C
  997 CONTINUE
C-----FILE-OPENING ERRORS
      WRITE(IOUT,9970)
 9970 FORMAT(1X,/1X,'PROGRAM STOPPED. PROBLEM OPENING RMNAME FILE.',
     1  '  BE SURE THAT THE',/,' RMNAME FILE EXISTS.')
      CALL GSTOP(' ')
C
!  998 WRITE(IOUT,9980)RMNAME
!      CALL GSTOP(' ')
! 9980 FORMAT(1X,/1X,'PROGRAM STOPPED. RESPONSE MATRIX FILE NAMED ',
!     1  /,A,/,
!     2  'ALREADY EXISTS.  GWM WILL NOT OVERWRITE RMNAME FILE')
C
  999 CONTINUE
C-----FILE-OPENING ERROR
      INQUIRE(FUNIT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9990)TRIM(FLNM),LOCAT,FILSTAT,FMTARG,ACCARG,FILACT
      WRITE(IOUT,9990)TRIM(FLNM),LOCAT,FILSTAT,FMTARG,ACCARG,
     &                 FILACT
 9990 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE STATUS: ',A,/
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1RMS3AR)')
      CALL GSTOP(' ')
C
      END SUBROUTINE GWM1RMS3AR
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3PL(IPERT,NPERT,FIRSTSIM,LASTSIM)
C***********************************************************************
C  VERSION: 21JAN2010
C  PURPOSE - SET THE TERMINAL VALUE OF THE PERTURBATION LOOP
C            AND READ THE RESPONSE MATRIX FROM FILE IF AVAILABLE 
C-----------------------------------------------------------------------
      USE GWM1RMS3, ONLY : SLPITCNT,SLPITMAX,SLPINFCNT,IBASE,IREF,
     &                     OBJOLD,MFCNVRG,DELTA,DELINC,NRMC,NCON,NDV,
     &                     IRM,NRESET
      USE GWM1DCV3, ONLY : NFVAR,NBVAR,NEVAR,FVBASE,FVMAX,FVDIR,FVINI
      USE GWM1HDC3, ONLY : HDCRHS,HDCSTATE0,HDCNAME,HDCDIR
      USE GWM1OBJ3, ONLY : SOLNTYP
      USE GWM1BAS3, ONLY : ZERO,BIGINF,RMFILE,GWMOUT
      INTEGER(I4B),INTENT(OUT)::IPERT,NPERT
      LOGICAL(LGT),INTENT(IN)::FIRSTSIM,LASTSIM
      INTEGER(I4B)::I,TNFVAR,TNEVAR,TNBVAR,TIBASE,RSTRT
! for debugging
!      integer :: isize      
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(FIRSTSIM .AND. .NOT. LASTSIM)THEN
C-------THIS IS THE FIRST TIME THROUGH THIS SUBROUTINE
        MFCNVRG(0) = .TRUE.                      ! INITIALIZE MF CONVERGENCE FOR VI
        IF(SOLNTYP.EQ.'SLP')THEN                 ! PERFORM SLP INITIALIZATIONS
           SLPITCNT=0                            ! SLP ITERATION COUNTER
           SLPINFCNT=0                           ! SLP FAILURE COUNTER
           NRESET=0                              ! BASE RESET COUNTER
           OBJOLD  = SQRT(BIGINF)                ! CONVERGENCE COMPARISON
        ENDIF
        IF(IRM.GE.1 .AND. IRM.LE.4)THEN          ! RESPONSE MATRIX NEEDED
          DO 100 I = 1,NFVAR                     ! CREATE INITIAL PERTURBATIONS
            DELINC(I) = DELTA*FVMAX(I)           ! SCALE EACH BY UPPER BOUND
            IF(FVDIR(I).EQ.2)DELINC(I)=-DELINC(I)! SWITCH FOR WITHDRAWAL VARIABLE
  100     ENDDO
          IF(IREF.EQ.1 .AND. IBASE.EQ.1)THEN     ! EXTRA SIMULATION NEEDED
            IPERT=-1                             ! SET INITIAL LOOP VALUE 
          ELSE                                   ! EITHER REFERENCE NOT NEEDED
            IPERT=0                              !   OR REFERENCE FLOWS
          ENDIF                                  !   SAME AS BASE FLOWS
          NPERT = NFVAR                          ! SET TERMINAL VALUE TO 
C                                                !   NUMBER OF VARIABLES
        ELSEIF(IRM.EQ.0 .OR. IRM.EQ.5)THEN       ! RESPONSE MATRIX READ FROM FILE
          CALL SGWM1RMS3PL
          IPERT = 0                              ! SET IPERT AND NPERT  
          NPERT = -1                             ! SO LOOP IS NOT EXECUTED
        ENDIF
      ELSEIF(.NOT.FIRSTSIM .AND. .NOT.LASTSIM)THEN
C-------THIS IS AN ITERATION OF THE SLP ALGORITHM
        IF(SLPINFCNT.EQ.0)THEN                   ! PRIOR ITERATION SUCCEEDED  
          CALL GWM1RMS3SLP_SETBASE               ! SET NEW BASE SOLUTION AND
          CALL GWM1RMS3SLP_SETDEL                ! SET NEW PERTUBATION VALUE
        ELSEIF(SLPINFCNT.GT.0)THEN               ! PRIOR ITERATION FAILED 
          CALL GWM1RMS3SLP_SETDEL                ! SET ONLY PERTUBATION VALUE
        ENDIF
        IPERT = 0                                ! SET INITIAL LOOP VALUE 
        NPERT = NFVAR                            ! SET TERMINAL VALUE 
C
      ELSEIF(LASTSIM)THEN 
C-------THIS IS LAST TIME THROUGH. ONLY ONE SIMULATION IS NEEDED
        IPERT = 0                                ! SET INITIAL LOOP VALUE 
        NPERT = 0                                ! SET TERMINAL VALUE TO ZERO
      ENDIF
C
! for debugging
!      write(gwmout,*)'at the bottom of GWM1RMS3PL, fvbase contains:'
!      isize = size(fvbase)
!      write(gwmout,987)(i,fvbase(i),i=1,isize)
!  987 format(i3,1x,g15.8)
      RETURN
C
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1RMS3PL
C***********************************************************************
C
C  PURPOSE: READ THE REFERENCE AND BASE STATES AND RESPONSE MATRIX FROM FILE 
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWM1BAS3PS,GWM1BAS3PF
      USE GWM1HDC3, ONLY : GWM1HDC3FPR 
      USE GWM1STC3, ONLY : GWM1STC3FPR
      USE GWM1STA3, ONLY : GWM1STA3FPR
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      READ(RMFILE,ERR=991)TNFVAR,TNEVAR,TNBVAR,TIBASE ! READ HEADER FOR TEST
      IF(TNFVAR.NE.NFVAR.OR.TNEVAR.NE.NEVAR.OR.
     &   TNBVAR.NE.NBVAR.OR.TIBASE.NE.IBASE)THEN 
        WRITE(GWMOUT,1000,ERR=990)TNFVAR,TNEVAR,TNBVAR,IBASE ! NOT A VALID MATRIX
        CALL GSTOP(' ')
      ENDIF  
      IF(IBASE.EQ.1)THEN                         ! BASE VALUES ARE AVAILABLE
        DO 100 I=1,NFVAR                         
          FVBASE(I)=FVINI(I)                     ! LOAD BASE FLOWS
  100   ENDDO                                    
      ENDIF
C
      CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
      CALL GWM1BAS3PF('               Solution Algorithm',0,ZERO)
      CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
      CALL GWM1BAS3PS('  Reading Response Matrix',0)
C
C-----FIRST, READ HEAD CONSTRAINT REFERENCE AND BASE STATE INTO FULL RHS ARRAYS     
      RSTRT = 1                                  ! START FROM THE FIRST ROW
      CALL GWM1HDC3FPR(RSTRT, 0, 0)        
C-----2ND, READ STREAM CONSTRAINT REFERENCE AND BASE STATE INTO FULL RHS ARRAYS
      CALL GWM1STC3FPR(RSTRT, 0, 0)              ! READ BASE STATE
C-----3RD, READ BASE STATE FOR STATE VARIABLES
      CALL GWM1STA3FPR(0,0)
C-----READ THE RESPONSE MATRIX INTO AMAT ARRAY, ONE FULL COLUMN AT A TIME
      DO 200 I=1,NFVAR                           ! LOOP OVER COLUMNS                                
        RSTRT = 1                                ! START FROM THE FIRST ROW
        CALL GWM1HDC3FPR(RSTRT, 1, I)            ! READ HEAD CON IN COLUMN
        CALL GWM1STC3FPR(RSTRT, 1, I)            ! READ STRM CON IN COLUMN
        CALL GWM1STA3FPR(1, I)                   ! READ STATE RESPONSE
  200 ENDDO
C
C-----WRITE CONSTRAINT STATUS
      CALL GWM1RMS3OT(1,0)
C
 1000 FORMAT('PROGRAM STOPPED: SPECIFIED RESPONSE MATRIX ',/,
     &       'DOES NOT MATCH THE CURRENT VARIABLE NUMBERS',4I5)
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (SGWM1RMS3PL)')
      CALL GSTOP(' ')
C
  991 CONTINUE
C-----FILE-READING ERROR
      INQUIRE(RMFILE,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9910)TRIM(FLNM),RMFILE,FMTARG,ACCARG,FILACT
 9910 FORMAT(/,1X,'*** ERROR READING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (SGWM1RMS3PL)')
      CALL GSTOP(' ')
C
      END SUBROUTINE SGWM1RMS3PL
C
C***********************************************************************
      SUBROUTINE GWM1RMS3SLP_SETDEL
C***********************************************************************
      USE GWM1RMS3, ONLY : DINIT,DMIN,DSC
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----CALCULATE A NEW VALUE OF THE PERTURBATION PARAMETER
      DELTA = DMIN + (DINIT-DMIN)/DSC**(SLPITCNT)
C
      DO 100 I = 1,NFVAR                     ! RECALULATE PERTURBATIONS
        DELINC(I) = DELTA*FVMAX(I)           ! SCALE EACH BY UPPER BOUND
        IF(FVDIR(I).EQ.2)DELINC(I)=-DELINC(I)! SWITCH FOR WITHDRAWAL VARIABLE
100   ENDDO
C
      RETURN
      END SUBROUTINE GWM1RMS3SLP_SETDEL
C
C***********************************************************************
      SUBROUTINE GWM1RMS3SLP_SETBASE 
C***********************************************************************
      USE GWM1RMS3, ONLY : CST,FVOLD
      USE GWM1DCV3, ONLY : NFVAR,FVBASE,NEVAR,EVBASE,NBVAR,BVBASE,
     &                     FVCURRENT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----STORE CURRENT VALUES OF FLOW VARIABLES IN FVBASE
      DO 100 I=1,NFVAR
        FVOLD(I)  = FVBASE(I)                    ! SAVE THE OLD SOLUTION
        FVBASE(I) = CST(I)
        FVCURRENT(I) = CST(I)
  100 ENDDO
C-----RECORD CURRENT VALUES OF EXTERNAL VARIABLES FOR CONVERGENCE TEST
      DO 200 I=1,NEVAR
        EVBASE(I) = CST(I+NFVAR)
  200 ENDDO
C-----RECORD CURRENT VALUES OF BINARY VARIABLES FOR CONVERGENCE TEST
      DO 300 I=1,NBVAR
        BVBASE(I) = CST(I+NFVAR+NEVAR)
  300 ENDDO
C
      RETURN
      END SUBROUTINE GWM1RMS3SLP_SETBASE
C
      END SUBROUTINE GWM1RMS3PL
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3PP(IOUT,IPERT,FIRSTSIM,LASTSIM,NPERT,GWMSTRG)
C***********************************************************************
C  VERSION: 30MAY2011
C  PURPOSE - SET OUTPUT FLAGS AT BEGINNING OF SIMULATION,
C            PERTURB THE PUMPING RATE FOR RESPONSE MATRIX GENERATION
C-----------------------------------------------------------------------
      USE GWM1RMS3, ONLY : NONLIN,NCON,NDV,FVOLD,AFACT,DEWATER,IPGNA,
     1                     DELINC,SLPITCNT,VARBASE,PGFACT,IBASE,IREF,
     2                     DEWATERQ,NPGNMX,MFCNVRG
      USE GWM1DCV3, ONLY : NFVAR,FVBASE,FVINI,FVMAX,FVNAME,FVCURRENT
      USE GWM1BAS3, ONLY : ZERO,GWMOUT
      USE GWM1OBJ3, ONLY : SOLNTYP,GWM1OBJ3OT2 
      USE GWM1BAS3, ONLY : GWM1BAS3PF,GWM1BAS3PS
      INTEGER(I4B),INTENT(IN)::IOUT,IPERT,NPERT
      LOGICAL(LGT),INTENT(IN):: FIRSTSIM,LASTSIM
      CHARACTER(LEN=200),INTENT(OUT)::GWMSTRG
      CHARACTER(LEN=18)::SPT
      INTEGER(I4B)::I
      ! for debugging
!      integer :: ii, np
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      GWMSTRG = ' '                              ! EMPTY STRING ARRAY
      IF(LASTSIM)THEN                            ! LAST SIMULATION
        CALL SGWM1RMS3PP(4)  
        DEWATER(IPERT) = .FALSE.                 ! INITIALIZE DEWATER FLAGS
        DEWATERQ(IPERT)= .FALSE.                 ! INITIALIZE DEWATERQ FLAGS
      ELSEIF(IPERT.LE.0)THEN  
C-------THIS IS A BASE OR REFERENCE SIMULATION 
        IF(FIRSTSIM .AND. .NOT.LASTSIM)THEN
          DO I=-1,NPERT+1
            DEWATER(I) = .FALSE.                 ! INITIALIZE DEWATER FLAGS
            DEWATERQ(I)= .FALSE.                 ! INITIALIZE DEWATERQ FLAGS
          ENDDO
          IF(SOLNTYP.NE.'FR')CALL SGWM1RMS3PP(1)
          IF(SOLNTYP.EQ.'FR')CALL SGWM1RMS3PP(0)
        ELSEIF(.NOT.FIRSTSIM .AND. .NOT.LASTSIM)THEN
          IF(MFCNVRG(IPERT).AND..NOT.DEWATER(IPERT)
     &               .AND..NOT.DEWATERQ(IPERT))THEN  ! FLOW PROCESS SUCCESSFUL
            CALL SGWM1RMS3PP(2)
          ELSE                                   ! FLOW PROCESS FAILED 
            IF(SLPITCNT.GT.0                     ! MUST BE AN SLP BASE SIMULATION
     &        .AND.NPGNMX.GT.0)THEN              ! AND RESETTING IS ALLOWED
	        CALL SGWM1RMS3PP(3)
              CALL GWM1RMS3PP_RESETBASE(IOUT)   
            ELSE                                 ! EITHER FIRST SIMULATION OR NO RESET
              IF(.NOT.MFCNVRG(IPERT))THEN
                CALL GWM1BAS3PF
     &           ('FLOW PROCESS FAILED TO CONVERGE ON FIRST SIMULATION'
     &             ,0,ZERO)
                CALL GSTOP
     &           ('FLOW PROCESS FAILED TO CONVERGE ON FIRST SIMULATION')
              ELSEIF(DEWATER(IPERT))THEN
                CALL GWM1BAS3PF
     &           ('FLOW PROCESS PRODUCED DEWATERED CONSTRAINED HEADS'
     &             ,0,ZERO)
                CALL GSTOP
     &           ('FLOW PROCESS PRODUCED DEWATERED CONSTRAINED HEADS')
              ELSEIF(DEWATERQ(IPERT))THEN
                CALL GWM1BAS3PF
     &           ('FLOW PROCESS PRODUCED DEWATERED ACTIVE WELL CELL'
     &             ,0,ZERO)
                CALL GSTOP
     &           ('FLOW PROCESS PRODUCED DEWATERED ACTIVE WELL CELL')
              ENDIF
            ENDIF 
            DO I=-1,NPERT+1
              DEWATER(I) = .FALSE.                    ! RESET DEWATER FLAGS
              DEWATERQ(I)= .FALSE.                    ! INITIALIZE DEWATERQ FLAGS
            ENDDO
          ENDIF
        ENDIF
        IF(IPERT.EQ.0  .AND.                     ! THIS IS A BASE SIMULATION
     &    (SOLNTYP.EQ.'LP'.OR.SLPITCNT.EQ.0)     ! THIS IS THE FIRST ONE 
     &    .AND. IBASE.EQ.1)THEN                  ! BASE VALUES ARE AVAILABLE
          DO 100 I=1,NFVAR                          
            FVBASE(I)=FVINI(I)                   ! LOAD BASE VALUES TO FVBASE 
  100     ENDDO                                  
        ENDIF                                    
      ELSEIF(IPERT.GT.0)THEN
C-------THIS IS A PERTURBATION SIMULATION
        IF(SOLNTYP.EQ.'FR')THEN
          CALL GWM1OBJ3OT2                       ! WRITE OBJECTIVE VALUE
          CALL GWM1BAS3PS('    ',0)
          CALL GWM1BAS3PS('    SOLNTYP is FR - GWM Run is Terminated',0)
          CALL GSTOP('  ')
        ENDIF
        IF(IPERT.EQ.1.AND.IPGNA(IPERT).EQ.0)THEN ! WRITE PERTURBATION HEADER
          CALL SGWM1RMS3PP(5)                    
        ENDIF
        IF(IPGNA(IPERT).GT.0)THEN                       ! PERTURBATION NOT SUCCESSFUL
          CALL GWM1RMS3PP_RESETDEL               ! RESET THE PERTURBATION
          DEWATER(IPERT) = .FALSE.               ! RESET DEWATER FLAG
          DEWATERQ(IPERT)= .FALSE.               ! INITIALIZE DEWATERQ FLAG
        ENDIF
        VARBASE(IPERT) = FVBASE(IPERT)                  ! SAVE BASE PUMPING RATE 
        FVBASE(IPERT) = VARBASE(IPERT) + DELINC(IPERT)  ! PERTURB PUMPING RATE
        CALL SGWM1RMS3PP(6)
      ENDIF
  !      add for debugging
!        if (ipert>0)then
!          write(gwmout,150)ipert,fvbase(ipert),dewater(ipert),
!     &    dewaterq(ipert)
!        else
!          np=size(fvbase)
!          do ii=1,np
!            write(gwmout,150)ii,fvbase(ii),dewater(ii),
!     &      dewaterq(ii)
!          enddo
!        endif
!  150 format('in GWM1RMS3PP, ipert fvbase dewater dewaterq = ',
!     &   i2,2x,g23.16,2x,l1,2x,l1)
C
      RETURN
C
      CONTAINS
C
C***********************************************************************
      SUBROUTINE SGWM1RMS3PP(IFLG)
C***********************************************************************
C
C  PURPOSE - WRITE OUTPUT FOR GWM1RMS3PP
      INTEGER(I4B)::IFLG
C-----------------------------------------------------------------------
C
      IF(IFLG.EQ.0)THEN
        CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
        CALL GWM1BAS3PF(
     &'   Running SOLNTYP=FR option (Groundwater Flow Process only):'
     &                    ,0,ZERO)
        CALL GWM1BAS3PF(
     &'  All External and Binary Variables are assigned a value of zero'
     &                    ,0,ZERO)
        CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
        IF(IPERT.LT.0)THEN                       ! THIS IS A REFERENCE RUN
          CALL GWM1BAS3PS('    Running Reference Simulation',0)
          GWMSTRG = 'REFERENCE FLOW PROCESS SIMULATION FOR GWM'
        ELSEIF(IPERT.EQ.0 .AND. IREF.EQ.1)THEN   ! THIS IS A BASE AND REF RUN
          CALL GWM1BAS3PS('    Running Flow Process Simulation',0)
          CALL GWM1BAS3PS('      for both Reference and Base ',0)
          GWMSTRG = 
     &     'FLOW PROCESS SIMULATION FOR GWM FOR BOTH REFERENCE AND BASE'
        ELSEIF(IPERT.EQ.0 .AND. IREF.EQ.0)THEN   ! THIS IS JUST A BASE RUN
          CALL GWM1BAS3PS('    Running Base Flow Process Simulation',0)
          GWMSTRG = 'BASE FLOW PROCESS SIMULATION FOR GWM'
        ENDIF
      ELSEIF(IFLG.EQ.1)THEN
        CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
        CALL GWM1BAS3PF('               Solution Algorithm',0,ZERO)
        CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
        CALL GWM1BAS3PS('  Begin Solution Algorithm',0)
        IF(IPERT.LT.0)THEN                       ! THIS IS A REFERENCE RUN
          CALL GWM1BAS3PS('    Running Reference Simulation',0)
          GWMSTRG = 'REFERENCE FLOW PROCESS SIMULATION FOR GWM'
        ELSEIF(IPERT.EQ.0 .AND. IREF.EQ.1)THEN   ! THIS IS A BASE AND REF RUN
          CALL GWM1BAS3PS('    Running Flow Process Simulation',0)
          CALL GWM1BAS3PS('      for both Reference and Base ',0)
          GWMSTRG = 
     &     'FLOW PROCESS SIMULATION FOR GWM FOR BOTH REFERENCE AND BASE'
        ELSEIF(IPERT.EQ.0 .AND. IREF.EQ.0)THEN   ! THIS IS JUST A BASE RUN
          CALL GWM1BAS3PS('    Running Base Flow Process Simulation',0)
          GWMSTRG = 'BASE FLOW PROCESS SIMULATION FOR GWM'
        ENDIF
      ELSEIF(IFLG.EQ.2)THEN
        IF(NONLIN)CALL GWM1BAS3PS(' ',0)
        IF(NONLIN)
     &    CALL GWM1BAS3PS('  SLP Algorithm: Begin Iteration ',
     &                  SLPITCNT+1)
        CALL GWM1BAS3PS('    Running Base Flow Process Simulation',0)
        GWMSTRG =  'BASE FLOW PROCESS SIMULATION FOR GWM'
      ELSEIF(IFLG.EQ.3)THEN
        CALL GWM1BAS3PS('      Flow Process Failed: Reset Base',0)
      ELSEIF(IFLG.EQ.4)THEN
        CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
        CALL GWM1BAS3PF('         Final Flow Process Simulation',0,ZERO)
        CALL GWM1BAS3PF(
     &'---------------------------------------------------------------'
     &                    ,0,ZERO)
        CALL GWM1BAS3PS('  Running Final Flow Process Simulation',0)
        CALL GWM1BAS3PS('    using Optimal Flow Variable Rates ',0)
        CALL GWM1BAS3PS('    ',0)
        GWMSTRG = 'FINAL FLOW PROCESS SIMULATION FOR GWM'
      ELSEIF(IFLG.EQ.5)THEN
        CALL GWM1BAS3PS('    ',0)
        CALL GWM1BAS3PS('    Calculating Response Matrix',0)
      ELSEIF(IFLG.EQ.6)THEN
        CALL GWM1BAS3PS('      Perturb Flow Variable',IPERT)
        CALL GWM1BAS3PF('       By Perturbation Value:',1,DELINC(IPERT))
        GWMSTRG = 'FLOW PROCESS SIMULATION FOR GWM PERTURBATION'
      ENDIF
      IF(FIRSTSIM)WRITE(IOUT,1005,ERR=990)GWMSTRG
C
 1005 FORMAT(/,5X,A,/)
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
     &2X,'-- STOP EXECUTION (GWM1RMS3PP)')
      CALL GSTOP(' ')
C
      END SUBROUTINE SGWM1RMS3PP
C
C***********************************************************************
      SUBROUTINE GWM1RMS3PP_RESETDEL
C***********************************************************************
C  VERSION: 16JULY2009
C  PURPOSE: CHANGE THE PERTURBATION VALUE BASED ON FAILURE MODE
C-----------------------------------------------------------------------
C
        IF(IPGNA(IPERT).EQ.1)THEN
C---------REDUCE PERTURBATION BUT RETAIN SIGN
          DELINC(IPERT) = DELINC(IPERT)*PGFACT
        ELSEIF(IPGNA(IPERT).EQ.2)THEN
C---------INCREASE PERTURBATION BUT RETAIN SIGN
          DELINC(IPERT) = DELINC(IPERT)/PGFACT
        ELSEIF(IPGNA(IPERT).EQ.3)THEN
C---------REDUCE PERTURBATION AND CHANGE SIGN
          DELINC(IPERT) = -DELINC(IPERT)*PGFACT
        ELSEIF(IPGNA(IPERT).EQ.4)THEN
C---------INCREASE PERTURBATION AND CHANGE SIGN
          DELINC(IPERT) = -DELINC(IPERT)/PGFACT
        ENDIF
C
      RETURN
      END SUBROUTINE GWM1RMS3PP_RESETDEL
C
C***********************************************************************
      SUBROUTINE GWM1RMS3PP_RESETBASE(IOUT)
C***********************************************************************
C  VERSION: 27AUG2013
C  PURPOSE: MOVE BASE FLOW RATES CLOSER TO MOST RECENT SUCCESSFUL SIMULATION 
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: IOUT  

      WRITE(IOUT,10)AFACT
      DO 100 I=1,NFVAR
        FVBASE(I) = (1.0 - AFACT)*FVBASE(I) + AFACT*FVOLD(I)
        FVCURRENT(I) = FVBASE(I)
        WRITE(IOUT,20)I,FVOLD(I),FVBASE(I) 
100   ENDDO
C
      RETURN
   10 FORMAT(/,'New base rates for flow variables (AFACT = ',
     & F5.2,'):',/,
     & 'Var. #',2x,'Rate last iter.',4x,'New rate')
   20 FORMAT(I6,2X,G14.7,5X,G14.7)
      END SUBROUTINE GWM1RMS3PP_RESETBASE
C
      END SUBROUTINE GWM1RMS3PP
C
C***********************************************************************
      SUBROUTINE GWM1RMS3OS(IGRID,KPER,IPERT,HDRY)
C***********************************************************************
C     VERSION: 21JAN2010
C     PURPOSE: CHECK THAT NO MANAGED WELL CELL HAS BEEN DEWATERED
C-----------------------------------------------------------------------
      USE GWM1RMS3, ONLY : DEWATERQ
      USE GWM1BAS3, ONLY : ZERO,GWMOUT
      USE GWM1DCV3, ONLY : NFVAR,FVNCELL,FVKLOC,FVILOC,FVJLOC,FVSP,
     &                     GRDLOCDCV,FVBASE,FVRATIO
      USE GWM1DCV3, ONLY : GWM1DCV3FVCPNT,FVMNW
      USE GWFMNW2MODULE, ONLY:MNW2,MNWMAX,WELLID
      USE GLOBAL,      ONLY: NCOL,NROW,NLAY,HNEW
      USE GWM_SUBS, ONLY: FUZZY_EQUALS, EPSQNET
      INTEGER(I4B),INTENT(IN)::IGRID,KPER,IPERT
      REAL(DP),INTENT(IN)::HDRY
C-----LOCAL VARIABLES
      REAL(DP)::STATE, QDES, QNET
      REAL(DP) :: EPS = 1.0D-5
      INTEGER(I4B)::I,K,IL,IR,IC,MNWID
      LOGICAL :: FAIL1, FAIL2
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----EXAMINE MANAGED WELLS
      DO 200 I=1,NFVAR                      ! LOOP OVER GWM FLOW VARIABLES
        IF(IGRID.EQ.GRDLOCDCV(I))THEN       ! CONSTRAINT ON GRID
        IF(FVSP(I,KPER)) THEN               ! FLOW VAR ACTIVE IN STRESS PERIOD
          IF(FVNCELL(I).GT.0)THEN           ! THIS IS WEL-TYPE FLOW VARIABLE
          CALL GWM1DCV3FVCPNT(I)            ! POINT TO CORRECT CELL INFO
          DO 210 K=1,FVNCELL(I)             ! LOOP OVER CELLS FOR THIS VARIABLE
            IL = FVKLOC(K)                  ! ASSIGN CELL LAYER
            IR = FVILOC(K)                  ! ASSIGN CELL ROW
            IC = FVJLOC(K)                  ! ASSIGN CELL COLUMN
            STATE = HNEW(IC,IR,IL)
            IF(FUZZY_EQUALS(STATE,HDRY,EPS))THEN  ! CELL HAS DEWATERED
              DEWATERQ(IPERT)=.TRUE. 
              WRITE(GWMOUT,1000)I,K,KPER,IR,IC,IL
            ENDIF               
  210     ENDDO
          ELSEIF(FVNCELL(I).LT.0)THEN       ! THIS IS MNW-TYPE FLOW VARIABLE
          CALL GWM1DCV3FVCPNT(I)            ! POINT TO CORRECT CELL INFO
          DO 220 K=1,ABS(FVNCELL(I))        ! LOOP OVER CELLS FOR THIS VARIABLE
            MNWID=FVMNW(I)%FVWTOMNW(K)      ! RETRIEVE LOCATION IN MNW ARRAYS
            QDES = FVBASE(I)*FVRATIO(K)     ! VALUE OF DESIRED FLOW RATE
            QNET = MNW2(18,MNWID)           ! VALUE OF ACTUAL FLOW AT MNW WELL
            ! Failed if Qdes/= stored value of Qnet
            FAIL1=.NOT.FUZZY_EQUALS(QNET,QDES,EPSQNET) 
            ! Failed if Qdes/=0 but well is deactivated at some time step
            FAIL2=(QDES.NE.ZERO) .AND. (MNW2(1,MNWID).NE.1)
            IF(FAIL1 .OR. FAIL2)THEN        ! CELL HAS DEWATERED
              DEWATERQ(IPERT)=.TRUE. 
              WRITE(GWMOUT,1010)I,K,KPER,WELLID(MNWID)
            ENDIF               
  220     ENDDO
          ENDIF
        ENDIF
        ENDIF
  200 ENDDO
C
 1000 FORMAT('GWM FINDS DRY CELL AT DECISION VARIABLE=',I5,
     & ' CELL=',I5,' STRESS PERIOD=',I5,
     & /,'     LOCATED AT ROW=',I8,' COL=',I8,' LAY=',I8)
 1010 FORMAT('GWM FINDS Qnet /= Qdes AT DECISION VARIABLE=',I5,
     & ' CELL=',I5,
     & /,'      STRESS PERIOD=',I5,' LOCATED AT MNW WELL ',A)
      RETURN
      END SUBROUTINE GWM1RMS3OS
C
C***********************************************************************
      SUBROUTINE GWM1RMS3FP(IPERT,NPERT,FIRSTSIM,LASTSIM)
C***********************************************************************
C  VERSION: 21JAN2010
C  PURPOSE - USE SIMULATION RESULTS TO COMPUTE RESPONSE MATRIX AND 
C            AUGMENTED RIGHT HAND SIDE; INCREMENT THE PERTURBATION INDEX
C-----------------------------------------------------------------------
      USE GWM1RMS3, ONLY : IPGNA,IRM,NCON,NDV,IBASE,DEWATER,DEWATERQ,
     1                     HCLOSEG,SLPITCNT,SLPITPRT,DELINC,NPGNMX,
     2                     VARBASE,NRESET,IREF,NPGNA,MFCNVRG
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,FVBASE
      USE GWM1STA3, ONLY : GWM1STA3FP  
      USE GWM1BAS3, ONLY : ZERO,RMFILE
      USE GWM1OBJ3, ONLY : SOLNTYP
      USE GWM1HDC3, ONLY : GWM1HDC3FP 
      USE GWM1STC3, ONLY : GWM1STC3FP 
      USE GWM1BAS3, ONLY : GWM1BAS3PF
      LOGICAL(LGT),INTENT(IN)::FIRSTSIM,LASTSIM
      INTEGER(I4B),INTENT(INOUT)::IPERT
      INTEGER(I4B),INTENT(IN)::NPERT 
C     LOCAL VARIABLES
      INTEGER(I4B)::RSTRT, I
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----IF LAST SIMULATION SKIP TESTS
      IF(LASTSIM)THEN  
        IF(MFCNVRG(IPERT).AND..NOT.DEWATER(IPERT)
     &            .AND..NOT.DEWATERQ(IPERT))THEN  ! FINAL FLOW PROCESS CONVERGED
          CALL GWM1RMS3OT(3,0)             ! WRITE CONSTRAINT STATUS
        ELSEIF(MFCNVRG(IPERT).AND.
     &        (DEWATER(IPERT).OR.DEWATERQ(IPERT)))THEN   ! FINAL FLOW PROCESS HAS DEWATERING
          CALL GWM1RMS3OT(4,0)             ! WRITE CONSTRAINTS AND WARNING 
        ELSEIF(.NOT.MFCNVRG(IPERT))THEN           ! FINAL FLOW PROCESS DID NOT CONVERGE
          CALL GWM1RMS3OT(5,0)             ! WRITE MESSAGE 
        ENDIF
        IPERT = IPERT+1                    ! INCREMENT TO TERMINATE LOOP
        RETURN
      ENDIF
C
C-----TEST FOR PROBLEMS IN THE SIMULATION RESULTS
      IF(IPERT.LE.0.AND.                   ! THIS IS A BASE OR REFERENCE
     &   (.NOT.MFCNVRG(IPERT)              ! SIMULATION AND FLOW PROCESS
     &     .OR.DEWATER(IPERT).OR.DEWATERQ(IPERT)))THEN   ! FAILED TO CONVERGE OR DEWATERED
        NRESET = NRESET + 1                ! INCREMENT BASE RESET COUNTER 
        IF(NRESET.GT.NPGNMX.AND.NPGNMX.GT.0)THEN ! TOO MANY ATTEMPTS 
          CALL GWM1BAS3PF('BASE RESET ATTEMPTS EXCEED MAXIMUM'
     &                     ,0,ZERO)
          CALL GSTOP('BASE RESET ATTEMPTS EXCEED MAXIMUM') ! STOP PROGRAM
        ENDIF
        RETURN                           
C
      ELSEIF(IPERT.EQ.0                    ! THIS IS A BASE SIMULATION
     &       .AND.MFCNVRG(IPERT) 
     &       .AND..NOT.DEWATER(IPERT)
     &       .AND..NOT.DEWATERQ(IPERT))THEN ! WHICH IS SUCCESSFUL, SO
        DO I=1,NPERT
          NPGNA(I) = 0                     ! INITIALIZE PERTURBATION COUNTER
        ENDDO
        NRESET = 0                         ! REZERO BASE RESET COUNTER
C
      ELSEIF(IPERT.GT.0)THEN               ! THIS IS A PERTURBATION SIMULATION
        IF(NPGNMX.EQ.0)THEN                ! DO NOT PERFORM FULL TESTS
          IF(.NOT.MFCNVRG(IPERT))THEN      ! JUST TEST MF CONVERGENCE
            CALL GWM1BAS3PF  
     &       ('FLOW PROCESS FAILED TO CONVERGE ON PERTURBATION RUN'
     &         ,0,ZERO)
            CALL GSTOP                     ! STOP THE RUN
     &       ('FLOW PROCESS FAILED TO CONVERGE ON PERTURBATION RUN')
          ENDIF
        ELSEIF(NPGNMX.GT.0)THEN            ! TEST FOR GOOD RESPONSE
          CALL SGWM1RMS3FP                
          IF(IPGNA(IPERT).EQ.0)THEN        ! PERTURBATION GENERATION SUCCEEDED
            NPGNA(IPERT) = 0               ! RESET PERTURBATION ATTEMPT COUNTER
C
          ELSEIF(IPGNA(IPERT).GT.0) THEN   ! PERTURBATION GENERATION FAILED
            NPGNA(IPERT) = NPGNA(IPERT) + 1  ! INCREMENT PERTURBATION COUNTER 
            IF(NPGNA(IPERT).GT.NPGNMX)THEN ! TOO MANY ATTEMPTS; STOP PROGRAM 
              CALL GWM1BAS3PF('PERTURBATION ATTEMPTS EXCEED MAXIMUM'
     &                       ,0,ZERO)
              CALL GSTOP('PERTURBATION ATTEMPTS EXCEED MAXIMUM')
            ENDIF
            FVBASE(IPERT) = VARBASE(IPERT) ! RESTORE ORIGINAL FLOW RATE 
            RETURN                         ! TRY PERTURBATION AGAIN
          ENDIF
        ENDIF
      ENDIF
C
C-----RESTORE ORIGINAL FLOW RATE FOR PERTURBATION SIMULATION 
      IF(IPERT.GT.0)THEN                   !THIS IS A PERTURBATION SIMULATION
        FVBASE(IPERT) = VARBASE(IPERT)        
      ENDIF
C
C-----WRITE RESPONSE MATRIX OUTPUT HEADER FOR FIRST SIMULATION
      IF((IRM.EQ.1.OR.IRM.EQ.4).AND.FIRSTSIM)THEN
        WRITE(RMFILE)NFVAR,NEVAR,NBVAR,IBASE ! WRITE HEADER FOR TEST
      ENDIF
C
C-----SET RSTRT, THE ROW LOCATION FOR NEXT SET OF CONSTRAINTS
      RSTRT = 1
C
C-----FORMULATE THE HEAD CONSTRAINT COEFFICIENTS FOR THE PERTURBED VARIABLE
      CALL GWM1HDC3FP(RSTRT,IPERT)
C
C-----FORMULATE THE STREAM CONSTRAINT COEFFICIENTS FOR THE PERTURBED VARIABLE
      CALL GWM1STC3FP(RSTRT,IPERT)
C    
C-----FORMULATE THE STATE RESPONSE COEFFICIENTS FOR THE PERTURBED VARIABLE
      CALL GWM1STA3FP(IPERT)
C
C-----WRITE OUTPUT FOR BASE SIMULATIONS
      IF(IPERT.EQ.0)THEN                   ! THIS IS A BASE SIMULATION
        IF(SOLNTYP.EQ.'LP' .OR. SLPITCNT.EQ.0 .OR. SLPITPRT.GE.1)THEN
          CALL GWM1RMS3OT(1,0)             ! WRITE CONSTRAINTS STATUS
        ENDIF
      ENDIF
C
C-----INCREMENT THE PERTURBATION INDEX
      IPERT = IPERT + 1
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(RMFILE,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),RMFILE,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1RMS3FP)')
      CALL GSTOP(' ')
C
      CONTAINS  
C***********************************************************************
      SUBROUTINE SGWM1RMS3FP
C***********************************************************************
C
C   PURPOSE - TEST RESULTS OF PERTURBED SIMULATION FOR GOOD RESPONSE
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO
      USE GWM1RMS3, ONLY : AMAT,DELINC,NRMC,NSIGDIG,DEWATER,DEWATERQ
      USE GWM1HDC3, ONLY : HDCNUM,HDCSTATE,HDCSTATE0
      USE GWM1STC3, ONLY : NSF,NSD
      USE GWM1STC3, ONLY : STCNUM,STCSTATE,STCSTATE0
      USE GWM1HDC3, ONLY : HDCRHS,NHB,NDD,NDF,NGD,HDCDIR
      USE GWM1BAS3, ONLY : GWM1BAS3PS
      ! use gwm1bas3, only : gwmout  ! add for debugging
C-----LOCAL VARIABLES
      REAL(DP)::XNUM
      INTEGER(I4B)::IROW,LNUM,LNUMMAX,COLZ
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  
!     debugging
!      write(gwmout,10)ipert,mfcnvrg(ipert),dewater(ipert),
!     &  dewaterq(ipert)
!   10 format('at top of SGWM1RMS3FP, ipert mfcnvrg dewater dewaterq = '
!     & i3,2x,3l3)
C-----CHECK IF MODFLOW FAILED TO CONVERGE
      IF(.NOT.MFCNVRG(IPERT))THEN
        CALL GWM1BAS3PS('        Perturbation Failed: ',0)
        CALL GWM1BAS3PS('          Flow Process Did Not Converge',0)
        IPGNA(IPERT) = 1                 ! SET PERTURBATION INSTRUCTION
        RETURN
      ENDIF
C
C-----CHECK IF RESPONSE COEFFICIENTS HAVE SUFFICIENT PRECISION
      IF(HDCNUM.GT.0 .OR. STCNUM.GT.0)THEN       ! CHECK RESPONSE PRECISION
        LNUMMAX = 0
        COLZ = 1                                 ! FLAG FOR PRESENCE OF COEFF.
        DO 100 IROW = 1,HDCNUM                   ! LOOP OVER HDC ROWS
          XNUM = HDCSTATE(IROW)-HDCSTATE0(IROW)  ! FIND STATE VALUE DIFFERENCE
          IF(XNUM.NE.ZERO)THEN                   ! DETERMINE THE MAGNITUDE 
            XNUM = XNUM/HCLOSEG                  ! RELATIVE TO MODFLOW CONVERGENCE
            LNUM = IDINT(DLOG10(DABS(XNUM))) + 1 ! FIND ORDER OF MAGNITUDE 
            LNUMMAX = MAX(LNUM,LNUMMAX)          ! COMPARE WITH CURRENT MAXIMUM
            COLZ = 0                             ! RESPONSE COEFF IS PRESENT
          ENDIF
  100   ENDDO
        DO 200 IROW = 1,STCNUM                   ! LOOP OVER STC ROWS
          XNUM = STCSTATE(IROW)-STCSTATE0(IROW)  ! FIND STATE VALUE DIFFERENCE
          IF(XNUM.NE.ZERO)THEN                   ! DETERMINE THE MAGNITUDE 
            XNUM = XNUM/HCLOSEG                  ! RELATIVE TO MODFLOW CONVERGENCE
            LNUM = IDINT(DLOG10(DABS(XNUM))) + 1 ! FIND ORDER OF MAGNITUDE 
            LNUMMAX = MAX(LNUM,LNUMMAX)          ! COMPARE WITH CURRENT MAXIMUM
            COLZ = 0                             ! RESPONSE COEFF IS PRESENT
          ENDIF
  200   ENDDO
        IF(COLZ.EQ.0.AND.LNUMMAX.LT.NSIGDIG)THEN ! INSUFFICIENT PRECISION
          CALL GWM1BAS3PS('        Perturbation Failed: ',0)
          CALL GWM1BAS3PS('          Response Precision Inadequate',0)
          IPGNA(IPERT) = 2                       ! SET PERTURBATION INSTRUCTION
          RETURN
        ENDIF
      ENDIF 
C
C-----CHECK IF ANY HEAD CELLS HAVE BECOME DEWATERED
      IF(DEWATER(IPERT))THEN                     ! A CONSTRAINTED CELL DEWATERED
        CALL GWM1BAS3PS('        Perturbation Failed: ',0)
        CALL GWM1BAS3PS(
     &             '          At Least One Constrained Head is Dry',0)
        IPGNA(IPERT) = 3                  ! SET PERTURBATION INSTRUCTION
        DEWATER(IPERT) = .FALSE.                        ! RESET DEWATER FLAG
        RETURN
      ENDIF
C
C-----CHECK IF ANY WELL CELLS HAVE BECOME DEWATERED
      IF(DEWATERQ(IPERT))THEN                    ! AN ACTIVE WELL DEWATERED
        CALL GWM1BAS3PS('        Perturbation Failed: ',0)
        CALL GWM1BAS3PS(
     &             '          At Least One Active Well is Dry',0)
        IPGNA(IPERT) = 3                  ! SET PERTURBATION INSTRUCTION
        DEWATERQ(IPERT) = .FALSE.                       ! RESET DEWATER FLAG
        RETURN
      ENDIF
C
C-----PASSED ALL CHECKS, PERTURBATION WAS SUCCESSFUL
      IPGNA(IPERT) = 0
C
      RETURN
      END SUBROUTINE SGWM1RMS3FP
C
      END SUBROUTINE GWM1RMS3FP
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3FM 
C***********************************************************************
C  VERSION: 11SEPT2009
C  PURPOSE - ASSEMBLE THE COEFFICIENT MATRIX, COST COEFFICIENT, UPPER BOUNDS
C              AND RIGHT HAND SIDE FOR THE LINEAR PROGRAM
C    INPUT:
C      AMAT    -ON INPUT, AMAT CONTAINS THE RESPONSE MATRIX 
C                 WITH THE FIRST ROWS AND NFVAR COLUMNS OCCUPIED
C    OUTPUT:
C      AMAT    -COEFFICIENT MATRIX OF THE LINEAR PROGRAM
C      CST     -COST COEFFICIENT ARRAY FOR THE LINEAR PROGRAM
C      BNDS    -UPPER BOUND ARRAY FOR THE LINEAR PROGRAM
C      RHS     -RIGHT HAND SIDE ARRAY FOR THE LINEAR PROGRAM
C
C  COLUMNS OF THE COEFFICIENT MATRIX ARE SET IN THIS ORDER:
C      COEFFICIENTS ON PUMPING RATES VARIABLES (NFVAR COLUMNS)
C      COEFFICIENTS ON PUMPING RATES VARIABLES (NEVAR COLUMNS)
C      COEFFICIENTS ON BINARY WELL VARIABLES (NBVAR COLUMNS)
C      COEFFICIENTS ON SLACK/SURPLUS VARIABLS (M COLUMS, HOWEVER, SINCE
C        EACH ROW HAS ONLY ONE SLACK/SURPLUS COEFFICIENT THEY ARE STORED IN 
C        AMAT IN A SINGLE COLUMN)
C
C  ROWS OF THE COEFFICIENT MATRIX CAN BE ASSIGNED IN ANY ORDER
C     EXCEPT THAT THE DECISION VARIABLE CONSTRAINTS MUST BE LAST
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO
      USE GWM1RMS3, ONLY : IRM,AMAT,DELINC,NRMC,NCON,NDV,HCLOSEG
      USE GWM1DCV3, ONLY : NFVAR,GWM1DCV3FM
      USE GWM1OBJ3, ONLY : GWM1OBJ3FM 
      USE GWM1HDC3, ONLY : GWM1HDC3FM 
      USE GWM1STC3, ONLY : GWM1STC3FM
      USE GWM1SMC3, ONLY : GWM1SMC3FM
      USE GWM1DCC3, ONLY : GWM1DCC3FM
C-----LOCAL VARIABLES
      INTEGER(I4B)::RSTRT,NSLK
      REAL(DP)::XNUM
      INTEGER(I4B)::ISIG,NSIG,IROW,IPERT,LNUM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----CALCULATE AND REPORT THE NUMBER OF SIGNIFICANT DIGITS IN RESMAT
      IF(IRM.NE.0)CALL SGWM1RMS3FM
C
C-----FORMULATE THE OBJECTIVE FUNCTION
      CALL GWM1OBJ3FM
C
C-----FORMULATE THE BOUNDS ON DECISION VARIABLES
      CALL GWM1DCV3FM
C
C-----SET NSLK, THE COLUMN OF THE LP MATRIX FOR SLACK/SURPLUS COEFFICIENTS
      NSLK = NDV
C
C-----SET RSTRT, THE ROW LOCATION FOR NEXT SET OF CONSTRAINTS
      RSTRT = 1
C
C-----FORMULATE THE HEAD CONSTRAINTS
      CALL GWM1HDC3FM(RSTRT,NSLK)
C
C-----FORMULATE THE STREAM CONSTRAINTS
      CALL GWM1STC3FM(RSTRT,NSLK)
C
C-----FORMULATE THE SUMMATION CONSTRAINTS
      CALL GWM1SMC3FM(RSTRT,NSLK)
C
C-----FORMULATE THE DECISION VARIABLE CONSTRAINTS
      CALL GWM1DCC3FM(RSTRT,NSLK)
C
      RETURN
C
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1RMS3FM
C***********************************************************************
C
C  PURPOSE - DETERMINE SIGNIFICANT DIGITS IN RESPONSE MATRIX ELEMENT
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWM1BAS3PF,RMFILEF
	USE GWM1STA3, ONLY : STATERES, STANUM
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----WRITE RESPONSE MATRIX TO FORMATTED FILE
      IF(IRM.GE.3)THEN
        IF(NRMC.GT.0)THEN 
          WRITE(RMFILEF,9000)
          CALL UCOLNO(1,NFVAR,5,5,14,RMFILEF)
          DO 300 IROW = 1,NRMC 
            WRITE(RMFILEF,9010) IROW,(AMAT(IROW,IPERT),IPERT=1,NFVAR)
  300     ENDDO
        ENDIF
        IF(STANUM.GT.0)THEN
          WRITE(RMFILEF,9005)
          CALL UCOLNO(1,NFVAR,5,5,14,RMFILEF)
          DO 400 IROW = 1,STANUM
            WRITE(RMFILEF,9010)IROW,(STATERES(IROW,IPERT),IPERT=1,NFVAR)
  400     ENDDO
        ENDIF
      ENDIF
C
      IF(IRM.EQ.5)RETURN ! THIS IS AN EXISTING RESPONSE MATRIX
      ISIG = 0
      NSIG = 0
C
C-----LOOP OVER ALL ELEMENTS OF THE RESPONSE MATRIX
      DO 100 IROW = 1,NRMC 
        DO 80 IPERT = 1,NFVAR
C---------RECOVER THE DIFFERENCE IN STATE VALUES
          XNUM = AMAT(IROW,IPERT)*DELINC(IPERT) 
          IF(XNUM.NE.ZERO)THEN                  
C-----------DETERMINE THE MAGNITUDE RELATIVE TO MODFLOW CONVERGENCE
            XNUM = XNUM/HCLOSEG
C-----------DETERMINE ORDER OF MAGNITUDE OF STATE DIFFERENCE
            LNUM = IDINT(DLOG10(DABS(XNUM))) + 1
            IF(LNUM .LT. 0) LNUM = 0
            ISIG = ISIG + LNUM                  ! ACCUMULATE STATISTICS
            NSIG = NSIG + 1
          ENDIF
   80   ENDDO
  100 ENDDO
C
C-----LOOP OVER ALL ELEMENTS OF THE STATE RESPONSE MATRIX
      DO 200 IROW = 1,STANUM 
        DO 180 IPERT = 1,NFVAR
C---------RECOVER THE DIFFERENCE IN STATE VALUES
          XNUM = STATERES(IROW,IPERT)*DELINC(IPERT)
          IF(XNUM.NE.ZERO)THEN
C-----------DETERMINE THE MAGNITUDE RELATIVE TO MODFLOW CONVERGENCE
            XNUM = XNUM/HCLOSEG
C-----------DETERMINE ORDER OF MAGNITUDE OF STATE DIFFERENCE
            LNUM = IDINT(DLOG10(DABS(XNUM))) + 1
            IF(LNUM .LT. 0) LNUM = 0
            ISIG = ISIG + LNUM                  ! ACCUMULATE STATISTICS
            NSIG = NSIG + 1
          ENDIF
  180   ENDDO
  200 ENDDO
C
      CALL GWM1BAS3PF(' ',0,ZERO)
	IF(ISIG.NE.0)THEN
	  CALL GWM1BAS3PF(
     &   '      Average Number of Significant Digits in Matrix',
     &                1,DBLE(ISIG)/DBLE(NSIG))
      ELSEIF(ISIG.EQ.0)THEN
	  CALL GWM1BAS3PF(
     &   '      Average Number of Sig. Digits in Matrix Less Than 1',
     &                0,ZERO)
      ENDIF

      RETURN
 9000 FORMAT('  GWM CONSTRAINT RESPONSE MATRIX ',/,
     &   '    ROWS IN READ-ORDER FOR SIMULATION-BASED CONSTRAINTS',/,
     &   '       HEAD CONSTRAINTS IN ORDER READ FROM HEDCON',/,
     &   '       STREAM CONSTRAINTS IN ORDER READ FROM STRMCON',/,
     &   '    COLUMNS IN READ-ORDER FOR FLOW RATE VARIABLES',/,
     &   '-------------------------------------------------------')
 9005 FORMAT(/,'  STATE VARIABLE RESPONSE MATRIX ',/,
     &   '    ROWS FOR STATE VARIABLES IN READ-ORDER FROM STAVAR',/,
     &   '    COLUMNS IN READ-ORDER FOR FLOW RATE VARIABLES',/,
     &   '-------------------------------------------------------')
 9010 FORMAT((I5,1X,5G14.6):/(6X,5G14.6))
      END SUBROUTINE SGWM1RMS3FM
C
      END SUBROUTINE GWM1RMS3FM   
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3AP(GWMCNVRG)
C***********************************************************************
C  VERSION: 16JULY2009
C  PURPOSE - CALL THE APPROPRIATE SOLVER FOR THE OPTIMIZATION PROBLEM
C              AND RETURN THE SOLUTION WITH SIGN GIVEN TO PUMPING RATES
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,BIGINF,SMALLEPS
      USE GWM1BAS3, ONLY : GWM1BAS3PS
      USE GWM1RMS3, ONLY : NRMC,NCON, NV, NDV,LPITMAX,AMAT,CST,BNDS,
     &                     RHS,OBJ,SLPITCNT,SLPITMAX,NCONF,NVF
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,FVBASE,FVDIR,FVNAME,
     &                     FVMIN,FVMAX,EVMIN,EVMAX,FVCURRENT,FVON
      USE GWM1OBJ3, ONLY : SOLNTYP,OBJTYP,OBJCNST
      USE GWM1RMS3, ONLY : LASTLP
      ! use gwm1bas3, only : gwmout    ! add for debugging
      LOGICAL(LGT),INTENT(OUT)::GWMCNVRG
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,IFLG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(SOLNTYP.EQ.'NS')THEN
        CALL GWM1BAS3PS('    ',0)
        CALL GWM1BAS3PS('    Response Matrix File Has Been Written',0)
        CALL GSTOP('  ')
      ELSEIF(SOLNTYP.EQ.'MPS')THEN
        CALL GWM1BAS3PS('    ',0)
        CALL GWM1BAS3PS('    Begin Writing MPS File',0)
        CALL GWM1RMS3MPS
        CALL GWM1BAS3PS('    ',0)
        CALL GWM1BAS3PS('    MPS File Has Been Written',0)
        CALL GSTOP('  ')
      ELSEIF(SOLNTYP.EQ.'LP')THEN
        LASTLP = .TRUE.
        CALL SGWM1RMS3AP                    ! SOLVE LP
        CALL GWM1RMS3LP_CHKSOL(IFLG,GWMCNVRG)! CHECK STATUS OF SOLUTION
      ELSEIF(SOLNTYP.EQ.'SLP')THEN
        SLPITCNT = SLPITCNT + 1             ! INCREMENT ITERATION COUNTER
        IF(SLPITCNT.GT.SLPITMAX)THEN        ! TERMINATE PROGRAM IF ITERATIONS EXCEEDED
          CALL GSTOP('PROGRAM STOPPED: SLP ITERATIONS EXCEEDED')
        ENDIF
        LASTLP = .TRUE.                     ! THIS MAY BE LAST - DON'T KNOW
        CALL SGWM1RMS3AP                    ! SOLVE LP AT THIS ITERATION
        CALL GWM1RMS3SLP(IFLG,OBJ,GWMCNVRG) ! CHECK STATUS OF SOLUTION
                                            ! & PREPARE FOR NEXT ITERATION 
      ENDIF

      IF(GWMCNVRG)THEN
C-----TRANSFER OPTIMAL FLOW VARIABLES TO FVBASE FOR FINAL SIMULATION
        DO 300 I=1,NFVAR          
          FVBASE(I) = CST(I)
          FVCURRENT(I) = CST(I)
  300   ENDDO
      ENDIF
C
!      add for debugging
!      write(gwmout,900)
!      write(gwmout,901)(i,cst(i),i=1,nfvar)
!  900 format('in GWM1RMS3AP, cst contains:')
!  901 format(i3,2x,g25.16)
      RETURN
C
      CONTAINS
C***********************************************************************
      SUBROUTINE SGWM1RMS3AP
C***********************************************************************
      USE GWM1RMS3LP, ONLY : GWM1SIMPLEX1,BSOLVE
      USE GWM1BAS3, ONLY : GWM1BAS3PS
      INTEGER(I4B)::J,II
      REAL(DP)::LHS
      REAL(DP), ALLOCATABLE :: RHS1(:)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----ALLOCATE MEMORY FOR USE BY BRANCH AND BOUND ALGORITHM
      ALLOCATE (RHS1(NCON),STAT=ISTAT)
      IF(ISTAT.NE.0)GOTO 992
C
C-----SET BOUNDS AND COSTS ON SLACK/SURPLUS VARIABLES
      DO 100 I=NDV,NV
        BNDS(I) = BIGINF
        CST(I) = ZERO
  100 ENDDO
C
C-----CHECK THAT A MATRIX HAS NO REDUNDANT ROWS
      DO 140 I=1,NCON                            ! LOOP OVER ALL ROWS
        DO 120 II=I+1,NCON                       ! COMPARE WITH ALL LOWER ROWS
          DO 110 J=1,NDV                         ! COMPARE EACH ELEMENT
            IF(AMAT(I,J).NE.AMAT(II,J))GOTO 120  ! ROWS ARE DIFFERENT
  110     ENDDO
C         ALL ELEMENTS ARE IDENTICAL
          CALL GWM1BAS3PS('    ',0)
          CALL GWM1BAS3PS('  WARNING:     CONSTRAINT NUMBER',I)
          CALL GWM1BAS3PS('  IDENTICAL TO CONSTRAINT NUMBER',II)
          CALL GWM1BAS3PS(
     &          '  REDUNDANT CONSTRAINTS MAY LEAD TO SIMPLEX FAILURE',0)
  120   ENDDO
        RHS1(I) = RHS(I)                         ! SAVE RHS FOR COMPARISON
  140 ENDDO
C
C-----SOLVERS TREAT EXTRACTION FLOWS AS POSITIVE.
C     SO, FOR EXTRACTION FLOWS, REVERSE THE COLUMN ON THE ASSOCIATE VARIABLE
      DO 210 J=1,NFVAR
        IF(FVDIR(J).EQ.2)THEN                    ! FLOW VARIABLE IS EXTRACTING
          DO 200 I=1,NRMC                        ! FOR ITS COLUMN IN RESPONSE MATRIX
            AMAT(I,J)= -AMAT(I,J)                ! REVERSE SIGN ON COEFFICIENTS
  200     ENDDO
        ENDIF
  210 ENDDO
C
C-----CALL SOLVER
      CALL GWM1BAS3PS('    ',0)
      CALL GWM1BAS3PS('    Solving Linear Program',0)
      IF(NBVAR.EQ.0) THEN                        ! CALL SIMPLEX LP SOLVER DIRECTLY 
        CALL GWM1SIMPLEX1(NCON,NV,NDV,AMAT,CST,BNDS,RHS,OBJ,
     &                    IFLG,LPITMAX,1) 
        NCONF = NCON
        NVF   = NDV-1
      ELSE                                       ! CALL BRANCH AND BOUND SOLVER 
        CALL BSOLVE(NCON,NV,NDV,AMAT,CST,BNDS,RHS,OBJ,IFLG)
      ENDIF
C
C-----SWITCH SIGN OF OBJ IF MAXIMIZATION
      IF(OBJTYP.EQ.'MAX')THEN                ! CONVERT TO MINIMIZATION
         OBJ = -OBJ
	  ENDIF
C
C-----ADD THE CONSTANT PART BACK INTO THE OBJECTIVE FUNCTION
      OBJ = OBJ + OBJCNST
C
C-----IF VARIABLE VALUES ARE WITHIN ROUNDOFF ERROR OF THEIR BOUNDS, THEN RESET
      DO 300 I=1,NFVAR
        IF(CST(I).GT.(FVMIN(I)-SMALLEPS) .AND. 
     &     CST(I).LT.(FVMIN(I)+SMALLEPS))CST(I)=FVMIN(I) 
        IF(CST(I).GT.(FVMAX(I)-SMALLEPS) .AND. 
     &     CST(I).LT.(FVMAX(I)+SMALLEPS))CST(I)=FVMAX(I) 
        ! When FSTAT=N, FVMAX is set to SMALLEPS; Here, CST reset to zero
        IF(FVON(I).LT.0)CST(I)=ZERO 
  300 ENDDO
C
      DO 400 J=1,NEVAR
	  I = J+NFVAR
        IF(CST(I).GT.(EVMIN(J)-SMALLEPS) .AND. 
     &     CST(I).LT.(EVMIN(J)+SMALLEPS))CST(I)=EVMIN(J) 
        IF(CST(I).GT.(EVMAX(J)-SMALLEPS) .AND. 
     &     CST(I).LT.(EVMAX(J)+SMALLEPS))CST(I)=EVMAX(J)
  400 ENDDO
C
C-----IF DUAL VALUES ARE WITHIN ROUNDOFF ERROR OF ZERO, RESET TO ZERO
      DO 500 I=1,NCON          
        IF(ABS(RHS(I)).LT.SMALLEPS) RHS(I)=ZERO
        IF(OBJTYP.EQ.'MAX')RHS(I)=-RHS(I)    ! CONVERT DUALS FOR MAXIMIZATION PROBLEM
  500 ENDDO
C
C-----FOR EACH CONSTRAINT DETERMINE IF IT IS BINDING AND IF DUAL IS ZERO
      DO 600 I=1,NCON
	  LHS = ZERO
        DO 590 J=1,NFVAR+NEVAR+NBVAR
          LHS = LHS + AMAT(I,J)*CST(J)           ! COMPUTE LHS OF CONSTRAINT
  590   ENDDO
        IF(ABS(LHS-RHS1(I)).LT.SMALLEPS)THEN     ! CONSTRAINT IS BINDING
          IF(RHS(I).EQ.ZERO)THEN                 ! BUT DUAL VALUE IS ZERO
            RHS(I) = BIGINF                      ! SO SET IT TO BIGINF AS FLAG
          ENDIF
        ENDIF
  600 ENDDO
      DEALLOCATE (RHS1,STAT=ISTAT)                ! DEALLOCATE STORED RHS
      IF(ISTAT.NE.0)GOTO 993
C
C-----MODFLOW TREATS EXTRACTION FLOWS AS NEGATIVE.
C     SO, FOR EXTRACTION FLOWS, REVERSE THE DIRECTION OF OPTIMAL FLOW RATE
      DO 700 I=1,NFVAR
        IF (FVDIR(I).EQ.2) CST(I) = -CST(I)
  700 ENDDO
C
      RETURN
C
  992 CONTINUE                                   ! ARRAY-ALLOCATING ERROR
      WRITE(*,9920)
 9920 FORMAT(/,'** ERROR ALLOCATING ARRAY RHS1-- STOP IN SGWM1RMS3AP')
      CALL GSTOP(' ')
C
  993 CONTINUE                                   ! ARRAY-DEALLOCATING ERROR
      WRITE(*,9930)
 9930 FORMAT(/,'** ERROR DEALLOCATING ARRAY RHS1-- STOP IN SGWM1RMS3AP')
      CALL GSTOP(' ')
C
      END SUBROUTINE SGWM1RMS3AP
C
      END SUBROUTINE GWM1RMS3AP
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3MPS
C***********************************************************************
C     VERSION: 21MAR2008
C     PURPOSE: WRITE A STANDARD MPS INPUT DECK 
C           MPS IS THE STANDARD INPUT FORMAT FOR MANY COMMERCIAL 
C           LINEAR PROGRAM SOLVERS
C
C  INPUT CONSISTS OF THE PROBLEM
C
C                  MINIMIZE CX
C                   
C               SUCH THAT AX = B
C
C                         0 <= X <= U
C
C  THIS VERSION ASSUMES THAT THE ORIGINAL PROBLEM IS FULLY INEQUALITY
C  AND THAT THE LAST COLUMN OF A CONTAINS ALL THE SLACK/SURPLUS VARIABLES
C---------------------------------------------------------------------------
      USE GWM1BAS3, ONLY : MPSFILE,ONE,ZERO
      USE GWM1RMS3, ONLY : NRMC,NCON,NV,NDV,AMAT,CST,BNDS,RHS
      USE GWM1DCV3, ONLY : NFVAR,NEVAR,NBVAR,FVNAME,EVNAME,BVNAME,
     &                     FVMIN,FVMAX,EVMIN,EVMAX,FVDIR
C-----LOCAL VARIABLES
      INTEGER(I4B)::I,J,K,NSLK,NPNT
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  
      WRITE(MPSFILE,1000,ERR=990)'GWM MPS'        ! NO NAME IS ASSIGNED
C
      WRITE(MPSFILE,1000,ERR=990)'ROWS   '
      NSLK = NDV                                  ! NDV IS THE SLACK COLUMN
      DO 100 I=1,NCON
        IF(AMAT(I,NSLK).GT.0)THEN
          WRITE(MPSFILE,2000,ERR=990)'L','ROW',I  ! ROWS ARE NAMED BY NUMBER
        ELSEIF(AMAT(I,NSLK).LT.0)THEN
          WRITE(MPSFILE,2000,ERR=990)'G','ROW',I  
        ELSEIF(AMAT(I,NSLK).EQ.0)THEN
          WRITE(MPSFILE,2000,ERR=990)'E','ROW',I
        ENDIF
  100 ENDDO
C
      WRITE(MPSFILE,2000,ERR=990)'N','OBJ'       ! WRITE OBJECTIVE AS LAST ROW
C
C-----SOLVERS TREAT EXTRACTION FLOWS AS POSITIVE.
C     SO, FOR EXTRACTION FLOWS, REVERSE THE COLUMN ON THE ASSOCIATE VARIABLE
      DO 200 J=1,NFVAR
        IF(FVDIR(J).EQ.2)THEN               ! FLOW VARIABLE IS EXTRACTING
          DO 210 I=1,NRMC                   ! FOR ITS COLUMN IN RESPONSE MATRIX
            AMAT(I,J)= -AMAT(I,J)           ! REVERSE SIGN ON COEFFICIENTS
  210     ENDDO
        ENDIF
  200 ENDDO

      WRITE(MPSFILE,1000,ERR=990)'COLUMNS'
      DO 300 J=1,NFVAR
        DO 310 I=1,NCON
          IF(AMAT(I,J).NE.ZERO)THEN           ! WRITE ALL NON-ZERO COEFFS.
            WRITE(MPSFILE,3000,ERR=990)FVNAME(J),'ROW',I,AMAT(I,J)
          ENDIF
  310   ENDDO
        WRITE(MPSFILE,4000,ERR=990)FVNAME(J),'OBJ',CST(J) ! WRITE OBJECTIVE COEFFS.
  300 ENDDO
      DO 400 J=NFVAR+1,NFVAR+NEVAR
        DO 410 I=1,NCON
          IF(AMAT(I,J).NE.0.0)THEN               ! WRITE ALL NON-ZERO COEFFS.
            WRITE(MPSFILE,3000,ERR=990)EVNAME(J-NFVAR),'ROW',I,AMAT(I,J)
          ENDIF
  410   ENDDO
        WRITE(MPSFILE,4000,ERR=990)EVNAME(J-NFVAR),'OBJ',CST(J) ! WRITE OBJECTIVE COEFFS.
  400 ENDDO
      IF(NBVAR.GT.0)THEN
      WRITE(MPSFILE,5000,ERR=990)'BINARY','MARKER','INTORG'
      NPNT = NFVAR+NEVAR
      DO 500 J=NPNT+1,NPNT+NBVAR
        DO 510 I=1,NCON
          IF(AMAT(I,J).NE.0.0)THEN               ! WRITE ALL NON-ZERO COEFFS.
            WRITE(MPSFILE,3000,ERR=990)BVNAME(J-NPNT),'ROW',I,AMAT(I,J)
          ENDIF
  510   ENDDO
        WRITE(MPSFILE,4000,ERR=990)BVNAME(J-NPNT),'OBJ',CST(J) ! WRITE OBJECTIVE COEFFS.
  500 ENDDO
      WRITE(MPSFILE,5000,ERR=990)'BINARY','MARKER','INTEND'
      ENDIF
C
      WRITE(MPSFILE,1000,ERR=990)'RHS    '
      DO 600 I=1,NCON
        WRITE(MPSFILE,6000,ERR=990)'RHS','ROW',I,RHS(I)
  600 ENDDO
C
      WRITE(MPSFILE,1000)'BOUNDS '
      DO 700 J=1,NFVAR
        WRITE(MPSFILE,7000,ERR=990)'UP','BOUND',FVNAME(J),FVMAX(J)
        IF(FVMIN(J).GT.ZERO)THEN 
          WRITE(MPSFILE,7000,ERR=990)'LO','BOUND',FVNAME(J),FVMIN(J)
        ENDIF
  700 CONTINUE
      DO 800 J=NFVAR+1,NFVAR+NEVAR
        K = J-NFVAR
        WRITE(MPSFILE,7000,ERR=990)'UP','BOUND',EVNAME(K),EVMAX(K)
        IF(EVMIN(K).GT.ZERO)THEN 
          WRITE(MPSFILE,7000,ERR=990)'LO','BOUND',EVNAME(K),EVMIN(K)
        ENDIF
  800 ENDDO
      DO 850 J=NFVAR+NEVAR+1,NFVAR+NEVAR+NBVAR
        WRITE(MPSFILE,7000,ERR=990)'UP','BOUND',BVNAME(J-NFVAR+NEVAR),
     &                             ONE
  850 ENDDO
C
      WRITE(MPSFILE,1000,ERR=990)'ENDATA '
      CLOSE(MPSFILE)
C
      RETURN
C  
 1000   FORMAT(A7)
 2000   FORMAT(T2,A1,T5,A3,I4)
 3000   FORMAT(      T5,A8, T15,A3,I4, T25,E12.4)
 4000   FORMAT(      T5,A8, T15,A3,    T25,E12.4)
 5000   FORMAT(      T5,A6, T15,A6,    T40,A6)
 6000   FORMAT(      T5,A3,    T15,A3,I4, T25,E12.4)
 7000   FORMAT(T2,A2,T5,A5,    T15,A8, T25,E12.4)
C  
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(MPSFILE,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),MPSFILE,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1RMS3MPS)')
      CALL GSTOP(' ')
C
      END SUBROUTINE GWM1RMS3MPS
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3SLP (IFLG,OBJ,GWMCNVRG)
C***********************************************************************
C  VERSION: 21JAN2010
C  PURPOSE: CHECK FOR CONVERGENCE OF SEQUENTIAL LINEAR PROGRAM.
C           IF CONVERGENCE NOT ACHIEVED, PREPARE FOR NEXT ITERATION
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : ZERO,ONE
      USE GWM1RMS3, ONLY : AMAT,CST,SLPVCRIT,SLPZCRIT,SLPINFCNT,
     1                     SLPITPRT,SLPITCNT,NINFMX,OBJOLD
      USE GWM1DCV3, ONLY : NFVAR,FVBASE,NEVAR,EVBASE,NBVAR,BVBASE
      USE GWM1BAS3, ONLY : GWM1BAS3PS,GWM1BAS3PF
      INTEGER(I4B),INTENT(IN)::IFLG
      LOGICAL(LGT),INTENT(OUT)::GWMCNVRG
      REAL(DP),INTENT(IN)::OBJ
      INTEGER(I4B)::I,J
      REAL(DP)::FVDMAX,EVDMAX,BVDMAX,FVAMAX,EVAMAX,BVAMAX,OBJCNG
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IFLG.EQ.0)THEN
C-------SUCCESSFUL SOLUTION OF LP: FIND MAXIMUM CHANGE IN VARIABLES
        FVDMAX = ZERO
        EVDMAX = ZERO
        BVDMAX = ZERO
        FVAMAX = ZERO
        EVAMAX = ZERO
        BVAMAX = ZERO
        DO 100 I=1,NFVAR                         ! CHECK ALL FLOW VARIABLES  
          FVDMAX = MAX(FVDMAX,ABS(CST(I)-FVBASE(I)))
          FVAMAX = MAX(FVAMAX,ABS(CST(I)))
  100   ENDDO
        FVDMAX = FVDMAX/(ONE+FVAMAX)
        DO 200 I=1,NEVAR                         ! CHECK ALL EXTERNAL VARIABLES  
          J=I+NFVAR
          EVDMAX = MAX(EVDMAX,ABS(CST(J)-EVBASE(I)))
          EVAMAX = MAX(EVAMAX,ABS(CST(J)))
  200   ENDDO
        EVDMAX = EVDMAX/(ONE+EVAMAX)
        DO 300 I=1,NBVAR                         ! CHECK ALL BINARY VARIABLES  
          J=I+NFVAR+NEVAR
          BVDMAX = MAX(BVDMAX,ABS(CST(J)-BVBASE(I)))
          BVAMAX = MAX(BVAMAX,ABS(CST(J)))
  300   ENDDO
        BVDMAX = BVDMAX/(ONE+BVAMAX)
C-------TEST FOR CONVERGENCE
        OBJCNG = ABS((OBJ-OBJOLD)/(ONE+ABS(OBJ)))
        IF(MAX(FVDMAX,EVDMAX,BVDMAX).GT.SLPVCRIT)THEN 
          GWMCNVRG = .FALSE.                     ! FAILS VARIABLE CRITERIA
        ELSEIF(OBJCNG.GT.SLPZCRIT)THEN 
          GWMCNVRG = .FALSE.                     ! FAILS OBJECTIVE CRITERIA
        ELSE                
          GWMCNVRG = .TRUE.                      ! CRITERIA PASSED                    
        ENDIF
        OBJOLD = OBJ                             ! RESET STORED OBJECTIVE VALUE
        IF(SLPITPRT.GE.1)THEN                    ! WRITE ITERATION STATUS
          CALL GWM1BAS3PS('    Optimal Solution Found',0)
          CALL GWM1BAS3PF('    Objective Value',1,OBJ)
          IF(SLPITCNT.GT.1)THEN
            CALL GWM1BAS3PF(
     &      '    Relative Change in Objective Value',0,0.0)
            CALL GWM1BAS3PF(
     &      '      Needs to be less than SLPZCRIT =',1,SLPZCRIT)
            CALL GWM1BAS3PF(
     &      '      Value at this iteration        =',1,OBJCNG)
            CALL GWM1BAS3PF(
     &      '    Maximum Relative Change in Variables',0,0.0)
            CALL GWM1BAS3PF(
     &      '      Needs to be less than SLPVCRIT =',1,SLPVCRIT)
            IF(NFVAR.GT.0)
     &      CALL GWM1BAS3PF(
     &      '      For Flow Variable Max Change   =',1,FVDMAX)
            IF(NEVAR.GT.0)
     &      CALL GWM1BAS3PF(
     &      '      For Ext Variables Max Change   =',1,EVDMAX)
            IF(NBVAR.GT.0)
     &      CALL GWM1BAS3PF(
     &      '      For Binary Variable Max Change =',1,BVDMAX)
          ENDIF
          CALL GWM1BAS3PS('  SLP Algorithm: End Iteration',SLPITCNT)
        ENDIF
        IF(GWMCNVRG)THEN                         ! CONVERGENCE ACHIEVED
          IF(SLPITPRT.GE.1)
     &    CALL GWM1BAS3PS('  SLP Iterations have converged',0)
        ELSE                                     ! CONVERGENCE NOT ACHIEVED
          IF(SLPITPRT.GE.1)
     &    CALL GWM1BAS3PS('  No SLP convergence at this iteration',0)
          SLPINFCNT = 0
        ENDIF
      ELSE
C-------LP SOLUTION FAILED
        GWMCNVRG = .FALSE.
        IF(SLPITPRT.GE.1) THEN                   ! WRITE ITERATION STATUS
          IF(IFLG.EQ.1)THEN
            CALL GWM1BAS3PS
     &           ('    Problem Is Infeasible At Iteration',SLPITCNT)
          ELSEIF(IFLG.EQ.2)THEN
            CALL GWM1BAS3PS
     &           ('    Problem Is Unbounded At Iteration',SLPITCNT)
          ENDIF
        ENDIF
        SLPINFCNT = SLPINFCNT + 1
        IF(SLPINFCNT.GT.NINFMX)THEN
C---------TERMINATE 
          CALL GWM1BAS3PS('  SLP has failed too many times',0)
          CALL GSTOP(' ')
        ENDIF
      ENDIF
C-----CLEAR COEFFICIENT ARRAY FOR NEXT SLP ITERATION
      AMAT = ZERO
C
      RETURN
      END SUBROUTINE GWM1RMS3SLP  
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3LP_CHKSOL(IFLG,GWMCNVRG)
C***********************************************************************
C  VERSION: 21MAR2008
C  PURPOSE: CHECK THE SOLUTION OF THE LINEAR PROGRAM
C-----------------------------------------------------------------------
      USE GWM1BAS3, ONLY : GWM1BAS3PS
      INTEGER(I4B),INTENT(IN)::IFLG
      LOGICAL(LGT),INTENT(OUT)::GWMCNVRG
      INTEGER(I4B)::I
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(IFLG.EQ.0)THEN
C-------SOLUTION HAS BEEN FOUND 
        CALL GWM1BAS3PS('  Optimal Solution Found',0)
        GWMCNVRG = .TRUE.                         ! SET CONVERGENCE VARIABLE
      ELSEIF(IFLG.EQ.1)THEN
C-------PROBLEM IS INFEASIBLE
        CALL GWM1BAS3PS('  No Feasible Solution Found',0)
        CALL GSTOP(' ')
      ELSEIF(IFLG.EQ.2)THEN
C-------SOLUTION HAS BEEN FOUND 
        CALL GWM1BAS3PS('  Problem is Unbounded',0)
        CALL GSTOP(' ')
      ENDIF
C
      RETURN
      END SUBROUTINE GWM1RMS3LP_CHKSOL
C
C
C***********************************************************************
      SUBROUTINE GWM1RMS3OT(IFLG,NGRIDS)
C***********************************************************************
C  VERSION: 09FEB2012
C  PURPOSE - WRITE OUTPUT FOR GWM SOLUTION  
C-----------------------------------------------------------------------
      USE GLOBAL,   ONLY : NPER
      USE GWM1BAS3, ONLY : ZERO,ONE,BIGINF,GWMOUT,GWMWFILE
      USE GWM1RMS3, ONLY : CST,RHS,BNDS,NRMC,NCON,NV,NDV,NONLIN,RHSIN,
     1                     RHSINF,RANGENAME,RANGENAMEF,RANGEFLG,NCONF,
     2                     NVF,RHSREL,RHSRLL,RHSREU,RHSRLU,CSTREL,
     3                     CSTRLL,CSTREU,CSTRLU,RHSRLB,RHSRUB,RHSROR,
     4                     CSTROR,CSTRLB,CSTRUB,CONTYP,SLPITPRT,SLPITCNT
      USE GWM1DCV3, ONLY : NFVAR,NBVAR,FVBASE,FVSP,FVNCELL,FVKLOC,
     1                     FVILOC,FVJLOC,FVRATIO,FVNAME,GRDLOCDCV
      USE GWM1DCV3, ONLY : GWM1DCV3FVCPNT
      USE GWM1OBJ3, ONLY : OBJTYP,GWM1OBJ3OT
      USE GWM1OBJ3, ONLY : SOLNTYP    
      USE GWM1BAS3, ONLY : GWM1BAS3PF
      USE GWM1HDC3, ONLY : GWM1HDC3OT
      USE GWM1STC3, ONLY : GWM1STC3OT
      USE GWM1SMC3, ONLY : GWM1SMC3OT
      USE GWM1STA3, ONLY : GWM1STA3OT,STANUM
      USE GWM1DCC3, ONLY : GWM1DCC3OT
      INTEGER(I4B),INTENT(IN)::IFLG,NGRIDS
C-----LOCAL VARIABLES
      CHARACTER(LEN=10)::NAME,ENTER,LEAVE
      REAL(DP)::RHSDIFF,SLACK,RHSUB,RHSLB
      REAL(SP)::Q
      INTEGER(I4B)::I,II,I2,K,KPER,RSTRT,INDEX,NCON2,OFLG,GRDNUM,G
      INTEGER(I4B)::LRHSREU,LRHSRLU,LRHSREL,LRHSRLL,ISTART,IEND,ILAST
      LOGICAL :: SLPON, GWMVION, GWMLAST
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C-----ADJUST VALUE OF IFLG AS NECESSARY
      IF(ABS(IFLG).GE.7)THEN
        GWMVION=.TRUE.                      ! THIS IS A CALL FROM GWM-VI
      ELSE
        GWMVION=.FALSE.                    
      ENDIF
      IF(IFLG.LT.0.AND.SLPITPRT.EQ.2)THEN
        SLPON=.TRUE.                        ! THIS IS A CALL FROM SLP 
      ELSE                                  ! AND SOLUTION AT EACH SLP REQUESTED
        SLPON=.FALSE.                    
      ENDIF
      IF(IFLG.EQ.0 .OR. IFLG.EQ.8)THEN
        GWMLAST =.TRUE.                     ! THIS IS CALL AFTER GWMCNVRG=TRUE
      ELSE                                  ! FROM GWM2005 OR GWMVI
        GWMLAST =.FALSE.                    
      ENDIF
      IF(SLPON.OR.GWMLAST)THEN  
        OFLG=2                              ! WRITE SOLUTION TO GLOBAL ONLY
 	  IF(SUM(GWMWFILE).GT.0               ! AT LEAST ONE 'WELL' FILE IS REQUESTED
     &     .AND..NOT.GWMVION)OFLG=6         ! WRITE TO GLOBAL AND 'WELL' FILE
      ELSEIF(IFLG.GT.0.AND.IFLG.LE.7)THEN
        OFLG=IFLG                           ! JUST SET TO ASSIGNED VALUE
	ELSE
        RETURN
	ENDIF
C
C-----WRITE CONSTRAINT STATUS FOR A FLOW PROCESS SIMULATION
      IF(OFLG.EQ.1)THEN     
        IF(STANUM.GT.0)THEN
          CALL GWM1BAS3PF('      Status of State Variable Values '
     &       ,0,ZERO)
          CALL GWM1BAS3PF(
     &'      State Variable Type    Name        Computed Value'
     &       ,0,ZERO)
          CALL GWM1BAS3PF(
     &'      -------------------    ----        --------------'
     &       ,0,ZERO)
          CALL GWM1STA3OT                        ! WRITE STATE VARIABLE STATUS
          CALL GWM1BAS3PF(' ',0,ZERO)              ! BLANK LINE
        ENDIF
        CALL GWM1BAS3PF('      Status of Simulation-Based Constraints '
     &       ,0,ZERO)
        CALL GWM1BAS3PF(
     &'      Constraint Type        Name       Status   Distance To RHS'
     &       ,0,ZERO)
        CALL GWM1BAS3PF(
     &'      ---------------        ----       ------   ---------------'
     &       ,0,ZERO)
C
        RSTRT = 0                           ! RSTRT IS A DUMMY ARGUMENT HERE
        CALL GWM1HDC3OT(RSTRT,1)            ! WRITE HEAD CONSTRAINT STATUS
        CALL GWM1STC3OT(RSTRT,1)            ! WRITE STREAM CONSTRAINT STATUS
        IF(SOLNTYP.NE.'SLP')THEN
          CALL GWM1SMC3OT(RSTRT,1)          ! WRITE SUMMATION CONSTRAINT STATUS
        ELSEIF(SOLNTYP.EQ.'SLP'.AND.SLPITCNT.EQ.0)THEN
          CALL GWM1SMC3OT(RSTRT,1)          ! WRITE SUMMATION CONSTRAINT STATUS
        ELSEIF(SOLNTYP.EQ.'SLP'.AND.SLPITCNT.GT.0)THEN ! EXTERNAL AND BINARY 
          CALL GWM1SMC3OT(RSTRT,3)          ! ARE ALSO WRITTEN IN STATUS
        ENDIF
        WRITE(GWMOUT,7005,ERR=990)
C
C-----WRITE LINEAR PROGRAM OUTPUT
      ELSEIF(OFLG.EQ.2.OR.OFLG.EQ.6)THEN    
        IF(SLPON)THEN                       ! THIS IS AN SLP ITERATION
	    WRITE(GWMOUT,1000,ERR=990)'Solution at This Iteration of SLP'
	  ELSE
	    WRITE(GWMOUT,1000,ERR=990)'Groundwater Management Solution'
          WRITE(GWMOUT,1010,ERR=990)
        ENDIF
        CALL GWM1OBJ3OT                     ! WRITE OBJECTIVE FUNCTION INFO
        RSTRT = 1                           ! SET LOCATION FOR NEXT CONSTRAINT
C
C-------WRITE CONSTRAINTS INFORMATION
        CALL GWM1BAS3PF('       BINDING CONSTRAINTS ',0,ZERO)
        CALL GWM1BAS3PF(
     &       'Constraint Type        Name     Status      Shadow Price'
     &       ,0,ZERO)
        CALL GWM1BAS3PF(
     &       '---------------        ----     ------      ------------'
     &       ,0,ZERO)
        CALL GWM1HDC3OT(RSTRT,2)            ! WRITE HEAD CONSTRAINT INFO
        CALL GWM1STC3OT(RSTRT,2)            ! WRITE STREAM CONSTRAINT INFO
        CALL GWM1SMC3OT(RSTRT,2)            ! WRITE SUMMATION CONSTRAINT INFO
        CALL GWM1DCC3OT(RSTRT)              ! WRITE VARIABLE CONSTRAINT INFO
C
C-------WRITE RANGE ANALYSIS INFORMATION
        IF(RANGEFLG.EQ.1)THEN
          CALL WRITE_RANGE
        ELSE
          WRITE(GWMOUT,5030,ERR=990)
          CALL GWM1BAS3PF(' ',0,ZERO)
          CALL GWM1BAS3PF('       Range Analysis Not Reported ',0,ZERO)
        ENDIF
C
C-----WRITE FINAL STATE VARIABLE AND CONSTRAINT STATUS
      ELSEIF(OFLG.EQ.3 .OR. OFLG.EQ.4)THEN          
        IF(STANUM.GT.0)THEN
          CALL GWM1BAS3PF('      Status of State Variable Values '
     &       ,0,ZERO)
          CALL GWM1BAS3PF(
     &       '        Using Optimal Flow Rate Variable Values'
     &       ,0,ZERO)
          CALL GWM1BAS3PF(
     &'      State Variable Type    Name        Computed Value'
     &       ,0,ZERO)
          CALL GWM1BAS3PF(
     &'      -------------------    ----        --------------'
     &       ,0,ZERO)
          CALL GWM1STA3OT                        ! WRITE STATE VARIABLE STATUS
	    WRITE(GWMOUT,7002,ERR=990)
          CALL GWM1BAS3PF(' ',0,ZERO)              ! BLANK LINE
        ENDIF
C
        CALL GWM1BAS3PF('      Status of Simulation-Based Constraints '
     &       ,0,ZERO)
        CALL GWM1BAS3PF(
     &       '        Using Optimal Flow Rate Variable Values'
     &       ,0,ZERO)
	  WRITE(GWMOUT,7007,ERR=990)
        RSTRT = 0                                 ! RSTRT IS A DUMMY ARGUMENT HERE
        CALL GWM1HDC3OT(RSTRT,4)                  ! WRITE HEAD CONSTRAINT STATUS
        CALL GWM1STC3OT(RSTRT,4)                  ! WRITE STREAM CONSTRAINT STATUS
        CALL GWM1SMC3OT(RSTRT,4)                  ! WRITE SUMMATION CONSTRAINT STATUS
	  WRITE(GWMOUT,7006,ERR=990)
	  WRITE(GWMOUT,7000,ERR=990)
        IF(OFLG.EQ.4)WRITE(GWMOUT,7010,ERR=990)   ! DEWATERING IN HEADS OR WELLS
C
C-----WRITE WARNING MESSAGE: FINAL FLOW PROCESS FAILED TO CONVERGE
      ELSEIF(OFLG.EQ.5)THEN          
        WRITE(GWMOUT,7020,ERR=990)
      ENDIF
C
C-----WRITE OPTIMAL SOLUTION WITH A FORMAT SIMILAR TO THE MODFLOW WELL FILE
      IF(OFLG.EQ.6.OR.OFLG.EQ.7)THEN 
        DO G=1,NGRIDS
	    IF(GWMWFILE(G).GT.0)REWIND(GWMWFILE(G))! IF ALREADY WRITTEN TO, REWIND
        ENDDO
        GRDNUM = 0 
        DO 630 KPER=1,NPER
          ISTART = 0
          I = 1                             
          DO 600 WHILE (I.LE.NFVAR)
C-----------LOOK FOR START AND END ON EACH GRID
            I2=I
            IEND=0
            DO 605 WHILE (IEND.EQ.0)                 
              IF(ISTART.EQ.0.AND.FVSP(I2,KPER))THEN ! FLOW VARIABLE ACTIVE
                GRDNUM = GRDLOCDCV(I2)          ! SET GRID NUMBER
                ISTART = I2                     ! MARK START OF GRID 
              ENDIF
              IF(I2.EQ.NFVAR)THEN               ! THIS IS END OF GRID LIST
                IEND = NFVAR                    ! SET IEND
              ELSEIF(GRDLOCDCV(I2+1).NE.GRDNUM.AND.FVSP(I2+1,KPER))THEN
                IEND = I2                       ! NEXT ACTIVE FV IS ON 
              ENDIF                             ! ANOTHER GRID; SET IEND
              I2 = I2 + 1
  605       ENDDO  
C-----------IF GRID IS FOUND, WRITE THE WELL INFORMATION     
            IF(ISTART.NE.0)THEN   
              I = IEND                          ! SET COUNTER TO IEND
              II = 0
              DO 610 I2=ISTART,IEND
                IF(FVSP(I2,KPER).AND.
     &              FVNCELL(I2).GT.0)II=II+FVNCELL(I2) ! COUNT NUMBER OF CELLS
  610         ENDDO
              IF(GWMWFILE(GRDNUM).NE.0)THEN
                WRITE(GWMWFILE(GRDNUM),6000)II,0,KPER ! WRITE THE NUMBER OF WELLS 
                DO 620 I2 = ISTART,IEND        ! LOOP OVER FV ON THIS GRID
                  IF(FVSP(I2,KPER).AND.        ! FLOW VARIABLE ACTIVE
     &               FVNCELL(I2).GT.0)THEN     ! WEL-TYPE DECISION VARIABLE
                    CALL GWM1DCV3FVCPNT(I2)    ! POINT TO CORRECT CELL INFO
                    DO 615 K=1,FVNCELL(I2)     ! LOOP OVER CELLS FOR THIS VARIABLE
                      Q=CST(I2)*FVRATIO(K)     ! ASSIGN FLOW RATE 
                      WRITE(GWMWFILE(GRDNUM),6010)FVKLOC(K),
     &                        FVILOC(K),FVJLOC(K),Q,
     &                        FVNAME(I2),FVRATIO(K)
  615               ENDDO
                  ENDIF
  620           ENDDO
              ENDIF
              ISTART = 0                       ! RESET GRID START FLAG
            ENDIF 
            I = I + 1                          ! INCREMENT FV COUNTER
  600     ENDDO     
  630   ENDDO
      ENDIF
      
C
 1000 FORMAT(/,70('-'),/,T16,A,/,70('-'))
 1010 FORMAT(1P,/,T8,'OPTIMAL SOLUTION FOUND ',/)
 5030 FORMAT(/,
     &  '  Binding constraint values are ',
     &  'determined from the linear program',/,'    and based on the',
     &  ' response matrix approximation of the flow process.')
 6000 FORMAT(2I10,T50,'Stress Period:',I5)
 6010 FORMAT(3I10,E15.7,T50,'Flow Variable:',A10,' RATIO=',F8.4)
 7000 FORMAT(
     &  '  Precision limitations and nonlinear ',
     &  'response may cause the ',/,'    values of the binding ',
     &  'constraints computed directly by the flow process ',/,
     &  '    to differ from those computed using ',
     &  'the linear program.  ')
 7002 FORMAT(
     &  '  Precision limitations and ',
     &  'nonlinear response may cause ',/,'    the state ',
     &  'variables computed directly by the flow process ',/,
     &  '    to differ from those computed using ',
     &  'the linear program.  ')
 7005 FORMAT(/,
     &'  Distance to RHS is the absolute value of the difference ',
     &'between the',/,'    the right hand side of the constraint and ',
     &'the left side of the',/,'    constraint evaluated using the ',
     &'current set of decision variable values.')    
 7006 FORMAT(/,
     &'  Difference is computed by subtracting right hand side ',
     &'of the constraint ',/,
     &'    from the left side of the constraint.')    
 7007 FORMAT(/,
     &'                                           Simulated    ',
     &'Specified',/,
     &'                                            By Flow        in',/,
     &'      Constraint Type        Name           Process    ',
     &'Constraints   Difference',/,
     &'      ---------------        ----          ----------   ',
     &'----------   ----------')
 7010 FORMAT(/,
     &  '    WARNING: At least one constrained head or active',
     &  ' well dewatered',/,
     &  '             in the final flow process simulation')
 7020 FORMAT(/,
     &  '  WARNING: Final flow process simulation failed',
     &  ' to converge ',/,'    when using the optimal',
     &  ' rates for flow variables. ')
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1RMS3OT)')
      CALL GSTOP(' ')
C
      CONTAINS
C***********************************************************************
      SUBROUTINE WRITE_RANGE
!
          WRITE(GWMOUT,2000,ERR=990)
          CALL GWM1BAS3PF(' ',0,ZERO)
          CALL GWM1BAS3PF('       RANGE ANALYSIS ',0,ZERO)
C---------WRITE CONSTRAINT RANGE ANALYSIS
          WRITE(GWMOUT,3000,ERR=990)
          K = NVF + NBVAR
          II = 0
		DO 100 I=1,NCONF
            NAME=RANGENAMEF(I+NDV-1)                ! NAME OF CONSTRAINT
            RHSDIFF = RHSINF(I)-RHSROR(I)           ! RELATIVE
            IF(CONTYP(I).EQ.0)THEN		          ! CONSTRAINT IS EQUALITY
              SLACK=ZERO
              IF(RHSRLB(I).GT.-BIGINF)THEN
                RHSLB = RHSRLB(I)+RHSDIFF
              ELSE
                RHSLB = -BIGINF
              ENDIF
              IF(RHSRUB(I).LT. BIGINF)THEN
                RHSUB = RHSRUB(I)+RHSDIFF
              ELSE
                RHSUB = BIGINF
              ENDIF
              LRHSREU = RHSREU(I)
              LRHSRLU = RHSRLU(I)
              LRHSREL = RHSREL(I)
              LRHSRLL = RHSRLL(I)
            ELSEIF(CONTYP(I).EQ.1)THEN              ! CONSTRAINT IS INEQUALITY
              II=II+1
              SLACK=CST(II+K)
              IF(RHSRLB(I).GT.-BIGINF)THEN
                RHSLB = RHSRLB(I)+RHSDIFF
              ELSE
                RHSLB = -BIGINF
              ENDIF
              IF(RHSRUB(I).LT. BIGINF)THEN
                RHSUB = RHSRUB(I)+RHSDIFF
              ELSE
                RHSUB = BIGINF
              ENDIF
              LRHSREU = RHSREU(I)
              LRHSRLU = RHSRLU(I)
              LRHSREL = RHSREL(I)
              LRHSRLL = RHSRLL(I)
            ELSEIF(CONTYP(I).EQ.2)THEN              ! CONSTRAINT IS TRANSFORMED 
              II=II+1                               ! INEQUALITY, SO NEED TO
              SLACK=CST(II+K)                       ! SWAP BOUNDS DIRECTION
              IF(RHSRLB(I).GT.-BIGINF)THEN          ! LOWER BD NOT -INFINITY
                RHSUB = RHSINF(I)+SLACK             ! SET UPPER BOUND
              ELSE                                  ! LOWER BOUND IS -INFINITY
                RHSUB = BIGINF                      ! SWAP TO UPPER BD=INF
              ENDIF
              IF(RHSRUB(I).LT. BIGINF)THEN          ! UPPER BD NOT INFINITY
                RHSLB = RHSINF(I)-SLACK
              ELSE                                  ! UPPER BOUND IS INFINITY
                RHSLB = -BIGINF                     ! SWAP TO LOWER BD=-INF
              ENDIF
              LRHSREU = RHSREL(I)                   ! SWAP ENTERING AND LEAVING
              LRHSRLU = RHSRLL(I)                   ! VARIABLES FOR UPPER AND
              LRHSREL = RHSREU(I)                   ! LOWER BOUNDS
              LRHSRLL = RHSRLU(I)
            ENDIF
C
C-----------WRITE THE LOWER BOUND RANGE INFORMATION
            IF(RHSLB.EQ.-BIGINF) THEN
              WRITE(GWMOUT,4000,ERR=990)NAME,SLACK,RHSINF(I)
            ELSEIF(LRHSREL.EQ.0)THEN                ! NO ENTERING VARIABLE
              WRITE(GWMOUT,4005,ERR=990)NAME,SLACK,RHSINF(I),RHSLB
            ELSE
              ENTER=GETNAME(LRHSREL)                ! NAME OF ENTERING VARIABLE
              LEAVE=GETNAME(LRHSRLL)                ! NAME OF LEAVING VARIABLE
              WRITE(GWMOUT,4010,ERR=990)NAME,SLACK,RHSINF(I),
     &                                  RHSLB,ENTER,LEAVE
            ENDIF
C
C-----------WRITE THE UPPER BOUND RANGE INFORMATION
            IF(RHSUB.EQ.BIGINF)THEN                                                
              WRITE(GWMOUT,4020,ERR=990)                                                  
            ELSEIF(LRHSREU.EQ.0)THEN                ! NO ENTERING VARIABLE
              WRITE(GWMOUT,4025,ERR=990)RHSUB
            ELSE  
              ENTER=GETNAME(LRHSREU)                ! NAME OF ENTERING VARIABLE
              LEAVE=GETNAME(LRHSRLU)                ! NAME OF LEAVING VARIABLE
              WRITE(GWMOUT,4030,ERR=990)RHSUB,ENTER,LEAVE
            ENDIF
  100     ENDDO
C
C---------WRITE COST RANGE ANALYSIS
          WRITE(GWMOUT,5000,ERR=990)
          DO 200 I=1,NVF
            NAME=RANGENAMEF(I)                      ! NAME OF CONSTRAINT
            IF(OBJTYP.EQ.'MAX')THEN
              CSTROR(I)=-CSTROR(I)  ! CONVERT BACK TO ORIGINAL COEFF.
              BNDS(I)=-BNDS(I)      ! REDUCED COSTS OF MAXIMIZATION PROBLEM
              CSTRLB(I)=-CSTRLB(I)
              CSTRUB(I)=-CSTRUB(I)
            ENDIF
            IF(OBJTYP.EQ.'MIN')THEN
C
C-------------WRITE THE LOWER BOUND RANGE INFORMATION
              IF(CSTRLB(I).LE. -BIGINF)THEN
                WRITE(GWMOUT,5010,ERR=990)NAME,BNDS(I),CSTROR(I)
              ELSEIF(CSTRLL(I).EQ.0) THEN
                WRITE(GWMOUT,5020,ERR=990)NAME,BNDS(I),CSTROR(I)
              ELSE
                ENTER=GETNAME(CSTREL(I))             ! NAME OF ENTERING VARIABLE
                LEAVE=GETNAME(CSTRLL(I))             ! NAME OF LEAVING VARIABLE
                WRITE(GWMOUT,4010,ERR=990)
     &                      NAME,BNDS(I),CSTROR(I),CSTRLB(I),ENTER,LEAVE
              ENDIF
C
C-------------WRITE THE UPPER BOUND RANGE INFORMATION
              IF(CSTRUB(I).GE.BIGINF)THEN
                WRITE(GWMOUT,4020,ERR=990)
              ELSE
                ENTER=GETNAME(CSTREU(I))             ! NAME OF ENTERING VARIABLE
                LEAVE=GETNAME(CSTRLU(I))             ! NAME OF LEAVING VARIABLE
                WRITE(GWMOUT,4030,ERR=990)CSTRUB(I),ENTER,LEAVE
              ENDIF
            ELSEIF(OBJTYP.EQ.'MAX')THEN ! FOR MAXIMIZATION PROBLEMS, LOWER BOUNDS
                                  ! ARE THE NEGATIVE UPPER BOUNDS DETERMINED
                                  ! FOR THE EQUIVALENT MINIMIZATION PROBLEM
C
C-------------WRITE THE LOWER BOUND RANGE INFORMATION
              IF(CSTRUB(I).LE.-BIGINF)THEN
                WRITE(GWMOUT,5010,ERR=990)NAME,BNDS(I),CSTROR(I)  !UNBOUNDED RANGE
              ELSE
                ENTER=GETNAME(CSTREU(I))             ! NAME OF ENTERING VARIABLE
                LEAVE=GETNAME(CSTRLU(I))             ! NAME OF LEAVING VARIABLE
                WRITE(GWMOUT,4010,ERR=990)
     &                      NAME,BNDS(I),CSTROR(I),CSTRUB(I),ENTER,LEAVE
              ENDIF
C
C-------------WRITE THE UPPER BOUND RANGE INFORMATION
              IF(CSTRLB(I).GE. BIGINF)THEN
                WRITE(GWMOUT,4020,ERR=990)
              ELSEIF(CSTRLL(I).EQ.0) THEN
                WRITE(GWMOUT,5020,ERR=990)NAME,BNDS(I),CSTROR(I)
              ELSE
                ENTER=GETNAME(CSTREL(I))             ! NAME OF ENTERING VARIABLE
                LEAVE=GETNAME(CSTRLL(I))             ! NAME OF LEAVING VARIABLE
                WRITE(GWMOUT,4030,ERR=990)CSTRLB(I),ENTER,LEAVE
              ENDIF
            ENDIF
 200      ENDDO
        RETURN
 2000 FORMAT(/,
     &'  Binding constraint and range analysis values are ',
     &'determined from the linear',/,'    program and based on the',
     &' response matrix approximation of the flow-process.')
 3000 FORMAT(/,/,T7,' Constraint Ranges' ,/,/,T3,
     &   'Lower/Upper Bound are the values of the RHS beyond which ',
     &    'basis will change.',/,
     &    T5,'Leaving is the variable which will leave the basis. ',/,
     &    T5,'Entering is the variable which will enter the basis.',/,
     &    T5,'If the entering or leaving variable is a constraint name,'
     &    ,/,T7,'then the constraint slack variable is active',/,/,
     &    T1,'Constraint',T27,'Original',T41,'Lower/Upper', /,
     &    T1,'Name',T13, 'Slack',T27,'RHS',T41,'Bound',T58,'Entering',
     &    T70,'Leaving',/,
     &    T1,10('-'),T13,10('-'),T27,10('-'),T41,10('-'),T58,10('-'),
     &    T70,10('-'))
 4000 FORMAT(1P,T1,A,T11,ES12.4,T24,ES13.4,T42,'-Infinity',
     &                                   T58,'----- No Change -----')
 4005 FORMAT(1P,T1,A,T11,ES12.4,T24,ES13.4,T38,ES13.4,
     &                                   T58,'-- Degenerate Basis -')
 4010 FORMAT(1P,T1,A,T11,ES12.4,T24,ES13.4,T38,ES13.4,T58,A,T70,A)
 4020 FORMAT(                            T42,' Infinity' 
     &                                   T58,'----- No Change -----',/)
 4025 FORMAT(1P,                         T38,ES13.4,
     &                                   T58,'-- Degenerate Basis -',/)
 4030 FORMAT(1P,                         T38,ES13.4,T58,A,T70,A,/)
 5000 FORMAT(/,/,T7,' Objective-Function Coefficient Ranges ',/,/,T3,
     &   'Lower/Upper Bound are the values of the coefficients',
     &    ' beyond which basis will change.',/,
     &    T5,'Leaving is the variable which will leave the basis. ',/,
     &    T5,'Entering is the variable which will enter the basis.',/,
     &    T5,'If the entering or leaving variable is a constraint name,'
     &    ,/,T7,'then the constraint slack variable is active',/,
     &    T5,'Basic variables are shown with zero reduced cost',/,/,
     &    T1,'Variable',T13,'Reduced',T27,'Original',
     &    T41,'Lower/Upper',/,
     &    T1,'Name',T13,'Cost',T27,'Coefficient',T41,'Bound',T58,
     &    'Entering',T70,'Leaving',/,
     &    T1,10('-'),T13,10('-'),T27,10('-'),T41,10('-'),T58,10('-'),
     &    T70,10('-'))
 5010 FORMAT(1P,T1,A,T11,ES12.4,T24,ES13.4,T42,' Infinity',
     &                                   T58,'----- No Change -----')
 5020 FORMAT(1P,T1,A,T11,ES12.4,T24,ES13.4,             T64,'UNBOUNDED')
 5030 FORMAT(/,
     &  '  Binding constraint values are ',
     &  'determined from the linear program',/,'    and based on the',
     &  ' response matrix approximation of the flow process.')
C
      RETURN
C
C-----ERROR HANDLING
  990 CONTINUE
C-----FILE-WRITING ERROR
      INQUIRE(GWMOUT,NAME=FLNM,FORM=FMTARG,ACCESS=ACCARG,ACTION=FILACT)
      WRITE(*,9900)TRIM(FLNM),GWMOUT,FMTARG,ACCARG,FILACT
 9900 FORMAT(/,1X,'*** ERROR WRITING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (GWM1RMS3OT)')
      CALL GSTOP(' ')

      END SUBROUTINE WRITE_RANGE
C***********************************************************************
      CHARACTER(LEN=10) FUNCTION GETNAME(INDEX)
C***********************************************************************
      INTEGER(I4B),INTENT(IN)::INDEX
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF(INDEX.EQ.0)THEN
        GETNAME = 'NA'
      ELSEIF(INDEX.LE.(NDV-1)-NBVAR)THEN
C-------INDEX NUMBER IS AN ACTUAL NON-BINARY VARIABLE
        GETNAME = RANGENAMEF(INDEX)
      ELSEIF(INDEX.LE.(NDV-1)-NBVAR+NCONF)THEN
C-------INDEX NUMBER IS A SLACK VARIABLE - ADJUST FOR BINARY LOCATIONS
        GETNAME = RANGENAMEF(INDEX+NBVAR)
	  ELSE
C-------INDEX NUMBER IS A BINARY-CONSTRAINED LOWER BOUND - NAME
        GETNAME = 'NA'
      ENDIF
C
      RETURN
      END FUNCTION GETNAME
C
      END SUBROUTINE GWM1RMS3OT
C
C
      END MODULE GWM1RMS3SUBS
C
C
C
C
C
C
C***********************************************************************
C      BEGIN LAPACK ROUTINES
C***********************************************************************
      SUBROUTINE DGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
     $                   EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR,
     $                   WORK, IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), C( * ), FERR( * ), R( * ),
     $                   WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DGESVX uses the LU factorization to compute the solution to a real
*  system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  Error bounds on the solution and a condition estimate are also
*  provided.
*
*  Description
*  ===========
*
*  The following steps are performed:
*
*  1. If FACT = 'E', real scaling factors are computed to equilibrate
*     the system:
*        TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
*        TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
*        TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
*     Whether or not the system will be equilibrated depends on the
*     scaling of the matrix A, but if equilibration is used, A is
*     overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
*     or diag(C)*B (if TRANS = 'T' or 'C').
*
*  2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
*     matrix A (after equilibration if FACT = 'E') as
*        A = P * L * U,
*     where P is a permutation matrix, L is a unit lower triangular
*     matrix, and U is upper triangular.
*
*  3. If some U(i,i)=0, so that U is exactly singular, then the routine
*     returns with INFO = i. Otherwise, the factored form of A is used
*     to estimate the condition number of the matrix A.  If the
*     reciprocal of the condition number is less than machine precision,
*     INFO = N+1 is returned as a warning, but the routine still goes on
*     to solve for X and compute error bounds as described below.
*
*  4. The system of equations is solved for X using the factored form
*     of A.
*
*  5. Iterative refinement is applied to improve the computed solution
*     matrix and calculate error bounds and backward error estimates
*     for it.
*
*  6. If equilibration was used, the matrix X is premultiplied by
*     diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
*     that it solves the original system before equilibration.
*
*  Arguments
*  =========
*
*  FACT    (input) CHARACTER*1
*          Specifies whether or not the factored form of the matrix A is
*          supplied on entry, and if not, whether the matrix A should be
*          equilibrated before it is factored.
*          = 'F':  On entry, AF and IPIV contain the factored form of A.
*                  If EQUED is not 'N', the matrix A has been
*                  equilibrated with scaling factors given by R and C.
*                  A, AF, and IPIV are not modified.
*          = 'N':  The matrix A will be copied to AF and factored.
*          = 'E':  The matrix A will be equilibrated if necessary, then
*                  copied to AF and factored.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Transpose)
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
*          not 'N', then A must have been equilibrated by the scaling
*          factors in R and/or C.  A is not modified if FACT = 'F' or
*          'N', or if FACT = 'E' and EQUED = 'N' on exit.
*
*          On exit, if EQUED .ne. 'N', A is scaled as follows:
*          EQUED = 'R':  A := diag(R) * A
*          EQUED = 'C':  A := A * diag(C)
*          EQUED = 'B':  A := diag(R) * A * diag(C).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  AF      (input or output) DOUBLE PRECISION array, dimension (LDAF,N)
*          If FACT = 'F', then AF is an input argument and on entry
*          contains the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.  If EQUED .ne. 'N', then
*          AF is the factored form of the equilibrated matrix A.
*
*          If FACT = 'N', then AF is an output argument and on exit
*          returns the factors L and U from the factorization A = P*L*U
*          of the original matrix A.
*
*          If FACT = 'E', then AF is an output argument and on exit
*          returns the factors L and U from the factorization A = P*L*U
*          of the equilibrated matrix A (see the description of A for
*          the form of the equilibrated matrix).
*
*  LDAF    (input) INTEGER
*          The leading dimension of the array AF.  LDAF >= max(1,N).
*
*  IPIV    (input or output) INTEGER array, dimension (N)
*          If FACT = 'F', then IPIV is an input argument and on entry
*          contains the pivot indices from the factorization A = P*L*U
*          as computed by DGETRF; row i of the matrix was interchanged
*          with row IPIV(i).
*
*          If FACT = 'N', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization A = P*L*U
*          of the original matrix A.
*
*          If FACT = 'E', then IPIV is an output argument and on exit
*          contains the pivot indices from the factorization A = P*L*U
*          of the equilibrated matrix A.
*
*  EQUED   (input or output) CHARACTER*1
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration (always true if FACT = 'N').
*          = 'R':  Row equilibration, i.e., A has been premultiplied by
*                  diag(R).
*          = 'C':  Column equilibration, i.e., A has been postmultiplied
*                  by diag(C).
*          = 'B':  Both row and column equilibration, i.e., A has been
*                  replaced by diag(R) * A * diag(C).
*          EQUED is an input argument if FACT = 'F'; otherwise, it is an
*          output argument.
*
*  R       (input or output) DOUBLE PRECISION array, dimension (N)
*          The row scale factors for A.  If EQUED = 'R' or 'B', A is
*          multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
*          is not accessed.  R is an input argument if FACT = 'F';
*          otherwise, R is an output argument.  If FACT = 'F' and
*          EQUED = 'R' or 'B', each element of R must be positive.
*
*  C       (input or output) DOUBLE PRECISION array, dimension (N)
*          The column scale factors for A.  If EQUED = 'C' or 'B', A is
*          multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
*          is not accessed.  C is an input argument if FACT = 'F';
*          otherwise, C is an output argument.  If FACT = 'F' and
*          EQUED = 'C' or 'B', each element of C must be positive.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit,
*          if EQUED = 'N', B is not modified;
*          if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
*          diag(R)*B;
*          if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
*          overwritten by diag(C)*B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X
*          to the original system of equations.  Note that A and B are
*          modified on exit if EQUED .ne. 'N', and the solution to the
*          equilibrated system is inv(diag(C))*X if TRANS = 'N' and
*          EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C'
*          and EQUED = 'R' or 'B'.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  RCOND   (output) DOUBLE PRECISION
*          The estimate of the reciprocal condition number of the matrix
*          A after equilibration (if done).  If RCOND is less than the
*          machine precision (in particular, if RCOND = 0), the matrix
*          is singular to working precision.  This condition is
*          indicated by a return code of INFO > 0.
*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (4*N)
*          On exit, WORK(1) contains the reciprocal pivot growth
*          factor norm(A)/norm(U). The "max absolute element" norm is
*          used. If WORK(1) is much less than 1, then the stability
*          of the LU factorization of the (equilibrated) matrix A
*          could be poor. This also means that the solution X, condition
*          estimator RCOND, and forward error bound FERR could be
*          unreliable. If factorization fails with 0<INFO<=N, then
*          WORK(1) contains the reciprocal pivot growth factor for the
*          leading INFO columns of A.
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, and i is
*                <= N:  U(i,i) is exactly zero.  The factorization has
*                       been completed, but the factor U is exactly
*                       singular, so the solution and error bounds
*                       could not be computed. RCOND = 0 is returned.
*                = N+1: U is nonsingular, but RCOND is less than machine
*                       precision, meaning that the matrix is singular
*                       to working precision.  Nevertheless, the
*                       solution and error bounds are computed because
*                       there are a number of situations where the
*                       computed solution can be more accurate than the
*                       value of RCOND would suggest.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU
      CHARACTER          NORM
      INTEGER            I, INFEQU, J
      DOUBLE PRECISION   AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN,
     $                   ROWCND, RPVGRW, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE, DLANTR
      EXTERNAL           LSAME, DLAMCH, DLANGE, DLANTR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGECON, DGEEQU, DGERFS, DGETRF, DGETRS, DLACPY,
     $                   DLAQGE, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      ELSE
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         SMLNUM = DLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      END IF
*
*     Test the input parameters.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) )
     $     THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT.
     $         ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -10
      ELSE
         IF( ROWEQU ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
   10       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -11
            ELSE IF( N.GT.0 ) THEN
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               ROWCND = ONE
            END IF
         END IF
         IF( COLEQU .AND. INFO.EQ.0 ) THEN
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
   20       CONTINUE
            IF( RCMIN.LE.ZERO ) THEN
               INFO = -12
            ELSE IF( N.GT.0 ) THEN
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            ELSE
               COLCND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -14
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -16
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGESVX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*        Compute row and column scalings to equilibrate the matrix A.
*
         CALL DGEEQU( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*           Equilibrate the matrix.
*
            CALL DLAQGE( N, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   EQUED )
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         END IF
      END IF
*
*     Scale the right hand side.
*
      IF( NOTRAN ) THEN
         IF( ROWEQU ) THEN
            DO 40 J = 1, NRHS
               DO 30 I = 1, N
                  B( I, J ) = R( I )*B( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( COLEQU ) THEN
         DO 60 J = 1, NRHS
            DO 50 I = 1, N
               B( I, J ) = C( I )*B( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the LU factorization of A.
*
         CALL DLACPY( 'Full', N, N, A, LDA, AF, LDAF )
         CALL DGETRF( N, N, AF, LDAF, IPIV, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.NE.0 ) THEN
            IF( INFO.GT.0 ) THEN
*
*              Compute the reciprocal pivot growth factor of the
*              leading rank-deficient INFO columns of A.
*
               RPVGRW = DLANTR( 'M', 'U', 'N', INFO, INFO, AF, LDAF,
     $                  WORK )
               IF( RPVGRW.EQ.ZERO ) THEN
                  RPVGRW = ONE
               ELSE
                  RPVGRW = DLANGE( 'M', N, INFO, A, LDA, WORK ) / RPVGRW
               END IF
               WORK( 1 ) = RPVGRW
               RCOND = ZERO
            END IF
            RETURN
         END IF
      END IF
*
*     Compute the norm of the matrix A and the
*     reciprocal pivot growth factor RPVGRW.
*
      IF( NOTRAN ) THEN
         NORM = '1'
      ELSE
         NORM = 'I'
      END IF
      ANORM = DLANGE( NORM, N, N, A, LDA, WORK )
      RPVGRW = DLANTR( 'M', 'U', 'N', N, N, AF, LDAF, WORK )
      IF( RPVGRW.EQ.ZERO ) THEN
         RPVGRW = ONE
      ELSE
         RPVGRW = DLANGE( 'M', N, N, A, LDA, WORK ) / RPVGRW
      END IF
*
*     Compute the reciprocal of the condition number of A.
*
      CALL DGECON( NORM, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO )
*
*     Set INFO = N+1 if the matrix is singular to working precision.
*
      IF( RCOND.LT.DLAMCH( 'Epsilon' ) )
     $   INFO = N + 1
*
*     Compute the solution matrix X.
*
      CALL DLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL DGETRS( TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X,
     $             LDX, FERR, BERR, WORK, IWORK, INFO )
*
*     Transform the solution matrix X to a solution of the original
*     system.
*
      IF( NOTRAN ) THEN
         IF( COLEQU ) THEN
            DO 80 J = 1, NRHS
               DO 70 I = 1, N
                  X( I, J ) = C( I )*X( I, J )
   70          CONTINUE
   80       CONTINUE
            DO 90 J = 1, NRHS
               FERR( J ) = FERR( J ) / COLCND
   90       CONTINUE
         END IF
      ELSE IF( ROWEQU ) THEN
         DO 110 J = 1, NRHS
            DO 100 I = 1, N
               X( I, J ) = R( I )*X( I, J )
  100       CONTINUE
  110    CONTINUE
         DO 120 J = 1, NRHS
            FERR( J ) = FERR( J ) / ROWCND
  120    CONTINUE
      END IF
*
      WORK( 1 ) = RPVGRW
      RETURN
*
*     End of DGESVX
*
      END
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $                   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $                   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*           Use SMALL plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END
      SUBROUTINE DGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEEQU computes row and column scalings intended to equilibrate an
*  M-by-N matrix A and reduce its condition number.  R returns the row
*  scale factors and C the column scale factors, chosen to try to make
*  the largest element in each row and column of the matrix B with
*  elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
*
*  R(i) and C(j) are restricted to be between SMLNUM = smallest safe
*  number and BIGNUM = largest safe number.  Use of these scaling
*  factors is not guaranteed to reduce the condition number of A but
*  works well in practice.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The M-by-N matrix whose equilibration factors are
*          to be computed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  R       (output) DOUBLE PRECISION array, dimension (M)
*          If INFO = 0 or INFO > M, R contains the row scale factors
*          for A.
*
*  C       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0,  C contains the column scale factors for A.
*
*  ROWCND  (output) DOUBLE PRECISION
*          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
*          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
*          AMAX is neither too large nor too small, it is not worth
*          scaling by R.
*
*  COLCND  (output) DOUBLE PRECISION
*          If INFO = 0, COLCND contains the ratio of the smallest
*          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
*          worth scaling by C.
*
*  AMAX    (output) DOUBLE PRECISION
*          Absolute value of largest matrix element.  If AMAX is very
*          close to overflow or very close to underflow, the matrix
*          should be scaled.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i,  and i is
*                <= M:  the i-th row of A is exactly zero
*                >  M:  the (i-M)-th column of A is exactly zero
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   BIGNUM, RCMAX, RCMIN, SMLNUM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
*     Get machine constants.
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
*     Compute row scale factors.
*
      DO 10 I = 1, M
         R( I ) = ZERO
   10 CONTINUE
*
*     Find the maximum element in each row.
*
      DO 30 J = 1, N
         DO 20 I = 1, M
            R( I ) = MAX( R( I ), ABS( A( I, J ) ) )
   20    CONTINUE
   30 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 50 I = 1, M
            IF( R( I ).EQ.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   50    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
*
*        Compute ROWCND = min(R(I)) / max(R(I))
*
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
*     Compute column scale factors
*
      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE
*
*     Find the maximum element in each column,
*     assuming the row scaling computed above.
*
      DO 90 J = 1, N
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), ABS( A( I, J ) )*R( I ) )
   80    CONTINUE
   90 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 110 J = 1, N
            IF( C( J ).EQ.ZERO ) THEN
               INFO = M + J
               RETURN
            END IF
  110    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
*
*        Compute COLCND = min(C(J)) / max(C(J))
*
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
      RETURN
*
*     End of DGEEQU
*
      END
      SUBROUTINE DLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   EQUED )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED
      INTEGER            LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( * ), R( * )
*     ..
*
*  Purpose
*  =======
*
*  DLAQGE equilibrates a general M by N matrix A using the row and
*  scaling factors in the vectors R and C.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M by N matrix A.
*          On exit, the equilibrated matrix.  See EQUED for the form of
*          the equilibrated matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  R       (input) DOUBLE PRECISION array, dimension (M)
*          The row scale factors for A.
*
*  C       (input) DOUBLE PRECISION array, dimension (N)
*          The column scale factors for A.
*
*  ROWCND  (input) DOUBLE PRECISION
*          Ratio of the smallest R(i) to the largest R(i).
*
*  COLCND  (input) DOUBLE PRECISION
*          Ratio of the smallest C(i) to the largest C(i).
*
*  AMAX    (input) DOUBLE PRECISION
*          Absolute value of largest matrix entry.
*
*  EQUED   (output) CHARACTER*1
*          Specifies the form of equilibration that was done.
*          = 'N':  No equilibration
*          = 'R':  Row equilibration, i.e., A has been premultiplied by
*                  diag(R).
*          = 'C':  Column equilibration, i.e., A has been postmultiplied
*                  by diag(C).
*          = 'B':  Both row and column equilibration, i.e., A has been
*                  replaced by diag(R) * A * diag(C).
*
*  Internal Parameters
*  ===================
*
*  THRESH is a threshold value used to decide if row or column scaling
*  should be done based on the ratio of the row or column scaling
*  factors.  If ROWCND < THRESH, row scaling is done, and if
*  COLCND < THRESH, column scaling is done.
*
*  LARGE and SMALL are threshold values used to decide if row scaling
*  should be done based on the absolute size of the largest matrix
*  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( ROWCND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE )
     $     THEN
*
*        No row scaling
*
         IF( COLCND.GE.THRESH ) THEN
*
*           No column scaling
*
            EQUED = 'N'
         ELSE
*
*           Column scaling
*
            DO 20 J = 1, N
               CJ = C( J )
               DO 10 I = 1, M
                  A( I, J ) = CJ*A( I, J )
   10          CONTINUE
   20       CONTINUE
            EQUED = 'C'
         END IF
      ELSE IF( COLCND.GE.THRESH ) THEN
*
*        Row scaling, no column scaling
*
         DO 40 J = 1, N
            DO 30 I = 1, M
               A( I, J ) = R( I )*A( I, J )
   30       CONTINUE
   40    CONTINUE
         EQUED = 'R'
      ELSE
*
*        Row and column scaling
*
         DO 60 J = 1, N
            CJ = C( J )
            DO 50 I = 1, M
               A( I, J ) = CJ*R( I )*A( I, J )
   50       CONTINUE
   60    CONTINUE
         EQUED = 'B'
      END IF
*
      RETURN
*
*     End of DLAQGE
*
      END
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of DLACPY
*
      END
      DOUBLE PRECISION FUNCTION DLANTR( NORM, UPLO, DIAG, M, N, A, LDA,
     $                 WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORM, UPLO
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANTR  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  trapezoidal or triangular matrix A.
*
*  Description
*  ===========
*
*  DLANTR returns the value
*
*     DLANTR = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANTR as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower trapezoidal.
*          = 'U':  Upper trapezoidal
*          = 'L':  Lower trapezoidal
*          Note that A is triangular instead of trapezoidal if M = N.
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A has unit diagonal.
*          = 'N':  Non-unit diagonal
*          = 'U':  Unit diagonal
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0, and if
*          UPLO = 'U', M <= N.  When M = 0, DLANTR is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0, and if
*          UPLO = 'L', N <= M.  When N = 0, DLANTR is set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The trapezoidal matrix A (A is triangular if M = N).
*          If UPLO = 'U', the leading m by n upper trapezoidal part of
*          the array A contains the upper trapezoidal matrix, and the
*          strictly lower triangular part of A is not referenced.
*          If UPLO = 'L', the leading m by n lower trapezoidal part of
*          the array A contains the lower trapezoidal matrix, and the
*          strictly upper triangular part of A is not referenced.  Note
*          that when DIAG = 'U', the diagonal elements of A are not
*          referenced and are assumed to be one.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UDIAG
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         IF( LSAME( DIAG, 'U' ) ) THEN
            VALUE = ONE
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 20 J = 1, N
                  DO 10 I = 1, MIN( M, J-1 )
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = J + 1, M
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            VALUE = ZERO
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 60 J = 1, N
                  DO 50 I = 1, MIN( M, J )
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  DO 70 I = J, M
                     VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         UDIAG = LSAME( DIAG, 'U' )
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 1, N
               IF( ( UDIAG ) .AND. ( J.LE.M ) ) THEN
                  SUM = ONE
                  DO 90 I = 1, J - 1
                     SUM = SUM + ABS( A( I, J ) )
   90             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 100 I = 1, MIN( M, J )
                     SUM = SUM + ABS( A( I, J ) )
  100             CONTINUE
               END IF
               VALUE = MAX( VALUE, SUM )
  110       CONTINUE
         ELSE
            DO 140 J = 1, N
               IF( UDIAG ) THEN
                  SUM = ONE
                  DO 120 I = J + 1, M
                     SUM = SUM + ABS( A( I, J ) )
  120             CONTINUE
               ELSE
                  SUM = ZERO
                  DO 130 I = J, M
                     SUM = SUM + ABS( A( I, J ) )
  130             CONTINUE
               END IF
               VALUE = MAX( VALUE, SUM )
  140       CONTINUE
         END IF
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 150 I = 1, M
                  WORK( I ) = ONE
  150          CONTINUE
               DO 170 J = 1, N
                  DO 160 I = 1, MIN( M, J-1 )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  160             CONTINUE
  170          CONTINUE
            ELSE
               DO 180 I = 1, M
                  WORK( I ) = ZERO
  180          CONTINUE
               DO 200 J = 1, N
                  DO 190 I = 1, MIN( M, J )
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  190             CONTINUE
  200          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               DO 210 I = 1, N
                  WORK( I ) = ONE
  210          CONTINUE
               DO 220 I = N + 1, M
                  WORK( I ) = ZERO
  220          CONTINUE
               DO 240 J = 1, N
                  DO 230 I = J + 1, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  230             CONTINUE
  240          CONTINUE
            ELSE
               DO 250 I = 1, M
                  WORK( I ) = ZERO
  250          CONTINUE
               DO 270 J = 1, N
                  DO 260 I = J, M
                     WORK( I ) = WORK( I ) + ABS( A( I, J ) )
  260             CONTINUE
  270          CONTINUE
            END IF
         END IF
         VALUE = ZERO
         DO 280 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
  280    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         IF( LSAME( UPLO, 'U' ) ) THEN
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = MIN( M, N )
               DO 290 J = 2, N
                  CALL DLASSQ( MIN( M, J-1 ), A( 1, J ), 1, SCALE, SUM )
  290          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 300 J = 1, N
                  CALL DLASSQ( MIN( M, J ), A( 1, J ), 1, SCALE, SUM )
  300          CONTINUE
            END IF
         ELSE
            IF( LSAME( DIAG, 'U' ) ) THEN
               SCALE = ONE
               SUM = MIN( M, N )
               DO 310 J = 1, N
                  CALL DLASSQ( M-J, A( MIN( M, J+1 ), J ), 1, SCALE,
     $                         SUM )
  310          CONTINUE
            ELSE
               SCALE = ZERO
               SUM = ONE
               DO 320 J = 1, N
                  CALL DLASSQ( M-J+1, A( J, J ), 1, SCALE, SUM )
  320          CONTINUE
            END IF
         END IF
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANTR = VALUE
      RETURN
*
*     End of DLANTR
*
      END
      DOUBLE PRECISION FUNCTION DLANGE( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix A.
*
*  Description
*  ===========
*
*  DLANGE returns the value
*
*     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANGE as described
*          above.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.  When M = 0,
*          DLANGE is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.  When N = 0,
*          DLANGE is set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASSQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL DLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      DLANGE = VALUE
      RETURN
*
*     End of DLANGE
*
      END
      SUBROUTINE DGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGECON estimates the reciprocal of the condition number of a general
*  real matrix A, in either the 1-norm or the infinity-norm, using
*  the LU factorization computed by DGETRF.
*
*  An estimate is obtained for norm(inv(A)), and the reciprocal of the
*  condition number is computed as
*     RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies whether the 1-norm condition number or the
*          infinity-norm condition number is required:
*          = '1' or 'O':  1-norm;
*          = 'I':         Infinity-norm.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  ANORM   (input) DOUBLE PRECISION
*          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*          If NORM = 'I', the infinity-norm of the original matrix A.
*
*  RCOND   (output) DOUBLE PRECISION
*          The reciprocal of the condition number of the matrix A,
*          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      DOUBLE PRECISION   AINVNM, SCALE, SL, SMLNUM, SU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACON, DLATRS, DRSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGECON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = DLAMCH( 'Safe minimum' )
*
*     Estimate the norm of inv(A).
*
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL DLACON( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(L).
*
            CALL DLATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N, A,
     $                   LDA, WORK, SL, WORK( 2*N+1 ), INFO )
*
*           Multiply by inv(U).
*
            CALL DLATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN, N,
     $                   A, LDA, WORK, SU, WORK( 3*N+1 ), INFO )
         ELSE
*
*           Multiply by inv(U').
*
            CALL DLATRS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N, A,
     $                   LDA, WORK, SU, WORK( 3*N+1 ), INFO )
*
*           Multiply by inv(L').
*
            CALL DLATRS( 'Lower', 'Transpose', 'Unit', NORMIN, N, A,
     $                   LDA, WORK, SL, WORK( 2*N+1 ), INFO )
         END IF
*
*        Divide X by 1/(SL*SU) if doing so will not cause overflow.
*
         SCALE = SL*SU
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IDAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL DRSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
      RETURN
*
*     End of DGECON
*
      END
      SUBROUTINE DGERFS( TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB,
     $                   X, LDX, FERR, BERR, WORK, IWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  DGERFS improves the computed solution to a system of linear
*  equations and provides error bounds and backward error estimates for
*  the solution.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B     (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices B and X.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The original N-by-N matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDAF    (input) INTEGER
*          The leading dimension of the array AF.  LDAF >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          The right hand side matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          On entry, the solution matrix X, as computed by DGETRS.
*          On exit, the improved solution matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  FERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The estimated forward error bound for each solution vector
*          X(j) (the j-th column of the solution matrix X).
*          If XTRUE is the true solution corresponding to X(j), FERR(j)
*          is an estimated upper bound for the magnitude of the largest
*          element in (X(j) - XTRUE) divided by the magnitude of the
*          largest element in X(j).  The estimate is as reliable as
*          the estimate for RCOND, and is almost always a slight
*          overestimate of the true error.
*
*  BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*          The componentwise relative backward error of each solution
*          vector X(j) (i.e., the smallest relative change in
*          any element of A or B that makes X(j) an exact solution).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Internal Parameters
*  ===================
*
*  ITMAX is the maximum number of steps of iterative refinement.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      DOUBLE PRECISION   TWO
      PARAMETER          ( TWO = 2.0D+0 )
      DOUBLE PRECISION   THREE
      PARAMETER          ( THREE = 3.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      CHARACTER          TRANST
      INTEGER            COUNT, I, J, K, KASE, NZ
      DOUBLE PRECISION   EPS, LSTRES, S, SAFE1, SAFE2, SAFMIN, XK
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DGETRS, DLACON, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGERFS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         DO 10 J = 1, NRHS
            FERR( J ) = ZERO
            BERR( J ) = ZERO
   10    CONTINUE
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
         TRANST = 'T'
      ELSE
         TRANST = 'N'
      END IF
*
*     NZ = maximum number of nonzero elements in each row of A, plus 1
*
      NZ = N + 1
      EPS = DLAMCH( 'Epsilon' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SAFE1 = NZ*SAFMIN
      SAFE2 = SAFE1 / EPS
*
*     Do for each right hand side
*
      DO 140 J = 1, NRHS
*
         COUNT = 1
         LSTRES = THREE
   20    CONTINUE
*
*        Loop until stopping criterion is satisfied.
*
*        Compute residual R = B - op(A) * X,
*        where op(A) = A, A**T, or A**H, depending on TRANS.
*
         CALL DCOPY( N, B( 1, J ), 1, WORK( N+1 ), 1 )
         CALL DGEMV( TRANS, N, N, -ONE, A, LDA, X( 1, J ), 1, ONE,
     $               WORK( N+1 ), 1 )
*
*        Compute componentwise relative backward error from formula
*
*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) )
*
*        where abs(Z) is the componentwise absolute value of the matrix
*        or vector Z.  If the i-th component of the denominator is less
*        than SAFE2, then SAFE1 is added to the i-th components of the
*        numerator and denominator before dividing.
*
         DO 30 I = 1, N
            WORK( I ) = ABS( B( I, J ) )
   30    CONTINUE
*
*        Compute abs(op(A))*abs(X) + abs(B).
*
         IF( NOTRAN ) THEN
            DO 50 K = 1, N
               XK = ABS( X( K, J ) )
               DO 40 I = 1, N
                  WORK( I ) = WORK( I ) + ABS( A( I, K ) )*XK
   40          CONTINUE
   50       CONTINUE
         ELSE
            DO 70 K = 1, N
               S = ZERO
               DO 60 I = 1, N
                  S = S + ABS( A( I, K ) )*ABS( X( I, J ) )
   60          CONTINUE
               WORK( K ) = WORK( K ) + S
   70       CONTINUE
         END IF
         S = ZERO
         DO 80 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               S = MAX( S, ABS( WORK( N+I ) ) / WORK( I ) )
            ELSE
               S = MAX( S, ( ABS( WORK( N+I ) )+SAFE1 ) /
     $             ( WORK( I )+SAFE1 ) )
            END IF
   80    CONTINUE
         BERR( J ) = S
*
*        Test stopping criterion. Continue iterating if
*           1) The residual BERR(J) is larger than machine epsilon, and
*           2) BERR(J) decreased by at least a factor of 2 during the
*              last iteration, and
*           3) At most ITMAX iterations tried.
*
         IF( BERR( J ).GT.EPS .AND. TWO*BERR( J ).LE.LSTRES .AND.
     $       COUNT.LE.ITMAX ) THEN
*
*           Update solution and try again.
*
            CALL DGETRS( TRANS, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,
     $                   INFO )
            CALL DAXPY( N, ONE, WORK( N+1 ), 1, X( 1, J ), 1 )
            LSTRES = BERR( J )
            COUNT = COUNT + 1
            GO TO 20
         END IF
*
*        Bound error from formula
*
*        norm(X - XTRUE) / norm(X) .le. FERR =
*        norm( abs(inv(op(A)))*
*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X)
*
*        where
*          norm(Z) is the magnitude of the largest component of Z
*          inv(op(A)) is the inverse of op(A)
*          abs(Z) is the componentwise absolute value of the matrix or
*             vector Z
*          NZ is the maximum number of nonzeros in any row of A, plus 1
*          EPS is machine epsilon
*
*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B))
*        is incremented by SAFE1 if the i-th component of
*        abs(op(A))*abs(X) + abs(B) is less than SAFE2.
*
*        Use DLACON to estimate the infinity-norm of the matrix
*           inv(op(A)) * diag(W),
*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) )))
*
         DO 90 I = 1, N
            IF( WORK( I ).GT.SAFE2 ) THEN
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I )
            ELSE
               WORK( I ) = ABS( WORK( N+I ) ) + NZ*EPS*WORK( I ) + SAFE1
            END IF
   90    CONTINUE
*
         KASE = 0
  100    CONTINUE
         CALL DLACON( N, WORK( 2*N+1 ), WORK( N+1 ), IWORK, FERR( J ),
     $                KASE )
         IF( KASE.NE.0 ) THEN
            IF( KASE.EQ.1 ) THEN
*
*              Multiply by diag(W)*inv(op(A)**T).
*
               CALL DGETRS( TRANST, N, 1, AF, LDAF, IPIV, WORK( N+1 ),
     $                      N, INFO )
               DO 110 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  110          CONTINUE
            ELSE
*
*              Multiply by inv(op(A))*diag(W).
*
               DO 120 I = 1, N
                  WORK( N+I ) = WORK( I )*WORK( N+I )
  120          CONTINUE
               CALL DGETRS( TRANS, N, 1, AF, LDAF, IPIV, WORK( N+1 ), N,
     $                      INFO )
            END IF
            GO TO 100
         END IF
*
*        Normalize error.
*
         LSTRES = ZERO
         DO 130 I = 1, N
            LSTRES = MAX( LSTRES, ABS( X( I, J ) ) )
  130    CONTINUE
         IF( LSTRES.NE.ZERO )
     $      FERR( J ) = FERR( J ) / LSTRES
*
  140 CONTINUE
*
      RETURN
*
*     End of DGERFS
*
      END
      SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   ABSXI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX ).NE.ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE.LT.ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of DLASSQ
*
      END
      SUBROUTINE DLACON( N, V, X, ISGN, EST, KASE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      INTEGER            KASE, N
      DOUBLE PRECISION   EST
*     ..
*     .. Array Arguments ..
      INTEGER            ISGN( * )
      DOUBLE PRECISION   V( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLACON estimates the 1-norm of a square, real matrix A.
*  Reverse communication is used for evaluating matrix-vector products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  V      (workspace) DOUBLE PRECISION array, dimension (N)
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and DLACON must be re-called with all the other parameters
*         unchanged.
*
*  ISGN   (workspace) INTEGER array, dimension (N)
*
*  EST    (output) DOUBLE PRECISION
*         An estimate (a lower bound) for norm(A).
*
*  KASE   (input/output) INTEGER
*         On the initial call to DLACON, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from DLACON, KASE will again be 0.
*
*  Further Details
*  ======= =======
*
*  Contributed by Nick Higham, University of Manchester.
*  Originally named SONEST, dated March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      DOUBLE PRECISION   ALTSGN, ESTOLD, TEMP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM
      EXTERNAL           IDAMAX, DASUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, NINT, SIGN
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
      IF( KASE.EQ.0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / DBLE( N )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
*
      GO TO ( 20, 40, 70, 110, 140 )JUMP
*
*     ................ ENTRY   (JUMP = 1)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
*
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
*        ... QUIT
         GO TO 150
      END IF
      EST = DASUM( N, X, 1 )
*
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
*
*     ................ ENTRY   (JUMP = 2)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
*
   40 CONTINUE
      J = IDAMAX( N, X, 1 )
      ITER = 2
*
*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
*
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( J ) = ONE
      KASE = 1
      JUMP = 3
      RETURN
*
*     ................ ENTRY   (JUMP = 3)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
   70 CONTINUE
      CALL DCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = DASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) ).NE.ISGN( I ) )
     $      GO TO 90
   80 CONTINUE
*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
*
   90 CONTINUE
*     TEST FOR CYCLING.
      IF( EST.LE.ESTOLD )
     $   GO TO 120
*
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
*
*     ................ ENTRY   (JUMP = 4)
*     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
*
  110 CONTINUE
      JLAST = J
      J = IDAMAX( N, X, 1 )
      IF( ( X( JLAST ).NE.ABS( X( J ) ) ) .AND. ( ITER.LT.ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
*
*     ITERATION COMPLETE.  FINAL STAGE.
*
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+DBLE( I-1 ) / DBLE( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
*
*     ................ ENTRY   (JUMP = 5)
*     X HAS BEEN OVERWRITTEN BY A*X.
*
  140 CONTINUE
      TEMP = TWO*( DASUM( N, X, 1 ) / DBLE( 3*N ) )
      IF( TEMP.GT.EST ) THEN
         CALL DCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
*
  150 CONTINUE
      KASE = 0
      RETURN
*
*     End of DLACON
*
      END
      SUBROUTINE DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X, SCALE,
     $                   CNORM, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, NORMIN, TRANS, UPLO
      INTEGER            INFO, LDA, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), CNORM( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DLATRS solves one of the triangular systems
*
*     A *x = s*b  or  A'*x = s*b
*
*  with scaling to prevent overflow.  Here A is an upper or lower
*  triangular matrix, A' denotes the transpose of A, x and b are
*  n-element vectors, and s is a scaling factor, usually less than
*  or equal to 1, chosen so that the components of x will be less than
*  the overflow threshold.  If the unscaled problem will not cause
*  overflow, the Level 2 BLAS routine DTRSV is called.  If the matrix A
*  is singular (A(j,j) = 0 for some j), then s is set to 0 and a
*  non-trivial solution to A*x = 0 is returned.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  TRANS   (input) CHARACTER*1
*          Specifies the operation applied to A.
*          = 'N':  Solve A * x = s*b  (No transpose)
*          = 'T':  Solve A'* x = s*b  (Transpose)
*          = 'C':  Solve A'* x = s*b  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  NORMIN  (input) CHARACTER*1
*          Specifies whether CNORM has been set or not.
*          = 'Y':  CNORM contains the column norms on entry
*          = 'N':  CNORM is not set on entry.  On exit, the norms will
*                  be computed and stored in CNORM.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of the array A contains the upper
*          triangular matrix, and the strictly lower triangular part of
*          A is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of the array A contains the lower triangular
*          matrix, and the strictly upper triangular part of A is not
*          referenced.  If DIAG = 'U', the diagonal elements of A are
*          also not referenced and are assumed to be 1.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max (1,N).
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the right hand side b of the triangular system.
*          On exit, X is overwritten by the solution vector x.
*
*  SCALE   (output) DOUBLE PRECISION
*          The scaling factor s for the triangular system
*             A * x = s*b  or  A'* x = s*b.
*          If SCALE = 0, the matrix A is singular or badly scaled, and
*          the vector x is an exact or approximate solution to A*x = 0.
*
*  CNORM   (input or output) DOUBLE PRECISION array, dimension (N)
*
*          If NORMIN = 'Y', CNORM is an input argument and CNORM(j)
*          contains the norm of the off-diagonal part of the j-th column
*          of A.  If TRANS = 'N', CNORM(j) must be greater than or equal
*          to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j)
*          must be greater than or equal to the 1-norm.
*
*          If NORMIN = 'N', CNORM is an output argument and CNORM(j)
*          returns the 1-norm of the offdiagonal part of the j-th column
*          of A.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*
*  Further Details
*  ======= =======
*
*  A rough bound on x is computed; if that is less than overflow, DTRSV
*  is called, otherwise, specific code is used which checks for possible
*  overflow or divide-by-zero at every operation.
*
*  A columnwise scheme is used for solving A*x = b.  The basic algorithm
*  if A is lower triangular is
*
*       x[1:n] := b[1:n]
*       for j = 1, ..., n
*            x(j) := x(j) / A(j,j)
*            x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
*       end
*
*  Define bounds on the components of x after j iterations of the loop:
*     M(j) = bound on x[1:j]
*     G(j) = bound on x[j+1:n]
*  Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}.
*
*  Then for iteration j+1 we have
*     M(j+1) <= G(j) / | A(j+1,j+1) |
*     G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
*            <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )
*
*  where CNORM(j+1) is greater than or equal to the infinity-norm of
*  column j+1 of A, not counting the diagonal.  Hence
*
*     G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
*                  1<=i<=j
*  and
*
*     |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
*                                   1<=i< j
*
*  Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTRSV if the
*  reciprocal of the largest M(j), j=1,..,n, is larger than
*  max(underflow, 1/overflow).
*
*  The bound on x(j) is also used to determine when a step in the
*  columnwise method can be performed without fear of overflow.  If
*  the computed bound is greater than a large constant, x is scaled to
*  prevent overflow, but if the bound overflows, x is set to 0, x(j) to
*  1, and scale to 0, and a non-trivial solution to A*x = 0 is found.
*
*  Similarly, a row-wise scheme is used to solve A'*x = b.  The basic
*  algorithm for A upper triangular is
*
*       for j = 1, ..., n
*            x(j) := ( b(j) - A[1:j-1,j]' * x[1:j-1] ) / A(j,j)
*       end
*
*  We simultaneously compute two bounds
*       G(j) = bound on ( b(i) - A[1:i-1,i]' * x[1:i-1] ), 1<=i<=j
*       M(j) = bound on x(i), 1<=i<=j
*
*  The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we
*  add the constraint G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1.
*  Then the bound on x(j) is
*
*       M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
*
*            <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
*                      1<=i<=j
*
*  and we can safely call DTRSV if 1/M(n) and 1/G(n) are both greater
*  than max(underflow, 1/overflow).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN, NOUNIT, UPPER
      INTEGER            I, IMAX, J, JFIRST, JINC, JLAST
      DOUBLE PRECISION   BIGNUM, GROW, REC, SMLNUM, SUMJ, TJJ, TJJS,
     $                   TMAX, TSCAL, USCAL, XBND, XJ, XMAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DASUM, DDOT, DLAMCH
      EXTERNAL           LSAME, IDAMAX, DASUM, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, DTRSV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Test the input parameters.
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $         LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT.
     $         LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      BIGNUM = ONE / SMLNUM
      SCALE = ONE
*
      IF( LSAME( NORMIN, 'N' ) ) THEN
*
*        Compute the 1-norm of each column, not including the diagonal.
*
         IF( UPPER ) THEN
*
*           A is upper triangular.
*
            DO 10 J = 1, N
               CNORM( J ) = DASUM( J-1, A( 1, J ), 1 )
   10       CONTINUE
         ELSE
*
*           A is lower triangular.
*
            DO 20 J = 1, N - 1
               CNORM( J ) = DASUM( N-J, A( J+1, J ), 1 )
   20       CONTINUE
            CNORM( N ) = ZERO
         END IF
      END IF
*
*     Scale the column norms by TSCAL if the maximum element in CNORM is
*     greater than BIGNUM.
*
      IMAX = IDAMAX( N, CNORM, 1 )
      TMAX = CNORM( IMAX )
      IF( TMAX.LE.BIGNUM ) THEN
         TSCAL = ONE
      ELSE
         TSCAL = ONE / ( SMLNUM*TMAX )
         CALL DSCAL( N, TSCAL, CNORM, 1 )
      END IF
*
*     Compute a bound on the computed solution vector to see if the
*     Level 2 BLAS routine DTRSV can be used.
*
      J = IDAMAX( N, X, 1 )
      XMAX = ABS( X( J ) )
      XBND = XMAX
      IF( NOTRAN ) THEN
*
*        Compute the growth in A * x = b.
*
         IF( UPPER ) THEN
            JFIRST = N
            JLAST = 1
            JINC = -1
         ELSE
            JFIRST = 1
            JLAST = N
            JINC = 1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 50
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, G(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 30 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              M(j) = G(j-1) / abs(A(j,j))
*
               TJJ = ABS( A( J, J ) )
               XBND = MIN( XBND, MIN( ONE, TJJ )*GROW )
               IF( TJJ+CNORM( J ).GE.SMLNUM ) THEN
*
*                 G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
*
                  GROW = GROW*( TJJ / ( TJJ+CNORM( J ) ) )
               ELSE
*
*                 G(j) could overflow, set GROW to 0.
*
                  GROW = ZERO
               END IF
   30       CONTINUE
            GROW = XBND
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 40 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 50
*
*              G(j) = G(j-1)*( 1 + CNORM(j) )
*
               GROW = GROW*( ONE / ( ONE+CNORM( J ) ) )
   40       CONTINUE
         END IF
   50    CONTINUE
*
      ELSE
*
*        Compute the growth in A' * x = b.
*
         IF( UPPER ) THEN
            JFIRST = 1
            JLAST = N
            JINC = 1
         ELSE
            JFIRST = N
            JLAST = 1
            JINC = -1
         END IF
*
         IF( TSCAL.NE.ONE ) THEN
            GROW = ZERO
            GO TO 80
         END IF
*
         IF( NOUNIT ) THEN
*
*           A is non-unit triangular.
*
*           Compute GROW = 1/G(j) and XBND = 1/M(j).
*           Initially, M(0) = max{x(i), i=1,...,n}.
*
            GROW = ONE / MAX( XBND, SMLNUM )
            XBND = GROW
            DO 60 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
*
               XJ = ONE + CNORM( J )
               GROW = MIN( GROW, XBND / XJ )
*
*              M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
*
               TJJ = ABS( A( J, J ) )
               IF( XJ.GT.TJJ )
     $            XBND = XBND*( TJJ / XJ )
   60       CONTINUE
            GROW = MIN( GROW, XBND )
         ELSE
*
*           A is unit triangular.
*
*           Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
*
            GROW = MIN( ONE, ONE / MAX( XBND, SMLNUM ) )
            DO 70 J = JFIRST, JLAST, JINC
*
*              Exit the loop if the growth factor is too small.
*
               IF( GROW.LE.SMLNUM )
     $            GO TO 80
*
*              G(j) = ( 1 + CNORM(j) )*G(j-1)
*
               XJ = ONE + CNORM( J )
               GROW = GROW / XJ
   70       CONTINUE
         END IF
   80    CONTINUE
      END IF
*
      IF( ( GROW*TSCAL ).GT.SMLNUM ) THEN
*
*        Use the Level 2 BLAS solve if the reciprocal of the bound on
*        elements of X is not too small.
*
         CALL DTRSV( UPLO, TRANS, DIAG, N, A, LDA, X, 1 )
      ELSE
*
*        Use a Level 1 BLAS solve, scaling intermediate results.
*
         IF( XMAX.GT.BIGNUM ) THEN
*
*           Scale X so that its components are less than or equal to
*           BIGNUM in absolute value.
*
            SCALE = BIGNUM / XMAX
            CALL DSCAL( N, SCALE, X, 1 )
            XMAX = BIGNUM
         END IF
*
         IF( NOTRAN ) THEN
*
*           Solve A * x = b
*
            DO 110 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) / A(j,j), scaling x if necessary.
*
               XJ = ABS( X( J ) )
               IF( NOUNIT ) THEN
                  TJJS = A( J, J )*TSCAL
               ELSE
                  TJJS = TSCAL
                  IF( TSCAL.EQ.ONE )
     $               GO TO 100
               END IF
               TJJ = ABS( TJJS )
               IF( TJJ.GT.SMLNUM ) THEN
*
*                    abs(A(j,j)) > SMLNUM:
*
                  IF( TJJ.LT.ONE ) THEN
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by 1/b(j).
*
                        REC = ONE / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE IF( TJJ.GT.ZERO ) THEN
*
*                    0 < abs(A(j,j)) <= SMLNUM:
*
                  IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                       Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
*                       to avoid overflow when dividing by A(j,j).
*
                     REC = ( TJJ*BIGNUM ) / XJ
                     IF( CNORM( J ).GT.ONE ) THEN
*
*                          Scale by 1/CNORM(j) to avoid overflow when
*                          multiplying x(j) times column j.
*
                        REC = REC / CNORM( J )
                     END IF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
                  X( J ) = X( J ) / TJJS
                  XJ = ABS( X( J ) )
               ELSE
*
*                    A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                    scale = 0, and compute a solution to A*x = 0.
*
                  DO 90 I = 1, N
                     X( I ) = ZERO
   90             CONTINUE
                  X( J ) = ONE
                  XJ = ONE
                  SCALE = ZERO
                  XMAX = ZERO
               END IF
  100          CONTINUE
*
*              Scale x if necessary to avoid overflow when adding a
*              multiple of column j of A.
*
               IF( XJ.GT.ONE ) THEN
                  REC = ONE / XJ
                  IF( CNORM( J ).GT.( BIGNUM-XMAX )*REC ) THEN
*
*                    Scale x by 1/(2*abs(x(j))).
*
                     REC = REC*HALF
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                  END IF
               ELSE IF( XJ*CNORM( J ).GT.( BIGNUM-XMAX ) ) THEN
*
*                 Scale x by 1/2.
*
                  CALL DSCAL( N, HALF, X, 1 )
                  SCALE = SCALE*HALF
               END IF
*
               IF( UPPER ) THEN
                  IF( J.GT.1 ) THEN
*
*                    Compute the update
*                       x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
*
                     CALL DAXPY( J-1, -X( J )*TSCAL, A( 1, J ), 1, X,
     $                           1 )
                     I = IDAMAX( J-1, X, 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               ELSE
                  IF( J.LT.N ) THEN
*
*                    Compute the update
*                       x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
*
                     CALL DAXPY( N-J, -X( J )*TSCAL, A( J+1, J ), 1,
     $                           X( J+1 ), 1 )
                     I = J + IDAMAX( N-J, X( J+1 ), 1 )
                     XMAX = ABS( X( I ) )
                  END IF
               END IF
  110       CONTINUE
*
         ELSE
*
*           Solve A' * x = b
*
            DO 160 J = JFIRST, JLAST, JINC
*
*              Compute x(j) = b(j) - sum A(k,j)*x(k).
*                                    k<>j
*
               XJ = ABS( X( J ) )
               USCAL = TSCAL
               REC = ONE / MAX( XMAX, ONE )
               IF( CNORM( J ).GT.( BIGNUM-XJ )*REC ) THEN
*
*                 If x(j) could overflow, scale x by 1/(2*XMAX).
*
                  REC = REC*HALF
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                  END IF
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.ONE ) THEN
*
*                       Divide by A(j,j) when scaling x if A(j,j) > 1.
*
                     REC = MIN( ONE, REC*TJJ )
                     USCAL = USCAL / TJJS
                  END IF
                  IF( REC.LT.ONE ) THEN
                     CALL DSCAL( N, REC, X, 1 )
                     SCALE = SCALE*REC
                     XMAX = XMAX*REC
                  END IF
               END IF
*
               SUMJ = ZERO
               IF( USCAL.EQ.ONE ) THEN
*
*                 If the scaling needed for A in the dot product is 1,
*                 call DDOT to perform the dot product.
*
                  IF( UPPER ) THEN
                     SUMJ = DDOT( J-1, A( 1, J ), 1, X, 1 )
                  ELSE IF( J.LT.N ) THEN
                     SUMJ = DDOT( N-J, A( J+1, J ), 1, X( J+1 ), 1 )
                  END IF
               ELSE
*
*                 Otherwise, use in-line code for the dot product.
*
                  IF( UPPER ) THEN
                     DO 120 I = 1, J - 1
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  120                CONTINUE
                  ELSE IF( J.LT.N ) THEN
                     DO 130 I = J + 1, N
                        SUMJ = SUMJ + ( A( I, J )*USCAL )*X( I )
  130                CONTINUE
                  END IF
               END IF
*
               IF( USCAL.EQ.TSCAL ) THEN
*
*                 Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
*                 was not used to scale the dotproduct.
*
                  X( J ) = X( J ) - SUMJ
                  XJ = ABS( X( J ) )
                  IF( NOUNIT ) THEN
                     TJJS = A( J, J )*TSCAL
                  ELSE
                     TJJS = TSCAL
                     IF( TSCAL.EQ.ONE )
     $                  GO TO 150
                  END IF
*
*                    Compute x(j) = x(j) / A(j,j), scaling if necessary.
*
                  TJJ = ABS( TJJS )
                  IF( TJJ.GT.SMLNUM ) THEN
*
*                       abs(A(j,j)) > SMLNUM:
*
                     IF( TJJ.LT.ONE ) THEN
                        IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                             Scale X by 1/abs(x(j)).
*
                           REC = ONE / XJ
                           CALL DSCAL( N, REC, X, 1 )
                           SCALE = SCALE*REC
                           XMAX = XMAX*REC
                        END IF
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE IF( TJJ.GT.ZERO ) THEN
*
*                       0 < abs(A(j,j)) <= SMLNUM:
*
                     IF( XJ.GT.TJJ*BIGNUM ) THEN
*
*                          Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
*
                        REC = ( TJJ*BIGNUM ) / XJ
                        CALL DSCAL( N, REC, X, 1 )
                        SCALE = SCALE*REC
                        XMAX = XMAX*REC
                     END IF
                     X( J ) = X( J ) / TJJS
                  ELSE
*
*                       A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
*                       scale = 0, and compute a solution to A'*x = 0.
*
                     DO 140 I = 1, N
                        X( I ) = ZERO
  140                CONTINUE
                     X( J ) = ONE
                     SCALE = ZERO
                     XMAX = ZERO
                  END IF
  150             CONTINUE
               ELSE
*
*                 Compute x(j) := x(j) / A(j,j)  - sumj if the dot
*                 product has already been divided by 1/A(j,j).
*
                  X( J ) = X( J ) / TJJS - SUMJ
               END IF
               XMAX = MAX( XMAX, ABS( X( J ) ) )
  160       CONTINUE
         END IF
         SCALE = SCALE / TSCAL
      END IF
*
*     Scale the column norms by 1/TSCAL for return.
*
      IF( TSCAL.NE.ONE ) THEN
         CALL DSCAL( N, ONE / TSCAL, CNORM, 1 )
      END IF
*
      RETURN
*
*     End of DLATRS
*
      END
      SUBROUTINE DRSCL( N, SA, SX, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SA
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
*     ..
*
*  Purpose
*  =======
*
*  DRSCL multiplies an n-element real vector x by the real scalar 1/a.
*  This is done without overflow or underflow as long as
*  the final result x/a does not overflow or underflow.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of components of the vector x.
*
*  SA      (input) DOUBLE PRECISION
*          The scalar a which is used to divide each component of x.
*          SA must be >= 0, or the subroutine will divide by zero.
*
*  SX      (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-1)*abs(INCX))
*          The n-element vector x.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector SX.
*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      DOUBLE PRECISION   BIGNUM, CDEN, CDEN1, CNUM, CNUM1, MUL, SMLNUM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Initialize the denominator to SA and the numerator to 1.
*
      CDEN = SA
      CNUM = ONE
*
   10 CONTINUE
      CDEN1 = CDEN*SMLNUM
      CNUM1 = CNUM / BIGNUM
      IF( ABS( CDEN1 ).GT.ABS( CNUM ) .AND. CNUM.NE.ZERO ) THEN
*
*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
*
         MUL = SMLNUM
         DONE = .FALSE.
         CDEN = CDEN1
      ELSE IF( ABS( CNUM1 ).GT.ABS( CDEN ) ) THEN
*
*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
*
         MUL = BIGNUM
         DONE = .FALSE.
         CNUM = CNUM1
      ELSE
*
*        Multiply X by CNUM / CDEN and return.
*
         MUL = CNUM / CDEN
         DONE = .TRUE.
      END IF
*
*     Scale the vector X by MUL
*
      CALL DSCAL( N, MUL, SX, INCX )
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of DRSCL
*
      END
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
      SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSV  solves one of the systems of equations
*
*     A*x = b,   or   A'*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular matrix.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A'*x = b.
*
*              TRANS = 'C' or 'c'   A'*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := inv( A )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := inv( A' )*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )
     $               TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSV .
*
      END
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      SUBROUTINE DLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  Purpose
*  =======
*
*  DLABAD takes as input the values computed by DLAMCH for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of LARGE is sufficiently large.  This subroutine is intended to
*  identify machines with a large exponent range, such as the Crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by DLAMCH.  This subroutine is needed because
*  DLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  Arguments
*  =========
*
*  SMALL   (input/output) DOUBLE PRECISION
*          On entry, the underflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (input/output) DOUBLE PRECISION
*          On entry, the overflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
*
      RETURN
*
*     End of DLABAD
*
      END

            SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      DOUBLE PRECISION   TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the LAPACK routines.
*  It is called by an LAPACK routine if an input parameter has an
*  invalid value.  A message is printed and execution stops.
*
*  Installers may consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Arguments
*  =========
*
*  SRNAME  (input) CHARACTER*6
*          The name of the routine which called XERBLA.
*
*  INFO    (input) INTEGER
*          The position of the invalid parameter in the parameter list
*          of the calling routine.
*
*
      WRITE( *, FMT = 9999 )SRNAME, INFO
*
      STOP
*
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',
     $      'an illegal value' )
*
*     End of XERBLA
*
      END
            subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xGELSD and xGESDD)
*          =10: ieee NaN arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. External Functions ..
      INTEGER            IEEECK
      EXTERNAL           IEEECK
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800, 900, 1000,
     $        1100 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 32
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 32
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE 
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
  900 CONTINUE
*
*     ISPEC = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xGELSD and xGESDD)
*
      ILAENV = 25
      RETURN
*
 1000 CONTINUE
*
*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 0, 0.0, 1.0 ) 
      END IF
      RETURN
*
 1100 CONTINUE
*
*     ISPEC = 11: infinity arithmetic can be trusted not to trap
*
C     ILAENV = 0
      ILAENV = 1
      IF( ILAENV.EQ.1 ) THEN
         ILAENV = IEEECK( 1, 0.0, 1.0 ) 
      END IF
      RETURN
*
*     End of ILAENV
*
      END
      INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1998
*
*     .. Scalar Arguments ..
      INTEGER            ISPEC
      REAL               ONE, ZERO
*     ..
*
*  Purpose
*  =======
*
*  IEEECK is called from the ILAENV to verify that Infinity and
*  possibly NaN arithmetic is safe (i.e. will not trap).
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and NaN arithmetic.
*          = 0: Verify infinity arithmetic only.
*          = 1: Verify infinity and NaN arithmetic.
*
*  ZERO    (input) REAL
*          Must contain the value 0.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  ONE     (input) REAL
*          Must contain the value 1.0
*          This is passed to prevent the compiler from optimizing
*          away this code.
*
*  RETURN VALUE:  INTEGER
*          = 0:  Arithmetic failed to produce the correct answers
*          = 1:  Arithmetic produced the correct answers
*
*     .. Local Scalars ..
      REAL               NAN1, NAN2, NAN3, NAN4, NAN5, NAN6, NEGINF,
     $                   NEGZRO, NEWZRO, POSINF
*     ..
*     .. Executable Statements ..
      IEEECK = 1
*
      POSINF = ONE / ZERO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = -ONE / ZERO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGZRO = ONE / ( NEGINF+ONE )
      IF( NEGZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = ONE / NEGZRO
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEWZRO = NEGZRO + ZERO
      IF( NEWZRO.NE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = ONE / NEWZRO
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      NEGINF = NEGINF*POSINF
      IF( NEGINF.GE.ZERO ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      POSINF = POSINF*POSINF
      IF( POSINF.LE.ONE ) THEN
         IEEECK = 0
         RETURN
      END IF
*
*
*
*
*     Return if we were only asked to check infinity arithmetic
*
      IF( ISPEC.EQ.0 )
     $   RETURN
*
      NAN1 = POSINF + NEGINF
*
      NAN2 = POSINF / NEGINF
*
      NAN3 = POSINF / POSINF
*
      NAN4 = POSINF*ZERO
*
      NAN5 = NEGINF*NEGZRO
*
      NAN6 = NAN5*0.0
*
      IF( NAN1.EQ.NAN1 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN2.EQ.NAN2 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN3.EQ.NAN3 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN4.EQ.NAN4 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN5.EQ.NAN5 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      IF( NAN6.EQ.NAN6 ) THEN
         IEEECK = 0
         RETURN
      END IF
*
      RETURN
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
* =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END

