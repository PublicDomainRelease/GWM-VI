PROGRAM MMPROC
! MODFLOW Management Pre- and Post-Processor.
! MMProc is designed to be used with GWM-VI, a version- 
! independent MODFLOW groundwater-management application.
! MMProc supports the following groundwater-modeling
! programs based on MODFLOW:
!   MODFLOW-2005
!   SEAWAT
!   CONDUIT FLOW PROCESS
!   MODFLOW-NWT 
!   FARM PROCESS
! MMProc may be usable with other MODFLOW-based programs 
! as well.
!
USE MMPROCMOD, ONLY: COMMAND, IFREFM, NCOL, NROW, NLAY, IOUT, &
                     NPER, NSTP, WELFILE, &
                     INITIALIZE_MMPMOD, MMP_INI, MMP_ALLOC, &
                     READ_DEPENDENTS, MMP_READ_DRN_FILE, &
                     READ_FILENAMES, READ_INPUT_1, &
                     MMP_READ_NAME_FILE, READ_WEL_DATA, &
                     SAVE_WEL_FILE, &
                     WRITE_WEL_FILE, RESTORE_WEL_FILE, &
                     EXTRACT_DRN_DEPS, EXTRACT_HEADS, &
                     EXTRACT_BUD_GRID_DEPS, EXTRACT_SFR_DEPS, &
                     WRITE_DEPENDENT_DATA, WRITE_STATUS, &
                     SFR_IDENTIFY, MMP_READ_MNW2_FILE, CHECK_MNW_QNET                
USE MMPROCMOD, ONLY: NDRNFDEP, NDRNVDEP, NHEDDEP, NMQWEL, &
                     NSFRFDEP, NSFRLDEP, NSTORDEP, DRNFILE, &
                     MNW2FILE
USE MMPWELMOD, ONLY: WELTYPE, WELTYPE_COMBINE, WELTYPE_DEALLOC, &
                     READ_WEL_FILE
USE MMPMNWMOD, ONLY: MNW2WELL, MNW2WELL_ALLOCATE, MNW2WELL_DEALLOC, &
                     READ_MNW2_DATA
USE GWM_UTLS, ONLY: APPEND_MESSAGE, OPEN_FILE
USE GWM_SUBS, ONLY: EPSQNET
USE UTILITIES, ONLY: UTL_STOP, UTL_SYSTEM_FUNC, UTL_WRITE_MESSAGE
USE MGT, ONLY: MGT_INI, MGT_SETSTATUS, MGT_WRITESTATUS
IMPLICIT NONE
! Local variables
INTEGER :: IUMAIN, NMQMNW
CHARACTER(LEN=20) :: MMPROCIN = 'MMProc.in', MMPROCOUT = 'MMProc.out'
CHARACTER(LEN=100) :: ERRMSG
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DHEADS
TYPE (WELTYPE) :: MANAGED_WELLS, ORIG_WELLS, ALL_WELLS
TYPE (MNW2WELL), ALLOCATABLE, DIMENSION(:) :: MNW2WELLS
LOGICAL :: OK, HEADFILEEMPTY
! Format statements
10 FORMAT(' Message output for MMProc',/)
!
IOUT = 7
IUMAIN = 8
NMQMNW = 0
HEADFILEEMPTY = .FALSE.
!
CALL MMP_INI()
CALL INITIALIZE_MMPMOD()
!
! Open output file for messages
OK = OPEN_FILE(IOUT,MMPROCOUT,'REPLACE')
IF (.NOT. OK) GOTO 9000
WRITE(IOUT,10)
!
! Open input file MMProc.in.  Data are read and stored as needed.
OK = OPEN_FILE(IUMAIN,MMPROCIN,'OLD')
IF (.NOT. OK) GOTO 9000
!
! Read items 1-14
CALL READ_INPUT_1(IUMAIN,ORIG_WELLS,MANAGED_WELLS,NMQMNW,OK)
IF (.NOT. OK) GOTO 9000
!
CALL MMP_READ_NAME_FILE(OK)
IF (.NOT. OK) GOTO 9000
!
IF (DRNFILE.NE.'') THEN
  ! Find IDRNCB, in case it's needed
  CALL MMP_READ_DRN_FILE(OK)
  ! Don't care at this point if reading of DRNFILE succeeds,
  ! because IDRNCB may not be needed.
  OK = .TRUE. 
ENDIF
!
IF (MNW2FILE.NE.'') THEN
  ! Find IMNWCB, in case it's needed
  CALL MMP_READ_MNW2_FILE(OK)
  ! Don't care at this point if reading of MNW2FILE succeeds,
  ! because IMNWCB may not be needed.
  OK = .TRUE. 
ENDIF
!
IF (NMQMNW > 0) THEN
  ! Read and store data for MNW2 wells (items 15)
  ALLOCATE(MNW2WELLS(NMQMNW))
  CALL READ_MNW2_DATA(IUMAIN,NMQMNW,NPER,MNW2WELLS,OK)
  IF (.NOT. OK) GOTO 9000
ENDIF
!
! Read items 17-18
CALL READ_DEPENDENTS(IUMAIN,OK)
IF (.NOT. OK) GOTO 9000
!
! Read items 19-20
CALL READ_FILENAMES(IUMAIN,OK)
IF (.NOT. OK) GOTO 9000
!
! Read item 21 (EPSQNET)
READ(IUMAIN,*,ERR=9000)EPSQNET
!
! If there are and SFR dependents, read SFR input and populate DEPS()%INDX(4) with 
! SEQNUM (index in list of reaches as written to binary budget file)
IF (NSFRFDEP+NSFRLDEP > 0) THEN
  CALL SFR_IDENTIFY(OK)
ENDIF
!
IF (NMQWEL > 0) THEN
  ! Revise WEL package data by adding managed wells to wells read from WEL file.
  CALL WELTYPE_COMBINE(ORIG_WELLS,MANAGED_WELLS,ALL_WELLS)
  ! Save original WEL file with a new file name
  OK = SAVE_WEL_FILE()
  IF (.NOT. OK) GOTO 9000
  ! If NMQWEL > 0, rewrite WEL Package input file.
  CALL WRITE_WEL_FILE(ALL_WELLS,OK)
  IF (.NOT. OK) GOTO 9000
  ! Any errors below this point should cause GOTO 8900 (or 8800???)
ENDIF
!  
! Deallocate unneeded memory before starting model.
IF (NMQWEL>0) THEN
  CALL WELTYPE_DEALLOC(ORIG_WELLS)
  CALL WELTYPE_DEALLOC(ALL_WELLS)
ENDIF
!
! for debugging
!write(*,*)'press Enter to continue mmproc'
!read(*,*)
! Execute model by invoking Command.
CALL MGT_INI()
OK = UTL_SYSTEM_FUNC(COMMAND,ERRMSG)
IF (.NOT. OK) THEN
  CALL APPEND_MESSAGE('Operating system returned error: ' // TRIM(ERRMSG))
  GOTO 8900
ENDIF
!
! Allocate arrays of MMPROCMOD module required for post processing
CALL MMP_ALLOC()
!  
! Extract head data for dependents from binary head file.
ALLOCATE(DHEADS(NCOL,NROW))
IF (NHEDDEP>0 .OR. MANAGED_WELLS%MXNPW>0 .OR. NMQMNW>0) THEN
  CALL EXTRACT_HEADS(NCOL,NROW,DHEADS,MANAGED_WELLS,NMQMNW,MNW2WELLS, &
                     OK,HEADFILEEMPTY)
  IF (.NOT. OK) GOTO 8800
  IF (HEADFILEEMPTY) GOTO 8800
ENDIF
!  
! Extract storage-change data for dependents from binary compact budget file.
IF (NSTORDEP>0) THEN
  CALL EXTRACT_BUD_GRID_DEPS(OK)
  IF (.NOT. OK) GOTO 8900
ENDIF
!  
! Extract SFR flow and sw-gw exchange data for dependents from binary compact 
! budget file(s).
IF (NSFRFDEP+NSFRLDEP>0) THEN
  CALL EXTRACT_SFR_DEPS(OK)
  IF (.NOT. OK) GOTO 8900
ENDIF
!  
! Extract drain data for flow and volume dependents from binary compact 
! budget file.
IF (NDRNFDEP+NDRNVDEP>0) THEN
  CALL EXTRACT_DRN_DEPS(OK)
  IF (.NOT. OK) GOTO 8900
ENDIF
!
! Check binary file containing MNW2 data to compare
! Qnet with Qdes for all managed MNW2 wells
IF (NMQMNW>0) THEN
  CALL CHECK_MNW_QNET(NMQMNW,MNW2WELLS,OK)
  IF (.NOT. OK) GOTO 8900
ENDIF
!  
8800 CONTINUE
! Write output to files SimulatedValues.out and modflow.status.
CALL WRITE_DEPENDENT_DATA(OK)
!IF (.NOT. OK) GOTO 8900
CALL WRITE_STATUS(MANAGED_WELLS,NMQMNW,MNW2WELLS,HEADFILEEMPTY,OK,NPER,NSTP)
CALL WELTYPE_DEALLOC(MANAGED_WELLS)
!
! Restore all modified MODFLOW input files
8900 CONTINUE
IF (NMQWEL > 0) THEN
  ! Restore original WEL file using original name
  CALL RESTORE_WEL_FILE(OK)
ENDIF
!
9000 CONTINUE
!
IF (OK) THEN
  IF (HEADFILEEMPTY) THEN
    WRITE(*,*)'MMProc terminates with warning: Head file empty'
  ELSE
    WRITE(*,*)'Normal termination of MMProc'
  ENDIF
ELSE
  CALL UTL_WRITE_MESSAGE(IOUT)
  WRITE(*,*)'Error encountered in MMProc'
ENDIF
CALL UTL_STOP()
END PROGRAM MMPROC
