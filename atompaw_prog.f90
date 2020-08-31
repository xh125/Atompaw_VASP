PROGRAM atompaw
  !*************************************************************
  !  Driver program for new version of atompaw program
  !      arguments can be SAVEAEATOM (save file name),  OR
  !                       LOADAEATOM (load file name)
  !      Note:  must have 0 or 2 arguments
  !************************************************************
  USE GlobalMath
  USE aeatom
  USE atomdata
  USE atompaw_report
  USE pseudo
  USE abinitinterface
  USE pwscfinterface
  USE xmlinterface
  USE libxc_mod
  USE pkginfo

  IMPLICIT NONE
  INTEGER :: i,j,iargc
  CHARACTER(120) :: inputfile,outputfile,token
  LOGICAL :: lotsofoutput=.false.    ! can be changed eventually
  LOGICAL :: saveaeatom=.false.,loadaeatom=.false.
  LOGICAL :: OK
  INTEGER :: ifinput=10,ifen=11

! First write out details on AtomPAW
  WRITE(6,*) atp_package,' v',atp_version
  WRITE(6,*) 'Compiled for ',atp_target
  WRITE(6,*)

  j=iargc()
  if (j==2) then
       call GetArg(1,token)
       call UpperCase(token)
       if (TRIM(token)=='SAVEAEATOM') then
           saveaeatom=.true.
           call GetArg(2,outputfile)
           write(6,*) 'AEATOM results will be saved to file ',TRIM(outputfile)
       else if (TRIM(token)=='LOADAEATOM') then
           loadaeatom=.true.
           call GetArg(2,inputfile)
           write(6,*) 'AEATOM results will be loaded from file ',TRIM(inputfile)
           Inquire(file=TRIM(inputfile),exist=OK)
           If (.not.OK) then
              write(6,*) 'Load file does not exist -- program will stop'
              stop
           endif
       else
          write(6,*) 'Argument form not recognized', token
          stop
       endif

  else if (j==0) then
       write(6,*) 'Input/output not saved/dumped in this run'
  else
       write(6,'(3a)') 'Argument form not recognized (', j,') !'
       stop
  endif
  write(6,*)

  exploremode=.false.

  OPEN(ifinput,file='dummy',form='formatted')

  CALL Init_GlobalConstants()

!  WRITE(6,*) 'loadaeatom', loadaeatom
  If (loadaeatom) then
     Call Load_AEatom(TRIM(inputfile),&
        Grid,AEOrbit,AEpot,AESCF,FCOrbit,FCpot,FCSCF,FC,ifinput)
     OPEN (ifen, file=TRIM(AEPot%sym), form='formatted')
     CALL Report_Atomres('AE',Grid,AEOrbit,AEPot,AESCF,ifen)
     CALL Report_Atomres('SC',Grid,FCOrbit,FCPot,FCSCF,ifen)
  else
     CALL SCFatom_Init(ifinput)
!     CALL SCFatom('AE',lotsofoutput)
!     OPEN (ifen, file=TRIM(AEPot%sym), form='formatted')
!     CALL Report_Atomres('AE',Grid,AEOrbit,AEPot,AESCF,ifen)
!     CALL SCFatom('SC',lotsofoutput,ifinput)
!     CALL Report_Atomres('SC',Grid,FCOrbit,FCPot,FCSCF,ifen)
!     if (saveaeatom) then
!        Call Dump_AEatom(TRIM(outputfile),&
!          Grid,AEOrbit,AEpot,AESCF,FCOrbit,FCpot,FCSCF,FC)
!     endif
  endif

END PROGRAM atompaw
