!**********************************************************************************************
!  File  snowdv22.f  
!  Driver program for the Utah Energy Balance
!  Snow Accumulation and Melt Model.
!  Version 3.0 (UEBGrid)  
!
!  Authors
!  David G. Tarboton (dtarb@usu.edu), Utah Water Research Laboratory, Utah State University
!  Avirup Sen Gupta, PhD student at USU 2012
!  Vinod Mahat, PhD student at USU completed 2011
!  Charlie Luce (cluce@fs.fed.us), USDA Forest Service, Rocky Mountain Research Station, Boise, ID
!  Jinsheng You, PhD student at USU completing 2004 
!  Tom Jackson, PhD student at USU completing, 1994 
!  Tanveer Chowdhury, MS student at USU completing 1993 
!
!  Last Change Sept 2012 during work to implement UEB for glacier snowmelt using NetCDF inputs and outputs (grid version)
!
!  Copyright (C) 2012  David Tarboton, Utah State University, dtarb@usu.edu.  http://www.engineering.usu.edu/dtarb
!
!  This file is part of UEB.
!
!    UEB is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    UEB is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    A copy of the GNU General Public License is included in the file gpl.txt.
!    This is also available at: http://www.gnu.org/licenses/.
!
!  If you wish to use or incorporate this program (or parts of it) into 
!  other software that does not meet the GNU General Public License 
!  conditions contact the author to request permission.
!  David G. Tarboton  
!  Utah State University 
!  8200 Old Main Hill 
!  Logan, UT 84322-8200 
!  USA 
!  http://www.engineering.usu.edu/dtarb/ 
!  email:  dtarb@usu.edu 
!
!**********************************************************************************************

      ! for netCDF declaration
      use netcdf
      PARAMETER(nsv=10,npar=32,nxv=6,niv=8,nov=14)
      integer :: OutCount,ioutv,totalgrid,numgrid
      integer, parameter :: NDIMS = 3
      integer:: irad, ireadalb, subtype, TotalModelTimeSteps, InputVTimeSteps(11), NumNCFILES(11),ISInputFROMNC(11)
      integer:: n, SVDT, reason, ilat, jlon, krec,dimlen2,dimlen1,dimlen3,MaxNumofFile, TimeSteps, ArrayElements
      integer ModelStartDate(3), ModelEndDate(3)
      REAL slope,azi,SingleArray(1), bca,bcc,subalb, ModelStartHour,ModelEndHour,Modeldt
      CHARACTER*200 svfile, vfile, Wfile, StateSiteFiles(16), xx, file_name,InputVName(11), InputTSFilename(11)
      CHARACTER*200, dimension(:,:), allocatable :: NCDFContainer
      !double precision, dimension(:), allocatable:: timeStepsnetCDF
      double precision, dimension(:,:,:), allocatable :: NCfiletimeDimesions
      REAL, dimension(:), allocatable :: OutputArray
      double precision, dimension(:), allocatable :: TempArray
      double precision, dimension(:,:),allocatable:: InputTime
      integer, dimension(:,:),allocatable:: OutPoint
      character*200, dimension(:),allocatable :: OutPointFiles
      Character (200):: WatershedVARID
      real, dimension(1):: WatershedGridValue
      integer:: NCVtimeStep(11), MaxInputVTimeStep, MaxNCVTimestep, timelength, FileOpenFlag(11)
      CHARACTER*200 MainHeading
      Character*50 varnameinncdf(11)
      CHARACTER*200 pfile,inputcon,StateSiteVName(32)
      CHARACTER* 200 Watershedfile,AggOutControl,AggdOutput
      CHARACTER*200 afile
      INTEGER YEAR,DAY,iflag(6)                        ! Pass a flag to control surface temperature modeling 
      REAL IO, LAT,Ts_last
      real INPT(niv,1)
      real sitev(nsv)
      real outv(nov)
      real statev(nxv)
      real param(npar)
      real dtbar(12)
      real mtime(4)                                     ! YJS pass to reflect the change of Snow (Year, month, Date, Hour)
      integer:: istep
      REAL:: IDNumber(1)
      REAL, dimension(:), allocatable :: DimValue1,DimValue2
      Double precision, allocatable :: TSV(:,:)
      integer, allocatable::StartEndNCDF(:,:)
      REAL, allocatable :: Allvalues(:,:)
      integer:: arrayx,NoofTS(11)
      integer::CURRENTARRAYPOS(11)
      real:: cumGM
!      REAL, dimension(:,:),Allocatable:: inputstorage
!      Double precision, dimension(:,:),Allocatable:: inputstorageJDT
      !CHANGES TO ACCOMODATE GLACIER
      real:: WGT ! WGT=WATER EQUIVALENT GLACIER THICKNESS
! Arrays to keep records of surface temperature and snowpace average 
! temperature. This is for the fourth model (Modified force restore approach) 

      real, allocatable :: tsprevday(:)
      real, allocatable :: taveprevday(:)
      Integer, dimension(:), allocatable :: NumTimeStepPerFile
      Character*200, dimension(:,:), allocatable :: OutputNCContainer
      Character*200, dimension(:,:), allocatable :: NCOutfileArr
      integer, dimension(:,:), allocatable :: NCIDARRAY
      REAL Values(1) !  one dimensional array to pass to netcdf
      Character*200:: FN,inputvarname(11) ! file name
! Add a control of ground heat input
      integer mQgOption
      integer:: m,t ! looping variable
      character (200) :: DimName1,DimName2
      character*(200) :: DimUnit1,DimUnit2
      Character* 200,dimension(:), allocatable :: OutFolder,outSampleFile,outputfolder
      integer,dimension(:),allocatable :: outvar  !  Variable to hold the index position of each output variable
      Character* 200:: OV
!  Arrays to manage the steping of time through the input files using the following rules
!  First value in a set of input files is the first time step value
!  If there are multiple netcdf file they are used in the order of their first time value
!  an input value persists until another value is encountered
      integer :: currentfile(11)
      integer :: currenttspos(11)  
      double precision :: timeofnextavailinput(11), FileNextDateJDT
      integer ::  fileofnextts(11)
      integer ::  posinfileofnextts(11),NumofFile,NumOutPoint
      integer, allocatable ::  NCfileNumtimesteps(:,:)   !  holds the number of time steps in each netcdf file
      double precision, allocatable :: NCStartTime(:,:)  ! holds the start time in each netcdf file
      real :: InpVals(11)
      integer :: CurrentNCFilePos(11),NextNCFilePos(11)
      integer :: NCFiletspos(11)
      integer :: FileCurrentDate(3,11),FileNextDate(3,11)
      double precision :: FileNextTime(11)
      double precision :: dhour, EJD, TJD,tol
      real :: FileNextVal(11)
      !CHARACTER(200), DIMENSION(63) :: outputfolder
      CHARACTER(200) :: OutControlFILE
      CHARACTER(200) ::Watershfile
      real,allocatable :: OutVarValue(:,:)
      REal:: OutArr(53),IniOutVals(16)
      integer:: timestepinfile,incfile
      INTEGER::  VARINDX(3),TIMEINDX(1)
      CHARACTER(200), DIMENSION(66) :: outSymbol
      logical towrite
      Double precision:: ReferenceHour,ReferenceTime,CTJD
      ! The start and count arrays will tell the netCDF library where to
      ! write our data.
      Double precision, allocatable:: FNDJDT(:)
      Double precision, allocatable:: TimeMaxPerFile(:,:)
      Double precision, allocatable:: TimeMinPerFile(:,:)
      REAL, allocatable::VarMissingValues(:,:),VarfILLValues(:,:)
      Double precision:: CurrentModelDT 
      ! FOR PROGRAM EXECUTION TIME 
      real, dimension(2) :: tarray
      real :: tresult, tio, tcomp,tout,tagg,tlast,toutnc
      REAL:: Hour
      Double precision:: SHOUR
      !aggregated output declaration
       integer:: AggOutNum,uniqueIDNumber,IDNum,ioutvar
       Character*200, dimension(:),Allocatable:: AggOutVar
       integer, dimension(:),Allocatable:: AggOutVarnum
       integer,dimension(:),allocatable:: UniqueIDArray
       REAL:: AggValues
       Real,dimension(:,:,:),Allocatable:: AggdWSVarVal ! Aggregated WaterShed variable Value
       integer, Allocatable:: yymmddarray(:,:)
       Real,dimension(:),Allocatable:: Timearray
       Character*200, dimension(:), allocatable :: AllNCDFfile
       integer:: totalNC
       !declares for check
       integer:: count,countinput,timerec
       integer:: AggUnit,ST
       REAL:: WsMissingValues,WsFillValues
      ! UTC time conversion declarations
       REAL::longitude,NHour,UTChour
       Integer:: MYEAR,MMONTH,MDAY
       Double precision::UTCJulDat,OHour,MHOUR,ModHour
       integer::modx,totaldayMOne,totalDay,StepInADay,modeldtInt,gg
       REAL:: TNCCreate,TInputCheck
      ! the symbol table element lengths have been expanded to match
      ! fixed width lengths.  This accomidates cross-compiler
      ! differences and conforms to later Fortran-2003 standards.
      outSymbol = (/ "ATF-BC   ","HRI      ","Eacl     ","Ema      ", &
         "Ta       ","P        ","V        ","RH       ","Qsi      ", &
         "Qli      ","Qnet     ","Cos      ","Ub       ","SWE      ", &
         "tausn    ","Prain    ","Psnow    ","Albedo   ","Qh       ", &
         "Qe       ","E        ","SWIT     ","Qm       ","Q        ", &
         "dMdt     ","Tave     ","Ts       ","CumP     ","CumE     ", &
         "CumMelt  ","NetRads  ","Smelt    ","refDep   ","totRefDep", &
         "Cf       ","Taufb    ","Taufd    ","Qsib     ","Qsid     ", &
         "Taub     ","Taud     ","Qsns     ","Qsnc     ","Qlns     ", &
         "Qlnc     ","Vz       ","Rkinsc   ","Rkinc    ","Inmax    ", &
         "int      ","ieff     ","Ur       ","SWEc     ","Tc       ", &
         "Tac      ","QHc      ","QEc      ","Ec       ","Qpc      ", &
         "Qmc      ","Mc       ","FMc      ","MassError","SWIGM    ", &
         "SWISM    ","SWIR     "/)
         
!  Constant data set
      data to /0.0/        !  Temperature of freezing (0 C)a
      data tk /273.15/     !  Temperature to convert C to K (273.15)
      data sbc /2.0747e-7/ !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data hf /333.5/      !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cw /4.18/       !  Water Heat Capacity (4.18 KJ/kg/C)
      data cs /2.09/       !  Ice heat capacity (2.09 KJ/kg/C)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data Rag /287.0/     !  Ideal Gas constant for dry air (287 J/kg/K) (name changed)
      data k  /0.4/        !  Von Karmans constant (0.4)
      data hff /3600.0/    !  Factor to convert /s into /hr (3600)
      data rhoi /917.0/    !  Density of Ice (917 kg/m^3)
      data rhow /1000.0/   !  Density of Water (1000 kg/m^3)
      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)
      data w1day /0.261799/!  Daily frequency (2pi/24 hr 0.261799 radians/hr) 
      data pi /3.141592654/!  Pi
      data Io/4914.0/      !  solar constant  Kj/m^2/hr
      tol=5./(60.*24.)     !  5 min tolerance
!  End of the constant data set declaration

!  Input Arguments
      narg=iargc()
      IF(narg .ge. 1)THEN
        CALL getarg(1,afile)
      ELSE
        WRITE(6,*) 'give filename of input files'
        READ(5,'(A200)') afile
      ENDIF


! Writing the Long file that includes all warnings
      OPEN(66,FILE='UEBWarning.log',STATUS='UNKNOWN')
      OPEN(636,FILE='UEBTiming.log',STATUS='UNKNOWN')
!  Open and Read Input File
      OPEN(1,FILE=afile,STATUS= 'OLD')
      READ(1,*)MainHeading
      READ(1,*)pfile,svfile,inputcon,OutControlFILE,Watershedfile,WatershedVARID,AggOutControl,AggdOutput
      CLOSE(1)

!  Read parameters
      CALL readvals(param,irad,ireadalb,bca,bcc,pfile)                 ! Parameter file

!     Flag to control radiation
!     0 is no measurements - radiation estimated from diurnal temperature range
!     1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
!     2 is incoming shortwave and longwave radiation read from file (measured)
!     3 is net radiation read from file (measured)

      
      !CALL nCDF2DArrayInfo(Watershedfile,dimlen2,dimlen1)
      CALL nCDF2DArrayInfo2(Watershedfile,dimlen2,dimlen1,WatershedVARID,WsMissingValues,WsFillValues)
!     Call function to go through input control file and find the maximum number of netcdf files for any time varying input variable     
      CAll InputMaxNCFiles(inputcon,MaxNumofFile,inputvarname,UTCOffSet)
      Allocate(NCDFContainer(MaxNumofFile,11))
      Allocate(NCfileNumtimesteps(MaxNumofFile,11))
      Allocate(VarMissingValues(MaxNumofFile,11))
      Allocate(VarfILLValues(MaxNumofFile,11))
!       Call function to read input files and determine the format of variable input (netcdf or text)
!       Outputs are 
!       IsInputFromNC is 0 for text and 1 for netcdf for corresponding variable
!       NumNCFiles is number of netcdf files.  Only filled for inputs from NetCDF
!       InputTSFilename contains file names for inputs that are from text files.  Only filled if var is text
!       NCDFContainer - Character array.  Rows (first dimension) contain file names of netCDF files for the variable.  
!       Columns correspond to each input variable.  Columns are only filled if var is netCDF

       CALL InputFiles(inputcon,MaxNumofFile,IsInputFromNC,NumNCFiles,InputTSFilename,NCDFContainer,&
        &ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,NCfileNumtimesteps,&
        &nrefyr,nrefmo,nrefday,varnameinncdf,arrayx,NoofTS,InpVals,VarMissingValues,VarfILLValues)
       Allocate(TimeMaxPerFile(MaxNumofFile,11))
       Allocate(TimeMinPerFile(MaxNumofFile,11))
       allocate(TSV(arrayx,11))
       allocate(AllValues(arrayx,11))
       CALL TimeSeriesAndTimeSteps(MaxNumofFile,NUMNCFILES,IsInputFromNC,InputTSFilename,NCDFContainer,&
        &ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,NCfileNumtimesteps,&
        &varnameinncdf,arrayx,NOofTS,TSV,Allvalues)
        
      ReferenceHour=0.00
      CALL JULDAT(nrefyr,nrefmo,nrefday,ReferenceHour,ReferenceTime)
      allocate(DimValue1(dimlen1))
      allocate(DimValue2(dimlen2))
      CALL SpatialCoordinate(Watershedfile,dimlen1,dimlen2,DimName1,DimName2,DimValue1,DimValue2,DimUnit1,DimUnit2)
      
      TInputCheck = ETIME(tarray)
      write(6,*)'Time to check inputs',TInputCheck
      write(636,*)'Check Inputs',TInputCheck, ' seconds'
! Output file creation
        CALL NumOutFiles(OutControlFILE,ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,&
                         &dimlen2,dimlen1,NumTimeStep,NumofFile,NumOutPoint,OutCount)                  
        Allocate(NumTimeStepPerFile(NumofFile)) ! This array will contains the number of time steps in the sequence of netcdf files
        Allocate(OutputNCContainer(NumofFile,OutCount))  !  This array will contain the file names for each output netcdf variable
        Allocate(NCOutfileArr(NumofFile,OutCount))  
        Allocate(OutPoint(NumOutPoint,2))
        Allocate(OutPointFiles(NumOutPoint))
        Allocate(outSampleFile(outcount))
        Allocate(OutVar(outcount))
        Allocate(outputfolder(outcount))
        Allocate(OutVarValue(NumTimeStep,66))
        CALL OutputFiles(OutControlFILE,NumTimeStep,Dimlen2,dimlen1,NumofFile,outSampleFile,NumTimeStepPerFile,OutVar,&
        &OutPoint,OutPointFiles,NumOutPoint,OutCount)
!  suggest concatenate outputfolder and outSamplefile into one string array
        Allocate(NCIDARRAY(NumofFile,outcount))
        Allocate(OutFolder(outcount))
        CALL DirectoryCreate(nrefyr,nrefmo,nrefday,dimlen1,dimlen2,DimName1,DimName2,DimValue1,DimValue2,DimUnit1,&
        &DimUnit2,NumofFile,outcount,Outvar,&
                        &NumTimeStepPerFile,outSampleFile,OutputNCContainer,NCIDARRAY)
        Allocate(StartEndNCDF(NumofFile,2))                
        CALL checks(svfile,MaxNumofFile,IsInputFromNC,NumNCFiles,totalNC,StateSiteVName)
        allocate(AllNCDFfile(totalNC))
        CALL  NCChecks(svfile,StateSiteVName,WatershedFile,MaxNumofFile,IsInputFromNC,NCDFContainer,NumNCFiles,totalNC,AllNCDFfile)
!   Work with aggregated outputs    
       CALL AggregatedOutNum(AggOutControl,outSymbol,AggOutNum)
       Allocate(AggOutVar(AggOutNum))
       Allocate(AggOutVarnum(AggOutNum))
       CALL  AggOutWSUniqueID(AggOutControl,outSymbol,AggOutNum,AggOutVar,Watershedfile,WatershedVARID,dimlen2,&
       &dimlen1,uniqueIDNumber,AggOutVarnum)
       Allocate(UniqueIDArray(uniqueIDNumber))   
       CALL WSUniqueArray(Watershedfile,WatershedVARID,dimlen1,dimlen2,uniqueIDNumber,UniqueIDArray)
       Allocate(AggdWSVarVal(NumTimeStep,uniqueIDNumber,AggOutNum))
       Allocate(yymmddarray(3,NumTimeStep))
       Allocate(Timearray(NumTimeStep))
       Allocate(FNDJDT(NumTimeStep))
       AggdWSVarVal=0
!  Suggest only use NCOutFileArr not OutputNCContainer                         
       AggUnit=887
       OPEN(AggUnit,FILE=AggdOutput,STATUS='unknown')
       write(Aggunit,*)'Year month day hour variable watershed value' ! write header
     
        St=0 
        ii=1
1056    If (ii .LE. NumofFile)THEN
            StartEndNCDF(ii,1)=ST+1
            StartEndNCDF(ii,2)=StartEndNCDF(ii,1)+NumTimeStepPerFile(ii)-1
            ST=StartEndNCDF(ii,2)
            ii=ii+1
            Go to 1056
        End if
        
       write(6,*)"Starting loop over grid cells"
      !  Initialize timing results
        TNCCreate= ETIME(tarray)
        write(6,*)'Time to create netCDFs ',(TNCCreate-TInputCheck)
        write(636,*)'Create netCDFs ',(TNCCreate-TInputCheck),' seconds'
        !write(6,*)'Initializing time tracking:',tarray(1),tarray(2),tresult
        tio=0.
        tcomp=0.
        tout=0.
        tagg=0.
        tlast=tresult
       totalgrid=dimlen1*dimlen2     
!      Start loop over space
       numgrid=0
       DO iycoord=1, dimlen1
       DO jxcoord=1, dimlen2

       iunit=119  !  unit for point output

       CALL nCDF2DRead(Watershedfile,WatershedVARID,IDNumber(1),jxcoord,iycoord)
        !  read site variables and initial conditions of state variables
       if((IDNumber(1) .ne. 0).or. (IDNumber(1) .ne. WsMissingValues) .or. (IDNumber(1) .ne. WsFillValues))then  !  Omit calculations if not in the watershed
         !  TODO change the above to also exclude netcdf no data values
       CALL readsv(param,statev,sitev,svfile,slope,azi,lat,subtype,dimlen2,dimlen1,iycoord,jxcoord,dtbar,Ts_last,longitude)
       CALL Values4VareachGrid(IsInputFromNC,MaxNumofFile,NUMNCFILES,NCDFContainer,varnameinncdf,iycoord,jxcoord,&
            &NCfileNumtimesteps,NOofTS,arrayx,TSV,Allvalues,VarMissingValues,VarfILLValues)

       StepInADay=24/modeldt     
       If (inputvarname(5) .NE. 'tmin')THEN
            NoofTS(5)=NoofTS(1)
            TSV(1:arrayx,5)=TSV(1:arrayx,1)
            modx=MOD(arrayx,StepInADay)
            if (modx .GT. 0)THEN
                totaldayMOne=(arrayx-modx)/StepInADay
                totalDDay=totalDayMOne+1
                Do gg=1,totaldayMOne
                    Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,5)=MINVAL(Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
                    Allvalues((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),5)&
                    &=MINVAL(Allvalues((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),1))
            Else
                totaldayMOne=(arrayx)/StepInADay
                totalDay=totalDayMOne
                Do gg=1,totalDay
                    Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,5)=MINVAL(Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
            End if
        END IF
        
        If (inputvarname(6) .NE. 'tmax')THEN
            NoofTS(6)=NoofTS(1)
            TSV(1:arrayx,6)=TSV(1:arrayx,1)
            modx=MOD(arrayx,StepInADay)
            if (modx .GT. 0)THEN
                totaldayMOne=(arrayx-modx)/StepInADay
                totalDDay=totalDayMOne+1
                Do gg=1,totaldayMOne
                    Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,6)=MINVAL(Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
                    Allvalues((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),6)=&
                    &MINVAL(Allvalues((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),1))
            Else
                totaldayMOne=(arrayx)/StepInADay
                totalDay=totalDayMOne
                Do gg=1,totalDay
                    Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,6)=MAXVAL(Allvalues(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
            End if
        END IF
        

                                        
!       OPEN(665,FILE='date.DAT',STATUS='UNKNOWN')
!       OPEN(666,FILE='Data.DAT',STATUS='UNKNOWN') 
!       OPEN(667,FILE='DataWrite.DAT',STATUS='UNKNOWN')
!       Do I = 1,arrayx
!            Write(665,37) TSV(i,1),TSV(i,2),TSV(i,3),TSV(i,4),TSV(i,5),TSV(i,6),&
!            &TSV(i,7),TSV(i,8),TSV(i,9),TSV(i,10),TSV(i,11)
!            Write(666,37) Allvalues(i,1),Allvalues(i,2),Allvalues(i,3),Allvalues(i,4),Allvalues(i,5),Allvalues(i,6),&
!            &Allvalues(i,7),Allvalues(i,8),Allvalues(i,9),Allvalues(i,10),Allvalues(i,11)
!37          format(f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5)
!38          format(i4,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5)
!       End do
!       Close(666)
!       Close(665)

!  Block of code to replicate variables from UEBVeg before time loop
       dt=Modeldt
      nstepday=24/dt

      ALLOCATE(Tsprevday(nstepday))
      ALLOCATE(Taveprevday(nstepday))
!  Initialize Tsbackup and TaveBackup
      DO 3 i = 1,nstepday
                Tsprevday(i)=-9999.
                Taveprevday(i)=-9999.0
 3    CONTINUE

      IF(ts_last .le. -9999.)THEN    
        Tsprevday(nstepday)=0  
      ELSE
        Tsprevday(nstepday)=ts_last                      !has measurements

      END IF 

!****************** To compute Ave.Temp (For Previous day temp.)  ******************************

          Us   = statev(1)                  ! Ub in UEB
      Ws   = statev(2)                          ! W in UEB
          Wc   = statev(4)
          Apr  = sitev(2)                                   ! Atm. Pressure  (PR in UEB)
      cg   = param(4)                                   ! Ground heat capacity (nominally 2.09 KJ/kg/C)
      rhog = param(8)                                   ! Soil Density (nominally 1700 kg/m^3)
      de   = param(11)                              ! Thermally active depth of soil (0.1 m)
      
      !  Glacier adjustment of ws
    IF(SITEV(10) .EQ. 0 .OR. SITEV(10) .EQ. 3)THEN
        WGT=0.0
    ELSE
        WGT=1.0
    END IF
    
          Tave = TAVG(Us,Ws+WGT,RHOW,CS,TO,RHOG,DE,CG,HF)        ! This call only
      Taveprevday(nstepday) = Tave
      
!   initialize variables for mass balance
      Ws1   = statev(2)
          Wc1   = statev(4)     
      cumP = 0.
      cumEs = 0.
      cumEc = 0.
      cumMr = 0.
      cumGM = 0.
!************************************************************************************************      

!  find out if this is an output point and if so open the point output file
      towrite=.false.
      do NumOP=1,NumOutPoint
       If (iycoord .eq. OutPoint(NumOP,1) .and. jxcoord .eq. OutPoint(NumOP,2))Then  
            OPEN(iUnit,FILE=OutPointFiles(NumOP),STATUS='unknown')
            towrite= .true.
            exit
        endif 
      enddo                     

!  Variables to keep track of which time step we are in and which netcdf output file we are in
        istep=0  ! time step
! map on to old UEB names
        YEAR=ModelStartDate(1)
        MONTH=ModelStartDate(2)
        DAY=ModelStartDate(3)
        Hour=ModelStartHour
        SHOUR=DBLE(ModelStartHour)
        call JULDAT(YEAR,MONTH,DAY,SHOUR,CurrentModelDT)
        dhour=dble(ModelEndHour)
        call JULDAT(ModelEndDate(1),ModelEndDate(2),ModelEndDate(3),dhour,EJD)
        !  This is the start of the main time loop      
        
1       istep=istep+1
        IF(sitev(10).NE. 3)THEN
        
        CALL InputVariableValue(INPUTVARNAME,IsInputFromNC,NoofTS,TSV,Allvalues,arrayx,ModelEndDate,&
        &ModelStartDate,ModelStartHour,ModelEndHour,&
        &nrefyr,nrefmo,nrefday,Year,Month,Day,Hour,ModelDt,InpVals,CurrentArrayPos,CurrentModelDT,istep)
         
!        Write(667,38)istep,InpVals(1),InpVals(2),InpVals(3),InpVals(4),InpVals(5),InpVals(6),&
!            &InpVals(7),InpVals(8),InpVals(9),InpVals(10),InpVals(11)

        tresult= ETIME(tarray)
        tio=tio+tresult-tlast
        tlast=tresult
        
 !      Map from wrapper input variables to UEB variables     
        TA=INPVals(1)
        P=INPVals(2)
        V=INPVals(3)
        RH=INPVals(4)
        Tmin=INPVals(5)
        Tmax=INPVals(6)
        trange=Tmax-Tmin
        QSIOBS=INPVals(7)
        qg=INPVals(8)
        qli=INPVals(9)
        QNETOB=INPVals(10)
        Snowalb=INPVals(11)
        
!  Below is code from point UEB 
        sitev(3)=qg   
        INPT(1,1)=TA
        INPT(2,1)=P
        INPT(3,1)=V
        INPT(4,1)=RH
        INPT(7,1)= QNETOB
        
!UTC to local time conversion
      UTCHour=Hour-UTCOffset
      NHour=UTCHour+longitude/15.0
      OHour=DBLE(NHour)
      CALL JulDat(YEAR,MONTH,DAY,OHour,UTCJulDat)
      CALL CALDat(UTCJulDat,MYEAR,MMONTH,MDAY,MHOUR)
      NHOUR=REAL(MHOUR)
!******************     Radiation Input Parameterization  ***************************************
       !CALL hyri(YEAR,MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT,HRI,COSZEN)
       CALL hyri(MYEAR,MMONTH,MDAY,NHOUR,DT,SLOPE,AZI,LAT,HRI,COSZEN)
       INPT(8,1)=COSZEN                
       IF(IRAD.LE.2)THEN               
            CALL atf(atff,trange,month,dtbar,bca,bcc)
            IF(IRAD.EQ.0) THEN                     
                INPT(5,1)=atff*IO*HRI
                          CALL cloud(param,atff,cf)   ! For cloudiness fraction
            ELSE 
               If(QSIOBS .lt. 0) then
                         write(66,*)"Warning! Negative incoming radiation: ",QSIOBS
                         write(66,*)"at date",year,month,day,hour
                         write(66,*)"was set to zero."
                 QSIOBS=0       
               Endif
!      Need to call HYRI for horizontal surface to perform horizontal
!      measurement adjustment
               !CALL hyri(YEAR,MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT,HRI,COSZEN)
               CALL hyri(MYEAR,MMONTH,MDAY,NHOUR,DT,0.0,AZI,LAT,HRI0,COSZEN)
!      If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
!      not it indicates a potential measurement problem. i.e. moonshine
                 if(HRI0 .GT. 0.0) then
                    atfimplied=min(qsiobs/(HRI0*IO),0.9) ! To avoid unreasonably large radiation when HRI0 is small
                    INPT(5,1)=atfimplied * HRI * IO
                 else
                   INPT(5,1)=QSIOBS
                     if(qsiobs .ne. 0.)then
                        write(66,*)"Warning ! you have nonzero nightime"
                        write(66,*)"incident radiation of",qsiobs
                        write(66,*)"at date",year,month,day,hour
                     endif
                 endif
               CALL cloud(param,atff,cf)   ! For cloudiness fraction  This is more theoretically correct
               !cf = 1.-atff/bca          ! Cf based on solar radiation measurement
            ENDIF
             
                  IF(irad .lt. 2)THEN    
                CALL qlif(TA,RH,TK,SBC,Ema,Eacl,cf,inpt(6,1) )
                  Else
                                inpt(6,1)=qli
          ENDIF 

          IRADFL=0                        ! Long wave or shortwave either measured and calculated

       ELSE 
              IRADFL=1                    ! This case is when given IRAD =3 (From Net Radiation)  
              INPT(7,1) = QNETOB          
       ENDIF
       
!************************************************************************************************
       if(towrite)WRITE(iunit,*)YEAR,MONTH,DAY,HOUR,atff,HRI,Eacl,Ema,(INPT(i,1),i=1,8) 
       IniOutVals(1)=YEAR
       IniOutVals(2)=MONTH
       IniOutVals(3)=DAY
       IniOutVals(4)=HOUR
       IniOutVals(5)=atff
       IniOutVals(6)=HRI
       IniOutVals(7)=Eacl
       IniOutVals(8)=Ema
       IniOutVals(9)=INPT(1,1)
       IniOutVals(10)=INPT(2,1)
       IniOutVals(11)=INPT(3,1)
       IniOutVals(12)=INPT(4,1)
       IniOutVals(13)=INPT(5,1)
       IniOutVals(14)=INPT(6,1)
       IniOutVals(15)=INPT(7,1)
       IniOutVals(16)=INPT(8,1)
       
!      set control flags
       iflag(1) = iradfl   ! radiation (0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7)
       !  In the code above radiation inputs were either computed or read from input files
       iflag(2) = 0        ! no 0 (/yes 1) printing
       iflag(3) = iunit        ! Output unit to which to print
       if(ireadalb .eq. 0)then
        iflag(4) = 1        ! Albedo Calculation (a value 1 means albedo is calculated, otherwise statev(3) is albedo
       else
        iflag(4)=0
        statev(3)=Snowalb
       endif 
!      iflag(5) = 4        ! yjs Surface temperature modeling method, read above from input file
!          iflag(5)  model option for surface temperature approximation (This is read as a model parameter stmflag
!              1) the Simple Gradient, almost the same as Original UEB,
!              2) Simple Gradient approach with shallow snow correction. 
!              3) The classical force-restore approach.
!              4) Modified force-restore approach.
!       iflag(6) = 2     !This is read as a model parameter forwsflag   ! Wind speed in forest canopy: 1 is observed and 2 is computed from above canopy observation
!      set modeling time
       mtime(1) = Year
       mtime(2) = month
       mtime(3) = day
       mtime(4) = hour

!**************************************************************************************************

        CALL SNOWUEB2(dt,1,inpt,sitev,statev,tsprevday, taveprevday, nstepday, param,iflag,&  
         &cump,cumes,cumEc,cummr,cumGM,outv,mtime,atff,cf,OutArr)
     
        tresult= ETIME(tarray)
        tcomp=tcomp+tresult-tlast
        tlast=tresult

        IF(sitev(10).EQ. 3)THEN  ! Substrate type is accumulation zone
           OutArr=0
        END IF
        
        if(towrite)WRITE(iunit,*)OutArr
        
        DStorage=statev(2)-Ws1+statev(4)-Wc1
        errmbal= cump-cumMr-cumEs-cumEc -DStorage+cumGM  
         
        IF(sitev(10).EQ. 3)THEN  ! Substrate type is accumulation zone
            ERRMBAL=0
        END IF
        if(towrite)WRITE(iunit,*)ERRMBAL
        
        !mapping to OutVarValue
        OutVarValue(istep,1:12)=IniOutVals(5:16)
        OutVarValue(istep,13:62)=OutArr(1:50)
        OutVarValue(istep,63)=ERRMBAL
        OutVarValue(istep,64:66)=OutArr(51:53)
        
        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(4) inside the loop below tells netCDF which
        ! timestep to write.)
    
        !  Here we rely on the even spread of time steps until the last file
         incfile=(istep-1)/NumTimeStepPerFile(1)+1  !  the netcdf file position
         timestepinfile=istep-(incfile-1)*NumTimeStepPerFile(1)  ! the time step in the file
         VARINDX = (/iycoord,jxcoord,timestepinfile/)
         TIMEINDX =(/timestepinfile/)
         ReferenceHour=DBLE(hour)
         call JULDAT(YEAR,MONTH,DAY,ReferenceHour,CTJD)  !  current julian date time
         FNDJDT(istep)=CTJD-ReferenceTime

     END IF
        

    IDNum=Int(IDNumber(1))
    do ioutvar=1,AggOutNum
        Do jUniqueID=1,uniqueIDNumber
            If(UniqueIDArray(jUniqueID) .eq. IDNum)then
                AggValues=OutVarValue(istep,AggOutVarnum(ioutvar))
                AggdWSVarVal(istep,jUniqueID,ioutvar)=AggdWSVarVal(istep,jUniqueID,ioutvar)+AggValues
            End if
        end do
    enddo

    yymmddarray(1,istep)=year
    yymmddarray(2,istep)=month
    yymmddarray(3,istep)=day
    Timearray(istep)=hour

    tresult= ETIME(tarray)
    tout=tout+tresult-tlast
    tlast=tresult

    CALL UPDATETIME(YEAR,MONTH,DAY,HOUR,DT)
!  End of time loop                     
!*************************************************************************************************   
        ModHour=DBLE(Hour)
        call JULDAT(YEAR,MONTH,DAY,ModHour,CurrentModelDT)
        
        If (EJD .GE. CurrentModelDT)Then      
            Go to 1
        End if
        deallocate(Tsprevday)
        deallocate(Taveprevday)
        Close(iunit)
        
!        close(667)
        

         do ioutv=1,outcount
           do incfile = 1,NumofFile
                CALL OutputnetCDF(NCIDARRAY,outvar,NumTimeStep,outcount,incfile,ioutv,jxcoord,iycoord,NumTimeStepPerFile,NumofFile,&
                &StartEndNCDF,OutVarValue)  
                CALL check(nf90_sync(NCIDARRAY(incfile,ioutv)))
          enddo
        enddo
    endif  !  this is the end of if we are in a watershed    
    tresult= ETIME(tarray)
    toutnc=toutnc+tresult-tlast
    tlast=tresult
    numgrid=numgrid+1
    write(6,FMT="(A1,A,t30,F6.2,A$)") achar(13), " Percent Grid completed: ", (real(numgrid)/real(totalgrid))*100.0, "%"
    
       END DO  !  These are the end of the space loop
     END DO
    
    do ioutv=1,outcount
       do incfile = 1,NumofFile
        CALL check(NF90_PUT_VAR(NCIDARRAY(incfile,ioutv),2,DimValue2))
        CALL check(NF90_PUT_VAR(NCIDARRAY(incfile,ioutv),3,DimValue1))
        CALL OutputTimenetCDF(NCIDARRAY,outvar,NumTimeStep,outcount,incfile,ioutv,NumTimeStepPerFile,NumofFile,StartEndNCDF,FNDJDT,jxcoord,iycoord)
        CALL check(nf90_sync(NCIDARRAY(incfile,ioutv)))
        CALL Check(nf90_close(NCIDARRAY(incfile,ioutv)))
      enddo
    enddo
    Write (6,FMT="(/A31/)") " now aggregation will be started" 
     do istep=1,NumTimestep
        do ivar=1,AggOutNum
            Do jUniqueID=1,uniqueIDNumber
                WRITE(Aggunit,47)yymmddarray(1,istep),yymmddarray(2,istep),yymmddarray(3,istep),Timearray(istep),outSymbol(AggOutVarnum(ivar)),&
                &UniqueIDArray(jUniqueID),AggdWSVarVal(istep,jUniqueID,ivar)
47             format(1x,i5,i3,i3,f8.3,1x,a11,i7,1x,g13.6)
            end do
        enddo
    end do
    Close(AggUnit)
    
    tresult= ETIME(tarray)
    tagg=tagg+tresult-tlast
    tlast=tresult

    tarray(1)=tarray(1)  !/60. convert seconds to minutes
    
    write(636,*) "Input time:",tio," Seconds"
    write(636,*) "Compute time:",tcomp," Seconds"
    write(636,*) "Out time:",tout," Seconds"
    write(636,*) "Out timeinNC:",toutnc," Seconds"
    write(636,*) "Aggregation time:",tagg," Seconds"
    write(636,*) "Complete runtime:",tarray(1)," Seconds"
    Close(636)
    Close(66)
    write(6,*) "Input time:",tio," Seconds"
    write(6,*) "Compute time:",tcomp," Seconds"
    write(6,*) "Out time:",tout," Seconds"
    write(6,*) "Out timeinNC:",toutnc," Seconds"
    write(6,*) "Aggregation time:",tagg," Seconds"
    write(6,*) "Complete runtime:",tarray(1)," Seconds"
    Write(6,*) "Your task is successfully performed! Plesae view the results in 'outputs' folder!"
end program
