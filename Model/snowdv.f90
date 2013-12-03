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
      Implicit None
      integer,PARAMETER:: LinesGridmodel=8  
      logical gridmodel
      ! Pointmodel = running UEB as a point model would have 4 lines in its control file
      ! Gridmodel = running UEB as a point model would have 8 lines in its control file
      integer:: snowdgtvariteflag
      integer::nsv,npar,nxv,niv,nov
      real::to,tk,sbc,hf,hneu,cw,cs,cp,rag,k,hff,rhoi,rhow,pi,narg,UTCOFFSET,Etime,g,w1day,dt
      integer::NUMtimeSTEP,NREFYR,NREFMO,NREFDAY,NumOP
      REAL::us,ws,wc,apr,cg,rhog,de,tave,ws1,wc1,cump,cumes,cumec,cummr,ta,p,v,rh,tmin,tmax,TRANGE,QSIOBS,QG,QLI
      real::tavg,qnetob,snowalb,coszen,atff,cf,hri,atfimplied,ema,eacl,dstorage,ERRMBAL,HRI0
      integer:: juniqueid,ivar,iradfl
      integer::month
      integer::ii,iycoord,jxcoord,iunit,totaldday,NSTEPDAY,i
      PARAMETER(nsv=10,npar=32,nxv=6,niv=8,nov=14)
      integer :: OutCount,ioutv,totalgrid,numgrid
      integer, parameter :: NDIMS = 3
      integer:: irad, ireadalb, subtype, NumNCFILES(11),ISInputFROMNC(11)
      integer:: dimlen2,dimlen1,MaxNumofFile
      integer ModelStartDate(3), ModelEndDate(3)
      REAL slope,azi, bca,bcc, ModelStartHour,ModelEndHour,Modeldt
      CHARACTER*200 svfile, InputTSFilename(11)
      CHARACTER*200, dimension(:,:), allocatable :: NCDFContainer
      integer, dimension(:,:),allocatable:: OutPoint
      character*200, dimension(:),allocatable :: OutPointFiles
      Character (200):: WatershedVARID
      CHARACTER*200 MainHeading
      Character*100 varnameinncdf(11)
      CHARACTER*200 pfile,inputcon,StateSiteVName(32)
      CHARACTER*200 Watershedfile,AggOutControl,AggdOutput
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
      Integer:: IDNumber(1)
      !REAL*8, dimension(:), allocatable :: DimValue1,DimValue2
      Double precision, allocatable :: TSV(:,:)
      integer, allocatable::StartEndNCDF(:,:)
      REAL, allocatable :: Allvalues(:,:)
      integer:: maxncfilents,NoofTS(11)
      integer::CURRENTARRAYPOS(11)
      real:: cumGM
      integer:: iycoord1
      CHARACTER(200) :: strwatershed
      CHARACTER(1) :: delimit1,delimit2
      integer:: nargs,nargs2,nargs3,nargs4,nargs5
      character(200),Allocatable:: wordsw(:)
      Character(200):: StateSiteFilesR, SitexcoordinateR, SiteycoordinateR, InputtcoordinateR
      Character(100):: VarNameinNCDFR
      integer:: ncidout, dimid1,dimid2
!      REAL, dimension(:,:),Allocatable:: inputstorage
!      Double precision, dimension(:,:),Allocatable:: inputstorageJDT
!      CHANGES TO ACCOMODATE GLACIER
      real:: WGT ! WGT=WATER EQUIVALENT GLACIER THICKNESS
      Double precision:: DBHour
      character(50)::daysstring
! Arrays to keep records of surface temperature and snowpace average 
! temperature. This is for the fourth model (Modified force restore approach) 
      real, allocatable :: tsprevday(:)
      real, allocatable :: taveprevday(:)
      Double precision, allocatable:: modelTimeJDT(:)
      integer, allocatable :: CurrentArrayPosRegrid(:,:)
      REAL, allocatable:: ReGriddedArray(:,:)
      Integer, dimension(:), allocatable :: NumtimeStepPerFile
      Character(200), dimension(:,:), allocatable :: OutputNCContainer
      Character(200), dimension(:,:), allocatable :: NCOutfileArr
      integer, dimension(:,:), allocatable :: NCIDARRAY
      Character(200):: inputvarname(11) ! file name
! Add a control of ground heat input
      character (200) :: DimName1,DimName2
      character*(200) :: DimUnit1,DimUnit2
      Character* 200,dimension(:), allocatable :: OutFolder,outSampleFile,outputfolder
      integer,dimension(:),allocatable :: outvar  !  Variable to hold the index position of each output variable

!  Arrays to manage the steping of time through the input files using the following rules
!  First value in a set of input files is the first time step value
!  If there are multiple netcdf file they are used in the order of their first time value
!  an input value persists until another value is encountered
      integer ::  NumofFile,NumOutPoint,ncidout1
      integer, allocatable ::  NCfileNumtimesteps(:,:)   !  holds the number of time steps in each netcdf file
      real :: InpVals(11)
      double precision :: dhour, EJD,tol
      CHARACTER(200) :: OutControlFILE
      real,allocatable :: OutVarValue(:,:)
      REal:: OutArr(53),IniOutVals(16)
      integer:: timestepinfile,incfile
      INTEGER::  VARINDX(3),timeINDX(1)
      CHARACTER(200), DIMENSION(66) :: outSymbol
      logical towrite
      Double precision:: ReferenceHour,Referencetime,CTJD
! The start and count arrays will tell the netCDF library where to
! write our data.
      Double precision, allocatable:: FNDJDT(:)
      Double precision, allocatable:: timeMaxPerFile(:,:)
      Double precision, allocatable:: timeMinPerFile(:,:)
      REAL, allocatable::VarMissingValues(:,:),VarfILLValues(:,:)
      Double precision:: CurrentModelDT 
! FOR PROGRAM EXECUTION time 
      real, dimension(2) :: tarray
      real :: tresult, tio, tcomp,tout,tagg,tlast,toutnc,taggre, tsitev
      REAL:: Hour
      Double precision:: SHOUR
!aggregated output declaration
       integer:: AggOutNum,uniqueIDNumber,IDNum,ioutvar
       Character*200, dimension(:),Allocatable:: AggOutVar
       integer, dimension(:),Allocatable:: AggOutVarnum
       integer,dimension(:),allocatable:: UniqueIDArray,UniqueIDArrayCount
       REAL:: AggValues
       Real,dimension(:,:,:),Allocatable:: AggdWSVarVal       ! Aggregated WaterShed variable Value
       Real,dimension(:,:,:),Allocatable:: AggdWSVarValAvg   ! Average over the grid of Aggregated values
       integer, Allocatable:: yymmddarray(:,:)
       Real,dimension(:),Allocatable:: timearray
       Character*200, dimension(:), allocatable :: AllNCDFfile
       integer:: totalNC
! declares for check
       integer:: AggUnit,ST
       REAL:: WsMissingValues,WsFillValues
! UTC time conversion declarations
       REAL::longitude,NHour,UTChour
       Integer:: MYEAR,MMONTH,MDAY
       Double precision::UTCJulDat,OHour,MHOUR,ModHour
       integer::modx,totaldayMOne,totalDay,StepInADay,gg
       REAL:: DimDiff1,DimDiff2
       Character*50 Inputxcoordinates(11),inputycoordinates(11),inputtcoordinates(11),Sitexcoordinates(32), Siteycoordinates(32)
       character*50 wsxcoord,wsycoord,wsxcoordinate,wsycoordinate
       integer:: nlines 
       Double precision:: MSTARTHOUR,MENDHOUR
       Double precision:: JMSD,JMED
       integer:: dimlen1OUT,dimlen2OUT
       REAL, allocatable:: InputvarMissVal(:,:), InputvarFillVal(:,:)
       !integer::WsMissingValuesInt,WsFillValuesInt
       Character (50):: DefaultDimNames(3)
       integer:: DefaultDimValues(3)
       REAL::RangeMin,RangeMax
       CHARACTER(1):: delimit3
       REAL:: InputVarRange(11,2),CumEG
       character(200), allocatable::Words7element(:)

       DefaultDimNames(1)='time'
       DefaultDimNames(2)='Y'
       DefaultDimNames(3)='X' 
! snowdgtvariteflag=1 means write all the warning messages
! snowdgtvariteflag=0 means do not write all the warning messages
       snowdgtvariteflag=0   ! TODO  add reading this as an input argument either on the command line or input file
       delimit1=';'
       delimit2=':'
       delimit3=','
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
     
! Opening the log files 
      OPEN(66,FILE='UEBWarning.log',STATUS='UNKNOWN')
      OPEN(636,FILE='UEBTiming.log',STATUS='UNKNOWN')
!  Open and Read Input File
      nlines = 0
      OPEN (875,file = afile,STATUS= 'OLD')
      DO
        READ (875,*, END=786)
        nlines = nlines + 1
      END DO
786   CLOSE (875)
      gridmodel=.true.
      if (nlines < linesGridModel)gridmodel = .false.
 
      if (GridModel)THEN
        OPEN(1,FILE=afile,STATUS= 'OLD')
        READ(1,*)MainHeading
        READ(1,*)pfile,svfile,inputcon,OutControlFILE
        READ(1,*)strwatershed
        CALL StringSep(strwatershed,delimit1,nargs)
        Allocate(wordsw(nargs))
        Allocate(Words7element(7))
        CALL StringSepWord(strwatershed,delimit1,nargs,wordsw)
        Words7element(1:nargs)=wordsw(1:nargs)
        CALL StringToVarName(nargs,Words7element,delimit2,StateSiteFilesR,SitexcoordinateR,SiteycoordinateR,VarNameinNCDFR,&
                             &InputtcoordinateR,DefaultDimValues,RangeMin,RangeMax,delimit3,.true.)
        deAllocate(wordsw)
        deAllocate(Words7element)
        Watershedfile=StateSiteFilesR
        ! DefaultDimValues(1:3) refers to time, y, x always
        ! DefaultDimValues is flagged as -9999 use the name of the dimension read from the input file and returned by StringToVarName
        If (DefaultDimValues(3) == -9999)THEN 
             wsxcoordinate=SitexcoordinateR
        else 
          CALL check(nf90_open(Watershedfile,NF90_NOWRITE, ncidout1))
          CALL check(nf90_inquire_dimension(ncidout1,DefaultDimValues(3),wsxcoordinate))  ! this reads wsxcoordinate dimension name from the file.  
          !  There is an assumption here valid for watershed files that the 2nd dimension is x
          CALL check(nf90_close(ncidout1))
        END IF
        If (DefaultDimValues(2) == -9999)THEN ! use the name of the dimension read from the input file and returned by StringToVarName
          wsycoordinate=SiteycoordinateR
        else
          CALL check(nf90_open(Watershedfile,NF90_NOWRITE, ncidout1))
          CALL check(nf90_inquire_dimension(NCIDout1,DefaultDimValues(2),wsycoordinate))
          CALL check(nf90_close(ncidout1))
        End if 
     
        WatershedVARID=VarNameinNCDFR
        read(1,*)AggOutControl,AggdOutput
        CLOSE(1)
      else
        OPEN(1,FILE=afile,STATUS= 'OLD')
        READ(1,*)MainHeading
        READ(1,*)pfile,svfile,inputcon
        CLOSE(1)
        dimlen1=1
        dimlen2=1
      END IF
      
!  Read parameters
      CALL readvals(param,irad,ireadalb,bca,bcc,pfile)                 ! pfile = Parameter file 
      
      if (GridModel)THEN
        CALL readsvcoordinate(svfile,Sitexcoordinates,Siteycoordinates)
  !     Flag to control radiation (irad)
  !     0 is no measurements - radiation estimated from diurnal temperature range
  !     1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
  !     2 is incoming shortwave and longwave radiation read from file (measured)
  !     3 is net radiation read from file (measured)
  !     Flag to control albedo (ireadalb)
  !     0 is no measurements - albedo estimated internally
  !     1 is net radiation read from file (provided: measured or obrained from another model)
      
        CALL nCDF2DArrayInfo2(Watershedfile,wsxcoordinate,wsycoordinate,dimlen2,dimlen1,WatershedVARID,WsMissingValues,WsFillValues)
      end if
      
      CAll InputMaxNCFiles(inputcon,MaxNumofFile,inputvarname,UTCOffSet)      
      Allocate(NCDFContainer(MaxNumofFile,11))
      Allocate(NCfileNumtimesteps(MaxNumofFile,11))
      Allocate(VarMissingValues(MaxNumofFile,11))
      Allocate(VarfILLValues(MaxNumofFile,11))
      Allocate(InputvarMissVal(MaxNumofFile,11))
      Allocate(InputvarFillVal(MaxNumofFile,11))
!       Call function to read input files and determine the format of variable input (netcdf or text)
!       Outputs are 
!       IsInputFromNC is 0 for text and 1 for netcdf for corresponding variable
!       NumNCFiles is number of netcdf files.  Only filled for inputs from NetCDF
!       InputTSFilename contains file names for inputs that are from text files.  Only filled if var is text
!       NCDFContainer - Character array.  Rows (first dimension) contain file names of netCDF files for the variable.  
!       Columns correspond to each input variable.  Columns are only filled if var is netCDF

      CALL InputFiles(inputcon,MaxNumofFile,IsInputFromNC,NumNCFiles,InputTSFilename,NCDFContainer,&
      &ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,NCfileNumtimesteps,&
      &nrefyr,nrefmo,nrefday,varnameinncdf,maxncfilents,NoofTS,InpVals,VarMissingValues,VarfILLValues,&
      &Inputxcoordinates,inputycoordinates,inputtcoordinates,InputVarRange,daysstring)
!  FIXME: what if the result is fractional
!  time steps must divide exactly in to a day because we use logic that requires the values from the same time
!  step on the previous day.  Consider in future making the specification of time step as number of time
!  steps in a day, not modeldt to ensure this 
!  modeldt is recalculated based on the integer timesteps in a day
!  assumption: number of model timesteps in a day must be an integer  
           
      StepInADay=int(24.0/modeldt+0.5)  ! closest rounding
      Modeldt=24.0/StepInADay
      nstepday=StepInADay
      Allocate(timeMaxPerFile(MaxNumofFile,11))
      Allocate(timeMinPerFile(MaxNumofFile,11))
      allocate(TSV(maxncfilents,11))
      allocate(AllValues(maxncfilents,11))
      CALL timeSeriesAndtimeSteps(MaxNumofFile,NUMNCFILES,IsInputFromNC,InputTSFilename,NCDFContainer,&
      maxncfilents,NOofTS,TSV,Allvalues,Inputtcoordinates,daysstring)
       
      ReferenceHour=0.00
      IF (MaxNumofFile .eq. 0)then
        nrefyr=ModelStartDate(1)
        nrefmo=ModelStartDate(2)
        nrefday=ModelStartDate(3)
      END IF
      CALL JULDAT(nrefyr,nrefmo,nrefday,ReferenceHour,Referencetime)

      
      ! Dimvalue1 and dimvalue2 are used to put the dimension values in output netCDFs
      ! this value have precision problem and therefore, ArcGIS canot read output netCDFs
!      if (nlines .eq. GridModel)THEN
!          allocate(DimValue1(dimlen1))
!          allocate(DimValue2(dimlen2))
!          !    SpatialCoordinate(File_name,dimlen1,dimlen2,DimName1,DimName2,DimValue1,DimValue2,DimUnit1,DimUnit2)
!          CALL SpatialCoordinate(Watershedfile,wsycoordinate,wsxcoordinate,DimValue1,DimValue2,DimUnit1,DimUnit2)
!          deallocate(DimValue1)
!          deallocate(DimValue2)
!      END IF
      
!      OPEN(665,FILE='date.DAT',STATUS='UNKNOWN')
!      Do I = 1,dimlen1
!            Write(665,37) DimValue1(i)
!      end do
!37    format(f17.10,1x)
!      close(665)

      tresult = Etime(tarray)
      tlast=tresult
      write(6,*)'time to check inputs: ',tresult
      write(636,*)'Check Inputs: ',tresult, ' seconds'
      
      if (GridModel)THEN
        ! Checking netCDFs starts here                    
        CALL checks(svfile,IsInputFromNC,NumNCFiles,totalNC,StateSiteVName)
        allocate(AllNCDFfile(totalNC))
        ! checking if all the netCDF files are provided in a desired format
        CALL  NCChecks(svfile,StateSiteVName,WatershedFile,MaxNumofFile,IsInputFromNC,NCDFContainer,NumNCFiles,totalNC,AllNCDFfile,&
                       Inputxcoordinates,Inputycoordinates,Sitexcoordinates,Siteycoordinates,wsxcoordinate,wsycoordinate)            
                    
  ! Checking netCDFs starts here 

  ! Output file creation starts here
        CALL NumOutFiles(OutControlFILE, ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,&
        NumtimeStep,NumofFile,NumOutPoint,OutCount,dimlen1,dimlen2)
                         
        Allocate(NumtimeStepPerFile(NumofFile))          ! This array will contains the number of time steps in the sequence of netcdf files
        Allocate(OutputNCContainer(NumofFile,OutCount))  !  This array will contain the file names for each output netcdf variable
        Allocate(NCOutfileArr(NumofFile,OutCount))  
        Allocate(OutPoint(NumOutPoint,2))
        Allocate(OutPointFiles(NumOutPoint))
        Allocate(outSampleFile(outcount))
        Allocate(OutVar(outcount))
        Allocate(outputfolder(outcount))
        Allocate(OutVarValue(NumtimeStep,66))
        CALL OutputFiles(OutControlFILE,NumtimeStep,NumofFile,outSampleFile,NumtimeStepPerFile,dimlen1,dimlen2,&
        OutVar,OutPoint,OutPointFiles,NumOutPoint,OutCount)
        Allocate(NCIDARRAY(NumofFile,outcount))
        Allocate(OutFolder(outcount))
        DimName1=wsycoordinate
        DimName2=wsxcoordinate
         
        ! dimlen1=340
        ! dimlen2=309
        CALL DirectoryCreate(nrefyr,nrefmo,nrefday,WATERSHEDfILE,wsycoordinate,&
        &wsxcoordinate,NumofFile,outcount,Outvar,NumtimeStepPerFile,&
        &outSampleFile,OutputNCContainer,NCIDARRAY)
     else
          MStartHour=dble(ModelStartHour)
          MEndHour=dble(ModelEndHour)
          Call JULDAT(ModelStartDate(1),ModelStartDate(2),ModelStartDate(3),MStartHour,JMSD) !JMSD= julian model start date-time
          Call JULDAT(ModelEndDate(1),ModelEndDate(2),ModelEndDate(3),MEndHour,JMED)         !JMED= julian model start date-time
          tol=5./(60.*24.)     !  5 min tolerance
          NumtimeStep=(((JMED+tol-JMSD)/(Modeldt/24.0)) + 1)   
      End if
        
      Allocate(ReGriddedArray(NumtimeStep,11))
      Allocate(modelTimeJDT(NumtimeStep))
      Allocate(CurrentArrayPosRegrid(NumtimeStep,11))
        
! Output file creation ends here
      CALL InputVariableValue(INPUTVARNAME,IsInputFromNC,NoofTS,TSV,Allvalues,maxncfilents,ModelStartDate,ModelStartHour,&
      ModelEndDate,ModelEndHour,nrefyr,nrefmo,nrefday,modeldt,NumtimeStep,CurrentArrayPosRegrid,modelTimeJDT)
       
       !*******************************************************************************           
!        OPEN(665,FILE='DataBeforeRegrid.DAT',STATUS='UNKNOWN')
!        Do I = 1,maxncfilents
!            Write(665,37) TSV(i,1),TSV(i,2),TSV(i,3),TSV(i,4),TSV(i,5),TSV(i,6),&
!            &TSV(i,7),TSV(i,8),TSV(i,9),TSV(i,10),TSV(i,11)
!37          format(f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5)
!        End do
!        Close(665)
!        OPEN(668,FILE='CurrentArrayPosRegrid.DAT',STATUS='UNKNOWN')
!        Do I = 1,NumtimeStep
!            Write(668,39) CurrentArrayPosRegrid(i,1),CurrentArrayPosRegrid(i,2),CurrentArrayPosRegrid(i,3),CurrentArrayPosRegrid(i,4),&
!            CurrentArrayPosRegrid(i,5),CurrentArrayPosRegrid(i,6),CurrentArrayPosRegrid(i,7),CurrentArrayPosRegrid(i,8),&
!            CurrentArrayPosRegrid(i,9),CurrentArrayPosRegrid(i,10),CurrentArrayPosRegrid(i,11)
!39          format(I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5,1x,I5)
!        END DO
!        Close(668)
      !*******************************************************************************         

      If (GridModel)THEN
        Allocate(StartEndNCDF(NumofFile,2))
  ! Work with aggregated outputs    
        CALL AggregatedOutNum(AggOutControl,outSymbol,AggOutNum)
        Allocate(AggOutVar(AggOutNum))
        Allocate(AggOutVarnum(AggOutNum))
        CALL  AggOutWSUniqueID(AggOutControl,outSymbol,AggOutNum,AggOutVar,Watershedfile,WatershedVARID,&
                              &uniqueIDNumber,AggOutVarnum,WsMissingValues,WsFillValues,wsxcoordinate,wsycoordinate)
        Allocate(UniqueIDArray(uniqueIDNumber))
        Allocate(UniqueIDArrayCount(uniqueIDNumber))
        CALL WSUniqueArray(Watershedfile,WatershedVARID,dimlen1,dimlen2,uniqueIDNumber,UniqueIDArray,WsMissingValues,WsFillValues,UniqueIDArrayCount,&
                          &wsxcoordinate,wsycoordinate)
        Allocate(AggdWSVarVal(NumtimeStep,uniqueIDNumber,AggOutNum))
        Allocate(AggdWSVarValAvg(NumtimeStep,uniqueIDNumber,AggOutNum))
        AggdWSVarVal=0                        
        AggUnit=887
        OPEN(AggUnit,FILE=AggdOutput,STATUS='unknown')
        write(Aggunit,*)'Year month day hour variable watershed value' ! write header
       
        St=0 
        ii=1
       
  1056   If (ii .LE. NumofFile)THEN
            StartEndNCDF(ii,1)=ST+1
            StartEndNCDF(ii,2)=StartEndNCDF(ii,1)+NumtimeStepPerFile(ii)-1
            ST=StartEndNCDF(ii,2)
            ii=ii+1
            Go to 1056
        End if
       
        ! required to calculate percent grid completed
        totalgrid=dimlen1*dimlen2     
        numgrid=0
      END IF
       
      Allocate(FNDJDT(NumtimeStep))
      Allocate(yymmddarray(3,NumtimeStep))
      Allocate(timearray(NumtimeStep))
      ! calculating model end date-time in julian date
      dhour=dble(ModelEndHour)
      call JULDAT(ModelEndDate(1),ModelEndDate(2),ModelEndDate(3),dhour,EJD)
       
      !  Initialize timing results
      tresult= Etime(tarray)
      write(6,*)'time to create netCDFs: ',(tresult-tlast)
      write(636,*)'Create netCDFs: ',(tresult-tlast),' seconds'  
      write(6,*)"Starting loop over grid cells"
       
      ! time tracking variables are inititated as 0.0
      tlast=tresult
      tio=0.0
      tcomp=0.0
      tout=0.0
      tagg=0.0
      taggre=0.0
      toutnc=0.0 
      tsitev=0.0
      dt=Modeldt
      !WsMissingValuesInt = int(WsMissingValues,8)
      !WsFillValuesInt = int(WsFillValues,8)
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Space loop starts here
      !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
      DO iycoord=1,dimlen1
      DO jxcoord=1,dimlen2
        iunit=119  !  unit for point output       
        if (GridModel)THEN
          CALL nCDF2DReadInteger(Watershedfile,WatershedVARID,IDNumber(1),iycoord,jxcoord,wsxcoordinate,wsycoordinate)
        ELSE
          IDNumber(1)=1
        END IF
        ! Perform calculations only if grid cell is in the watershed
        if((IDNumber(1) .ne. -9999) .and. (IDNumber(1) .ne. WsMissingValues) .and. (IDNumber(1) .ne. WsFillValues))then  
!  TODO change the above to also exclude netcdf no data values

          !  read site variables and initial conditions of state variables
          CALL readsv(param,statev,sitev,svfile,slope,azi,lat,subtype,iycoord,jxcoord,dtbar,Ts_last,longitude,Sitexcoordinates,Siteycoordinates)
          ! Set water equivalent thickness of glacier based on substrate type
          ! subtype = sitev(10)
          !   0 = Ground/Non Glacier, 
          !   1 = Clean Ice/glacier, 
          !   2 = Debris covered ice/glacier, 
          !   3 = Glacier snow accumulation zone
          ! 
          IF(SITEV(10) .EQ. 0 .OR. SITEV(10) .EQ. 3)THEN
              WGT=0.0
          ELSE
              WGT=1.0
          END IF
    
          tresult= Etime(tarray)
          tsitev=tsitev+tresult-tlast
          tlast=tresult
       
          ALLOCATE(Tsprevday(nstepday))
          ALLOCATE(Taveprevday(nstepday))
       
          if(sitev(10) .ne. 3) then   !  Only do this work for non accumulation cells where model is run 
            ! TODO check logic and consolidate with to perform calculations only on grid
            CALL Values4VareachGrid(inputvarname,IsInputFromNC,MaxNumofFile,NUMNCFILES,NCDFContainer,varnameinncdf,iycoord,jxcoord,&
            &NCfileNumtimesteps,NOofTS,maxncfilents,Allvalues,VarMissingValues,VarfILLValues,StepInADay,NumtimeStep,CurrentArrayPosRegrid,ReGriddedArray,&
            &Inputxcoordinates,inputycoordinates,inputtcoordinates,InputVarRange,InpVals)
       
      !*******************************************************************************              
!       OPEN(666,FILE='DataAfterRegrid.DAT',STATUS='UNKNOWN')
!       Do I = 1,NumtimeStep
!            Write(666,40) ReGriddedArray(i,1),ReGriddedArray(i,2),ReGriddedArray(i,3),ReGriddedArray(i,4),ReGriddedArray(i,5),ReGriddedArray(i,6)
!40          format(f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5,1x,f17.5)
!       End do
!       Close(666)
      !*******************************************************************************  

      ! FIXME: what if the result is fractional
      !  should not be fractional (except numerically)
            nstepday=StepInADay  ! number of time steps/day


            !  Initialize Tsbackup and TaveBackup
            DO 3 i = 1,nstepday
                Tsprevday(i)=-9999.
                Taveprevday(i)=-9999.0
3            CONTINUE

            ! Take surface temperature as 0 where it is unknown the previous time step
            ! This is for first day of the model to get the force restore going
            IF(ts_last .le. -9999.)THEN    
              Tsprevday(nstepday)=0  
            ELSE
              Tsprevday(nstepday)=ts_last                      !has measurements
            END IF 

    ! compute Ave.Temp for previous day 
            Us   = statev(1)                                  ! Ub in UEB
            Ws   = statev(2)                                  ! W in UEB
            Wc   = statev(4)                                  ! Canopy SWE
            Apr  = sitev(2)                                   ! Atm. Pressure  (PR in UEB)
            cg   = param(4)                                   ! Ground heat capacity (nominally 2.09 KJ/kg/C)
            rhog = param(8)                                   ! Soil Density (nominally 1700 kg/m^3)
            de   = param(11)                                  ! Thermally active depth of soil (0.1 m)
 
            Tave = TAVG(Us,Ws+WGT,RHOW,CS,TO,RHOG,DE,CG,HF)     ! This call only
            Taveprevday(nstepday) = Tave
      
        !  initialize variables for mass balance
            Ws1   = statev(2)
            Wc1   = statev(4)     
            cumP = 0.
            cumEs = 0.
            cumEc = 0.
            cumMr = 0.
            cumGM = 0.
            cumEG = 0.
    !************************************************************************************************      

        !  find out if this is an output point and if so open the point output file
            IF(GridModel)THEN
              towrite=.false.
              do NumOP=1,NumOutPoint
                If (iycoord .eq. OutPoint(NumOP,1) .and. jxcoord .eq. OutPoint(NumOP,2))Then  
                    OPEN(iUnit,FILE=OutPointFiles(NumOP),STATUS='unknown')
                    towrite= .true.
                    exit
                endif 
              enddo
            ELSE ! For point mode always write output
                towrite=.TRUE.
                OPEN(iUnit,FILE='PointOutput.DAT',STATUS='unknown')
            END IF                     
          endif  !  end the skip block done only for accumulation cells

          ! Variables to keep track of which time step we are in and which netcdf output file we are in
          istep=0  ! time step initiated as 0
        
          ! map on to old UEB names
          YEAR=ModelStartDate(1)
          MONTH=ModelStartDate(2)
          DAY=ModelStartDate(3)
          Hour=ModelStartHour
          SHOUR=DBLE(ModelStartHour)
          call JULDAT(YEAR,MONTH,DAY,SHOUR,CurrentModelDT)

          tresult= Etime(tarray)
          tio=tio+tresult-tlast
          tlast=tresult


          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! This is the start of the main time loop 
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++              
          do while(EJD >= CurrentModelDT)
            istep=istep+1
            yymmddarray(1,istep)=year
            yymmddarray(2,istep)=month
            yymmddarray(3,istep)=day
            timearray(istep)=hour
        
            IF(sitev(10).NE. 3)THEN
            Do i= 1,11
                InpVals(i)=ReGriddedArray(istep,i)
            End do
                     
    !        Write(667,38)istep,InpVals(1),InpVals(2),InpVals(3),InpVals(4),InpVals(5),InpVals(6),&
    !            &InpVals(7),InpVals(8),InpVals(9),InpVals(10),InpVals(11)

      !      Map from wrapper input variables to UEB variables     
            TA=INPVals(1)
            P=INPVals(2)
            V=INPVals(3)
            RH=INPVals(4)
            Tmin=INPVals(5)
            Tmax=INPVals(6)
            trange=Tmax-Tmin
            if (trange .LE. 0)THEN
                If (snowdgtvariteflag .EQ. 1)then
                    write(6,*) "Input Diurnal temperature range is less than or equal to 0 which is unrealistic "
                    write(6,*) "Diurnal temperature range is assumed as eight degree celsius "
                    write(66,*)"on ",year,month,day
                End  if
                trange=8.0
            End if
            QSIOBS=INPVals(7)
            qg=INPVals(8)
            qli=INPVals(9)
            QNETOB=INPVals(10)
            Snowalb=INPVals(11)
            if((ireadalb==1 .and. Snowalb .LT. 0.0) .or. (ireadalb==1 .and. Snowalb .GT. 1.0) .or. ireadalb==0)then
               ireadalb=0
            else
               ireadalb=1   
            end if
    !  Below is code from point UEB 
            sitev(3)=qg   
            INPT(1,1)=TA
            INPT(2,1)=P
            INPT(3,1)=V
            INPT(4,1)=RH
            INPT(7,1)= QNETOB
        
    ! UTC to local time conversion
            CALL CALDAT(modelTimeJDT(istep),YEAR,MONTH,DAY,DBHour)
            Hour=REAL(DBHour)
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
!      
! We found that Model reanalysis and dowscaled data may produce some unreasonably negative solar radiation.
! this is simply bad data and it is generally better policy to try to give a model good data. 
! If this is not possible, then the UEB checks will avoid the model doing anything too bad, 
! it handles negative solar radiation in following way:

! "no data in radiation would be to switch to the temperature method just for time steps 
! when radiation is negative." 
                IF(IRAD.EQ.0 .or. (QSIOBS .lt. 0.)) THEN    !  For cases where input is strictly negative we calculate QSI from HRI and air temp range.  This covers the case of missing data being flagged with negative number, i.e. -9999.                 
                    INPT(5,1)=atff*IO*HRI
                    CALL cloud(param,atff,cf)   ! For cloudiness fraction
                ELSE   ! Here incoming solar is input
    !      Need to call HYRI for horizontal surface to perform horizontal
    !      measurement adjustment
                    CALL hyri(MYEAR,MMONTH,MDAY,NHOUR,DT,0.0,AZI,LAT,HRI0,COSZEN)
    !      If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
    !      not it indicates a potential measurement problem. i.e. moonshine
                    if(HRI0 .GT. 0.0) then
                        atfimplied=min(qsiobs/(HRI0*IO),0.9) ! To avoid unreasonably large radiation when HRI0 is small
                        INPT(5,1)=atfimplied * HRI * IO
                    else
                        INPT(5,1)=QSIOBS
                          if(qsiobs .ne. 0.)then
                            If (snowdgtvariteflag .EQ. 1)then
                                write(66,*)"Warning ! you have nonzero nightime"
                                write(66,*)"incident radiation of",qsiobs
                                write(66,*)"at date",year,month,day,hour
                            End if
                          endif
                      endif
                      CALL cloud(param,atff,cf)   ! For cloudiness fraction  This is more theoretically correct
                  ENDIF        
              IF(irad .lt. 2)THEN    
                CALL qlif(TA,RH,TK,SBC,Ema,Eacl,cf,inpt(6,1) )
              Else
                Ema=-99  !  These values are not evaluated but may need to be written out so are assigned for completeness
                Eacl=-99
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
            if(day .eq. 29)then
              mtime(3) = day
            endif
            CALL SNOWUEB2(dt,1,inpt,sitev,statev,tsprevday, taveprevday, nstepday, param,iflag,&  
              &cump,cumes,cumEc,cummr,cumGM,outv,mtime,atff,cf,OutArr,CumEG)
     
            tresult= Etime(tarray)
            tcomp=tcomp+tresult-tlast
            tlast=tresult

    !        IF(sitev(10).EQ. 3)THEN  ! Substrate type is accumulation zone
    !           OutArr=0
    !        END IF
        
            if(towrite)WRITE(iunit,*)OutArr       
            DStorage=statev(2)-Ws1+statev(4)-Wc1
            errmbal= cump-cumMr-cumEs-cumEc-DStorage+cumGM+cumEG

            if(towrite)WRITE(iunit,*)ERRMBAL
            if (GridModel)THEN
            !mapping to OutVarValue
            OutVarValue(istep,1:12)=IniOutVals(5:16)
            OutVarValue(istep,13:62)=OutArr(1:50)
            OutVarValue(istep,63)=ERRMBAL
            OutVarValue(istep,64:66)=OutArr(51:53)
            END IF
            ! These settings tell netcdf to write one timestep of data. (The
            ! setting of start(4) inside the loop below tells netCDF which
            ! timestep to write.)
  
            ELSE  ! this block entered only if sitev(10)= 3
                OutArr=0.0
                ERRMBAL=0.0
                if(towrite)WRITE(iunit,*)OutArr 
                if(towrite)WRITE(iunit,*)ERRMBAL 
                OutVarValue(istep,1:66)=0.00
            END IF
        
            tresult= Etime(tarray)
            tout=tout+tresult-tlast
            tlast=tresult
        
            if (GridModel)THEN
              ! Here we rely on the even spread of time steps until the last file
              incfile=(istep-1)/NumtimeStepPerFile(1)+1  !  the netcdf file position
            END IF
        
            ReferenceHour=DBLE(hour)
            call JULDAT(YEAR,MONTH,DAY,ReferenceHour,CTJD)  !  current julian date time
            FNDJDT(istep)=DBLE(CTJD-Referencetime)
        
            IF (GridModel)THEN      
              IDNum=Int(IDNumber(1))
              if((IDNum .ne. -9999) .and. (IDNum .ne. WsMissingValues) .and. (IDNum .ne. WsFillValues))then
                  Do jUniqueID=1,uniqueIDNumber
                      If(UniqueIDArray(jUniqueID) .eq. IDNum)then
                          do ioutvar=1,AggOutNum
                              AggValues=OutVarValue(istep,AggOutVarnum(ioutvar))
                              AggdWSVarVal(istep,jUniqueID,ioutvar)=AggdWSVarVal(istep,jUniqueID,ioutvar)+AggValues
                              ! Averaging required
                              AggdWSVarValAvg(istep,jUniqueID,ioutvar) = AggdWSVarVal(istep,jUniqueID,ioutvar)/UniqueIDArrayCount(jUniqueID)
                          END DO
                          Exit
                      End if
                  end do
              end if
            END IF


            CALL UPDATEtime(YEAR,MONTH,DAY,HOUR,DT)   
            ModHour=DBLE(Hour)
            call JULDAT(YEAR,MONTH,DAY,ModHour,CurrentModelDT)
        
            tresult= Etime(tarray)
            taggre=taggre+tresult-tlast
            tlast=tresult
          end do   
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! This is the End of the main time loop 
          !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
      !If (EJD .GE. CurrentModelDT)Then      
      !    Go to 1
      !End if
      !
       
          deallocate(Tsprevday)
          deallocate(Taveprevday)
          Close(iunit)
        
          IF (GridModel)THEN 
            do ioutv=1,outcount
              do incfile = 1,NumofFile
                  CALL OutputnetCDF(NCIDARRAY,outvar,NumtimeStep,outcount,incfile,ioutv,jxcoord,iycoord,NumtimeStepPerFile,NumofFile,&
                  &StartEndNCDF,OutVarValue)  
                  CALL check(nf90_sync(NCIDARRAY(incfile,ioutv)))
              enddo
            enddo
          END IF
        endif  !  this is the end of if we are in a watershed    
        
        tresult= Etime(tarray)
        toutnc=toutnc+tresult-tlast
        tlast=tresult
    
        ! grid count progress bar is calculated and written here
        numgrid=numgrid+1
        write(6,FMT="(A1,A,t30,F6.2,A$)") achar(13), " Percent Grid completed: ", (real(numgrid)/real(totalgrid))*100.0, "%"
    
      END DO  !  These are the end of the space loop
      END DO    
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Space loop ends here
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
      IF (GridModel)THEN 
         ! Putting all the dimension values inside the netcdf files
        do ioutv=1,outcount
          do incfile = 1,NumofFile
            CALL OutputtimenetCDF(NCIDARRAY,NumtimeStep,outcount,incfile,ioutv,NumtimeStepPerFile,NumofFile,StartEndNCDF,FNDJDT) ! time is dimension 1
            CALL check(nf90_sync(NCIDARRAY(incfile,ioutv)))     !syncing the netcdf
            CALL Check(nf90_close(NCIDARRAY(incfile,ioutv)))    !closing the netcdf
          enddo
        enddo
        Write (6,FMT="(/A32)") " Now aggregation will be started"
        Write (6,FMT="(A35/)")  " Time taken for various operations:"
    
        do istep=1,Numtimestep
          do ivar=1,AggOutNum
            Do jUniqueID=1,uniqueIDNumber
              if((UniqueIDArray(jUniqueID) .ne. -9999) .and. (UniqueIDArray(jUniqueID) .ne. WsMissingValues) .and.&
                &(UniqueIDArray(jUniqueID) .ne. WsFillValues))then
                  WRITE(Aggunit,47)yymmddarray(1,istep),yymmddarray(2,istep),yymmddarray(3,istep),timearray(istep),outSymbol(AggOutVarnum(ivar)),&
                  &UniqueIDArray(jUniqueID),AggdWSVarValAvg(istep,jUniqueID,ivar)
47                format(1x,i5,i3,i3,f8.3,1x,a11,i7,1x,g13.6)
              end if
            end do
          enddo
        end do
        Close(AggUnit)
      END IF
    
      tresult= Etime(tarray)
      tagg=tagg+tresult-tlast
      tlast=tresult

      tarray(1)=tarray(1)
      write(636,*) "Reading SiteVariable: ",tsitev," Seconds"
      write(636,*) "Input time: ",tio," Seconds"
      write(636,*) "Compute time: ",tcomp," Seconds"
      write(636,*) "Out time: ",tout," Seconds"
      write(636,*) "Out timeinNC: ",toutnc," Seconds"
      write(636,*) "Aggregated StorageArray: ",taggre," Seconds"
      write(636,*) "Aggregation time: ",tagg," Seconds"
      write(636,*) "Complete runtime: ",tarray(1)," Seconds"
      Close(636)
    
      write(6,*) "Reading SiteVariable: ",tsitev," Seconds"
      write(6,*) "Input time:",tio," Seconds"
      write(6,*) "Compute time:",tcomp," Seconds"
      write(6,*) "Out time:",tout," Seconds"
      write(6,*) "Out timeinNC:",toutnc," Seconds"
      write(6,*) "Aggregated StorageArray:",taggre," Seconds"
      write(6,*) "Aggregation time:",tagg," Seconds"
      write(6,*) "Complete runtime:",tarray(1)," Seconds"
      Write(6,*) "Your task is successfully performed! Plesae view the results in 'outputs' folder!"
    
end program