!**********************************************************************************************
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

! ============================== This subroutine reads the values from input control file ========================================

 subroutine InputMaxNCFiles(inputcon,MaxNumofFile,inputvarname,UTCOffSet)
 
        ! inputcon (input) is name of control file
        ! MaxNumofFile (input) is the maximum number of NC files for any variable 
        ! inputvarname (output) Name of the variables that are provided inside inputcontrol.dat file
        ! UTCOffSet(output) UTC offset for the region/wateshed we are dealing with
        
        PARAMETER(n=11)
        ! TA,P,V,RH,trange,QSIOBS,QNETOB,qg,qli,inputcon, iloop
        integer:: IPDT,reason,ilatp,dimlen2,dimlen1, IsVarFromNC(n),InumOfFile(n),IsInputFromNC(n)
        integer:: Syear, Smonth, sdate,  Eyear, Emonth, Edate !s for tartng date-time and e for ending date-time
        real::  shour, Ehour, dt, EndingDate,VARVALUES(n),UTCOffSet     !dt=time increment in hours   
        CHARACTER*200 svfile, inputHeading, inputcode, inputname, inputcon, inputVName(n), InputNCFilename(n), InputTSFilename(n),Readfile
        CHARACTER*200 NCfileContain, NCfileContainer,varnameinncdf(n),inputvarname(n)
        !InputNCFilename(n)=here inpdex file names will be stored for those variables are both spatially-temporally variable (SVTV)
        !InputNCFilename(n)=here time series text file names will be stored for those variables are both temporally variable but spatially constant(SCTV)
        InputVName= (/ "Ta     ","Prec   ","v      ","RH     ", &
             "Tmin   ","Tmax   ","Qsi    ","Qg     ","Qli    ", &
             "Qnet   ","Snowalb"  /)
        IsInputFromNC=(/ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)
          OPEN(19,FILE=inputcon,STATUS='OLD')
          Read (19,*)  inputHeading
          Read (19,*)  Syear, Smonth, sdate, shour
          Read (19,*)  Eyear, Emonth, Edate, Ehour
          Read (19,*)  dt
          Read (19,*)  UTCOffSet
500    	  Read(19,*,iostat=reason, end=600)inputcode
    	    if(reason .eq. 0) then
            do i=1,n,1
                inputname=inputcode(1:(SCAN (inputcode, ':')-1))
                CALL lowercase(inputname,inputname)
                CALL lowercase(InputVName(i),InputVName(i))
                if(inputname .eq. trim(inputVName(i))) then
                    inputvarname(i)=inputname
                    xx=1
                    inputvarname(i)=inputname
                    read(19,*)IsInputFromNC(i)                                                                   
                    if(IsInputFromNC(i) .eq. 0) then
                        read(19,*)InputTSFilename(i)
                    else if(IsInputFromNC(i) .eq. 1) then
                        read(19,*)InputNCFilename(i)
                        read(19,*)varnameinncdf(i) 
                    else
                        read(19,*)VARVALUES(i) 
                    end if
                end if 
            end do
            if(xx .ne. 1)write(6,*)'Input variable code not matched:  ',inputname
        end if
        xx=0
        go to 500
600     CLOSE(19)
        

        Do i =1,n
        Count=0
        if (IsInputFromNC(i) .eq. 1) then
            NCfileContainer=InputNCFilename(i) !TSFile= time series file
            OPEN(49,FILE=NCfileContainer,STATUS='OLD')
1500        Read(49,*,end=1600)NCfileContain
            count=1+count
            go to 1500
1600        CLOSE(49)
            InumOfFile(i)=count
        end if
        end do
            
        MaxNumofFile = MAXVAL(InumOfFile)
        end subroutine InputMaxNCFiles
!================================================================================================================================
        Subroutine InputFiles(inputcon,MaxNumofFile,IsInputFromNC,NumNCFiles,InputTSFilename,NCDFContainer,&
        &ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,NCfileNumtimesteps,&
        &nrefyr,nrefmo,nrefday,varnameinncdf,arrayx,NoofTS,InpVals,VarMissingValues)
        ! inputcon (input) is name of control file
        ! MaxNumofFile (input) is the maximum number of NC files for any variable 
        ! IsInputFromNC(n) (Output) Array indicating whether variable is from NC (0 for TS, 1 for NC, 2 for value, 3 for not provided)
        ! NumNCFiles(n) (Output)  Array giving the number of NC files for NC variables
        ! InputTSFilename(n) (Output) Array giving file names for variables with time series input
        ! NCDFContainer(MaxNF, n) (Output)  Array giving file names for variables with NC files
        ! ModelStartDate(3) (Output) Array giving start year, month, day
        ! ModelStartHour (output) Start hour
        ! ModelEndDate(3) (output)  Array giving end year, month, day
        ! ModelEndHour (output)  end hour
        ! Modeldt (output)  time step
        ! NCfileNumtimesteps(MaxNF,n) (output) array giving the number of time steps in each NC file
        ! nrefyr,nrefmo,nrefday,(output) Reference year, month, day from netcdf time units specified as "days from y/m/d"
        !  Our convention is that our start hour must be 0
        ! varnameinncdf(n) (output)  variable name in control file specifying the netcdf variable to use for reading from NC
        ! arrayx  (output).  Maximum across variables of the sum of number of time steps in all NC files for that variable
        ! NoofTS(n) (output).  The number of combined NC time steps for each variable
        ! InpVals(n) (output).  Variable to hold the current value of each input variable.  This subroutine fills this for inputs that are constant
        
        Use netCDF
        
        parameter(n=11)                                         !n is  loop variable
        integer:: InumOfFile(n), IsInputFromNC(n), reason, NumNCFiles(n),nrefyr,nrefmo,nrefday,count
        integer::syear,smonth,sday,first
        integer:: ModelStartDate(3),  ModelEndDate(3)  
        real:: ModelStartHour, ModelEndHour,Modeldt,InpVals(n)  
        integer:: NCfileNumtimesteps(MaxNumofFile,n)
        integer:: NOofTS(n)
        double precision NCStartTime(MaxNumofFile,n) 
        real:: UTCOffSet
        character*200:: NCDFContainer(MaxNumofFile,n)
        !  Arrays to hold temporary information before sort
        character*200, allocatable :: tempfilelist(:)
        integer, allocatable :: tempfilesteps(:)
        double precision, allocatable :: tempstarttime(:)
        integer,allocatable :: iy(:)           
        CHARACTER*200 inputHeading, inputcode, inputname, inputcon, inputVName(n), InputNCFilename(n)
        CHARACTER*200 InputTSFilename(n),Readfile, NCfileContainer
        CHARACTER*200 NCfileContain 
        Character*50::varnameinncdf(n)
        Character*200:: Rec_name,File_name
        integer:: VarID, ArrayStart,ArrayEnd
        Integer:: FileNextT(3)
        REAL:: FileNextH, FileNextVal
        Double precision:: FileNextHR
        integer:: FileOpenFlag(n)
        Double precision:: FileNextDateJDT
        Character*200:: TSFile, CurrentInputVariable(n)
        Double precision:: TimeMaxPerFile(MaxNumofFile,n)
        Double precision:: TimeMinPerFile(MaxNumofFile,n)
        integer:: TSCounts
        Double precision, Allocatable :: AllTimeSteps(:)
        Double precision:: ModelStartJDT, ModelEndJDT,NCRefernceJDT,RefHour
        Integer:: TotalTS,StartY,StartM,StartD
        Real:: TVal
        integer::arrayx
        Character*200:: file_nameM
        character (len = *), parameter :: missing_value = "missing_value"
        real:: VarMissingValues(MaxNumofFile,n),VarMissval
        Character*200, Allocatable:: AllFiles(:)
        allocate(tempfilelist(MaxNumofFile))
        allocate(tempfilesteps(MaxNumofFile))
        allocate(tempstarttime(MaxNumofFile))
        allocate(iy(MaxNumofFile))

        !Double precision:: NCfiletimeDimesions(MaxNumofFile,n,3)
        !InputNCFilename(9)=here inpdex file names will be stored for those variables are both spatially-temporally variable (SVTV)
        !InputNCFilename(9)=here time series text file names will be stored for those variables are both temporally variable but spatially constant(SCTV)
        InputVName= (/ "ta     ","prec   ","v      ","rh     ", &
             "tmin   ","tmax   ","qsi    ","qg     ","qli    ", &
             "qnet   ","snowalb"  /)
        IsInputFromNC=(/ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 /)  !  3 is used to indicate this variable is not input
          OPEN(99,FILE=inputcon,STATUS='OLD')
          Read (99,*)  inputHeading
          Read (99,*)  ModelStartDate(1),ModelStartDate(2),ModelStartDate(3),ModelStartHour
          Read (99,*)  ModelEndDate(1),ModelEndDate(2),ModelEndDate(3),ModelEndHour
          Read (99,*)  Modeldt
          Read (99,*)  UTCOffSet
1300      Read(99,*,iostat=reason, end=1400)inputcode
    	    if(reason .eq. 0) then
    	     inputname=inputcode(1:(SCAN (inputcode, ':')-1))
             CALL lowercase(inputname,inputname)
             do i=1,n,1
               if(inputname .eq. trim(inputVName(i))) then
                    xx=1
                    read(99,*)IsInputFromNC(i)                                                                   
                    if(IsInputFromNC(i) .eq. 0) then
                        read(99,*)InputTSFilename(i)
                    elseif(IsInputFromNC(i) .eq. 1)then
                        read(99,*)InputNCFilename(i)  !  This is file that contains list of NC files
                         read(99,*)varnameinncdf(i)
                    elseif(IsInputFromNC(i) .eq. 2)then
                        read(99,*)InpVals(i)
                    elseif(IsInputFromNC(i) .ne. 3)then
                        write(6,*)'Input variable flag not matched:  ',inputname, IsInputFromNC(i)
                    end if
                end if 
            end do
            if(xx .ne. 1)write(6,*)'Input variable code not matched:  ',inputname
        end if
        xx=0
        go to 1300
1400    CLOSE(99)
!  At this point we have IsInputFromNC, InputTSFilename, InputNCFilename, varnameinncdf, InpVals populated
!  Done with control file
        first=1
        do i=1,n
        count=0
        if (IsInputFromNC(i) .eq. 1) then
   !     NCfileContainer=InputNCFilename(i) !NCfileContainer= contains the name of netCDF files for 1 variable
        OPEN(59,FILE=InputNCFilename(i),STATUS='OLD')
1700    Read(59,*,end=1800)NCfileContain
        count=count+1
        tempfilelist(count)=NCfileContain 
        NumNCFiles(i)=count
        iy(count)=count
        CALL nCDF3DTimeRead (NCfileContain,1,tempstarttime(count),tempfilesteps(count),syear,smonth,sday)
        if (first .eq. 1)then
           first=0
           nrefyr=syear
           nrefmo=smonth
           nrefday=sday
        else
           if((nrefyr .ne. syear) .or. (nrefmo .ne. smonth) .or. (nrefday .ne. sday))then
              write(6,*)"Inconsistent NetCDF file time units encountered at file"
              write(6,*)NCfileContain
              write(6,*)"First file days since ",nrefyr, nrefmo,nrefday
              write(6,*)"current file days since ",syear,smonth,sday
            endif
          endif
        go to 1700
1800    CLOSE(59)
!  Sort NCfiles into order of start time
        call ssort(tempstarttime,iy,count,kflag)
        do ii=1,count
          NCDFContainer(ii,i)=tempfilelist(iy(count-ii+1))
          NCfileNumtimesteps(ii,i)=tempfilesteps(iy(count-ii+1))
        enddo  
        end if
        end do
        deallocate(tempfilelist)
        deallocate(tempfilesteps)
        deallocate(tempstarttime)
        deallocate(iy)
        NOofTS=0
        Varid=3  !  We require that time is the 3rd dimension
        Do i = 1,n 
            if (IsInputFromNC(i) .eq. 1)then
                Do k=1,NumNCFiles(i)
                    File_name=NCDFContainer(k,i)
                    call check(nf90_open(File_name, nf90_nowrite, ncidout))                         ! open the netcdf file
                    call check(nf90_inquire_dimension(ncidout, Varid, Rec_name,NumTimeStepEachNC))  ! information about dimensionID 3
                    Call check(nf90_inq_varid(ncidout,varnameinncdf(i),InputVarId))
                    CALL check(nf90_get_att(ncidout,InputVarId,missing_value,VarMissingValues(k,i)))  
                    NOofTS(i)=NOofTS(i)+NumTimeStepEachNC
                end do
            end if
        end do
        
        Do i=1,n
            FileOpenFlag(i)=999+i
            if(IsInputFromNC(i) .eq. 0)then                                                         !  This is from a text file
                TSFile=InputTSFilename(i)
                OPEN((FileOpenFlag(i)),FILE=TSFile,STATUS='OLD')
                READ ((FileOpenFlag(i)),*) CurrentInputVariable(i)
2749            READ ((FileOpenFlag(i)),*,END=2750) FileNextT(1), FileNextT(2), FileNextT(3), FileNextH, FileNextVal
                NOofTS(i)=NOofTS(i)+1
                GO to 2749
2750            close(FileOpenFlag(i))
            end if
        End do
        arrayx=MAXVAL(NOofTS)
        End subroutine
        
        subroutine TimeSeriesAndTimeSteps(MaxNumofFile,NUMNCFILES,IsInputFromNC,InputTSFilename,NCDFContainer,&
        &ModelStartDate,ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,NCfileNumtimesteps,&
        &varnameinncdf,arrayx,NOofTS,TSV,Allvalues)
        
        ! MaxNumofFile (input) is the maximum number of NC files for any variable 
        ! NumNCFiles(n) (input)  Array giving the number of NC files for NC variables
        ! IsInputFromNC(n) (inputt) Array indicating whether variable is from NC (0 for TS, 1 for NC, 2 for value, 3 for not provided)
        ! InputTSFilename(n) (inputt) Array giving file names for variables with time series input
        ! NCDFContainer(MaxNF, n) (input)  Array giving file names for variables with NC files
        ! ModelStartDate(3) (input) Array giving start year, month, day
        ! ModelStartHour (input) Start hour
        ! ModelEndDate(3) (input)  Array giving end year, month, day
        ! ModelEndHour (input)  end hour
        ! Modeldt (input)  time step
        ! NCfileNumtimesteps (MaxNF,n) (input) array giving the number of time steps in each NC file
        ! nrefyr,nrefmo,nrefday,(input) Reference year, month, day from netcdf time units specified as "days from y/m/d"
        !                                Our convention is that our start hour must be 0
        ! varnameinncdf(n) (input)  variable name in control file specifying the netcdf variable to use for reading from NC
        ! arrayx  (input).  Maximum across variables of the sum of number of time steps in all NC files for that variable
        ! NoofTS(n) (input).  The number of combined NC time steps for each variable
        ! TSV (arrayx,11) (output) holds all the timsteps bot from NC and time series (TS) text files. Time from TS files are stored as
        !                          julian and Time from NC files are stored as day/hour from the reference date (time unit in NC).
        ! Allvalues (arrayx,11) (output) holds values of enite timeseries for each variable for a particular grid  point
        
        use netcdf
        Parameter(n=11)
        integer:: arrayx
        integer:: VarID,ncidout,NCfileNumtimesteps(MaxNumofFile,n)
        Character:: Rec_Name
        integer:: MaxNumofFile
        integer:: TScounts,FileOpenFlag(n),IsInputFromNC(n),NUMNCFILES(n)
        character*200:: TSFile
        character*200:: InputTSFilename(n)
        character*200:: NCDFContainer(MaxNumofFile,n),File_nameM
        Integer:: NOofTS(n)
        integer:: FileNextT(3)
        real:: FileNextH
        Real:: FileNextVal  
        Double precision:: FileNextHR
        Double precision:: FileNextDateJDT
        REAL:: AllValues(arrayx,11)
        Double precision:: TSV(arrayx,11)
        Double precision,allocatable ::AllTimeSteps(:)
        character*200:: CurrentInputVariable(n)
        integer::TScount,StartY,StartM,StartD
        Double precision:: NCRfeferenceJDT
        Double precision:: arraymin,arraymax
        REAL:: MODELSTARTHOUR,MODELENDHOUR,MODELDT
        integer:: MODELSTARTDATE(3),MODELENDDATE(3)
        Character*50:: VARNAMEINNCDF(n)
        double precision:: RefHour
        Integer::TotalTS
        Double precision:: TVal(1)
        integer:: ArrayStart,ArrayEnd
        integer:: InputVarId
        
        Varid=3  !  We require that time is the 3rd dimension
        RefHour=0.00
        FileCount=1
        
        Do i=1,n
            ArrayEnd=0
            TScounts=1
            FileOpenFlag(i)=999+i
            if(IsInputFromNC(i) .eq. 0)then   !  This is from a text file
                TSFile=InputTSFilename(i)
                OPEN((FileOpenFlag(i)),FILE=TSFile,STATUS='OLD')
                READ ((FileOpenFlag(i)),*) CurrentInputVariable(i)
                If (TScounts .LE. NOofTS(i))Then   !  Make in to do loop - good practice
2849                READ ((FileOpenFlag(i)),*,END=2850) FileNextT(1), FileNextT(2), FileNextT(3), FileNextH, FileNextVal
                    FileNextHR=DBLE(FileNextH)
                    CALL JULDAT(FileNextT(1), FileNextT(2), FileNextT(3),FileNextHR,FileNextDateJDT)
                    TSV(TScounts,i)=FileNextDateJDT
                    AllValues(TScounts,i)=FileNextVal
                    TScounts=TScounts+1
                    GO to 2849
2850                close(FileOpenFlag(i))
                End if
            end if

            if (IsInputFromNC(i) .eq. 1)then
                Do k=1,NumNCFiles(i)
                    File_nameM=NCDFContainer(k,i)
                    FileCount=FileCount+1
                    call check(nf90_open(File_nameM, nf90_nowrite, ncidout))                         ! open the netcdf file
                    call check(nf90_inquire_dimension(ncidout,Varid,Rec_name,NumTimeStepEachNC))  ! information about dimensionID 3  
                    CALL nCDF3DTimeRead(file_nameM,1,TVal,TotalTS,StartY,StartM,StartD)
                    CALL JULDAT(StartY,StartM,StartD,RefHour,NCRfeferenceJDT)
                    Allocate(AllTimeSteps(numTimeStepEachNC))                                    
                    call check(nf90_get_var(ncidout,VarId,AllTimeSteps))                            ! Read the first time value from the file
                    ArrayStart=ArrayEnd+1                           
                    ArrayEnd=NumTimeStepEachNC+ArrayStart-1
                    TSV(ArrayStart:ArrayEnd,i)=AllTimeSteps
                    deallocate(AllTimeSteps)
                end do
            End if
        End do
        return
        end subroutine
        

        
        Subroutine Values4VareachGrid(IsInputFromNC,MaxNumofFile,NUMNCFILES,NCDFContainer,varnameinncdf,iycoord,jxcoord,&
            &NCfileNumtimesteps,NOofTS,arrayx,TSV,Allvalues,VarMissingValues)
            
        ! IsInputFromNC(n) (inputt) Array indicating whether variable is from NC (0 for TS, 1 for NC, 2 for value, 3 for not provided)    
        ! MaxNumofFile (input) is the maximum number of NC files for any variable 
        ! NumNCFiles(n) (input)  Array giving the number of NC files for NC variables
        ! NCDFContainer(MaxNF, n) (input)  Array giving file names for variables with NC files
        ! varnameinncdf(n) (input)  variable name in control file specifying the netcdf variable to use for reading from NC
        ! iycoord (input) Y-position of a grid for which values for entire time-series is obtianed
        ! jxcoord (input) X-position of a grid for which values for entire time-series is obtianed
        ! NCfileNumtimesteps(MaxNF,n) (input) array giving the number of time steps in each NC file
        ! arrayx  (input).  Maximum across variables of the sum of number of time steps in all NC files for that variable
        ! NoofTS(n) (oinput).  The number of combined NC time steps for each variable
        ! TSV (arrayx,11) (output) holds all the timsteps bot from NC and time series (TS) text files. Time from TS files are stored as
        !                          julian and Time from NC files are stored as day/hour from the reference date (time unit in NC).
        ! Allvalues (arrayx,11) (output) holds values of enite timeseries for each variable for a particular grid  point
        
        Parameter(n=11)
        integer:: MaxNumofFile,NOofTS(n),NUMNCFILES(n)
        Character*200:: NCDFContainer(MaxNumofFile,n)
        Character*50:: varnameinncdf(n)
        integer::iycoord,jxcoord,NCfileNumtimesteps(MaxNumofFile,n)
        Double precision:: TimeMinPerFile(MaxNumofFile,n),TimeMaxPerFile(MaxNumofFile,n)
        character*50:: var_name
        Character*200::File_nameM
        integer:: arrayx,Dt
        integer:: ArrayStart,ArrayEnd,IsInputFromNC(n),rec
        Real:: Allvalues(arrayx,n)
        Real, allocatable :: AllVal(:)
        Double precision:: TSV(arrayx,n)
        REAL:: VarMissingValues(MaxNumofFile,n)
        Do i = 1,n 
            ArrayEnd=0
            if (IsInputFromNC(i) .eq. 1)then
                Do k=1,NumNCFiles(i)
                    File_nameM=NCDFContainer(k,i)
                    var_name=varnameinncdf(i)
                    Allocate(AllVal(NCfileNumtimesteps(k,i)))
                    rec=NCfileNumtimesteps(k,i) 
                    CALL NCDFReadAllTS(file_nameM,Var_name,AllVal,iycoord,jxcoord,rec)
                    ArrayStart=ArrayEnd+1                      
                    ArrayEnd=rec+ArrayStart-1
                    Allvalues(ArrayStart:ArrayEnd,i)=AllVal(1:rec)
                    Deallocate(AllVal)
                    Do jj=1,NoofTS(i)
                        if (Allvalues(jj,i) .EQ. VarMissingValues(k,i))Then
                            Allvalues(jj,i)=Allvalues((jj-1),i)
                        End if
                    End do
                End do
            End if
        End do
        
        End Subroutine
        
        Subroutine InputVariableValue(INPUTVARNAME,IsInputFromNC,NoofTS,TSV,Allvalues,arrayx,ModelEndDate,ModelStartDate,ModelStartHour,&
        &ModelEndHour,nrefyr,nrefmo,nrefday,Year,Month,Day,Hour,ModelDt,InpVals,CurrentArrayPos,CurrentModelDT,istep)
        
        ! inputvarname (n) (input) Name of the variables that are provided inside inputcontrol.dat file
        ! IsInputFromNC(n) (inputt) Array indicating whether variable is from NC (0 for TS, 1 for NC, 2 for value, 3 for not provided)   
        ! NoofTS(n) (input).  The number of combined NC time steps for each variable
        ! TSV (arrayx,11) (input) holds all the timsteps bot from NC and time series (TS) text files. Time from TS files are stored as
        !                         julian and Time from NC files are stored as day/hour from the reference date (time unit in NC).
        ! Allvalues (arrayx,11) (input) holds values of enite timeseries for each variable for a particular grid  point
        ! arrayx  (input).  Maximum across variables of the sum of number of time steps in all NC files for that variable
        ! ModelEndDate(3) (input)  Array giving end year, month, day
        ! ModelStartDate(3) (input) Array giving start year, month, day
        ! ModelStartHour (input) Start hour
        ! ModelEndHour (input)  end hour
        ! nrefyr,nrefmo,nrefday,(input) Reference year, month, day from netcdf time units specified as "days from y/m/d"
        !                                Our convention is that our start hour must be 0
        ! Year,Month,Day,Hour (input) year, month, day and hour after updating model time at each time step
        ! ModelDt (input) model timestep/resolution
        ! InpVals(n) (output).  Variable to hold the current value of each input variable. 
        ! CurrentArrayPos (n) (outnput) current postion in the array of TSV and Allvalues array
        ! CurrentModelDT (input) current model timestep in julian
        ! istep (input) Current model time position

        Parameter(n=11)
        Character*200:: INPUTVARNAME(n)
        integer::arrayx,ModelEndDate(3),ModelStartDate(3),Year,Month,Day,ArrayPos(n),istep
        REAL:: Allvalues(Arrayx,n),ModelStartHour,ModelEndHour,modeldt,InpVals(n),dt
        Double precision:: TSV(Arrayx,n)
        REAL:: Hour
        Double precision:: CurrentModelDT,RefJD
        Double precision:: FileCurrentDT(n),FileNextDT(n),SJD,EJD,RefHour,SHOUR
        Integer:: CurrentArrayPos(n),NextArrayPos(n),IsInputFromNC(n),nrefyr,nrefmo,nrefday
        integer:: NOOFTS(n)
        Double precision:: Tol
        RefHour=0.00
        SHOUR=dble(ModelStartHour)
        TOl=1./(60.*24.) ! 1 minute tolerance
        call JULDAT(ModelStartDate(1),ModelStartDate(2),ModelStartDate(3),SHour,SJD)
        call JULDAT(nrefyr,nrefmo,nrefday,RefHour,RefJD)
        If (inputvarname(5) .NE. 'tmin')Then
            IsInputFromNC(5)=1
        End if
        If (inputvarname(6) .NE. 'Tmax')Then
            IsInputFromNC(6)=1
        End if
        Do i=1,n
            if(istep==1)Then
                CurrentArrayPos(i)=1
                NextArrayPos(i)=CurrentArrayPos(i)+1
            Else
                If(CurrentArrayPos(i) .LT. NoofTS(i))Then
                    CurrentArrayPos(i)=CurrentArrayPos(i)
                    NextArrayPos(i)=CurrentArrayPos(i)+1
                Else
                    CurrentArrayPos(i)=NoofTS(i)-1
                    NextArrayPos(i)=CurrentArrayPos(i)+1
                End If
            End if
            If(IsInputFromNC(i) .Eq. 1)then
                FileCurrentDT(i)=dble(REfJd+TSV(CurrentArrayPos(i),i))
                FileNextDT(i)=dble(RefJD+TSV(NextArrayPos(i),i))
            End if
            If(IsInputFromNC(i) .Eq. 0)then
                FileCurrentDT(i)=TSV(CurrentArrayPos(i),i)
                FileNextDT(i)=TSV(NextArrayPos(i),i)
            End if
        End Do
        
        Do i=1,n
            IF(IsInputFromNC(i) .LT. 2)then
                If((FileCurrentDT(i) .GT. CurrentModelDT) .AND. (CurrentArrayPos(i)==1))THEN
                    InpVals(i)=Allvalues(1,i)
                    CurrentArrayPos(i)=1
                End IF
                
1234            If((FileCurrentDT(i)+Tol) .LT. CurrentModelDT)THEN 
                    FileCurrentDT(i)=FileNextDT(i)
                    CurrentArrayPos(i)=NextArrayPos(i)
                    If(CurrentArrayPos(i) .LT. NoofTS(i))Then
                        NextArrayPos(i)=NextArrayPos(i)+1
                        If(IsInputFromNC(i) .Eq. 1)then
                            FileNextDT(i)=dble(RefJD+TSV(NextArrayPos(i),i))
                        End if
                        If(IsInputFromNC(i) .Eq. 0)then
                            FileNextDT(i)=TSV(NextArrayPos(i),i)
                        End if
                        Go to 1234
                    Else
                        CurrentArrayPos(i)=NoofTS(i)
                        InpVals(i)=Allvalues(CurrentArrayPos(i),i)   
                    End if
                    
                Else if(((FileCurrentDT(i)-Tol) .LT. CurrentModelDT) .AND. ((FileCurrentDT(i)+Tol) .GT. CurrentModelDT))THEN
                    If(CurrentArrayPos(i) .LT. NoofTS(i))Then
                        InpVals(i)=Allvalues(CurrentArrayPos(i),i)
                        CurrentArrayPos(i)=CurrentArrayPos(i)+1
                    Else
                        CurrentArrayPos(i)=NoofTS(i)
                        InpVals(i)=Allvalues(CurrentArrayPos(i),i)
                    End If
                    
                Else If (((FileCurrentDT(i)-Tol) .GT. CurrentModelDT) .AND. (CurrentArrayPos(i).NE. 1))THEN
                    InpVals(i)=Allvalues(CurrentArrayPos(i),i)
                End If   
                
            End If
        End Do
        
        If (inputvarname(5) .NE. 'tmin')Then
            IsInputFromNC(5)=0
        End if
        If (inputvarname(6) .NE. 'Tmax')Then
            IsInputFromNC(6)=0
        End if
        End subroutine