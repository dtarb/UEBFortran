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
   implicit none
 
        integer :: n, i, xx, MaxNumofFile, count
        ! inputcon (input) is name of control file
        ! MaxNumofFile (input) is the maximum number of NC files for any variable 
        ! inputvarname (output) Name of the variables that are provided inside inputcontrol.dat file
        ! UTCOffSet(output) UTC offset for the region/wateshed we are dealing with
        
        PARAMETER(n=11)
        ! TA,P,V,RH,trange,QSIOBS,QNETOB,qg,qli,inputcon, iloop
        integer:: reason,InumOfFile(n),IsInputFromNC(n)
        integer:: Syear, Smonth, sdate,  Eyear, Emonth, Edate !s for tartng date-time and e for ending date-time
        real::  shour, Ehour, dt,VARVALUES(n),UTCOffSet     !dt=time increment in hours   
        CHARACTER*200 inputHeading, inputcode, inputname, inputcon, inputVName(n), InputNCFilename(n), InputTSFilename(n)
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
500       Read(19,*,iostat=reason, end=600)inputcode
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
        &nrefyr,nrefmo,nrefday,varnameinncdf,arrayx,NoofTS,InpVals,VarMissingValues,VarfILLValues)
        ! inputcon (input) is name of control file
        ! MaxNumofFile (input) is the maximum number of NC files for any variable 
        ! IsInputFromNC(n) (Output) Array indicating whether variable is from NC (0 for TS, 1 for NC, 2 for value, 3 for not provided)
        ! NumNCFiles(n) (Output)  Array giving the number of NC files for NC variables
        ! InputTSFilename(n) (Output) Array giving file names for variables with time series input
        ! NCDFContainer(MaxNC, n) (Output)  Array giving file names for variables with NC files
        ! ModelStartDate(3) (Output) Array giving start year, month, day
        ! ModelStartHour (output) Start hour
        ! ModelEndDate(3) (output)  Array giving end year, month, day
        ! ModelEndHour (output)  end hour
        ! Modeldt (output)  time step
        ! NCfileNumtimesteps(MaxNC,n) (output) array giving the number of time steps in each NC file
        ! nrefyr,nrefmo,nrefday,(output) Reference year, month, day from netcdf time units specified as "days from y/m/d"
        !  Our convention is that our start hour must be 0
        ! varnameinncdf(n) (output)  variable name in control file specifying the netcdf variable to use for reading from NC
        ! arrayx  (output).  Maximum across variables of the sum of number of time steps in all NC files for that variable
        ! NoofTS(n) (output).  The number of combined NC time steps for each variable
        ! InpVals(n) (output).  Variable to hold the current value of each input variable.  This subroutine fills this for inputs that are constant
        ! VarMissingValues (MaxNC,n) (out) contains the missing values in each netCDF
        ! VarfILLValues (MaxNC,n) (out) contains the filling values in each netCDF
        Use netCDF
        Implicit None
        integer :: n, i, ii, iii, xx, k, kflag, ncidout, NumtimeStepEachNC, MaxNumofFile, InputVarId
        parameter(n=11)                                         !n is  loop variable
        integer:: IsInputFromNC(n), reason, NumNCFiles(n),nrefyr,nrefmo,nrefday,count
        integer::syear,smonth,sday,first
        integer:: ModelStartDate(3),  ModelEndDate(3)  
        real:: ModelStartHour, ModelEndHour,Modeldt,InpVals(n)  
        integer:: NCfileNumtimesteps(MaxNumofFile,n)
        integer:: NOofTS(n)
        real:: UTCOffSet
        character*200:: NCDFContainer(MaxNumofFile,n)
        !  Arrays to hold temporary information before sort
        character*200, allocatable :: tempfilelist(:)
        integer, allocatable :: tempfilesteps(:)
        double precision, allocatable :: tempstarttime(:)
        integer,allocatable :: iy(:)           
        CHARACTER*200 inputHeading, inputcode, inputname, inputcon, inputVName(n), InputNCFilename(n)
        CHARACTER*200 InputTSFilename(n)
        CHARACTER*200 NCfileContain 
        Character*50::varnameinncdf(n)
        Character*200:: Rec_name,File_name
        integer:: VarID
        Integer:: FileNextT(3)
        REAL:: FileNextH, FileNextVal
        integer:: FileOpenFlag(n)
        Character*200:: TSFile, CurrentInputVariable(n)
        integer::arrayx
        character (len = *), parameter :: missing_value = "missing_value",fillvalue='_FillValue'
        real:: VarMissingValues(MaxNumofFile,n),VarfILLValues(MaxNumofFile,n)
        integer::  numAtts
        character (len = 50):: AttName

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
        CALL nCDF3DtimeRead (NCfileContain,1,tempstarttime(count),tempfilesteps(count),syear,smonth,sday)
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
        call ssort(tempstarttime,iy,count)
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
                    call check(nf90_inquire_dimension(ncidout, Varid, Rec_name,NumtimeStepEachNC))  ! information about dimensionID 3
                    Call check(nf90_inq_varid(ncidout,varnameinncdf(i),InputVarId))
                    CALL check(nf90_inquire_variable(ncidout,InputVarId,natts = numAtts))
                    DO iii=1,numAtts
                        CALL check(nf90_inq_attname(ncidout,InputVarId,iii,AttName))
                        if(AttName .eq. missing_value)THEN
                            CALL check(nf90_get_att(ncidout,InputVarId,missing_value,VarMissingValues(k,i))) 
                        ELSE
                            VarMissingValues(k,i)=-9999
                        end IF
                        if(AttName .eq. fillvalue)THEN
                            CALL check(nf90_get_att(ncidout,InputVarId,fillvalue,VarfILLValues(k,i)))  
                        ELSE
                            VarfILLValues(k,i)=-9999
                        end IF
                    END DO
                    NOofTS(i)=NOofTS(i)+NumtimeStepEachNC
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
        
        subroutine timeSeriesAndtimeSteps(MaxNumofFile,NUMNCFILES,IsInputFromNC,InputTSFilename,NCDFContainer,&
        arrayx,NOofTS,TSV,Allvalues)
        
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
        ! TSV (arrayx,11) (output) holds all the timsteps bot from NC and time series (TS) text files. time from TS files are stored as
        !                          julian and time from NC files are stored as day/hour from the reference date (time unit in NC).
        ! Allvalues (arrayx,11) (output) holds values of enite timeseries for each variable for a particular grid  point
        use netcdf
        Implicit None
        integer :: n, i, k, FileCount, NumtimeStepEachNC
        Parameter(n=11)
        integer:: arrayx
        integer:: VarID,ncidout
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
        Double precision,allocatable ::AlltimeSteps(:)
        character*200:: CurrentInputVariable(n)
        integer::StartY,StartM,StartD
        Double precision:: NCRfeferenceJDT
        double precision:: RefHour
        Integer::TotalTS
        Double precision:: TVal(1)
        integer:: ArrayStart,ArrayEnd
        
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
                    call check(nf90_inquire_dimension(ncidout,Varid,Rec_name,NumtimeStepEachNC))  ! information about dimensionID 3  
                    CALL nCDF3DtimeRead(file_nameM,1,TVal,TotalTS,StartY,StartM,StartD)
                    CALL JULDAT(StartY,StartM,StartD,RefHour,NCRfeferenceJDT)
                    Allocate(AlltimeSteps(numtimeStepEachNC))                                    
                    call check(nf90_get_var(ncidout,VarId,AlltimeSteps))                            ! Read the first time value from the file
                    ArrayStart=ArrayEnd+1                           
                    ArrayEnd=NumtimeStepEachNC+ArrayStart-1
                    TSV(ArrayStart:ArrayEnd,i)=AlltimeSteps
                    deallocate(AlltimeSteps)
                end do
            End if
        End do
        return
        end subroutine
        

        
        Subroutine Values4VareachGrid(inputvarname,IsInputFromNC,MaxNumofFile,NUMNCFILES,NCDFContainer,varnameinncdf,iycoord,jxcoord,&
        &NCfileNumtimesteps,NOofTS,arrayx,Allvalues,VarMissingValues,VarfILLValues,StepInADay,NumtimeStep,CurrentArrayPosRegrid,ReGriddedArray)
            
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
        ! TSV (arrayx,11) (output) holds all the timsteps bot from NC and time series (TS) text files. time from TS files are stored as
        !                          julian and time from NC files are stored as day/hour from the reference date (time unit in NC).
        ! Allvalues (arrayx,11) (output) holds values of enite timeseries for each variable for a particular grid  point
        ! VarMissingValues (MaxNC,n) (out) contains the missing values in each netCDF
        ! VarfILLValues (MaxNC,n) (out) contains the filling values in each netCDF
        ! ReGriddedArray(arrayx,n) (out) stored all regridded values
        ! NumtimeStep (in)
        ! CurrentArrayPosRegrid(NumtimeStep,n) (in)
        
        Implicit None
        integer :: n,i,j,jj,k
        Parameter(n=11)
        integer:: MaxNumofFile,NOofTS(n),NUMNCFILES(n),NumtimeStep,StepInADay
        integer:: CurrentArrayPosRegrid(NumtimeStep,n)
        Character*200:: NCDFContainer(MaxNumofFile,n),inputvarname(n)
        Character*50:: varnameinncdf(n)
        integer::iycoord,jxcoord,NCfileNumtimesteps(MaxNumofFile,n)
        character*50:: var_name
        Character*200::File_nameM
        integer:: arrayx
        integer:: ArrayStart,ArrayEnd,IsInputFromNC(n),rec
        Real:: Allvalues(arrayx,n)
        REAL:: ReGriddedArray(NumtimeStep,n)
        Real, allocatable :: AllVal(:)
        REAL:: VarMissingValues(MaxNumofFile,n),VarfILLValues(MaxNumofFile,n)
        Integer:: Modx,TOTALDAYMONE,TOTALDAY,GG
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
                        if (Allvalues(jj,i) .EQ. VarFillValues(k,i))Then
                            Allvalues(jj,i)=Allvalues((jj-1),i)
                        End if
                    End do
                End do
            End if
        End do
        
        Do i = 1,n 
            if (IsInputFromNC(i) .LT. 2)then
                Do j = 1,NumtimeStep
                    ReGriddedArray(j,i)=Allvalues(CurrentArrayPosRegrid(j,i),i)
                end do
            END IF
        end do

        If (inputvarname(5) .NE. 'tmin')THEN
            CurrentArrayPosRegrid(1:NumtimeStep,5)=CurrentArrayPosRegrid(1:NumtimeStep,1)
            modx=MOD(NumtimeStep,StepInADay)
            if (modx .GT. 0)THEN
                totalDayMOne=(NumtimeStep-modx)/StepInADay
                totalDay=totalDayMOne+1
                Do gg=1,totaldayMOne
                    ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,5)=MINVAL(ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
                    ReGriddedArray((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),5)&
                    &=MINVAL(ReGriddedArray((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),1))
            Else
                totaldayMOne=(NumtimeStep)/StepInADay
                totalDay=totalDayMOne
                Do gg=1,totalDay
                    ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,5)=MINVAL(ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
            End if
        END IF
       
        If (inputvarname(6) .NE. 'tmax')THEN
            modx=MOD(NumtimeStep,StepInADay)
            if (modx .GT. 0)THEN
                totaldayMOne=(NumtimeStep-modx)/StepInADay
                totalDay=totalDayMOne+1
                Do gg=1,totaldayMOne
                    ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,6)=MAXVAL(ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
                    ReGriddedArray((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),6)=&
                    &MAXVAL(ReGriddedArray((totalDayMOne*StepInADay+1):((totalDayMOne*StepInADay)+modx),1))
            Else
                totaldayMOne=(NumtimeStep)/StepInADay
                totalDay=totalDayMOne
                Do gg=1,totalDay
                    ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,6)=MAXVAL(ReGriddedArray(((gg-1)*StepInADay+1):gg*StepInADay,1))
                END DO
            End if
        END IF
        End Subroutine
        
        Subroutine InputVariableValue(INPUTVARNAME,IsInputFromNC,NoofTS,TSV,Allvalues,arrayx,ModelStartDate,ModelStartHour,&
        ModelEndDate,ModelEndHour,nrefyr,nrefmo,nrefday,modeldt,NumtimeStep,CurrentArrayPosRegrid,modelTimeJDT)
        Implicit None
        
        ! inputvarname (n) (input) Name of the variables that are provided inside inputcontrol.dat file
        ! IsInputFromNC(n) (inputt) Array indicating whether variable is from NC (0 for TS, 1 for NC, 2 for value, 3 for not provided)   
        ! NoofTS(n) (input).  The number of combined NC time steps for each variable
        ! TSV (arrayx,11) (input) holds all the timsteps bot from NC and time series (TS) text files. time from TS files are stored as
        !                         julian and time from NC files are stored as day/hour from the reference date (time unit in NC).
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
        ! CurrentArrayPosRegrid (n) (outnput) Indexing array postions by comparing with current model time of TSV and Allvalues array
        ! CurrentModelDT (input) current model timestep in julian
        ! istep (input) Current model time position
        
        integer :: n,i
        Parameter(n=11)
        integer:: NumtimeStep
        integer:: CurrentArrayPosRegrid(NumtimeStep,n)
        Character*200:: INPUTVARNAME(n)
        integer::arrayx,ModelStartDate(3),istep,ModelEndDate(3)
        REAL:: Allvalues(Arrayx,n),ModelStartHour,ModelEndHour
        Double precision:: TSV(Arrayx,n)
        Double precision:: CurrentModelDT,RefJD
        Double precision:: FileCurrentDT(n),FileNextDT(n),SJD,RefHour,SHOUR
        Integer:: CurrentArrayPos(n),NextArrayPos(n),IsInputFromNC(n),nrefyr,nrefmo,nrefday
        Integer:: YEAR,MONTH,DAY
        REAL:: HOUR,MODELDT,DT,tol
        integer:: NOOFTS(n)
        Double precision:: dhour,DBLEHOUR,EJD,modelTimeJDT(NumtimeStep)
        RefHour=0.00
        SHOUR=dble(ModelStartHour)
        TOl=5./(60.*24.) ! 1 minute tolerance
        call JULDAT(ModelStartDate(1),ModelStartDate(2),ModelStartDate(3),SHour,SJD)
        call JULDAT(ModelStartDate(1),ModelStartDate(2),ModelStartDate(3),SHour,CurrentModelDT)
        call JULDAT(nrefyr,nrefmo,nrefday,RefHour,RefJD)
        dhour=dble(ModelEndHour)
        call JULDAT(ModelEndDate(1),ModelEndDate(2),ModelEndDate(3),dhour,EJD)
        tol=5.0/(24*60)
        dt=Modeldt
        YEAR=ModelStartDate(1)
        MONTH=ModelStartDate(2)
        DAY=ModelStartDate(3)
        Hour=REAL(SHour)
        istep=0
        
8887   istep=istep+1 
        Do i=1,n
            if(istep==1)Then
                CurrentArrayPos(i)=1
                NextArrayPos(i)=CurrentArrayPos(i)+1
            Else
                If(CurrentArrayPos(i) .LT. NoofTS(i))Then
                    CurrentArrayPos(i)=CurrentArrayPos(i)
                    CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
                    NextArrayPos(i)=CurrentArrayPos(i)+1
                Else
                    CurrentArrayPos(i)=NoofTS(i)-1
                    CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
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
                    CurrentArrayPos(i)=1
                    CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
                End IF
                
1234            If((FileCurrentDT(i)+Tol) .LT. CurrentModelDT)THEN 
                    FileCurrentDT(i)=FileNextDT(i)
                    CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
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
                        CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
                    End if
                    
                Else if(((FileCurrentDT(i)-Tol) .LT. CurrentModelDT) .AND. ((FileCurrentDT(i)+Tol) .GT. CurrentModelDT))THEN
                    If(CurrentArrayPos(i) .LT. NoofTS(i))Then
                        CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
                        CurrentArrayPos(i)=CurrentArrayPos(i)+1
                        If (CurrentArrayPos(i)== NoofTS(i))THEN
                            CurrentArrayPos(i)=CurrentArrayPos(i)
                        END IF
                    Else
                        CurrentArrayPos(i)=NoofTS(i)
                        CurrentArrayPosRegrid(istep,i)=CurrentArrayPos(i)
                    End If
                End If   
            End If
        End Do
        
        CALL UPDATEtime(YEAR,MONTH,DAY,HOUR,DT)
        DBLEHOUR=DBLE(hour)
        call JULDAT(YEAR,MONTH,DAY,DBLEHOUR,CurrentModelDT)
        If (EJD .GE. CurrentModelDT)Then  
            modelTimeJDT(istep)=CurrentModelDT    
            Go to 8887 
        End if

        End subroutine