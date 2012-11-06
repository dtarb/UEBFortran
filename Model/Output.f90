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

!   This subroutine evaluates based on the grid dimensions and begin and end date the 
!   total number of time steps in the model run and the number of netcdf files to be 
!   created for each output variable based on an assumed maximum number of variables per netcdf file
!   It also parses the Output control file to identify how many point detail outputs there are so that 
!   arrays can be allocated
        subroutine NumOutFiles(OutControlFILE, ModelStartDate, &
        ModelStartHour,ModelEndDate,ModelEndHour,Modeldt, &
        dimlen2,dimlen1,NumtimeStep,NumofFile,NumOutPoint,OutCount)

            
        ! NumtimeStep is the total number of time steps output
        ! NumofFile is the number of output files per variable output to netCDF
        parameter (n=66)
        integer:: dimlen2,dimlen1,i
        Integer:: ModelStartDate(3),ModelEndDate(3)
        Real:: ModelStartHour,ModelEndHour,Modeldt
        Double precision:: JMSD,JMED,MStartHour,MEndHour,tol
        integer:: ntsperfile
        real xx
        integer:: NumtimeStep, NumofFile,NumOutPoint,reason,OutCount
        character*200 trimcode,filecode,xxs,OutControlFILE
        CHARACTER(200), DIMENSION(n) :: outSymbol
        integer:: netCDFDataValPerfile 
        
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

        MStartHour=dble(ModelStartHour)
        MEndHour=dble(ModelEndHour)
        Call JULDAT(ModelStartDate(1),ModelStartDate(2), &
        ModelStartDate(3),MStartHour,JMSD) !JMSD= julian model start date-time
        Call JULDAT(ModelEndDate(1),ModelEndDate(2),ModelEndDate(3) &
        ,MEndHour,JMED) !JMED= julian model start date-time
     
        netCDFDataValPerfile=int(1.5*1000*1000*1000/5.5) !  Determined experimentally
        !netCDFDataValPerfile=1.5*1000*50/5.5   !   Small for debugging
        tol=5./(60.*24.)     !  5 min tolerance
        NumtimeStep=int(((JMED+tol-JMSD)/(Modeldt/24)) + 1)
        ntsperfile=max(netCDFDataValPerfile/(dimlen1*dimlen2),1)
        NumofFile=NumtimeStep/ntsperfile+1  ! The +1 is to round up integer calculation
        NumOutPoint=0
        OutCount=0
        OPEN(109,FILE=OutControlFile,STATUS='OLD', ACTION='READ')
        READ(109,*) xxs
        ! Read until end of file
3200      Read(109,*,iostat=reason, end=3300)filecode
            If(reason .eq. 0) then
                trimcode=ADJUSTL(filecode(1:(SCAN(filecode, ':')-1)))
                Call lowercase(trimcode,trimcode)
                if(trimcode .eq. 'pointdetail')then
                    read(109,*,end=3300)xx
                    read(109,*,end=3300)xxs
                    NumOutPoint=NumOutPoint+1
                else 
                    do i=1,n
                    Call lowercase(outSymbol(i),outSymbol(i))
                    if(trimcode .eq. trim(outSymbol(i)))then
                       OutCount=OutCount+1
                       exit
                    endif
                    enddo
                    read(109,*,end=3300)xxs
                endif
                go to 3200
            endif
 3300     close(109)       
        return
        end subroutine NumOutFiles
!===End of reading the values from outconrol.dat file ==

!==This subroutine reads the outputcontrol file and returns
!  The folder path to each variable
!  The file name to each variable that is outpus
!  An array indicating whether each variable is output or not 

        subroutine OutputFiles(OutControlFILE,NumtimeStep,Dimlen2, &
       dimlen1,NumofFile,outSampleFile,NumtimeStepPerFile, &
       OutVar,OutPoint,OutPointFiles,NumOutPoint,OutCount)
        ! OutNCfiles output a array that contains all the putput NC files
        ! outputvarfile input a file that lists all the variables that a user wants as outputs
        implicit none
        integer n,i,netCDFDataValPerfile,ntsperfile,lastfilents
        parameter(n = 66)
        integer:: OutCount
        integer:: dimlen2,dimlen1,NumofFile
        CHARACTER*200 outname,OutControlFILE,outHeading,filecode
        CHARACTER(200), DIMENSION(n) :: outSymbol ! unused: ,OutUnits
        CHARACTER(200),DIMENSION(OutCount)::outSampleFile
        integer OutVar(Outcount)
        CHARACTER(200) :: outFolderName
        ! unused: CHARACTER(200) ::string
        Integer:: reason
        integer NumOutPoint
        integer::NumtimeStep
        Integer:: NumtimeStepPerFile(NumofFile)
        integer OutPoint(NumOutpoint,2),numOutPtRead
        logical matched
        integer:: OutNum
        character*200 OutPointFiles(NumOutPoint)

        ! FIXME: this is duplicated in several places and can simply
        ! be included from a globas or def file.

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
     
        if( NumofFile .gt. 1)then
            netCDFDataValPerfile=int(1.5*1000*1000*1000/5.5) !  Determined experimentally
            !netCDFDataValPerfile=1.5*1000*50/5.5   !   Small for debugging
            ntsperfile=max(netCDFDataValPerfile/(dimlen1*dimlen2),1)
            lastfilents=NumtimeStep-(NumofFile-1)*ntsperfile
            if(lastfilents .le. 0) &
             write(6,*)'Logic error in computing time steps per file'
            NumtimeStepPerFile(1:(NumofFile-1))=ntsperfile
            NumtimeStepPerFile(NumofFile)=lastfilents
        else
            NumtimeStepPerFile(1)=NumtimeStep
        endif
        numOutPtRead=0
        OPEN(109,FILE=OutControlFile,STATUS='OLD', ACTION='READ')
        READ(109,*) outHeading
        OutNum=0
        ! Read until end of file
4500      Read(109,*,iostat=reason, end=4600)filecode
            If(reason .eq. 0) then
            matched = .false.
            OUTname=ADJUSTL(filecode(1:(SCAN(filecode, ':')-1)))
            CALL lowercase(OutName,OutName)
             do i=1,n,1
                CALL lowercase(outSymbol(i),outSymbol(i))
                 If(OutNum .le. OutCount)then
               if(OutName .eq. trim(outSymbol(i)))then
                    OutNum=OutNum+1
                    OutVar(OutNum)=i  !OutName
                    
                    Read(109,*)outFolderName
                    outSampleFile(OutNum)=ADJUSTL(outFolderName)
                    matched = .true.
                    exit
                    
                else if(OutName .eq. 'pointdetail')then
                    numOutPtRead=numOutPtRead+1
                    read(109,*)OutPoint(numOutPtRead,1),  &
                     OutPoint(numOutPtRead,2)
                    read(109,*)OutPointFiles(numOutPtRead)
                    matched = .true.
                    exit
                end if
                end if
            end do
            if(.not. matched) &
             write(6,*)'Output code not matched:  ',OUTname
        end if
        
        go to 4500
4600    CLOSE(109)
        return
        end subroutine OutputFiles
        !===End of reading the values from outconrol.dat file ==

 Subroutine DirectoryCreate(nrefyr,nrefmo,nrefday,dimlen1,dimlen2,DimName1,DimName2,DimUnit1,&
        &DimUnit2,NumofFile,outcount,Outvar,&
        &NumtimeStepPerFile,outSampleFile,OutputNCContainer,NCIDARRAY)
        use netCDF
        implicit none
        integer n,i,j,LengthOutFile
        parameter(n=66)
        integer:: outcount
        ! unused: character*200::outputfolder(outcount)
        ! unused: character*200::DIRECTORY
        character*200:: outSymbol(n)
        integer::outvar(outcount)
        character*200::OutUnits(n)
        ! unused: character*50:: makedirectory ,removedirectory
        ! unused: LOGICAL:: exist3
        integer:: NumofFile
        Character(200)::OutputNCContainer(NumofFile,outcount)
        ! unused: integer:: lengthout
        character(len=10) :: IntToChar
        character(len=50) ::outputfile,NextFileName
        integer, parameter :: NDIMS=3  ! x-coord, y-coord and time-coord
        integer:: nrefyr,nrefmo,nrefday
        character (4):: refyear
        character (2)::refmonth,refday
        
        Integer:: dimlen1,dimlen2
        ! unused: character(len=500):: netCDFfileName(NumofFile,outcount)
        ! unused: character(len=500):: NCOutfileArr(NumofFile,outcount)
        Integer:: NumtimeStepPerFile(NumofFile)
        character (20) :: DimName1,DimName2
        character (100) :: DimUnit1,DimUnit2
        character (200) :: time_unit,NCOutfileName

        Integer:: ncid
        integer :: dimids(NDIMS)
        integer::x_dimid,y_dimid,time_dimid,x_varid,y_varid,time_varid
        integer:: Varid4
        ! unused: character (50) :: watershedFile
        ! unused: integer:: openunit
        character (len = *), parameter :: UNITS = "units"
        character (len = *), parameter :: missing_value = "missing_value"
        integer:: NCIDARRAY(NumofFile,outcount)
        character (200) :: outSampleFile(outcount)
        Real:: MissingValues
        write(refyear,'(i4)') nrefyr
        write(refmonth,'(i2.2)') nrefmo
        write(refday,'(i2.2)') nrefday
        MissingValues=-9999.000
        time_unit='days since'//' '//trim(refyear)//'-'//trim(refmonth) &
       //'-'//trim(refday)//'T'//'00:00'
       
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

        OutUnits = (/ "unitless","unitless","unitless","unitless", &
           "deg C   ","m/hr    ","m/s     ","unitless","kJ/m2/hr", &
           "kJ/m2/hr","kJ/m2/hr","Zen(deg)","kJ/m2   ","m       ", &
           "unitless","m/hr    ","m/hr    ","unitless","kJ/m2/hr", &
           "kJ/m2/hr","m       ","m/hr    ","kJ/m2/hr","kJ/m2/hr", &
           "m/hr    ","deg C   ","deg C   ","m       ","m       ", &
           "m       ","kJ/m2/hr","m/hr    ","m       ","m       ", &
           "unitless","deg C   ","deg C   ","kJ/m2/hr","kJ/m2/hr", &
           "deg C   ","deg C   ","kJ/m2/hr","kJ/m2/hr","kJ/m2/hr", &
           "kJ/m2/hr","m/s     ","unitless","m       ","m       ", &
           "m/hr    ","m/hr    ","m       ","m       ","deg C   ", &
           "deg C   ","kJ/m2/hr","kJ/m2/hr","m/hr    ","kJ/m2/hr", &
           "kJ/m2/hr","m/hr    ","m/hr    ","m       ","m/hr    ", &
           "m/hr    ","m/hr    "/)
      

        do i=1,outcount
            CALL lowercase(outSymbol(i),outSymbol(i))
            CALL lowercase(OutUnits(i),OutUnits(i))
            OutputNCContainer(1,i)=outSampleFile(i)
            do j=2,NumofFile
            !assumption=given netCDF file names cannot ave a "." withn filename.
            !           for eaxmple: it cannot be Qsi.0001.nc
            !                        it can be Qsi_001.nc or Qsi_0001
            
                LengthOutFile = SCAN(outSampleFile(i),'.nc',BACK=.TRUE.)
                If (LengthOutFile .ge. 0)then
                    NextFileName=outSampleFile(i)(1: &
                   (SCAN(outSampleFile(i),'.',BACK=.TRUE.)-1))
                    write(IntToChar,'(i4.4)')j
                    NextFileName=NextFileName
                    outputfile=trim(NextFileName)//trim(IntToChar)//'.nc'
                    OutputNCContainer(j,i)=outputfile
                else
                    outputfile=trim(NextFileName)//trim(IntToChar)
                    OutputNCContainer(j,i)=outputfile
                end if
            end do
            
            do j=1,NumofFile
                NCOutfileName=trim(OutputNCContainer(j,i))
                ! Create the file.
                call check(nf90_create(NCOutfileName,nf90_clobber,ncid))
               ! Define the dimensions. The record dimension is defined to have
               ! unlimited length - it can grow as needed. In this example it is
               ! the time dimension.
               !               x_varid=1
               !               y_varid=2
               !               x_dimid=1
               !               y_dimid=2
               !               time_dimid=3
               call check(nf90_def_dim(ncid,'time',NumtimeStepPerFile(j),time_dimid))
               call check(nf90_def_dim(ncid,DimName2,dimlen2,x_dimid))
               call check(nf90_def_dim(ncid,DimName1,dimlen1,y_dimid))
               
               
               ! Define the coordinate variables. We will only define coordinate
               ! variables for lat and lon.  Ordinarily we would need to provide
               ! an array of dimension IDs for each variable's dimensions, but
               ! since coordinate variables only have one dimension, we can
               ! simply provide the address of that dimension ID (lat_dimid) and
               ! similarly for (lon_dimid).
               call check(nf90_def_var(ncid,"time",NF90_REAL,time_dimid, &
                    time_varid))
               call check(nf90_def_var(ncid,DimName2,NF90_REAL,x_dimid, &
                    x_varid))
               call check(nf90_def_var(ncid,DimName1,NF90_REAL,y_dimid, &
                    y_varid))


!                ! Assign units attributes to coordinate variables.
               call check(nf90_put_att(ncid,y_varid,UNITS,DimUnit1))
               call check(nf90_put_att(ncid,x_varid,UNITS,DimUnit2))
               call check(nf90_put_att(ncid,time_varid,UNITS,time_unit))
                ! The dimids array is used to pass the dimids of the dimensions of
                ! the netCDF variables. Both of the netCDF variables we are creating
                ! share the same four dimensions. In Fortran, the unlimited
                ! dimension must come last on the list of dimids.
               dimids = (/time_dimid,x_dimid,y_dimid/)               
               ! Define the netCDF variables for the pressure and temperature data.
               call check( nf90_def_var(ncid,trim(OutSymbol(outvar(i))), &
                    NF90_REAL, &
                    dimids,Varid4))
                ! Assign units attributes to the netCDF variables.
               call check( nf90_put_att(ncid,Varid4,UNITS, &
                    trim(OutUnits(outvar(i)))))
               call check(nf90_put_att(ncid,Varid4,missing_value, &
                    MissingValues)) 
              call check(nf90_put_att(ncid,Varid4,'_FillValue', -9999.00))
               ! End define mode.
               call check(nf90_enddef(ncid))
               CALL check(nf90_sync(ncid))
               NCIDARRAY(j,i)=ncid
            end do
        end do
        end Subroutine DirectoryCreate
        
    Subroutine OutputnetCDF(NCIDARRAY,outvar,NumtimeStep,outcount,incfile,ioutv,jxcoord,iycoord,NumtimeStepPerFile,NumofFile,StartEndNCDF,OutVarValue)
    use NETCDF
    
    integer, parameter :: NDIMS = 3
    integer::incfile,ioutv,outcount,NumtimeStep
    integer :: start(NDIMS), count(NDIMS)
    integer:: iycoord,jxcoord, timerec
    integer::NumtimeStepPerFile(NumofFile),StartEndNCDF(NumofFile,2)
    integer::NCIDARRAY(NumofFile,outcount),OutVar(outcount)
    REAL:: OutVarValue(NumtimeStep,64)
    
    timerec=NumtimeStepPerFile(incfile)
!    count = (/ 1, 1, timerec /)
!    start = (/ 1, 1, 1 /)
    count = (/ timerec, 1, 1 /)
    start = (/ 1, 1, 1 /)
    start(2) = jxcoord
    start(3) = iycoord 
    call check(nf90_put_var(NCIDARRAY(incfile,ioutv),4,OutVarValue(StartEndNCDF(incfile,1):StartEndNCDF(incfile,2),outvar(ioutv)),start, count))
    End Subroutine 
   

    Subroutine OutputtimenetCDF(NCIDARRAY,NumtimeStep,outcount,incfile,ioutv,NumtimeStepPerFile,NumofFile,StartEndNCDF,FNDJDT)
        use NETCDF
        implicit none
        integer, parameter :: NDIMS = 3
        integer::outcount,NumtimeStep,NumofFile
        integer::incfile,ioutv,OutVar(outcount)
        character (50) :: FILE_NAME
        integer :: start(NDIMS), count(NDIMS),VarId,recid
        integer:: iycoord,jxcoord, timerec,NCIDARRAY(NumofFile,outcount)
        integer::NumtimeStepPerFile(NumofFile),StartEndNCDF(NumofFile,2)
        Double precision::FNDJDT(NumtimeStep)
        timerec=NumtimeStepPerFile(incfile)
        call Check(NF90_PUT_VAR(NCIDARRAY(incfile,ioutv),1,FNDJDT(StartEndNCDF(incfile,1):StartEndNCDF(incfile,2))))
    End Subroutine 