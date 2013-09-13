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
        subroutine NumOutFiles(OutControlFILE, ModelStartDate,&
        ModelStartHour,ModelEndDate,ModelEndHour,Modeldt,&
        NumtimeStep,NumofFile,NumOutPoint,OutCount,dimlen1,dimlen2)
        
        ! OutControlFILE (in) input control file
        ! ModelStartDate(3) (input) Array giving start year, month, day
        ! ModelEndDate(3) (input)  Array giving end year, month, day
        ! ModelStartHour (input) Start hour
        ! ModelEndHour (input)  end hour
        ! Modeldt (input) model time resolution
        ! Dimlen2 (out) length of y-coordinate
        ! dimlen1 (out) length of x-coordinate
        ! NumtimeStep (out) the total number of time steps output
        ! NumofFile (out) the number of output files per variable output to netCDF
        ! NumOutPoint (out) Number of points for point detail output
        ! OutCount (out) Number of variables for utputting in NC files
        
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
        integer:: netCDFDataValPerfile, LineNumberCount
        
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
        NumtimeStep=int(((JMED+tol-JMSD)/(Modeldt/24.0)) + 1)
        ntsperfile=max(netCDFDataValPerfile/(dimlen1*dimlen2),1)
        NumofFile=NumtimeStep/ntsperfile+1  ! The +1 is to round up integer calculation
        NumOutPoint=0
        OutCount=0
        LineNumberCount=0
        
        OPEN(109,FILE=OutControlFile,STATUS='OLD', ACTION='READ')
3316    Read(109,*,iostat=reason, end=3315)filecode
        LineNumberCount=LineNumberCount+1
        go to 3316
3315    close(109)  
        
        OPEN(109,FILE=OutControlFile,STATUS='OLD', ACTION='READ')
        READ(109,*) xxs
        ! Read until end of file
        If (LineNumberCount .GT. 1)THEN
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
        END IF
 3300     close(109)       
        return
        end subroutine NumOutFiles
!===End of reading the values from outconrol.dat file ==

!==This subroutine reads the outputcontrol file and returns
!  The folder path to each variable
!  The file name to each variable that is outpus
!  An array indicating whether each variable is output or not 

        subroutine OutputFiles(OutControlFILE,NumtimeStep,NumofFile,outSampleFile,NumtimeStepPerFile,dimlen1,dimlen2,&
        OutVar,OutPoint,OutPointFiles,NumOutPoint,OutCount)
        
        ! OutControlFILE
        ! NumtimeStep
        ! Dimlen2
        ! dimlen1
        ! NumofFile
        ! outSampleFile
        ! NumtimeStepPerFile
        ! OutVar
        ! OutPoint
        ! OutPointFiles
        ! NumOutPoint
        ! OutCount

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
        integer:: OutNum,LINENUMBERCOUNT
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
        LineNumberCount=0
        OPEN(109,FILE=OutControlFile,STATUS='OLD', ACTION='READ')
3316    Read(109,*,iostat=reason, end=3315)filecode
        LineNumberCount=LineNumberCount+1
        go to 3316
3315    close(109) 

        OPEN(109,FILE=OutControlFile,STATUS='OLD', ACTION='READ')
        READ(109,*) outHeading
        OutNum=0
        If (LineNumberCount .GT. 1)THEN
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
        End if
4600    CLOSE(109)
        return
        end subroutine OutputFiles
        !===End of reading the values from outconrol.dat file ==

        Subroutine DirectoryCreate(nrefyr,nrefmo,nrefday,WATERSHEDfILE,wsycoordinate,&
        &wsxcoordinate,NumofFile,outcount,Outvar,NumtimeStepPerFile,&
        &outSampleFile,OutputNCContainer,NCIDARRAY)
        
        ! nrefyr
        ! nrefmo
        ! nrefday
        ! dimlen1
        ! dimlen2
        ! DimName1
        ! DimName2
        ! DimUnit1
        ! DimUnit2
        ! NumofFile
        ! outcount
        ! Outvar
        ! NumtimeStepPerFile
        ! outSampleFile
        ! OutputNCContainer
        ! NCIDARRAY
        
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
        character*200:: WATERSHEDFILE
        character*50:: wsycoordinate,wsxcoordinate
        integer:: ncidout, dimid1, dimid2, vartype, varid
        Integer:: dimlen1,dimlen2,dimlen1OUT,dimlen2OUT
        byte, allocatable:: inbyte1(:)
        character, allocatable:: inchar1(:)
        integer*2, allocatable:: inshort1(:)
        integer*4, allocatable:: ininteger1(:)
        real*4, allocatable:: inreal1(:)
        REAL*8, allocatable:: indouble1(:)

        byte, allocatable:: inbyte2(:)
        character, allocatable:: inchar2(:)
        integer*2, allocatable:: inshort2(:)
        integer*4, allocatable:: ininteger2(:)
        real*4, allocatable:: inreal2(:)
        real*8, allocatable:: indouble2(:)

        Real*8, allocatable:: dimvalue1(:),dimvalue2(:)
        Real*8, allocatable:: DimValue1Org(:),DimValue2Org(:)
        
        integer, parameter:: NF90_BYTEs1 = 1, NF90_CHARs1 = 2, NF90_SHORTs1 = 3, NF90_INTs1 = 4, NF90_FLOATs1 = 5, NF90_DOUBLEs1 = 6
        integer, parameter:: NF90_BYTEs2 = 1, NF90_CHARs2 = 2, NF90_SHORTs2 = 3, NF90_INTs2 = 4, NF90_FLOATs2 = 5, NF90_DOUBLEs2 = 6
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
        Real:: MissingValues,Fillvalues

        double precision:: DimValueX
        
        
        call check(nf90_open(WATERSHEDfILE,nf90_nowrite, ncidout))                 ! open the netcdf file
        ! get dimension IDs
        call check(nf90_inq_dimid(ncidout,wsycoordinate,DimID1))
        call check(nf90_inq_dimid(ncidout,wsxcoordinate,DimID2))
        call check(nf90_inquire_dimension(ncidout,DimID1,len=dimlen1))       ! Information about dimensionID 1
        call check(nf90_inquire_dimension(ncidout,DimID2,len=dimlen2))       ! information about dimensionID 2  
        
        allocate(dimvalue1(dimlen1))
        allocate(dimvalue2(dimlen2))
        
        Allocate(inbyte1(dimlen1))
        Allocate(inchar1(dimlen1))
        Allocate(inshort1(dimlen1))
        Allocate(ininteger1(dimlen1))
        Allocate(inreal1(dimlen1))
        Allocate(indouble1(dimlen1))
        Allocate(DimValue1Org(dimlen1))

        Allocate(inbyte2(dimlen2))
        Allocate(inchar2(dimlen2))
        Allocate(inshort2(dimlen2))
        Allocate(ininteger2(dimlen2))
        Allocate(inreal2(dimlen2))
        Allocate(indouble2(dimlen2))
        Allocate(DimValue2Org(dimlen2))
                
        call check(nf90_inq_varid(ncidout,wsycoordinate, VarId))                    ! information about variableID for a given VariableName
        call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
        call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit1)) 

        if(vartype .eq. NF90_BYTEs1)then
            call check(nf90_get_var(ncidout,Varid,inbyte1)) 
            DimValue1Org=REAL(inbyte1,8)
        elseif(vartype .eq. NF90_CHARs1)then
            call check(nf90_get_var(ncidout,Varid,inchar1)) 
            Write (6, *) "Error: A site variable can't be a character in ", WATERSHEDfILE
        elseif(vartype .eq. NF90_SHORTs1)then
            call check(nf90_get_var(ncidout,Varid,inshort1)) 
            DimValue1Org=REAL(inshort1,8)
        elseif(vartype .eq. NF90_INTs1)then
            call check(nf90_get_var(ncidout,Varid,ininteger1)) 
            DimValue1Org=REAL(ininteger1,8)
        elseif(vartype .eq. NF90_FLOATs1)then 
            call check(nf90_get_var(ncidout,Varid,inshort1)) 
            DimValue1Org=REAL(inreal1,8) 
        elseif(vartype .eq. NF90_DOUBLEs1)then
            call check(nf90_get_var(ncidout,Varid,indouble1))  
            DimValue1Org=REAL(indouble1,8) 
        END IF 
        DimValue1=DimValue1Org

        call check(nf90_inq_varid(ncidout,wsxcoordinate, VarId))                    ! information about variableID for a given VariableName 
        call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
        call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit2)) 

        if(vartype .eq. NF90_BYTEs2)then
            call check(nf90_get_var(ncidout,Varid,inbyte2))
            DimValue2Org=REAL(inbyte2,8)
        elseif(vartype .eq. NF90_CHARs2)then
            call check(nf90_get_var(ncidout,Varid,inchar2))
            Write (6, *) "Error: A site variable can't be a character in ", WATERSHEDfILE
        elseif(vartype .eq. NF90_SHORTs2)then
            call check(nf90_get_var(ncidout,Varid,inshort2))
            DimValue2Org=REAL(inshort2,8)
        elseif(vartype .eq. NF90_INTs2)then
            call check(nf90_get_var(ncidout,Varid,ininteger2))
            DimValue2Org=REAL(ininteger2,8)
        elseif(vartype .eq. NF90_FLOATs2)then
            call check(nf90_get_var(ncidout,Varid,inreal2)) 
            DimValue2Org=REAL(inreal2,8) 
        elseif(vartype .eq. NF90_DOUBLEs2)then
            call check(nf90_get_var(ncidout,Varid,indouble2)) 
            DimValue2Org=REAL(indouble2,8) 
        END IF 
        DimValue2=DimValue2Org
        
        Deallocate(inbyte1)
        Deallocate(inchar1)
        Deallocate(inshort1)
        Deallocate(ininteger1)
        Deallocate(inreal1)
        Deallocate(indouble1)
        Deallocate(DimValue1Org)
        
        Deallocate(inbyte2)
        Deallocate(inchar2)
        Deallocate(inshort2)
        Deallocate(ininteger2)
        Deallocate(inreal2)
        Deallocate(indouble2)
        Deallocate(DimValue2Org)
        
        write(refyear,'(i4)') nrefyr
        write(refmonth,'(i2.2)') nrefmo
        write(refday,'(i2.2)') nrefday
        MissingValues=-9999.000
        Fillvalues=-9999.0
        time_unit='days since'//' '//trim(refyear)//'-'//trim(refmonth) &
       //'-'//trim(refday)//' '//'00:00:00'
       
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
            If(NumofFile .GE. 2)then
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
            End If
            
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
               call check(nf90_def_dim(ncid,wsycoordinate,dimlen1,y_dimid))
               call check(nf90_def_dim(ncid,wsxcoordinate,dimlen2,x_dimid))

               ! Define the coordinate variables. We will only define coordinate
               ! variables for lat and lon.  Ordinarily we would need to provide
               ! an array of dimension IDs for each variable's dimensions, but
               ! since coordinate variables only have one dimension, we can
               ! simply provide the address of that dimension ID (lat_dimid) and
               ! similarly for (lon_dimid).
               call check(nf90_def_var(ncid,"time",NF90_REAL,time_dimid, &
                    time_varid))
               call check(nf90_def_var(ncid,wsycoordinate,NF90_DOUBLE,y_dimid, &
                    y_varid))
                call check(nf90_def_var(ncid,wsxcoordinate,NF90_DOUBLE,x_dimid, &
                    x_varid))

               ! Assign units attributes to coordinate variables.
               call check(nf90_put_att(ncid,y_varid,UNITS,DimUnit1))
               call check(nf90_put_att(ncid,x_varid,UNITS,DimUnit2))
               call check(nf90_put_att(ncid,time_varid,UNITS,time_unit))
               ! The dimids array is used to pass the dimids of the dimensions of
               ! the netCDF variables. Both of the netCDF variables we are creating
               ! share the same four dimensions. In Fortran, the unlimited
               ! dimension must come last on the list of dimids.
               dimids = (/time_dimid,y_dimid,x_dimid/)               
               ! Define the netCDF variables for the pressure and temperature data.
               call check( nf90_def_var(ncid,trim(OutSymbol(outvar(i))), &
                    NF90_REAL, &
                    dimids,Varid4))
               ! Assign units attributes to the netCDF variables.
               call check( nf90_put_att(ncid,Varid4,UNITS, &
                    trim(OutUnits(outvar(i)))))
               call check(nf90_put_att(ncid,Varid4,missing_value, MissingValues)) 
               !call check(nf90_put_att(ncid,Varid4,'_FillValue', NF90_FILL_DOUBLE))
               call check(nf90_put_att(ncid,Varid4,'_FillValue', FillValues))
               ! End define mode.
               call check(nf90_enddef(ncid))
               CALL check(nf90_sync(ncid))
               CALL Check(NF90_PUT_VAR(ncid,y_varid,DimValue1)) ! latitude/y is dimension 2
               CALL Check(NF90_PUT_VAR(ncid,x_varid,DimValue2)) ! longtude/x is dimension 3
              !CALL NF_PUT_VAR_DOUBLE(ncid,x_varid,DimValue1) ! longtude/x is dimension 3
              !CALL NF_PUT_VAR_DOUBLE(ncid,y_varid,DimValue2) ! longtude/y is dimension 2
              NCIDARRAY(j,i)=ncid
            end do
        end do
        end Subroutine DirectoryCreate
        
    Subroutine OutputnetCDF(NCIDARRAY,outvar,NumtimeStep,outcount,incfile,ioutv,jxcoord,iycoord,&
    NumtimeStepPerFile,NumofFile,StartEndNCDF,OutVarValue)
    use NETCDF
    
    ! NCIDARRAY
    ! outvar
    ! NumtimeStep
    ! outcount
    ! incfile
    ! ioutv
    ! jxcoord
    ! iycoord
    ! NumtimeStepPerFile
    ! NumofFile
    ! StartEndNCDF
    ! OutVarValue
    
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
    start(2) = iycoord
    start(3) = jxcoord 
    call check(nf90_put_var(NCIDARRAY(incfile,ioutv),4,OutVarValue(StartEndNCDF(incfile,1):StartEndNCDF(incfile,2),outvar(ioutv)),start, count))
    End Subroutine 
   

    Subroutine OutputtimenetCDF(NCIDARRAY,NumtimeStep,outcount,incfile,ioutv,NumtimeStepPerFile,NumofFile,StartEndNCDF,FNDJDT)
    
    ! NCIDARRAY
    ! NumtimeStep
    ! outcount
    ! incfile 
    ! ioutv
    ! NumtimeStepPerFile
    ! NumofFile
    ! StartEndNCDF
    ! FNDJDT
    
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