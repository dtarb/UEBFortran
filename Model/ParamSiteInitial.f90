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

!==This subroutine reads the values from parameter file ==
        subroutine readvals(param,irad,ireadalb,bca,bcc,vfile)
!param (output) an array that holds all the paramter values 
!irad (out) a falg for radiation calculation
!ireadalb Albedo reading control flag (0=albedo is computed internally, 1 albedo is read)
!bca (out) A in Bristow-Campbell formula for atmospheric transmittance
!bcc (out) C in Bristow-Campbell formula for atmospheric transmittance
!vfile (in0 file that stores paramter values
        Implicit None
        Character (100) ParamHeading
        CHARACTER*15 filecode, ParamName
        integer::n
        parameter(n = 32)
        CHARACTER(10), DIMENSION(n) :: ParamSymbol
        Real ParamValue(n), param(32), bca, bcc
        integer:: irad, ireadalb
        CHARACTER*200 vfile
        ! define loop parameter
        integer reason
        integer matched,i
        ParamSymbol = (/ "irad    ","ireadalb","tr      ","ts      ", &
              "ems     ","cg      ","z       ","zo      ","rho     ", &
              "rhog    ","lc      ","ks      ","de      ","avo     ", &
              "anir0   ","lans    ","lang    ","wlf     ","rd1     ", &
              "dnews   ","emc     ","alpha   ", &
              "alphal  ","g       ","uc      ","as      ","Bs      ", &
              "lambda  ","rimax   ","wcoeff  ","a       ","c       " /)
        OPEN(11,FILE=vfile,STATUS='OLD', ACTION='READ')
        READ(11,*) ParamHeading
        ! Read until end of file
100     READ(11,*,iostat=reason, end=200)filecode
        CALL lowercase(filecode,filecode)
        if (reason .eq. 0)then
            do i=1,n,1
                CALL lowercase(ParamSymbol(i),ParamSymbol(i))
                paramname=ADJUSTL(filecode(1:(SCAN (filecode, ':')-1)))
                CALL lowercase(paramname,paramname)
                if(paramname .eq. trim(ParamSymbol(i))) then
                    matched=1
                    !paramcode=1 !used to check if all the paramters are provided with correction notation
                    read(11,*)ParamValue(i)
                    exit
                end if
            end do
        if(matched .ne. 1) write(6,*) 'Parameter: ', paramname
        end if
        matched=0
        go to 100
200     CLOSE(11)
!  Mapping from parameters read to UEB internal interpretation which follows UEBVeg scheme
        irad=int(paramValue(1))
        ireadalb=int(paramValue(2))
        Param(1:11)=ParamValue(3:13)
        Param(13:18)=ParamValue(14:19)
        Param(19)=-9999
        Param(20)=-9999
        Param(21)=ParamValue(20)
        Param(23:32)=ParamValue(21:30)
        bca = ParamValue(31)
        bcc = ParamValue(32)           
        return
        end subroutine readvals
!===End of reading the values from parameter file ==

! =================== for one point within the grid that the model is looping over 
      subroutine readsv(param,statev,sitev,svfile,slope,azi,lat,subtype, &
      ilat,jlon,dtbar,ts_last,longitude,Sitexcoordinates,Siteycoordinates)
    ! param (input and output) an array that holds all the paramter values
    ! statev (output) is state variable array that returns the initial conditions read from input files
    ! sitev (output)is site variables array that returns the site variables read from input files
    ! svfile (input)is name of control file for site variables and initial conditions
    ! slope (output)- site variable slope which for historical consistency is not part of sitev array
    ! azi (output)- azimuth site variable also not part of sitev array
    ! lat (output)- latitude site variable also not part of sitev array
    ! dimlen2 (input) - longitude dimension (number of columns) in netcdf grid covering the model domain
    ! dimlen1 (input) - latitude dimension (number of rows) in netcdf grid covering the model domain
    ! ilat (input) - latitude (row) index of the site where data is to be read and returned (from outside loop over space)
    ! jlat (input) - longitude (col) index of the site where data is to be read and returned (from outside loop over space)
      use netCDF 
      Implicit None
      integer:: n,i
      PARAMETER(n=32)
      integer:: reason,matched,isVarFromNC(n),subtype            
      integer:: ilat,jlon
      real param(32)
      REAL statev(6),sitev(10),slope,azi,lat, StateSiteValue(n), &
           SingleArray(1), vardefaults(n), dtbar(12),ts_last,longitude
      !REAL, DIMENSION(16):: vardefaults
      CHARACTER*50 SiteHeading, SVcode, SVname, StateSiteVName(n)
      CHARACTER*50 SiteVarNameinNCDF(n)
      Character*50 Sitexcoordinate, Siteycoordinate
      Character*50 Sitexcoordinates(n), Siteycoordinates(n)
      ! File names are allowed to be up to 512 characters
      character*50 svfile, StateSiteFiles(n), file_name, varname
      Character*200,dimension(:),allocatable :: args  !  Variable to hold the index position of each output variable
      integer:: nargs1
      CHARACTER(200) :: str 
      CHARACTER(1) :: delimit1,delimit2,delimit3
      integer:: nargs,nargs2,nargs3,nargs4,nargs5
      character(200),Allocatable:: words(:)
      Character(200):: StateSiteFilesR, SitexcoordinateR, SiteycoordinateR, InputtcoordinateR
      Character(200):: VarNameinNCDFR
      Character(200),Allocatable:: words1(:),words2(:),words3(:),words4(:),words5(:)
      integer:: DefaultDimValues(3),SiteDefDimval(n,3),ncidout
      REAL:: RangeMin,RangeMax
      character(200),allocatable::Words7element(:)
      delimit1=';'
      delimit2=':'
      delimit3=','
      StateSiteVName= (/ "USic     ","WSis     ","Tic      ","WCic     ", &
             "df       ","apr      ","Aep      ","cc       ","hcan     ", &
             "lai      ","Sbar     ","ycage    ","slope    ","aspect   ", &
             "latitude ","subalb   ","subtype  ","gsurf    ","b01      ", &
             "b02      ","b03      ","b04      ","b05      ","b06      ", &
             "b07      ","b08      ","b09      ","b10      ","b11      ", &
             "b12      ","ts_last  ","longitude"/)

      vardefaults= (/ 0.0, 0.0, 0.0, 0.0, 1.0, 100000.0, 0.1, 0.0, 0.0, &
           0.0, 6.6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.25, 0.98, 5.712903, &
           4.350000, 6.890322, 8.660001, 8.938710, 10.010000, 9.541936, &
           9.038710, 7.160001, 8.106450, 5.923332, 5.058064, -9999.0, 111.00 /)
      do i=1,n  
      SiteDefDimval=-9999
      isVarFromNC(i)= -1  ! use the coding 0 to indicate SCTC from file, 
                            ! 1 to indicate from raster, -1 not yet read
      enddo
      OPEN(8,FILE=svfile,STATUS='OLD')
      Read (8,*)  SiteHeading
      ! Here we start a loop that reads 3 lines at a time and determines the 
300       Read(8,*,iostat=reason, end=400)SVcode
          if(reason .eq. 0) then  ! anything other than 0 means some sort of 
                                  ! quirk in the read - we are not currently checking for this
            SVname=SVcode(1:(SCAN (SVcode, ':')-1))
            CALL lowercase(SVcode,SVcode)
            matched=0
        do i=1,n,1
            CALL lowercase(StateSiteVName(i),StateSiteVName(i))
            CALL lowercase(SVname,SVname)
            if(SVname .eq. trim(StateSiteVName(i))) then
                matched=1
                read(8,*)isVarFromNC(i)                                 !SVDT=site variable data type ditector
                if(isVarFromNC(i).eq. 0) then
                    read(8,*)StateSiteValue(i)
                else
                    read(8,fmt='(A)')str
                    CALL StringSep(str,delimit1,nargs)
                    Allocate(words(nargs))
                    Allocate(Words7element(7))
                    CALL StringSepWord(str,delimit1,nargs,words)
                    Words7element(1:nargs)=words(1:nargs)
                    CALL StringToVarName(nargs,words7element,delimit2,StateSiteFilesR,SitexcoordinateR,SiteycoordinateR,&
                                         &VarNameinNCDFR,InputtcoordinateR,DefaultDimValues,RangeMin,RangeMax,delimit3,.true.)
                    SiteDefDimval(i,1)=DefaultDimValues(1) !2-D netCDF file mapping (there is no time in a 2-D file)
                    SiteDefDimval(i,2)=DefaultDimValues(2) !2-D netCDF file mapping (get y coordinate)
                    SiteDefDimval(i,3)=DefaultDimValues(3) !2-D netCDF file mapping (get x coordinate)
                    StateSiteFiles(i)=StateSiteFilesR
                    Sitexcoordinates(i)=SitexcoordinateR
                    Siteycoordinates(i)=SiteycoordinateR                                        
                    SiteVarNameinNCDF(i)=VarNameinNCDFR
                    Deallocate(words)
                    deAllocate(Words7element)
                end if
                exit
            end if 
        end do
        if(matched .ne. 1)write(6,*)'site or initial condition variable code not matched:  ', &
             SVname
      end if
      go to 300
400   CLOSE(8)

       Do i=1,n
        if(isVarFromNC(i).eq. 1) then
            If(SiteDefDimval(i,2) .NE. -9999)THEN ! in a 2-D netCDF file first dimensionn is y and x is second
                                                  ! therefore, to get y-dim we need to ask for second element of 
                                                  ! DefaultDimValues
                CALL check(nf90_open(StateSiteFiles(i),NF90_NOWRITE, ncidout))
                CALL check(nf90_inquire_dimension(ncidout,SiteDefDimval(i,2),Siteycoordinates(i)))
                CALL check(nf90_close(ncidout))
            End if
            If(SiteDefDimval(i,3) .NE. -9999)THEN ! in a 2-D netCDF file first dimensionn is y and x is second
                                                  ! therefore, to get x-dim we need to ask for second element of 
                                                  ! DefaultDimValues
                CALL check(nf90_open(StateSiteFiles(i),NF90_NOWRITE, ncidout))
                CALL check(nf90_inquire_dimension(NCIDout,SiteDefDimval(i,3),Sitexcoordinates(i)))
                CALL check(nf90_close(ncidout))
            End if 
        end if   
      End do
      
    ! At this point all information has been read from the input file
    ! Now loop over each variable and if it is a netcdf file read the value out of the right place
    ! also check to make sure we got a value for every variable
      Do i=1,n,1                                                        !ivariables is a looping variable
        if (isVarFromNC(i) .eq. 1) then ! here 
            file_name=StateSiteFiles(i)
            varname=trim(StateSiteVName(i))
            CALL lowercase(StateSiteVName(i),StateSiteVName(i))
            CALL nCDF2DReadReal(file_name,SiteVarNameinNCDF(i),SingleArray,ilat,jlon,Sitexcoordinates(i),Siteycoordinates(i))
            StateSiteValue(i)=SingleArray(1)
        end if
        if (isVarFromNC(i) .eq. -1)then
          write(6,*)'No input for site or initial condition variable: ', &
               trim(StateSiteVName(i))
          write(6,*)'Assuming default value ',trim(StateSiteVName(i)),' = ', &
               vardefaults(i)
          StateSiteValue(i)=vardefaults(i)
        endif
      end do
        statev(1:4)=StateSiteValue(1:4)
        sitev(1:2)=StateSiteValue(5:6)
        sitev(4:9)=StateSiteValue(7:12)
        slope=StateSiteValue(13)
        azi=StateSiteValue(14)
        lat=StateSiteValue(15)
        param(12)=StateSiteValue(16)
        !subalb=StateSiteValue(15)
        sitev(10)=StateSiteValue(17)
        subtype=int(StateSiteValue(17))
        param(22)=StateSiteValue(18)
        !gsurf = StateSiteValue(17)
        dtbar(1:12) = StateSiteValue(19:30)
        ts_last= StateSiteValue(31)
        longitude=StateSiteValue(32)
      return
      end subroutine readsv
! ====== End of reading the values from site variable file===

! =================== for one point within the grid that the model is looping over 
      subroutine readsvcoordinate(svfile,Sitexcoordinates,Siteycoordinates)
    ! param (input and output) an array that holds all the paramter values
    ! statev (output) is state variable array that returns the initial conditions read from input files
    ! sitev (output)is site variables array that returns the site variables read from input files
    ! svfile (input)is name of control file for site variables and initial conditions
    ! slope (output)- site variable slope which for historical consistency is not part of sitev array
    ! azi (output)- azimuth site variable also not part of sitev array
    ! lat (output)- latitude site variable also not part of sitev array
    ! dimlen2 (input) - longitude dimension (number of columns) in netcdf grid covering the model domain
    ! dimlen1 (input) - latitude dimension (number of rows) in netcdf grid covering the model domain
    ! ilat (input) - latitude (row) index of the site where data is to be read and returned (from outside loop over space)
    ! jlat (input) - longitude (col) index of the site where data is to be read and returned (from outside loop over space)
      use netCDF
      Implicit None
      integer:: n,i
      PARAMETER(n=32)
      integer:: reason,matched,isVarFromNC(n)         
      !REAL, DIMENSION(16):: vardefaults
      CHARACTER*50 SiteHeading, SVcode, SVname, StateSiteVName(n)
      CHARACTER*50 SiteVarNameinNCDF(n)
      REAL:: StateSiteValue(n)
      Character*50 Sitexcoordinate, Siteycoordinate
      Character*50 Sitexcoordinates(n), Siteycoordinates(n)
      ! File names are allowed to be up to 512 characters
      character*50 svfile, StateSiteFiles(n), file_name, varname
      CHARACTER(200) :: str 
      CHARACTER(1) :: delimit1,delimit2
      integer:: nargs,nargs2,nargs3,nargs4,nargs5
      character(200),Allocatable:: words(:)
      Character(200):: StateSiteFilesR, SitexcoordinateR, SiteycoordinateR, InputtcoordinateR
      Character(200):: VarNameinNCDFR
      Character(200),Allocatable:: words1(:),words2(:),words3(:),words4(:),words5(:)
      integer:: DefaultDimValues(3)
      integer:: SiteDefDimval(n,3),NCIDout
      character(1)::delimit3
      REAL:: RangeMin,RangeMax
      character(200),Allocatable::Words7element(:)
      delimit1=';'
      delimit2=':'
      delimit3=','
      SiteDefDimval=-9999
      StateSiteVName= (/ "USic     ","WSis     ","Tic      ","WCic     ", &
             "df       ","apr      ","Aep      ","cc       ","hcan     ", &
             "lai      ","Sbar     ","ycage    ","slope    ","aspect   ", &
             "latitude ","subalb   ","subtype  ","gsurf    ","b01      ", &
             "b02      ","b03      ","b04      ","b05      ","b06      ", &
             "b07      ","b08      ","b09      ","b10      ","b11      ", &
             "b12      ","ts_last  ","longitude"/)
             
      do i=1,n  
        isVarFromNC(i)= -1  ! use the coding 0 to indicate SCTC from file, 
                            ! 1 to indicate from raster, -1 not yet read
      enddo
      OPEN(888,FILE=svfile,STATUS='OLD')
      Read(888,*)  SiteHeading
      ! Here we start a loop that reads 3 lines at a time and determines the 
300       Read(888,*,iostat=reason, end=400)SVcode
          if(reason .eq. 0) then  ! anything other than 0 means some sort of 
                                  ! quirk in the read - we are not currently checking for this
            SVname=SVcode(1:(SCAN (SVcode, ':')-1))
            CALL lowercase(SVcode,SVcode)
            matched=0
        do i=1,n,1
            CALL lowercase(StateSiteVName(i),StateSiteVName(i))
            CALL lowercase(SVname,SVname)
            if(SVname .eq. trim(StateSiteVName(i))) then
                matched=1
                read(888,*)isVarFromNC(i)                                 !SVDT=site variable data type ditector
                if(isVarFromNC(i).eq. 0) then
                    read(888,*)StateSiteValue(i)
                else
                    read(888,fmt='(A)')str
                    CALL StringSep(str,delimit1,nargs)
                    Allocate(words(nargs))
                    Allocate(Words7element(7))
                    CALL StringSepWord(str,delimit1,nargs,words)
                    words7element(1:nargs)=words(1:nargs)
                    CALL StringToVarName(nargs,words7element,delimit2,StateSiteFilesR,SitexcoordinateR,SiteycoordinateR,&
                                         &VarNameinNCDFR,InputtcoordinateR,DefaultDimValues,RangeMin,RangeMax,delimit3,.true.)
                    SiteDefDimval(i,1)=DefaultDimValues(1) !2-D netCDF file mapping (get y coordinate)
                    SiteDefDimval(i,2)=DefaultDimValues(2) !2-D netCDF file mapping (get x coordinate)
                    SiteDefDimval(i,3)=DefaultDimValues(3) !2-D netCDF file mapping (there is no time in a 2-D file)
                    StateSiteFiles(i)=StateSiteFilesR
                    Sitexcoordinates(i)=SitexcoordinateR
                    Siteycoordinates(i)=SiteycoordinateR                                        
                    SiteVarNameinNCDF(i)=VarNameinNCDFR
                    Deallocate(words)
                    deAllocate(Words7element)
                end if
                exit
            end if 
        end do
        if(matched .ne. 1)write(6,*)'site or initial condition variable code not matched:  ', &
             SVname
      end if
      go to 300
400   CLOSE(8)

       Do i=1,n
        if(isVarFromNC(i).eq. 1) then
            If(SiteDefDimval(i,2) .NE. -9999)THEN ! in a 2-D netCDF file first dimensionn is y and x is second
                                                  ! therefore, to get y-dim we need to ask for second element of 
                                                  ! DefaultDimValues
                CALL check(nf90_open(StateSiteFiles(i),NF90_NOWRITE, ncidout))
                CALL check(nf90_inquire_dimension(ncidout,SiteDefDimval(i,2),Siteycoordinates(i)))
                CALL check(nf90_close(ncidout))
            End if
            If(SiteDefDimval(i,3) .NE. -9999)THEN ! in a 2-D netCDF file first dimensionn is y and x is second
                                                  ! therefore, to get x-dim we need to ask for second element of 
                                                  ! DefaultDimValues
                CALL check(nf90_open(StateSiteFiles(i),NF90_NOWRITE, ncidout))
                CALL check(nf90_inquire_dimension(NCIDout,SiteDefDimval(i,3),Sitexcoordinates(i)))
                CALL check(nf90_close(ncidout))
            End if 
        end if   
      End do
      
     return
     end subroutine readsvcoordinate
! ====== End of reading the values from site variable file===