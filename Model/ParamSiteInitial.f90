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
        Character (100) ParamHeading
        CHARACTER*15 filecode, ParamName
        parameter(n = 34)
        CHARACTER(10), DIMENSION(n) :: ParamSymbol
        Real ParamValue(n), param(*), bca, bcc
        integer:: irad, ireadalb
        CHARACTER*200 vfile
        ! define loop parameter
        integer reason, paramcode
        integer matched
        ParamSymbol = (/ "irad    ","ireadalb","tr      ","ts      ", &
              "ems     ","cg      ","z       ","zo      ","rho     ", &
              "rhog    ","lc      ","ks      ","de      ","avo     ", &
              "anir0   ","lans    ","lang    ","wlf     ","rd1     ", &
              "fstab   ","tref    ","dnews   ","emc     ","alpha   ", &
              "alphal  ","g       ","uc      ","as      ","Bs      ", &
              "lambda  ","rimax   ","wcoeff  ","a       ","c       " /)
        OPEN(11,FILE=vfile,STATUS='OLD', ACTION='READ')
        READ(11,*) ParamHeading
        ! Read until end of file
100    	READ(11,*,iostat=reason, end=200)filecode
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
        irad=paramValue(1)
        ireadalb=paramValue(2)
        Param(1:11)=ParamValue(3:13)
        Param(13:21)=ParamValue(14:22)
        Param(23:32)=ParamValue(23:32)
        bca = ParamValue(33)
        bcc = ParamValue(34)           
        return
        end subroutine readvals
!===End of reading the values from parameter file ==

! =================== for one point within the grid that the model is looping over 
      subroutine readsv(param,statev,sitev,svfile,slope,azi,lat,subtype, &
      dimlen2,dimlen1,ilat,jlon,dtbar,ts_last,longitude)
    ! param (input and output) is 
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
      PARAMETER(n=32)
      integer:: SVDT,reason, krecp,dimlen2,dimlen1, &
           matched,isVarFromNC(n),subtype            
      integer:: ilat,jlon
      real param(32)
      REAL statev(6),sitev(10),slope,azi,lat, StateSiteValue(n), &
           SingleArray(1), vardefaults(n), dtbar(12),gsurf, subalb,ts_last,longitude
      !REAL, DIMENSION(16):: vardefaults
      CHARACTER*50 SiteHeading, SVcode, SVname, StateSiteVName(n)
      CHARACTER*50 SiteVarNameinNCDF(n)
      ! File names are allowed to be up to 512 characters
      character*50 svfile, StateSiteFiles(n), file_name, varname

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
        isVarFromNC(i)= -1  ! use the coding 0 to indicate SCTC from file, 
                            ! 1 to indicate from raster, -1 not yet read
      enddo
      OPEN(8,FILE=svfile,STATUS='OLD')
      Read (8,*)  SiteHeading
      ! Here we start a loop that reads 3 lines at a time and determines the 
300    	  Read(8,*,iostat=reason, end=400)SVcode
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
                    read(8,*)StateSiteFiles(i)
                    read(8,*)SiteVarNameinNCDF(i)
                end if
                exit
            end if 
        end do
        if(matched .ne. 1)write(6,*)'site or initial condition variable code not matched:  ', &
             SVname
      end if
      go to 300
400     CLOSE(8)
    ! At this point all information has been read from the input file
    ! Now loop over each variable and if it is a netcdf file read the value out of the right place
    ! also check to make sure we got a value for every variable
      Do i=1,n,1                                                        !ivariables is a looping variable
        if (isVarFromNC(i) .eq. 1) then ! here 
            file_name=StateSiteFiles(i)
            varname=trim(StateSiteVName(i))
            CALL lowercase(StateSiteVName(i),StateSiteVName(i))
            CALL nCDF2DRead(file_name,SiteVarNameinNCDF(i),SingleArray, &
                 jlon,ilat)
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
        subtype=StateSiteValue(17)
        param(22)=StateSiteValue(18)
        !gsurf = StateSiteValue(17)
        dtbar(1:12) = StateSiteValue(19:30)
        ts_last= StateSiteValue(31)
        longitude=StateSiteValue(32)
      return
      end subroutine readsv
! ====== End of reading the values from site variable file===
