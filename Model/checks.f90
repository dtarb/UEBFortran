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

subroutine checks(svfile,IsInputFromNC,NumNCFiles,totalNC,StateSiteVName)
Implicit None
integer::n,m,i
parameter(n=32,m=11)
Character*200:: svfile
Character*200:: SVname,SiteHeading,SVcode,StateSiteVName(n),StateSiteNCFile,StateSiteV
Real:: StateSiteValue
integer:: reason,IsInputFromNC(m),countinput,count,NumNCFiles(m),isSiteVarFromNC  
integer:: matched,x,totalNC


        StateSiteVName=   (/ "usic     ","wsis     ","tic      ","wcic     ", &
                 "df       ","apr      ","aep      ","cc       ","hcan     ", &
                 "lai      ","sbar     ","ycage    ","slope    ","aspect   ", &
                 "latitude ","subalb   ","subtype  ","gsurf    ","b01      ", &
                 "b02      ","b03      ","b04      ","b05      ","b06      ", &
                 "b07      ","b08      ","b09      ","b10      ","b11      ", &
                 "b12      ","ts_last  ","longitude"/)
                  
          count=0
          countinput=0
          OPEN(81,FILE=svfile,STATUS='OLD')
          Read (81,*) SiteHeading
          ! Here we start a loop that reads 3 lines at a time and determines the 
3000      Read(81,*,iostat=reason, end=4000)SVcode
          if(reason .eq. 0) then  ! anything other than 0 means some sort of 
                                  ! quirk in the read - we are not currently checking for this
            SVname=SVcode(1:(SCAN (SVcode, ':')-1))
            CALL lowercase(SVname,SVname)
            matched=0
            do i=1,n
                if(SVname .eq. trim(StateSiteVName(i))) then
                    matched=1
                    read(81,*)isSiteVarFromNC                                 !SVDT=site variable data type ditector
                    if(isSiteVarFromNC .eq. 0) then
                        read(81,*)StateSiteValue
                    else
                        read(81,*)StateSiteNCFile 
                        read(81,*)StateSiteV
                        count=count+1
                    end if
                    exit
                end if 
            end do
            if(matched .ne. 1)write(6,*)'Site variable code not matched:  ', SVname             
          end if
          go to 3000
4000      CLOSE(81)

        do i=1,11
            if (IsInputFromNC(i) .eq. 1)Then
                x=NumNCFiles(i)
                countinput=countinput+x
            end if
        end do
        totalNC=count + countinput
        
  End Subroutine checks
  
  subroutine NCChecks(svfile,StateSiteVName,WatershedFile,MaxNumofFile,IsInputFromNC,NCDFContainer,NumNCFiles,totalNC,AllNCDFfile)
        Use netCDF
        Implicit None
        integer::n,m,i,j,ii
        parameter(n=32,m=11)
        Character*200:: StateSiteVName(n)
        integer:: MaxNumofFile,dimlen2,dimlen1,dimlen3,dimlen4
        Character*200:: svfile,NCDFContainer(MaxNumofFile,n)
        Character*200:: SVname,SiteHeading,SVcode,StateSiteV
        Character*200:: inputVName(m),WatershedFile
        Real:: StateSiteValue
        integer:: reason,IsInputFromNC(m),NumNCFiles(m),isVarFromNC(n)  
        integer:: matched,totalNC,NCNumber
        Character*200:: AllNCDFfile(totalNC),NCFile
        Real, dimension(:), allocatable:: DimValue1,DimValue2,DimValue3,DimValue4
        Character*200:: DimUnit1,DimUnit2,DimUnit3,DimUnit4,DimName1,DimName2,DimName3,DimName4
        
        InputVName= (/ "Ta     ","Prec   ","v      ","RH     ", &
             "Tmin   ","Tmax   ","Qsi    ","Qg     ","Qli    ", &
             "Qnet   ","Snowalb"  /)
                  
        NCNumber=1
        OPEN(8,FILE=svfile,STATUS='OLD')
          Read (8,*) SiteHeading
3010      Read(8,*,iostat=reason,end=4010)SVcode
          if(reason .eq. 0) then  ! anything other than 0 means some sort of 
                                  ! quirk in the read - we are not currently checking for this
            SVname=SVcode(1:(SCAN (SVcode, ':')-1))
            CALL lowercase(SVname,SVname)
            matched=0
            do i=1,n,1
                CALL lowercase(StateSiteVName(i),StateSiteVName(i))
                if(SVname .eq. trim(StateSiteVName(i))) then
                    matched=1
                    read(8,*)isVarFromNC(i)                                 !SVDT=site variable data type ditector
                    if(isVarFromNC(i).eq. 0) then
                        read(8,*)StateSiteValue
                    else
                        read(8,*)NCfile
                        read(8,*)StateSiteV
                        AllNCDFfile(NCNumber)= NCfile
                        NCNumber=NCNumber+1
                    end if
                    exit
                end if 
            end do
            if(matched .ne. 1)write(6,*)'Input variable code not matched:  ', SVname
          end if
          go to 3010
4010     CLOSE(8)
        
        do i=1,m
            if(IsInputFromNC(i) .eq. 1)then
                do j=1,NumNCFiles(i)
                    AllNCDFfile(NCNumber) = NCDFContainer(j,i)
                    NCNumber=NCNumber+1
                end do
             end if
        end do
        
        CALL nCDF2DArrayInfo(WatershedFile,dimlen2,dimlen1)
        II=1
501     If (II .le. (NCNumber-1))Then
        CALL nCDF2DArrayInfo(AllNCDFfile(II),dimlen4,dimlen3)
            If (dimlen2 .eq. dimlen4 .and. dimlen1 .eq. dimlen3)Then
                II=II+1
                GO To 501
            else
                Write(6,*) "Spatial coordinates of",AllNCDFfile(ii),"mismatch with watershed file's aspatial coordinates"
                Write(6,*) "Inconsistent number of grid cells in ", AllNCDFfile(ii), " netCDF file"
                GO To 501
            End if
        End if
        
        allocate(DimValue1(dimlen1))
        allocate(DimValue2(dimlen2))

        CALL SpatialCoordinate(WatershedFile,dimlen1,dimlen2,DimName1,DimName2,DimValue1,DimValue2,DimUnit1,DimUnit2)
!        CALL lowercase(inputname,inputname)
!        CALL lowercase(inputname,inputname)
        II=1
601     If (II .le. (NCNumber-1))Then
        allocate(DimValue3(dimlen3))
        allocate(DimValue4(dimlen4))
        CALL SpatialCoordinate(AllNCDFfile(II),dimlen1,dimlen2,DimName3,DimName4,DimValue3,DimValue4,DimUnit3,DimUnit4)
!            CALL lowercase(inputname,inputname)
!            CALL lowercase(inputname,inputname)
            If (((MAXVAL(DimValue1) - MAXVAL(DimValue3)) .LE. 1E-5) .and. ((MINVAL(DimValue2) - MINVAL(DimValue4)) .LE. 1E-5)&
            &.and. DimUnit1 .eq. DimUnit3 .and. DimUnit2 .eq. DimUnit4)Then
                II=II+1
                Deallocate(DimValue3)
                Deallocate(DimValue4)
                GO To 601
            else
                Write(6,*) "Spatial coordinates of",AllNCDFfile(ii),"mismatch with watershed file's aspatial coordinates"
                Write(6,*) "Spatial bounrary of",AllNCDFfile(ii), "doesn't match with watershed netCDF file"
                II=II+1
                Deallocate(DimValue3)
                Deallocate(DimValue4)
                GO To 601
            End if
        End if
        
end subroutine NCChecks
