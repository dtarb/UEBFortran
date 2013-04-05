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
Character*50 Sitexcoordinates(n), Siteycoordinates(n)
Character*200:: svfile
Character*200:: SVname,SiteHeading,SVcode,StateSiteVName(n),StateSiteNCFile,StateSiteV
Real:: StateSiteValue
integer:: reason,IsInputFromNC(m),countinput,count,NumNCFiles(m),isSiteVarFromNC  
integer:: matched,x,totalNC

CHARACTER(200) :: str 
CHARACTER(1) :: delimit1,delimit2,delimit3
integer:: nargs,nargs2,nargs3,nargs4,nargs5
character(200),Allocatable:: words(:)
character(200)::Word7element(7)
Character(200):: SitexcoordinateR, SiteycoordinateR, InputtcoordinateR
Character(100):: VarNameinNCDFR
Character(200),Allocatable:: words1(:),words2(:),words3(:),words4(:),words5(:)
Character(200):: StateSiteFilesR
CHARACTER*50 SiteVarNameinNCDF(n)
Character*50 Sitexcoordinate, Siteycoordinate
Integer:: DefaultDimValues(3)
REAL::RangeMin,RangeMax

delimit1=';'
delimit2=':'
delimit3=','      
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
                        count=count+1
                        read(81,fmt='(A)')str
                        CALL StringSep(str,delimit1,nargs)
                        Allocate(words(nargs))
                        CALL StringSepWord(str,delimit1,nargs,words)
                        Word7element(1:nargs)=words(1:nargs)
                        CALL StringToVarName(nargs,Word7element,delimit2,StateSiteFilesR,SitexcoordinateR,SiteycoordinateR,&
                                             &VarNameinNCDFR,InputtcoordinateR,DEFAULTDIMVALUES,RangeMin,RangeMax,delimit3,.true.)
                        Sitexcoordinates(i)=SitexcoordinateR
                        Siteycoordinates(i)=SiteycoordinateR                                        
                        SiteVarNameinNCDF(i)=VarNameinNCDFR
                        Deallocate(words)
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
  
  subroutine NCChecks(svfile,StateSiteVName,WatershedFile,MaxNumofFile,IsInputFromNC,NCDFContainer,NumNCFiles,totalNC,AllNCDFfile,&
                      Inputxcoordinates,Inputycoordinates,Sitexcoordinates,Siteycoordinates,wsxcoordinate,wsycoordinate)
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
        Character*200:: ALLXcoord(totalNC), ALLYcoord(totalNC)
        integer:: InputDefDimval(totalNC,3)
        REAL*8, dimension(:), allocatable:: DimValue1,DimValue2,DimValue3,DimValue4
        Character*200:: DimUnit1,DimUnit2,DimUnit3,DimUnit4,DimName1,DimName2,DimName3,DimName4
        character(50):: Inputxcoordinates(m),Inputycoordinates(m)
        character(50):: Sitexcoordinates(n),Siteycoordinates(n)
        character(50):: WSXCOORDINATE,WSYCOORDINATE,Sitexcoordinate,Siteycoordinate
        character(200):: STATESITEFILESR,INPUTTCOORDINATER
        CHARACTER(200) :: str 
        CHARACTER(1) :: delimit1,delimit2,delimit3
        integer:: nargs,nargs2,nargs3,nargs4,nargs5
        character(200),Allocatable:: words(:),Words7element(:)
        Character(200):: SitexcoordinateR, SiteycoordinateR !, InputtcoordinateR,StateSiteFilesR
        Character(200):: VarNameinNCDFR
        Character(200),Allocatable:: words1(:),words2(:),words3(:),words4(:),words5(:)
        integer:: DefaultDimValues(3), NCIDOUT
        REAL::RangeMin,RangeMax
        delimit1=';'
        delimit2=':'
        delimit3=','
        InputVName= (/ "Ta     ","Prec   ","v      ","RH     ", &
             "Tmin   ","Tmax   ","Qsi    ","Qg     ","Qli    ", &
             "Qnet   ","Snowalb"  /)
        InputDefDimval=-9999          
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
                        read(8,fmt='(A)')str
                        CALL StringSep(str,delimit1,nargs)
                        Allocate(words(nargs))
                        Allocate(Words7element(7))
                        CALL StringSepWord(str,delimit1,nargs,words)
                        Words7element(1:nargs)=words(1:nargs)                      
                        CALL StringToVarName(nargs,Words7element,delimit2,StateSiteFilesR,SitexcoordinateR,SiteycoordinateR,VarNameinNCDFR,&
                                             &InputtcoordinateR,DefaultDimValues,RangeMin,RangeMax,delimit3,.false.)
                        InputDefDimval(NCNumber,1)=DefaultDimValues(1)
                        InputDefDimval(NCNumber,2)=DefaultDimValues(2)
                        InputDefDimval(NCNumber,3)=DefaultDimValues(3)  
                        AllNCDFfile(NCNumber)=StateSiteFilesR
                        ALLXcoord(NCNumber)=SitexcoordinateR
                        ALLYcoord(NCNumber)=SiteycoordinateR
                        Deallocate(words)
                        Deallocate(Words7element)
                        NCNumber=NCNumber+1
                    end if
                    exit
                end if 
            end do
            if(matched .ne. 1)write(6,*)'Input variable code not matched:  ', SVname
          end if
          go to 3010
4010     CLOSE(8)

        
        If (NCNumber .GE. 1)THEN
            Do i=1,NCNumber
                If(InputDefDimval(i,2) .NE. -9999)THEN ! Get the name of Y-coordinate
                    CALL check(nf90_open(AllNCDFfile(i),NF90_NOWRITE, ncidout))
                    CALL check(nf90_inquire_dimension(ncidout,InputDefDimval(i,1),ALLYcoord(i)))
                    CALL check(nf90_close(ncidout))
                End if
                If(InputDefDimval(i,3) .NE. -9999)THEN ! Get the name of X-coordinate
                    CALL check(nf90_open(AllNCDFfile(i),NF90_NOWRITE, ncidout))
                    CALL check(nf90_inquire_dimension(NCIDout,InputDefDimval(i,2),ALLXcoord(i)))
                    CALL check(nf90_close(ncidout))
                End if  
           End do
        end if
      
        do i=1,m
            if(IsInputFromNC(i) .eq. 1)then
                do j=1,NumNCFiles(i)
                    AllNCDFfile(NCNumber) = NCDFContainer(j,i)
                    ALLXcoord(NCNumber) = Inputxcoordinates(i)
                    ALLYcoord(NCNumber) = Inputycoordinates(i)
                    NCNumber=NCNumber+1
                end do
             end if
        end do
        
        II=1
        !CALL nCDF2DArrayInfo(WatershedFile,dimlen2,dimlen1)
        CALL nCDF2DArrayInfo(WatershedFile,wsxcoordinate,wsycoordinate,dimlen2,dimlen1)

501     If (II .le. (NCNumber-1))Then
        !CALL nCDF2DArrayInfo(AllNCDFfile(II),dimlen4,dimlen3)
        CALL nCDF2DArrayInfo(AllNCDFfile(II),ALLXcoord(II),ALLYcoord(II),dimlen4,dimlen3)
            If (dimlen2 .eq. dimlen4 .and. dimlen1 .eq. dimlen3)Then
                II=II+1
                GO To 501
            else
                Write(6,*) "Spatial coordinates of",AllNCDFfile(ii),"mismatch with watershed file's aspatial coordinates"
                Write(6,*) "Inconsistent number of grid cells in ", AllNCDFfile(ii), " netCDF file"
                II=II+1
                GO To 501
            End if
        End if
        
        allocate(DimValue1(dimlen1))
        allocate(DimValue2(dimlen2))

        CALL FindDimensionLength(WatershedFile,wsxcoordinate,wsycoordinate,Dimlen1,DimLen2)
        !Subroutine SpatialCoordinate('LangtangKholaWatershed.nc','latitude','longitude',DimValue1,DimValue2,DimUnit1,DimUnit2)
        CALL SpatialCoordinate(WatershedFile,wsxcoordinate,wsycoordinate,DimValue1,DimValue2,DimUnit1,DimUnit2,dimlen1,dimlen2)
        II=1
601     If (II .le. (NCNumber-1))Then
        allocate(DimValue3(dimlen3))
        allocate(DimValue4(dimlen4))
        CALL FindDimensionLength(AllNCDFfile(II),ALLXcoord(II),ALLYcoord(II),dimlen3,dimlen4)
        !Subroutine SpatialCoordinate('LangtangKholaWatershed.nc','latitude','longitude',DimValue1,DimValue2,DimUnit1,DimUnit2)
        CALL SpatialCoordinate(AllNCDFfile(II),ALLXcoord(II),ALLYcoord(II),DimValue3,DimValue4,DimUnit3,DimUnit4,dimlen3,dimlen4)
            If (((MAXVAL(DimValue1) - MAXVAL(DimValue3)) .LE. 1E-5) .and. ((MINVAL(DimValue2) - MINVAL(DimValue4)) .LE. 1E-5)&
            &.and. DimUnit1 .eq. DimUnit3 .and. DimUnit2 .eq. DimUnit4)Then
                II=II+1
                Deallocate(DimValue3)
                Deallocate(DimValue4)
                GO To 601
            else
                Write(6,*) "Spatial coordinates of",AllNCDFfile(ii),"mismatch with watershed file's spatial coordinates"
                Write(6,*) "Spatial bounrary of",AllNCDFfile(ii), "doesn't match with watershed netCDF file"
                II=II+1
                Deallocate(DimValue3)
                Deallocate(DimValue4)
                GO To 601
            End if
        End if
        deallocate(DimValue1)
        deallocate(DimValue2)
end subroutine NCChecks