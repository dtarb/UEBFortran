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
        Subroutine AggregatedOutNum(AggOutControl,outSymbol,AggOutNum)
        integer, parameter:: n=66
        character* 200:: AggOutControl, outSymbol(n),AggoutHeading,filecode,OUTname
        Integer:: AggOutNum,reason
        logical:: matched
        AggOutNum=0

        OPEN(1099,FILE=AggOutControl,STATUS='OLD', ACTION='READ')
        READ(1099,*,iostat=reason)AggoutHeading
45000   Read(1099,*,iostat=reason, end=46000)filecode
        
            If(reason .eq. 0) then
                matched = .false.
                OUTname=ADJUSTL(filecode(1:(SCAN(filecode, ':')-1)))
                CALL lowercase(OUTname,OUTname)
                do i=1,n
                   CALL lowercase(outSymbol(i),outSymbol(i))
               if(OutName .eq. outSymbol(i))then
                    AggOutNum=AggOutNum+1
                    matched = .true.
                    exit
                end if
            end do
            if(.not. matched)then
                write(6,*)'Output code not matched:  ',OUTname
            end if
        End if
        go to 45000
46000   CLOSE(1099)
        End Subroutine AggregatedOutNum
        
        
        
        Subroutine AggOutWSUniqueID(AggOutControl,outSymbol,AggOutNum,AggOutVar,Watershedfile,WatershedVARID,dimlen2,dimlen1,uniqueIDNumber,AggOutVarnum)
        use netCDF
        
        integer, parameter:: n=66
        character*200:: AggOutControl, outSymbol(n),AggoutHeading,OutName,filecode,Watershedfile,WatershedVARID
        integer::AggOutNum,OutNum,Reason,dimlen2,dimlen1,AOutNum
        Logical:: matched,found
        Character*200:: AggOutVar(AggOutNum)
        integer:: AggOutVarnum(AggOutNum)
        integer:: IDVALUEDim,uniqueIDNumber,ncidout,WSVarId
        integer,dimension(:,:),allocatable:: IDVALUE
        integer,dimension(:),allocatable:: IDVALUEAllocate
        IDVALUEDim=dimlen2*dimlen1
        allocate(IDVALUEAllocate(IDVALUEDim))
        allocate(IDVALUE(dimlen1,dimlen2))
        OPEN(10999,FILE=AggOutControl,STATUS='OLD', ACTION='READ')
        READ(10999,*,iostat=reason)AggoutHeading
        OutNum=1
        AOutNum=0
        ! Read until end of file
35000      Read(10999,*,iostat=reason, end=36000)filecode
           If(reason .eq. 0) then
           matched = .false.
           OUTname=ADJUSTL(filecode(1:(SCAN(filecode, ':')-1)))
           CALL lowercase(OUTname,OUTname)
           do i=1,n,1
                CALL lowercase(outSymbol(i),outSymbol(i))
                If(OutNum .le. AggOutNum)then
                    if(OutName .eq. trim(outSymbol(i)))then
                        AOutNum=AOutNum+1
                        AggOutVarnum(AOutNum)=i  !OutName
                        AggOutVar(OutNum)=OutName
                        OutNum=OutNum+1
                        matched = .true.
                        exit
                    end if
                end if
            end do
            if(.not. matched) &
             write(6,*)'Output code not matched:  ',OUTname
        end if
        go to 35000
36000   CLOSE(10999)

        call check(nf90_open(Watershedfile,nf90_nowrite,ncidout))
        call check(nf90_inq_varid(ncidout,WatershedVARID,WSVarId))
        CALL check(nf90_get_var(ncidout,WSVarId,IDValue))
        call check(nf90_close(ncidout)) 

            Do i=1,dimlen1
                Do j=1,dimlen2
                    IDVALUEAllocate((i-1)*dimlen2+j)=IDValue(i,j)
                end do
            end do

        deallocate(IDVALUE)
        uniqueIDNumber=0
        Do i=1,IDVALUEDim
            found=.true.
            Do j=(i+1),IDVALUEDim
                if (IDVALUEAllocate(i) .eq. IDVALUEAllocate(j))then
                    found=.false.
                    exit
                end if
            End do
            if(found .and. IDVALUEAllocate(i) .ne. 0)then
                uniqueIDNumber=uniqueIDNumber+1
            end if
        end do

        deallocate(IDVALUEAllocate)
        End Subroutine AggOutWSUniqueID
         
        Subroutine WSUniqueArray(Watershedfile,WatershedVARID,dimlen1,dimlen2,uniqueIDNumber,UniqueIDArray)
        use netcdf
        
        integer:: dimlen1,dimlen2
        integer, parameter:: n=66
        character*200:: Watershedfile,WatershedVARID
        integer:: IDVALUEDim,uniqueIDNumber,ncidout,WSVarId
        logical:: found
        integer:: UniqueIDArray(uniqueIDNumber)
        integer,dimension(:,:),allocatable:: IDVALUE
        integer,dimension(:),allocatable:: IDVALUEAllocate
        IDVALUEDim=dimlen2*dimlen1
        allocate(IDVALUEAllocate(IDVALUEDim))
        allocate(IDVALUE(dimlen1,dimlen2))
        
        call check(nf90_open(Watershedfile,nf90_nowrite,ncidout))
        call check(nf90_inq_varid(ncidout,WatershedVARID,WSVarId))
        CALL check(nf90_get_var(ncidout,WSVarId,IDValue))
        call check(nf90_close(ncidout)) 

        Do i=1,dimlen1
            Do j=1,dimlen2
                IDVALUEAllocate((i-1)*dimlen2+j)=IDValue(i,j)
            end do
        end do
        
        uniqueIDNumber=0
        Do i=1,IDVALUEDim
            found=.true.
            Do j=(i+1),IDVALUEDim
                if (IDVALUEAllocate(i) .eq. IDVALUEAllocate(j))then
                    found=.false.
                    exit
                end if
            End do
            if(found .and. IDVALUEAllocate(i) .ne. 0)then
                uniqueIDNumber=uniqueIDNumber+1
                UniqueIDArray(uniqueIDNumber)=IDVALUEAllocate(i)
            end if
        end do
    deallocate(IDVALUEAllocate)
    deallocate(IDVALUE)
    End Subroutine  WSUniqueArray
