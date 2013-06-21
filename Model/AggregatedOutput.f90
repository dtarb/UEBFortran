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
        Implicit None
        integer, parameter:: n=66
        character* 200:: AggOutControl, outSymbol(n),AggoutHeading,filecode,OUTname
        Integer:: AggOutNum,reason, i
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
        
        
        
        Subroutine AggOutWSUniqueID(AggOutControl,outSymbol,AggOutNum,AggOutVar,Watershedfile,WatershedVARID,&
        &uniqueIDNumber,AggOutVarnum,WsMissingValues,WsFillValues,wsxcoordinate,wsycoordinate)
        use netCDF
        Implicit None
        integer, parameter:: NF90_BYTEs = 1, NF90_CHARs = 2, NF90_SHORTs = 3, NF90_INTs = 4, NF90_FLOATs = 5, NF90_DOUBLEs = 6
        integer, parameter:: n=66
        character*200:: AggOutControl, outSymbol(n),AggoutHeading,OutName,filecode,Watershedfile,WatershedVARID
        integer::AggOutNum,OutNum,Reason,dimlen2,dimlen1,AOutNum
        Logical:: matched,found
        Character*200:: AggOutVar(AggOutNum)
        integer:: AggOutVarnum(AggOutNum),i,j
        integer:: IDVALUEDim,uniqueIDNumber,ncidout,WSVarId
        integer,dimension(:,:),allocatable:: IDVALUE
        integer,dimension(:),allocatable:: IDVALUEAllocate
        REAL:: WsMissingValues,WsFillValues
        integer:: IDnumber
        byte,dimension(:,:),allocatable:: inbyte
        character,dimension(:,:),allocatable:: inchar
        integer*2,dimension(:,:),allocatable:: inshort
        integer*4,dimension(:,:),allocatable:: ininteger
        real*4,dimension(:,:),allocatable:: inreal
        real*8,dimension(:,:),allocatable:: indouble
        integer*4,dimension(:,:),allocatable:: arrayint
        integer:: vartype
        character*50:: wsxcoordinate,wsycoordinate
        integer:: dimid1,dimid2
        
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
        call check(NF90_inquire_variable(ncidout,WSVarId,xtype=vartype))
        !CALL check(nf90_get_var(ncidout,WSVarId,IDValue))
        call check(nf90_inq_dimid(ncidout,wsycoordinate,dimid1))
        call check(nf90_inq_dimid(ncidout,wsxcoordinate,dimid2))
        call check(nf90_inquire_dimension(ncidout, DimID1,len=dimlen1))       ! Information about dimensionID 1
        call check(nf90_inquire_dimension(ncidout, DimID2,len=dimlen2))       ! information about dimensionID 2  
        
        if (dimid1==1 .and. dimid2==2)THEN
            allocate(inbyte(dimlen1,dimlen2))
            allocate(inchar(dimlen1,dimlen2))
            allocate(inshort(dimlen1,dimlen2))
            allocate(ininteger(dimlen1,dimlen2))
            allocate(inreal(dimlen1,dimlen2))
            allocate(indouble(dimlen1,dimlen2))
            allocate(arrayint(dimlen1,dimlen2))
            allocate(IDVALUE(dimlen1,dimlen2))
        end if
        
        if (dimid1==2 .and. dimid2==1)THEN
            allocate(inbyte(dimlen2,dimlen1))
            allocate(inchar(dimlen2,dimlen1))
            allocate(inshort(dimlen2,dimlen1))
            allocate(ininteger(dimlen2,dimlen1))
            allocate(inreal(dimlen2,dimlen1))
            allocate(indouble(dimlen2,dimlen1))
            allocate(arrayint(dimlen2,dimlen1))
            allocate(IDVALUE(dimlen2,dimlen1))
        end if

        
        if(vartype .eq. NF90_BYTE)then
            call check(nf90_get_var(ncidout,WSVarId,inbyte))        
            arrayint=INT(inbyte)
        elseif(vartype .eq. NF90_CHARs)then
            call check(nf90_get_var(ncidout,WSVarId, inchar))
            Write (6, *) "Error: A site variable can't be a character in ", Watershedfile
        elseif(vartype .eq. NF90_SHORTs)then
            call check(nf90_get_var(ncidout,WSVarId, inshort))
            arrayint=INT(inshort)
        elseif(vartype .eq. NF90_INTs)then
            call check(nf90_get_var(ncidout,WSVarId, ininteger))
            arrayint=INT(ininteger)
        elseif(vartype .eq. NF90_FLOATs)then 
            call check(nf90_get_var(ncidout,WSVarId, inreal))
            arrayint=INT(inreal) 
        elseif(vartype .eq. NF90_DOUBLEs)then 
            call check(nf90_get_var(ncidout,WSVarId, indouble))
            arrayint=INT(indouble)
        END IF  
        IDValue=arrayint
        call check(nf90_close(ncidout)) 
        
        Deallocate(inbyte)
        Deallocate(inchar)
        Deallocate(inshort)
        Deallocate(ininteger)
        Deallocate(inreal)
        Deallocate(indouble)
        IDVALUEDim=dimlen2*dimlen1
        allocate(IDVALUEAllocate(IDVALUEDim))
        
        if (dimid1==2 .and. dimid2==1)THEN
            Do i=1,dimlen2
                Do j=1,dimlen1
                    IDVALUEAllocate((i-1)*dimlen1+j)=IDValue(i,j)
                end do
            end do
        end if
        
        if (dimid1==1 .and. dimid2==2)THEN
            Do i=1,dimlen1
                Do j=1,dimlen2
                    IDVALUEAllocate((i-1)*dimlen2+j)=IDValue(i,j)
                end do
            end do
        end if
        
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
            if(found .and. ((IDVALUEAllocate(i) .ne. -9999) .and. (IDVALUEAllocate(i) .ne. WsMissingValues) .and. (IDVALUEAllocate(i) .ne. WsFillValues)))then
                uniqueIDNumber=uniqueIDNumber+1
            end if
        end do

        deallocate(IDVALUEAllocate)
        End Subroutine AggOutWSUniqueID
         
        Subroutine WSUniqueArray(Watershedfile,WatershedVARID,dimlen1,dimlen2,uniqueIDNumber,UniqueIDArray,WsMissingValues,WsFillValues,UniqueIDArrayCount,&
                          &wsxcoordinate,wsycoordinate)
        use netcdf
        Implicit None
        integer, parameter:: NF90_BYTEs = 1, NF90_CHARs = 2, NF90_SHORTs = 3, NF90_INTs = 4, NF90_FLOATs = 5, NF90_DOUBLEs = 6
        integer:: dimlen1,dimlen2
        character*200:: Watershedfile,WatershedVARID
        integer:: IDVALUEDim,uniqueIDNumber,ncidout,WSVarId
        logical:: found
        integer:: UniqueIDArray(uniqueIDNumber), i, j, kk, ll, UniqueIDArrayCount(uniqueIDNumber),counts
        integer,dimension(:,:),allocatable:: IDVALUE
        integer,dimension(:),allocatable:: IDVALUEAllocate
        REAL:: WsMissingValues,WsFillValues
        character*50:: wsxcoordinate,wsycoordinate
        integer:: dimid1,dimid2
        
        byte,dimension(:,:),allocatable:: inbyte
        character,dimension(:,:),allocatable:: inchar
        integer*2,dimension(:,:),allocatable:: inshort
        integer*4,dimension(:,:),allocatable:: ininteger
        real*4,dimension(:,:),allocatable:: inreal
        real*8,dimension(:,:),allocatable:: indouble
        integer*4,dimension(:,:),allocatable:: arrayint
        integer:: vartype

        call check(nf90_open(Watershedfile,nf90_nowrite,ncidout))
        call check(nf90_inq_varid(ncidout,WatershedVARID,WSVarId))
        call check(NF90_inquire_variable(ncidout,WSVarId,xtype=vartype))
        
        call check(nf90_inq_dimid(ncidout,wsycoordinate,dimid1))
        call check(nf90_inq_dimid(ncidout,wsxcoordinate,dimid2))
        call check(nf90_inquire_dimension(ncidout, DimID1,len=dimlen1))       ! Information about dimensionID 1
        call check(nf90_inquire_dimension(ncidout, DimID2,len=dimlen2))       ! information about dimensionID 2  
        
        if (dimid1==1 .and. dimid2==2)THEN
            allocate(inbyte(dimlen1,dimlen2))
            allocate(inchar(dimlen1,dimlen2))
            allocate(inshort(dimlen1,dimlen2))
            allocate(ininteger(dimlen1,dimlen2))
            allocate(inreal(dimlen1,dimlen2))
            allocate(indouble(dimlen1,dimlen2))
            allocate(arrayint(dimlen1,dimlen2))
            allocate(IDVALUE(dimlen1,dimlen2))
        end if
        
        if (dimid1==2 .and. dimid2==1)THEN
            allocate(inbyte(dimlen2,dimlen1))
            allocate(inchar(dimlen2,dimlen1))
            allocate(inshort(dimlen2,dimlen1))
            allocate(ininteger(dimlen2,dimlen1))
            allocate(inreal(dimlen2,dimlen1))
            allocate(indouble(dimlen2,dimlen1))
            allocate(arrayint(dimlen2,dimlen1))
            allocate(IDVALUE(dimlen2,dimlen1))
        end if

        
        if(vartype .eq. NF90_BYTE)then
            call check(nf90_get_var(ncidout,WSVarId,inbyte))        
            arrayint=INT(inbyte)
        elseif(vartype .eq. NF90_CHARs)then
            call check(nf90_get_var(ncidout,WSVarId, inchar))
            Write (6, *) "Error: A site variable can't be a character in ", Watershedfile
        elseif(vartype .eq. NF90_SHORTs)then
            call check(nf90_get_var(ncidout,WSVarId, inshort))
            arrayint=INT(inshort)
        elseif(vartype .eq. NF90_INTs)then
            call check(nf90_get_var(ncidout,WSVarId, ininteger))
            arrayint=INT(ininteger)
        elseif(vartype .eq. NF90_FLOATs)then 
            call check(nf90_get_var(ncidout,WSVarId, inreal))
            arrayint=INT(inreal) 
        elseif(vartype .eq. NF90_DOUBLEs)then 
            call check(nf90_get_var(ncidout,WSVarId, indouble))
            arrayint=INT(indouble)
        END IF  
        IDValue=arrayint
        call check(nf90_close(ncidout)) 
        
        Deallocate(inbyte)
        Deallocate(inchar)
        Deallocate(inshort)
        Deallocate(ininteger)
        Deallocate(inreal)
        Deallocate(indouble)
        
        IDVALUEDim=dimlen2*dimlen1
        allocate(IDVALUEAllocate(IDVALUEDim))
        
        if (dimid1==2 .and. dimid2==1)THEN
            Do j=1,dimlen1
                Do i=1,dimlen2
                    IDVALUEAllocate((j-1)*dimlen2+i)=IDValue(i,j)
                end do
            end do
        end if
        
        if (dimid1==1 .and. dimid2==2)THEN
            Do i=1,dimlen1
                Do j=1,dimlen2
                    IDVALUEAllocate((i-1)*dimlen2+j)=IDValue(i,j)
                end do
            end do
        end if

        uniqueIDNumber=0
        Do i=1,IDVALUEDim
            found=.true.
            Do j=(i+1),IDVALUEDim
                if (IDVALUEAllocate(i) .eq. IDVALUEAllocate(j))then
                    found=.false.
                    exit
                end if
            End do
            if(found .and. ((IDVALUEAllocate(i) .ne. -9999) .and. (IDVALUEAllocate(i) .ne. WsMissingValues) .and. (IDVALUEAllocate(i) .ne. WsFillValues)))then
                uniqueIDNumber=uniqueIDNumber+1
                UniqueIDArray(uniqueIDNumber)=IDVALUEAllocate(i)
            end if
        end do
        counts=0
        Do kk=1,uniqueIDNumber
            do ll=1,IDVALUEDim
                if (UniqueIDArray(kk) .eq. IDVALUEAllocate(ll))then
                    counts=counts+1
                end if
            end do
            UniqueIDArrayCount(kk)=counts
            counts=0
        end do
    !write (6,*) UniqueIDArray, UniqueIDArrayCount
    deallocate(IDVALUEAllocate)
    deallocate(IDVALUE)
    End Subroutine  WSUniqueArray