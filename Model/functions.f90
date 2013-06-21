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

!==================3-D netcdf file reading starts and reads its diemension,variables etc.===========================================
Subroutine nCDF3DArrayInfo (FILE_NAME,dimname1,dimname2,dimname3,dimlen2,dimlen1,dimlen3)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 3-D netccdf file
!dimname1 (in) name of dimension-x
!dimname2 (in) name of dimension-y
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
!dimlen3 (out) length of time steps
use netcdf
Implicit None
character*50 dimname1,dimname2,dimname3
integer, parameter :: NDIMS = 3
integer:: nDimensions, nVariables, nAttributes, dimlen1, dimlen2, dimlen3, ncidout
integer:: DimID1,DimID2,DimID3
character (50) :: FILE_NAME

!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute  
! get dimension IDs
call check(nf90_inq_dimid(ncidout,dimname1, DimID1))
call check(nf90_inq_dimid(ncidout,dimname2, DimID2))
call check(nf90_inq_dimid(ncidout,dimname3, DimID3))
! get information based on IDs
call check(nf90_inquire_dimension(ncidout, DimID1,len=dimlen1))       ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout, DimID2,len=dimlen2))       ! information about dimensionID 2         
call check(nf90_inquire_dimension(ncidout, DimID3,len=dimlen3))         
CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF3DArrayInfo
!==================3-D netcdf file reading ends  ==================================================================================

!!==================2-D netcdf file reading starts and reads a single value at a time==============================================
Subroutine nCDF2DReadInteger(file_name,varname,SingleArray,j,i,wsxcoordinate,wsycoordinate)
!Task: provides the value of a variable for a particular x- and y-coordinate
!file_name (in) 2-D netccdf file
!varname (in) variable name 
!SingleArray (output) array that holds the value
!i (in) partcular y-coordinate
!j (in) partcular x-coordinate
use netcdf
Implicit None
integer, parameter:: NF90_BYTEs = 1, NF90_CHARs = 2, NF90_SHORTs = 3, NF90_INTs = 4, NF90_FLOATs = 5, NF90_DOUBLEs = 6
integer, parameter :: NDIMS = 2 
integer:: nDimensions, nVariables, nAttributes, VarID, ncidout
character (50) :: varname
character (50) :: FILE_NAME
integer :: start(NDIMS), count(NDIMS)
integer:: i, j
Integer, dimension(1):: SingleArray,arrayint
byte inbyte(1)
character inchar(1)
integer*2 inshort(1) 
integer*4 ininteger(1)
real*4 inreal(1)
real*8 indouble(1)
integer:: vartype
character*50:: wsxcoordinate,wsycoordinate
integer:: dimid1,dimid2,varndim
integer,allocatable::vardmids(:)

count = (/ 1, 1 /)
start = (/ 1, 1/)
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute

call check(nf90_inq_varid(ncidout, varname, VarID))
!call check(nf90_inquire_variable(ncidout, 3, varname1))                 ! information about variableID 3
!  Inquire the type of the data
call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,ndims=varndim))
Allocate(vardmids(varndim))
call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,dimids=vardmids))
call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
call check(nf90_inq_dimid(ncidout,wsycoordinate,dimid1))
call check(nf90_inq_dimid(ncidout,wsxcoordinate,dimid2))
! iycoord=j
! jxcoord=i
! tcoord=1, ycoord=2, and xcoord=3 (DEFAULT of UEBGrid)
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html (section 2.4)
if (vardmids(1)==dimid1 .and. vardmids(2)==dimid2)THEN
    start(1) = j
    start(2) = i        
END IF

! tcoord=1, ycoord=3, and xcoord=2
if (vardmids(2)==dimid1 .and. vardmids(1)==dimid2)THEN
    start(1) = i
    start(2) = j
END IF

if(vartype .eq. NF90_BYTE)then
    call check(nf90_get_var(ncidout, VarID, inbyte, start, count))        
    arrayint(1)=INT(inbyte(1))
elseif(vartype .eq. NF90_CHARs)then
    call check(nf90_get_var(ncidout, VarID, inchar, start, count))
    Write (6, *) "Error: A site variable can't be a character in ", file_name
elseif(vartype .eq. NF90_SHORTs)then
    call check(nf90_get_var(ncidout, VarID, inshort, start, count))
    arrayint(1)=INT(inshort(1))
elseif(vartype .eq. NF90_INTs)then
    call check(nf90_get_var(ncidout, VarID, ininteger, start, count))
    arrayint(1)=INT(ininteger(1))
elseif(vartype .eq. NF90_FLOATs)then 
    call check(nf90_get_var(ncidout, VarID, inreal, start, count))
    arrayint(1)=INT(inreal(1)) 
elseif(vartype .eq. NF90_DOUBLEs)then 
    call check(nf90_get_var(ncidout, VarID, indouble, start, count))
    arrayint(1)=INT(indouble(1))
END IF   

SingleArray(1)=arrayint(1)  

CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))   
return                                      ! Closing the netcdf file
end subroutine nCDF2DReadInteger
!!==================2-D netcdf file reading ends  ==================================================================================

!CALL nCDF2DReadReal(file_name,SiteVarNameinNCDF(i),SingleArray, jlon,ilat,Sitexcoordinates(i),Siteycoordinates(i))
Subroutine nCDF2DReadReal(file_name,varname,SingleArray,j,i,Sitexcoord,Siteycoord)
!Task: provides the value of a variable for a particular x- and y-coordinate
!file_name (in) 2-D netccdf file
!varname (in) variable name 
!SingleArray (output) array that holds the value
!i (in) partcular y-coordinate
!j (in) partcular x-coordinate
use netcdf
Implicit None
integer, parameter:: NF90_BYTEs = 1, NF90_CHARs = 2, NF90_SHORTs = 3, NF90_INTs = 4, NF90_FLOATs = 5, NF90_DOUBLEs = 6
integer, parameter :: NDIMS = 2 
integer:: nDimensions, nVariables, nAttributes, VarID, ncidout
character (50) :: varname
character (50) :: FILE_NAME
integer :: start(NDIMS), count(NDIMS)
integer:: i, j
REAL, dimension(1):: SingleArray, arrayreal
byte inbyte(1)
character inchar(1)
integer*2 inshort(1) 
integer*4 ininteger(1)
real*4 inreal(1)
real*8 indouble(1)
integer:: vartype
character*50:: Sitexcoord,Siteycoord
integer:: dimid1,dimid2,varndim
integer,allocatable::vardmids(:)
count = (/ 1, 1 /)
start = (/ 1, 1/)
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute

call check(nf90_inq_varid(ncidout, varname, VarID))
!call check(nf90_inquire_variable(ncidout, 3, varname1))                 ! information about variableID 3
!  Inquire the type of the data
call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,ndims=varndim))
Allocate(vardmids(varndim))
call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,dimids=vardmids))
call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
call check(nf90_inq_dimid(ncidout,Siteycoord,dimid1))
call check(nf90_inq_dimid(ncidout,Sitexcoord,dimid2))
!
! tcoord=1, ycoord=2, and xcoord=3 (DEFAULT of UEBGrid)
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html (section 2.4)
if (vardmids(1)==dimid1 .and. vardmids(2)==dimid2)THEN
    start(1) = j
    start(2) = i        
END IF

! tcoord=1, ycoord=3, and xcoord=2
if (vardmids(2)==dimid1 .and. vardmids(1)==dimid2)THEN
    start(1) = i
    start(2) = j
END IF

call check(nf90_get_var(ncidout,VarID, arrayreal, start, count))      
if(vartype .eq. NF90_BYTE)then
    call check(nf90_get_var(ncidout,VarID, inbyte, start, count)) 
    arrayreal(1)=real(inbyte(1))
elseif(vartype .eq. NF90_CHARs)then
    call check(nf90_get_var(ncidout,VarID, inchar, start, count)) 
    Write (6, *) "Error: A site variable can't be a character in ", file_name
elseif(vartype .eq. NF90_SHORTs)then
    call check(nf90_get_var(ncidout,VarID, inshort, start, count)) 
    arrayreal(1)=real(inshort(1))
elseif(vartype .eq. NF90_INTs)then
    call check(nf90_get_var(ncidout,VarID, ininteger, start, count)) 
    arrayreal(1)=real(ininteger(1))
elseif(vartype .eq. NF90_FLOATs)then 
    call check(nf90_get_var(ncidout,VarID, inreal, start, count)) 
    arrayreal(1)=real(inreal(1)) 
elseif(vartype .eq. NF90_DOUBLEs)then 
    call check(nf90_get_var(ncidout,VarID, indouble, start, count)) 
    arrayreal(1)=real(indouble(1)) 
END IF   

SingleArray(1)=arrayreal(1)  !  loop over all values

CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))                                         ! Closing the netcdf file

end subroutine nCDF2DReadReal
!!==================2-D netcdf file reading ends  ==================================================================================


!!==================3-D netcdf file reading starts and reads timeSTEPS=============================================
Subroutine nCDF3DtimeRead(file_name,rec_name,time_pos,time_val,numtimeStep,syear,smonth,sday,daysstring)
!!  file_name  (in)  Netcdf file to read from
!!  time_pos  (in)  NetCDF time position to read value from
!!  time_val  (out)  time value read from netcdf file
!!  numtimeStep  (out)  Number of time steps in netcdf file
!!  syear   (out)  year from the "days since YYYY..." netcdf time units attribute
!!  smonth  (out)  month from the "days since YYYY-MM..."  time units attribute
!!  sday    (out)  day from the "days since YYYY-MM-DD..." time units attribute
!!  rec_name (in) name of time-dimension
!!  Note - the first call of this should have time_pos as 1 to determine numtimeStep. 
!!  Thereafter it is responsibility of calling program to not input a time_pos greater than numtimeStep

use netcdf
Implicit None
integer, parameter:: NF90_BYTEs = 1, NF90_CHARs = 2, NF90_SHORTs = 3, NF90_INTs = 4, NF90_FLOATs = 5, NF90_DOUBLEs = 6
! integer, parameter :: NDIMS = 3 
integer:: numtimeStep
character (50) :: FILE_NAME, Rec_name
! integer :: start(NDIMS), count(NDIMS)
integer :: VarId, isep1, isep2, isep3, ncidout
integer start1(1),count1(1),time_pos
character (len = *), parameter :: UNITS = "units"
DOUBLE precision:: time_val(1),Time_real(1)
integer, parameter :: MAX_ATT_LEN = 100
character*(MAX_ATT_LEN) :: time_units_in
integer :: att_len
character*50:: daysstring,sincestring,datestring
integer syear,smonth,sday
LOGICAL:: is_numeric
byte inbyte(1)
character inchar(1)
integer*2 inshort(1) 
integer*4 ininteger(1)
real*4 inreal(1)
real*8 indouble(1)
integer:: vartype,dimID
!Open the file and see whats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                     
call check(nf90_inq_varid(ncidout, Rec_name, VarId))                    ! information about variableID for a given VariableName
call check(nf90_get_att(ncidout, Varid, UNITS, time_units_in))  
call check(nf90_inquire_attribute(ncidout, Varid, UNITS, len = att_len))
!call check(nf90_inq_dimid(ncidout, REC_NAME, time_dimid))
! TODO  fix up so that the rec_name comes from output of an appropriate function
!  we are following the logic that the 3rd dimension no matter what its name is is the time dimension
call check(nf90_inq_dimid(ncidout, Rec_name, DimID))
call check(nf90_inquire_dimension(ncidout, DimID, Rec_name, numtimeStep))       ! information about dimensionID 3
!days since 2010-03-10T00:00
!allocate(time_in(dimlen3))
CALL check(nf90_sync(ncidout))
count1(1) = 1
start1(1) = time_pos
if(time_pos .gt. numtimeStep)then  !  here asked for a time step not in ncdf file
  write(6,*)"Warning - in nCDF3DtimeRead a time value was requested greater than contents of file"
  start1(1) = numtimeStep
endif

call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
  
if(vartype .eq. NF90_BYTE)then
    call check(nf90_get_var(ncidout,VarId,inbyte,start1,count1))   
    Time_real(1)=DBLE(inbyte(1))
elseif(vartype .eq. NF90_CHARs)then
    call check(nf90_get_var(ncidout,VarId,inchar,start1,count1))   
    Write (6, *) "Error: A site variable can't be a character in ", file_name
elseif(vartype .eq. NF90_SHORTs)then
    call check(nf90_get_var(ncidout,VarId,inshort,start1,count1))
    Time_real(1)=DBLE(inshort(1))
elseif(vartype .eq. NF90_INTs)then
    call check(nf90_get_var(ncidout,VarId,ininteger,start1,count1))
    Time_real(1)=DBLE(ininteger(1))
elseif(vartype .eq. NF90_FLOATs)then 
    call check(nf90_get_var(ncidout,VarId,inreal,start1,count1))
    Time_real(1)=DBLE(inreal(1)) 
elseif(vartype .eq. NF90_DOUBLEs)then 
    call check(nf90_get_var(ncidout,VarId,indouble,start1,count1))
    Time_real(1)=DBLE(indouble(1))
END IF   
    time_val(1)=Time_real(1) !  loop over all values
   
call check(nf90_close(ncidout))                                         ! Closing the netcdf file

read(time_units_in,*)daysstring,sincestring,datestring
isep1=1
23  if(is_numeric(datestring(isep1:isep1)))then
      isep1=isep1+1
    go to 23
    endif
    isep2=isep1+1
24  if(is_numeric(datestring(isep2:isep2)))then
      isep2=isep2+1
    go to 24
    endif
    isep3=isep2+1
25  if(is_numeric(datestring(isep3:isep3)))then
      isep3=isep3+1
    go to 25
    endif
!  Here we have found positions that separate numeric year, month, day
if(isep1+1 .ge. isep2 .or. isep2+1 .ge. isep3)then
  write(6,*)"Malformed year month date in netcdf file time unit in ", File_name
endif
if ((daysstring .eq. 'hours') .or. (daysstring .eq. 'hour'))then
    time_val=time_val/24
endif
if ((daysstring .eq. 'days') .or. (daysstring .eq. 'day'))then
    time_val=time_val
endif
read(datestring(1:(isep1-1)),*)syear
read(datestring((isep1+1):(isep2-1)),*)smonth
read(datestring((isep2+1):(isep3-1)),*)sday


! TODO Had trouble getting the code below to work - fix it
!call lowercase(daysstring,instring)
!if(instring .ne. "days")then
!write(6,*)"NetCDF units have to be 'days since'.  The file has them as ",time_units_in
!endif
!call lowercase(sincestring,instring)
!if(instring .ne. "since")write(6,*)"NetCDF units have to be 'days since'.  The file has them as ",time_units_in

end subroutine nCDF3DtimeRead
!!==================3-D netcdf file reading ends  ==================================================================================

!==================2-D netcdf file reading starts and reads its diemension,variables etc.===========================================
Subroutine nCDF2DArrayInfo(FILE_NAME,dimname1,dimname2,dimlen2,dimlen1)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 2-D netccdf file
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
!dimname1 (in) name of dimension-x
!dimname2 (in) name of dimension-y
use netcdf
Implicit None
!character (len = *), parameter :: LAT_NAME = "ycoord"
!character (len = *), parameter :: LON_NAME = "xcoord"
integer, parameter :: NDIMS = 2
integer:: dimlen1, dimlen2, ncidout, DimID1, DimID2
character (200) :: FILE_NAME
Character (50) :: dimname1,dimname2
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file
! get dimension IDs
call check(nf90_inq_dimid(ncidout,dimname1, DimID1))
call check(nf90_inq_dimid(ncidout,dimname2, DimID2))
call check(nf90_inquire_dimension(ncidout, DimID1,len=dimlen1))       ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout, DimID2,len=dimlen2))       ! information about dimensionID 2                    
!call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute  
!call check(nf90_inq_dimid(ncidout, LAT_NAME, y_dimid))
!call check(nf90_inquire_dimension(ncidout, y_dimid, y_NAME, dimlen1))       ! information about dimensionID 3
!call check(nf90_inq_dimid(ncidout, LON_NAME, x_dimid))
!call check(nf90_inquire_dimension(ncidout, x_dimid, x_NAME, dimlen2))       ! information about dimensionID
CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF2DArrayInfo
!==================2-D netcdf file reading ends  ==================================================================================
!==================2-D netcdf file reading starts and reads its diemension,variables etc.===========================================
Subroutine nCDF2DArrayInfo2(FILE_NAME,dimname2,dimname1,dimlen2,dimlen1,WatershedVARID,WsMissingValues,WsFillValues)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 2-D netccdf file
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
!WatershedVARID (in) vaiable ID (in this case watershed variable ID)
!WsMissingValues (out) missing value attribute for a variable in a netCDF 
!WsFillValues (out) missing value attribute in a netCDF 
!dimname1 (in) name of dimension-x
!dimname2 (in) name of dimension-y
use netcdf
Implicit None
integer, parameter :: NDIMS = 2
integer ::  iii
integer:: dimlen1, dimlen2,WSVarId
character (50) :: dimname1,  dimname2
character (200) :: FILE_NAME
character (200) :: WatershedVARID
character (len = *), parameter :: missing_value = "missing_value",fillvalue="_FillValue"
integer::numAtts, ncidout
character (len = 50):: AttName
real::WSMissingValues,wsfillvalues
integer:: DimID1,DImID2
!Open the file and see what's inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file
! get dimension IDs
! TODO - the logic below seems unnecessarily complex.  
call check(nf90_inq_dimid(ncidout,dimname1,DimID1))
call check(nf90_inq_dimid(ncidout,dimname2,DimID2))
if (DimID1==1 .and. dimID2==2)THEN
    call check(nf90_inquire_dimension(ncidout,DimID1,len=dimlen1))       ! Information about dimensionID 1
    call check(nf90_inquire_dimension(ncidout,DimID2,len=dimlen2))       ! information about dimensionID 2  
END IF
if (DimID1==2 .and. dimID2==1)THEN
    call check(nf90_inquire_dimension(ncidout,DimID1,len=dimlen2))       ! Information about dimensionID 1
    call check(nf90_inquire_dimension(ncidout,DimID2,len=dimlen1))       ! information about dimensionID 2  
END IF
!  DGT suggests making dimname and dimlen arrays and replacing the 10 lines above by the following
! do idim=1,2
! call check(nf90_inq_dimid(ncidout,dimname(i),DimID))
! call check(nf90_inquire_dimension(ncidout,DimID,len=dimlen(DimID)))
! enddo
call check(nf90_inq_varid(ncidout,WatershedVARID,WSVarId))
CALL check(nf90_inquire_variable(ncidout,WSVarId,natts=numAtts))

! The logic below will work if missing or fill value is not the last attribute in the list.
!DO iii=1,numAtts
!CALL check(nf90_inq_attname(ncidout,WSVarId,iii,AttName))
!    if(AttName .eq. missing_value)THEN
!        CALL check(nf90_get_att(ncidout,WSVarId,missing_value,WsMissingValues)) 
!    ELSE
!        WsMissingValues=0
!    end IF
!    if(AttName .eq. fillvalue)THEN
!        CALL check(nf90_get_att(ncidout,WSVarId,fillvalue,WsFillValues))  
!    ELSE
!        WsFillValues=0
!    end IF
!END DO
WsMissingValues=-9999
DO iii=1,numAtts
  CALL check(nf90_inq_attname(ncidout,WSVarId,iii,AttName))
  if(AttName .eq. missing_value)then
    CALL check(nf90_get_att(ncidout,WSVarId,missing_value,WsMissingValues))
    exit
  endif
enddo
WsFillValues=-9999
DO iii=1,numAtts
  CALL check(nf90_inq_attname(ncidout,WSVarId,iii,AttName))
  if(AttName .eq. fillvalue)then
    CALL check(nf90_get_att(ncidout,WSVarId,fillvalue,WsFillValues))
    exit
  endif
enddo
CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF2DArrayInfo2
!==================2-D netcdf file reading ends  ==================================================================================

!==================Start checking each of the netCDf command functions  ============================================================
subroutine check(status)
!Task: converts year month date and hours (yyyy mm dd hh.hh) into a simgle double precision julian date
use netcdf
implicit none
integer, intent ( in) :: status
if(status /= nf90_noerr) then
print *, trim(nf90_strerror(status))
READ *
stop 2
end if
end subroutine check
!================== end checking eavh of the netCDf command functions  ====

!================================================================
SUBROUTINE JULDAT (I,M,K,H,TJD)
!THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND
!time.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT time VALUE
!CAN BE IN ANY UT-LIKE time SCALE (UTC, UT1, TT, ETC.) - OUTPUT
!JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
!VAN FLANDERN.
!SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
!I = YEAR (IN)
!M = MONTH NUMBER (IN)
!K = DAY OF MONTH (IN)
!H = UT HOURS (IN)
!TJD = JULIAN DATE (OUT)
Implicit None
DOUBLE PRECISION H,TJD,JD
Integer:: I,M,K
!JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
JD = K-32075+1461*(I+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12-3*((I+4900+(M-14)/12)/100)/4
TJD = JD - 0.5D0 + H/24.D0
RETURN
END
!================================================================
!================================================================
SUBROUTINE CALDAT (TJD,I,M,K,H)
!THIS SUBROUTINE COMPUTES CALENDAR DATE AND time, GIVEN JULIAN
!DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE time SCALE
!(UTC, UT1, TT, ETC.) - OUTPUT time VALUE WILL HAVE SAME BASIS.
!OUTPUT CALENDAR DATE WILL BE GREGORIAN.  ALGORITHM BY FLIEGEL AND
!VAN FLANDERN.
!SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
!
!TJD = JULIAN DATE (IN)
!I = YEAR (OUT)
!M = MONTH NUMBER (OUT)
!K = DAY OF MONTH (OUT)
!H = UT HOURS (OUT)
Implicit None
DOUBLE PRECISION TJD,H,DJD,DMOD,JD
Integer L,N,I,M,K
DJD = TJD + 0.5D0
JD = DJD
H = DMOD (DJD,1.D0)*24 ! 24.D0
!JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
L = JD + 68569
N = 4*L/146097
L = L - (146097*N+3)/4
!I=YEAR, M=MONTH, K=DAY
I = 4000*(L+1)/1461001
L = L - 1461*I/4 + 31
M = 80*L/2447
K = L - 2447*M/80
L = M / 11
M = M + 2 - 12*L
I = 100*(N-49) + I + L
RETURN
END 
!================================================================

!====================Sorrting function===========================
SUBROUTINE SSORT (X, IY, N)
Implicit None
!Source: http://www.personal.psu.edu/jhm/f90/lectures/28.html

!      Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      IY - array to be carried with X (all swaps of X elements are
!          matched in IY .  After the sort IY(J) contains the original
!          postition of the value X(J) in the unsorted X array.
!      N - number of values in array X to be sorted
!      KFLAG - Not used in this implementation

      INTEGER N

      double precision X(1:N)
      INTEGER IY(N)

      double precision TEMP
      INTEGER I, ISWAP(1), ITEMP, ISWAP1
      INTRINSIC MAXLOC

      DO 47 I=1,N-1

     ISWAP=MAXLOC(X(I:N))
     ISWAP1=ISWAP(1)+I-1
     IF(ISWAP1.NE.I) THEN
        TEMP=X(I)
        X(I)=X(ISWAP1)
        X(ISWAP1)=TEMP
        ITEMP=IY(I)
        IY(I)=IY(ISWAP1)
        IY(ISWAP1)=ITEMP
         ENDIF
47  CONTINUE
      RETURN
      END
      
!==================Sorting function end here====================================

!!==================3-D netcdf file reading ends  ==================================================================================

! function to convert string to lowercase obtained from http://www.gbenthien.net/strings/index.html 4/20/12
! No copyright indicated on web page so presumed public domain
subroutine lowercase(str,lcstr)
!Task: converts a upper case (capital)/mixed string to a lower case string
!Str (in) string that needs to be converted to loer class
!lcstr (out) string after lower class conversion
Implicit None

character (len=*):: str
character (len=len_trim(str)):: lcstr
integer :: i, iav, ilen, ioffset, iqc, iquote

ilen=len_trim(str)
ioffset=iachar('A')-iachar('a')
iquote=0
lcstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('A') .and. iav <= iachar('Z')) then
    lcstr(i:i)=achar(iav-ioffset)
  else
    lcstr(i:i)=str(i:i)
  end if
end do
return
end

FUNCTION is_numeric(string)
  IMPLICIT NONE
  CHARACTER(len=*), INTENT(IN) :: string
  LOGICAL :: is_numeric
  REAL :: x
  INTEGER :: e
  READ(string,*,IOSTAT=e) x
  is_numeric = e == 0
END FUNCTION is_numeric

!!!================== Spatial coordinates of netCDF file  =============================================
Subroutine FindDimensionLength(File_name,DimName1,DimName2,Dimlen1,Dimlen2)
use netcdf
Implicit None
character (50) :: FILE_NAME
character (50):: DimName1,DimName2
integer:: ncidout,dimlen1,dimlen2,dimID1,DImID2

!Open the file and see whats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                ! open the netcdf file    
call check(nf90_inq_dimid(ncidout,dimname1, DimID1))
call check(nf90_inq_dimid(ncidout,dimname2, DimID2))
call check(nf90_inquire_dimension(ncidout, DimID1,len=dimlen1))        ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout, DimID2,len=dimlen2))        ! information about dimensionID 2  

End Subroutine FindDimensionLength

!!!================== Spatial coordinates of netCDF file  =============================================
Subroutine SpatialCoordinate(File_name,DimName1,DimName2,DimValue1,DimValue2,DimUnit1,DimUnit2,Dimlen1,DimLen2)
!Subroutine SpatialCoordinate('LangtangKholaWatershed.nc','latitude','longitude',DimValue1,DimValue2,DimUnit1,DimUnit2)
! this will give us spatial coornate names, their unit and value
! file_name  (in)  Netcdf file to read from
! dimlen1 (out) length of dimension (lon) 1
! dimlen2(out) length of dimension (lat) 2
! DimName1(in) name of dimension (lon) 1
! DimName2(in) name of dimension (lat) 2
! DimValue1(out) values in dimension (lon) 1
! DimValue2(out) values in dimension (lat) 1
! DimUnit1(out) unit of dimension (lon) 1
! DimUnit2(out) unit of dimension (lat) 2
use netcdf
Implicit None
integer, parameter:: NF90_BYTEs1 = 1, NF90_CHARs1 = 2, NF90_SHORTs1 = 3, NF90_INTs1 = 4, NF90_FLOATs1 = 5, NF90_DOUBLEs1 = 6
integer, parameter:: NF90_BYTEs2 = 1, NF90_CHARs2 = 2, NF90_SHORTs2 = 3, NF90_INTs2 = 4, NF90_FLOATs2 = 5, NF90_DOUBLEs2 = 6
character (50) :: FILE_NAME
character (50):: DimName1,DimName2
integer, parameter :: MAX_ATT_LEN = 100
character(100) :: DimUnit1,DimUnit2
integer:: ncidout,dimlen1,dimlen2
character (len = *), parameter :: UNITS = "units",missing_value = "missing_value"
integer :: varid
integer:: DimID1,DimID2
integer:: vartype
byte, allocatable:: inbyte1(:)
character, allocatable:: inchar1(:)
integer*2, allocatable:: inshort1(:)
integer*4, allocatable:: ininteger1(:)
real*4, allocatable:: inreal1(:)
real*8, allocatable:: indouble1(:)

byte, allocatable:: inbyte2(:)
character, allocatable:: inchar2(:)
integer*2, allocatable:: inshort2(:)
integer*4, allocatable:: ininteger2(:)
real*4, allocatable:: inreal2(:)
real*8, allocatable:: indouble2(:)
Real*8, allocatable:: DimValue1Org(:),DimValue2Org(:)

Real*8:: dimvalue1(dimlen1),dimvalue2(dimlen2)

!Open the file and see whats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                ! open the netcdf file  

! names are important to create output NetCDFs. So plz donnot delete following two lines  
call check(nf90_inq_dimid(ncidout,dimname1, DimID1))
call check(nf90_inq_dimid(ncidout,dimname2, DimID2))
!call check(nf90_inquire_dimension(ncidout, DimID1,len=dimlen1))        ! Information about dimensionID 1
!call check(nf90_inquire_dimension(ncidout, DimID2,len=dimlen2))        ! information about dimensionID 2  
 
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
               
call check(nf90_inq_varid(ncidout,dimname1, VarId))                    ! information about variableID for a given VariableName
call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit1)) 

if(vartype .eq. NF90_BYTEs1)then
    call check(nf90_get_var(ncidout,Varid,inbyte1)) 
    DimValue1Org=REAL(inbyte1,8)
elseif(vartype .eq. NF90_CHARs1)then
    call check(nf90_get_var(ncidout,Varid,inchar1)) 
    Write (6, *) "Error: A site variable can't be a character in ", file_name
elseif(vartype .eq. NF90_SHORTs1)then
    call check(nf90_get_var(ncidout,Varid,inshort1)) 
    DimValue1Org=REAL(inshort1,8)
elseif(vartype .eq. NF90_INTs1)then
    call check(nf90_get_var(ncidout,Varid,ininteger1)) 
    DimValue1Org=REAL(ininteger1,8)
elseif(vartype .eq. NF90_FLOATs1)then 
    call check(nf90_get_var(ncidout,Varid,inreal1)) 
    DimValue1Org=REAL(inreal1,8) 
elseif(vartype .eq. NF90_DOUBLEs1)then
    call check(nf90_get_var(ncidout,Varid,indouble1))  
    DimValue1Org=REAL(indouble1,8) 
END IF 
DimValue1=DimValue1Org

call check(nf90_inq_varid(ncidout,dimname2, VarId))                    ! information about variableID for a given VariableName 
call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit2)) 

if(vartype .eq. NF90_BYTEs2)then
    call check(nf90_get_var(ncidout,Varid,inbyte2))
    DimValue2Org=REAL(inbyte2,8)
elseif(vartype .eq. NF90_CHARs2)then
    call check(nf90_get_var(ncidout,Varid,inchar2))
    Write (6, *) "Error: A site variable can't be a character in ", file_name
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

!Deallocate(inbyte)
!Deallocate(inchar)
!Deallocate(inshort)
!Deallocate(ininteger)
!Deallocate(inreal)
!Deallocate(indouble)

Deallocate(inbyte1)
Deallocate(inchar1)
Deallocate(inshort1)
Deallocate(ininteger1)
Deallocate(inreal1)
Deallocate(indouble1)

Deallocate(inbyte2)
Deallocate(inchar2)
Deallocate(inshort2)
Deallocate(ininteger2)
Deallocate(inreal2)
Deallocate(indouble2)

CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout)) 
End Subroutine SpatialCoordinate


!!!==================3-D netcdf file reading starts and reads timeSTEPS=============================================
SubRoutine NCDFReadAllTS(file_name,Var_name,AllVal,iycoord,jxcoord,rec,&
                         &Inputxcoord,inputycoord,inputtcoord)
use netcdf
Implicit None
integer, parameter :: NDIMS = 3 
integer, parameter:: NF90_BYTEs = 1, NF90_CHARs = 2, NF90_SHORTs = 3, NF90_INTs = 4, NF90_FLOATs = 5, NF90_DOUBLEs = 6
character (50) :: FILE_NAME, Var_name
integer :: start(NDIMS), count(NDIMS),VarId
integer:: iycoord,jxcoord, rec, ncidout
real, dimension(rec):: AllVal,AllVals  
byte, allocatable:: inbyte(:)
character, allocatable:: inchar(:)
integer*2, allocatable:: inshort(:)
integer*4, allocatable:: ininteger(:)
real*4, allocatable:: inreal(:)
real*8, allocatable:: indouble(:)
integer:: vartype
Character*50:: Inputxcoord,inputycoord,inputtcoord
INTEGER:: DIMID1,DIMID2,DIMID3,ncid,varndim
integer, allocatable:: vardmids(:)

Allocate(inbyte(rec))
Allocate(inchar(rec))
Allocate(inshort(rec))
Allocate(ininteger(rec))
Allocate(inreal(rec))
Allocate(indouble(rec))


!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inq_varid(ncidout, Var_name, VarId))                    ! information about variableID for a given VariableName
call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,ndims=varndim))
Allocate(vardmids(varndim))
call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,dimids=vardmids))
call check(NF90_inquire_variable(ncidout,VarID,xtype=vartype))
call check(nf90_inq_dimid(ncidout,Inputtcoord,dimid1))
call check(nf90_inq_dimid(ncidout,Inputycoord,dimid2))
call check(nf90_inq_dimid(ncidout,Inputxcoord,dimid3))
!call check(nf90_inquire_variable(ncid=ncidout,varid=VarId,dimids=vardmids))

! tcoord=1, ycoord=2, and xcoord=3 (DEFAULT of UEBGrid)
! http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html (section 2.4)
if (dimid1==vardmids(1) .and. dimid2==vardmids(2) .and. dimid3==vardmids(3))THEN
    count = (/ rec, 1, 1 /)
    start = (/ 1, 1, 1 /)
    start(2) = iycoord
    start(3) = jxcoord         
END IF

! tcoord=1, ycoord=3, and xcoord=2
if (dimid1==vardmids(1) .and. dimid2==vardmids(3) .and. dimid3==vardmids(2))THEN
    count = (/ rec, 1, 1 /)
    start = (/ 1, 1, 1 /)
    start(2) = jxcoord
    start(3) = iycoord         
END IF

! tcoord=2, ycoord=1, and xcoord=3
if (dimid1==vardmids(2) .and. dimid2==vardmids(1) .and. dimid3==vardmids(3))THEN
    count = (/ 1, rec, 1 /)
    start = (/ 1, 1, 1 /)
    start(1) = iycoord
    start(3) = jxcoord         
END IF

! tcoord=2, ycoord=3, and xcoord=1
if (dimid1==vardmids(2) .and. dimid2==vardmids(3) .and. dimid3==vardmids(1))THEN
    count = (/ 1, rec, 1 /)
    start = (/ 1, 1, 1 /)
    start(1) = jxcoord 
    start(3) = iycoord        
END IF

! tcoord=3, ycoord=1, and xcoord=2
if (dimid1==vardmids(3) .and. dimid2==vardmids(1) .and. dimid3==vardmids(2))THEN
    count = (/ 1, 1, rec /)
    start = (/ 1, 1, 1 /)
    start(1) = iycoord
    start(2) = jxcoord         
END IF

! tcoord=3, ycoord=2, and xcoord=1
if (dimid1==vardmids(3) .and. dimid2==vardmids(2) .and. dimid3==vardmids(1))THEN
    count = (/ 1, 1, rec /)
    start = (/ 1, 1, 1 /)
    start(1) = jxcoord
    start(2) = iycoord         
END IF

if(vartype .eq. NF90_BYTEs)then
    call check(nf90_get_var(ncidout,VarId, inbyte, start, count))
    AllVals=REAL(inbyte)
elseif(vartype .eq. NF90_CHARs)then
    call check(nf90_get_var(ncidout,VarId, inchar, start, count))
    Write (6, *) "Error: A site variable can't be a character in ", file_name
elseif(vartype .eq. NF90_SHORTs)then
    call check(nf90_get_var(ncidout,VarId, inshort, start, count))
    AllVals=REAL(inshort)
elseif(vartype .eq. NF90_INTs)then
    call check(nf90_get_var(ncidout,VarId, ininteger, start, count))
    AllVals=REAL(ininteger)
elseif(vartype .eq. NF90_FLOATs)then 
    call check(nf90_get_var(ncidout,VarId, inreal, start, count))
    AllVals=REAL(inreal) 
elseif(vartype .eq. NF90_DOUBLEs)then 
    call check(nf90_get_var(ncidout,VarId, indouble, start, count))
    AllVals=REAL(indouble) 
END IF 
AllVal=AllVals

CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))                                         ! Closing the netcdf file

Deallocate(inbyte)
Deallocate(inchar)
Deallocate(inshort)
Deallocate(ininteger)
Deallocate(inreal)
Deallocate(indouble)

End SubRoutine NCDFReadAllTS
!**********************************************************************
! Count of words separated by a delimiter
Subroutine StringSep(str,delims,nargs)
Implicit none
CHARACTER(200) :: str 
CHARACTER(1) :: delims
integer:: nargs
INTEGER :: pos1, pos2, n, i
pos1=1
n=0
DO
    pos2 = INDEX(str(pos1:), delims)
    IF (pos2 == 0) THEN
       n = n + 1
       EXIT
    END IF
    n = n + 1
    pos1 = pos2+pos1
END DO
nargs=n
return
End Subroutine StringSep
!**********************************************

! Obtain words separated by a delimiter
Subroutine StringSepWord(str,delims,nargs,word)
Implicit none
CHARACTER(200) :: str 
CHARACTER(1) :: delims
integer:: nargs
character(200), dimension(nargs) :: word
INTEGER :: pos1, pos2, n, i
pos1 = 1
n = 0

DO i=1,nargs
    pos2 = INDEX(str(pos1:), delims)
    IF (pos2 == 0) THEN
       word(i) = str(pos1:)
       EXIT
    END IF
    word(i) = str(pos1:pos1+pos2-2)
    pos1 = pos2+pos1
END DO
return
END Subroutine StringSepWord 
!**********************************************

! Obtain name of the NC/Index files, Xcoordinate, Ycoordinate, Time and Variable names
Subroutine StringToVarName(nargs,wordgroup,delimit,StateSiteFilesR,SitexcoordinateR,SiteycoordinateR,VarNameinNCDFR,&
                &InputtcoordinateR,DefaultDimValues,RangeMin,RangeMax,delimit3,is2d)
                         
Implicit none
logical:: is2d
integer:: nargs,nargs2,i,nargsrange
Character(200):: SitexcoordinateR, SiteycoordinateR, InputtcoordinateR
Character(100):: VarNameinNCDFR
Character(200),Allocatable:: words(:),wordsrange(:)
character(200):: wordgroup(nargs)
Character(1):: delimit,delimit3
Character(200):: StateSiteFilesR
Integer:: DefaultDimValues(3)
real:: RangeMin, RangeMax

StateSiteFilesR = wordgroup(1)  ! first argument is always the file name
if(is2d)then  ! default dimension values for 2d files (not space time varying)
  DefaultDimValues(1)=-9999
  DefaultDimValues(2)=1
  DefaultDimValues(3)=2
else ! default dimension values for 3d files (space time varying)
  DefaultDimValues(1)=1
  DefaultDimValues(2)=2
  DefaultDimValues(3)=3
endif  
! default range values
RangeMin=-9999
RangeMax=-9999

do i=2,nargs
  CALL StringSep(wordgroup(i),delimit,nargs2)  !  separates the input word group into its parts
  Allocate(words(nargs2))
  CALL StringSepWord(wordgroup(i),delimit,nargs2,words)
  CALL lowercase(words(1),words(1))
  if (words(1) == 'x')Then
      SitexcoordinateR=words(2)
      DefaultDimValues(3)=-9999
  elseif (words(1) == 'y')Then
    SiteycoordinateR=words(2)
    DefaultDimValues(2)=-9999
  elseif (words(1) == 'time')Then
    InputtcoordinateR=words(2)
    DefaultDimValues(1)=-9999
  elseif (words(1) == 'd')Then
    VarNameinNCDFR=words(2)
  elseif (words(1) == 'range')Then
    CALL StringSep(words(2),delimit3,nargsrange)
    Allocate(wordsrange(nargsrange))
    CALL StringSepWord(words(2),delimit3,nargsrange,wordsrange) 
    READ(wordsrange(1),*)RangeMin
    READ(wordsrange(2),*)RangeMax
    deallocate(wordsrange)
  else
    write(6,*) "variable name must be writen along with file name", StateSiteFilesR
  end if
  deallocate(words)
enddo
End Subroutine StringToVarName

!==================Obtain Input Variable Missing and Filling value ===========================================
Subroutine VarMissFill(FILE_NAME,Varname,VarMissingValues,varFillValues)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 2-D netccdf file
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
!WatershedVARID (in) vaiable ID (in this case watershed variable ID)
!WsMissingValues (out) missing value attribute for a variable in a netCDF 
!WsFillValues (out) missing value attribute in a netCDF 
!dimname1 (in) name of dimension-x
!dimname2 (in) name of dimension-y
use netcdf
Implicit None
integer ::  iii
integer:: VarId
character (200) :: FILE_NAME
character (50) :: Varname
character (len = *), parameter :: missing_value = "missing_value",fillvalue="_FillValue"
integer::numAtts, ncidout
character (len = 50):: AttName
real::VarMissingValues,varFillValues

!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file
call check(nf90_inq_varid(ncidout,VarName,VarId))
CALL check(nf90_inquire_variable(ncidout,VarId,natts=numAtts))

DO iii=1,numAtts
CALL check(nf90_inq_attname(ncidout,VarId,iii,AttName))
    if(AttName .eq. missing_value)THEN
        CALL check(nf90_get_att(ncidout,VarId,missing_value,VarMissingValues)) 
    ELSE
        VarMissingValues=-9999
    end IF
    if(AttName .eq. fillvalue)THEN
        CALL check(nf90_get_att(ncidout,VarId,fillvalue,varFillValues))  
    ELSE
        VarMissingValues=-9999
    end IF
END DO
CALL check(nf90_sync(ncidout))
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end Subroutine VarMissFill
!==================Obtain Input Variable Missing and Filling value ===========================================