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

!!==================3-D netcdf file reading starts and reads a single value ata time=============================================
Subroutine nCDF3DRead (file_name,Var_name,SingleArray,i,j,rec)

!Task: provides the value of a variable for a particular x- and y-coordinate
!file_name (in) 2-D netccdf file name
!varname (in) variable name 
!SingleArray (output) array that holds the value
!i (in) partcular y-coordinate
!j (in) partcular x-coordinate
Implicit None
use netcdf

integer, parameter :: NDIMS = 3 
integer:: nDimensions, nVariables, nAttributes
character (50) :: FILE_NAME, Var_name
integer :: start(NDIMS), count(NDIMS),VarId
integer:: i, j, rec
real, dimension(1):: SingleArray  

count = (/ 1, 1, 1 /)
start = (/ 1, 1, 1 /)
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute  
call check(nf90_inq_varid(ncidout, Var_name, VarId))                       ! information about variableID for a given VariableName
!call check(nf90_inquire_variable(ncidout, 3, varname1))                ! information about variableID 3

start(1) = j
start(2) = i
start(3) = rec
call check(nf90_get_var(ncidout,VarId, SingleArray, start, count))          ! Read the surface  Elevation Data from the file
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF3DRead
!!==================3-D netcdf file reading ends  ==================================================================================



!==================3-D netcdf file reading starts and reads its diemension,variables etc.===========================================
Subroutine nCDF3DArrayInfo (FILE_NAME,dimlen2,dimlen1,dimlen3)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 3-D netccdf file
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
!dimlen3 (out) length of time steps
Implicit None
use netcdf
character (len = *), parameter :: LAT_NAME = "ycoord"
character (len = *), parameter :: LON_NAME = "xcoord"
character (len = *), parameter :: REC_NAME = "time"
integer, parameter :: NDIMS = 3
integer:: nDimensions, nVariables, nAttributes, dimlen1, dimlen2, dimlen3
character (20) :: dimname1,  dimname2, dimname3
character (50) :: FILE_NAME

!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute  
call check(nf90_inquire_dimension(ncidout, 1, dimname1, dimlen1))       ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout, 2, dimname2, dimlen2))       ! information about dimensionID 2
call check(nf90_inquire_dimension(ncidout, 3, dimname3, dimlen3))       ! information about dimensionID 3  
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF3DArrayInfo
!==================3-D netcdf file reading ends  ==================================================================================

!!==================2-D netcdf file reading starts and reads a single value at a time==============================================
Subroutine nCDF2DRead(file_name,varname,SingleArray,i,j)
!Task: provides the value of a variable for a particular x- and y-coordinate
!file_name (in) 2-D netccdf file
!varname (in) variable name 
!SingleArray (output) array that holds the value
!i (in) partcular y-coordinate
!j (in) partcular x-coordinate
Implicit None
use netcdf

integer, parameter :: NDIMS = 2 
integer:: nDimensions, nVariables, nAttributes, VarID
character (50) :: varname
character (50) :: FILE_NAME
integer :: start(NDIMS), count(NDIMS)
integer:: i, j
REAL, dimension(1):: SingleArray  

count = (/ 1, 1 /)
start = (/ 1, 1/)
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute  
call check(nf90_inq_varid(ncidout, varname, VarID))
!call check(nf90_inquire_variable(ncidout, 3, varname1))                 ! information about variableID 3

start(1) = j
start(2) = i
call check(nf90_get_var(ncidout,VarID, SingleArray, start, count))      ! Read the surface  Elevation Data from the file
call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF2DRead
!!==================2-D netcdf file reading ends  ==================================================================================


!!==================3-D netcdf file reading starts and reads TIMESTEPS=============================================
Subroutine nCDF3DTimeRead(file_name,time_pos,time_val,numTimeStep,syear,smonth,sday)
!!  file_name  (in)  Netcdf file to read from
!!  time_pos  (in)  NetCDF time position to read value from
!!  time_val  (out)  Time value read from netcdf file
!!  numTimeStep  (out)  Number of time steps in netcdf file
!!  syear   (out)  year from the "days since YYYY..." netcdf time units attribute
!!  smonth  (out)  month from the "days since YYYY-MM..."  time units attribute
!!  sday    (out)  day from the "days since YYYY-MM-DD..." time units attribute

!!  Note - the first call of this should have time_pos as 1 to determine numTimeStep. 
!!  Thereafter it is responsibility of calling program to not input a time_pos greater than numTimeStep

Implicit None
use netcdf

! integer, parameter :: NDIMS = 3 
integer:: numTimeStep
character (50) :: FILE_NAME, Rec_name="time"
! integer :: start(NDIMS), count(NDIMS)
integer :: VarId
integer start1(1),count1(1),time_pos
character (len = *), parameter :: UNITS = "units"
Double precision:: time_val(1)
integer, parameter :: MAX_ATT_LEN = 100
character*(MAX_ATT_LEN) :: time_units_in
integer :: att_len
character*50:: daysstring,sincestring,datestring
integer syear,smonth,sday
LOGICAL:: is_numeric

!Open the file and see whats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                     
!call check(nf90_inq_varid(ncidout, Var_name, VarId))                       ! information about variableID for a given VariableName
Varid=3  !  We require that time is the 3rd dimension
call check(nf90_get_att(ncidout, Varid, UNITS, time_units_in))  
call check(nf90_inquire_attribute(ncidout, Varid, UNITS, len = att_len))
!call check(nf90_inq_dimid(ncidout, REC_NAME, time_dimid))
! TODO  fix up so that the rec_name comes from output of an appropriate function
!  we are following the logic that the 3rd dimension no matter what its name is is the time dimension
call check(nf90_inquire_dimension(ncidout, Varid, Rec_name, numTimeStep))       ! information about dimensionID 3
!days since 2010-03-10T00:00
!allocate(time_in(dimlen3))

count1(1) = 1
start1(1) = time_pos
if(time_pos .gt. numTimeStep)then  !  here asked for a time step not in ncdf file
  write(6,*)"Warning - in nCDF3DTimeRead a time value was requested greater than contents of file"
  start1(1) = numTimeStep
endif

call check(nf90_get_var(ncidout,VarId,time_val,start1,count1))          ! Read the first time value from the file
call check(nf90_close(ncidout))                          ! Closing the netcdf file

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
    time_val=Time_val/24
endif
if ((daysstring .eq. 'days') .or. (daysstring .eq. 'day'))then
    time_val=Time_val
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

end subroutine nCDF3DTimeRead
!!==================3-D netcdf file reading ends  ==================================================================================

!==================2-D netcdf file reading starts and reads its diemension,variables etc.===========================================
Subroutine nCDF2DArrayInfo(FILE_NAME,dimlen2,dimlen1)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 2-D netccdf file
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
Implicit None
use netcdf
!character (len = *), parameter :: LAT_NAME = "ycoord"
!character (len = *), parameter :: LON_NAME = "xcoord"
integer, parameter :: NDIMS = 2
integer:: dimlen1, dimlen2
character (20) :: dimname1,  dimname2
character (200) :: FILE_NAME

!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file
call check(nf90_inquire_dimension(ncidout, 1, dimname1, dimlen1))       ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout, 2, dimname2, dimlen2))       ! information about dimensionID 2                    
!call check(nf90_inquire(ncidout, nDimensions, nVariables, nAttributes)) ! returns no od dimension, variable and global attribute  
!call check(nf90_inq_dimid(ncidout, LAT_NAME, y_dimid))
!call check(nf90_inquire_dimension(ncidout, y_dimid, y_NAME, dimlen1))       ! information about dimensionID 3
!call check(nf90_inq_dimid(ncidout, LON_NAME, x_dimid))
!call check(nf90_inquire_dimension(ncidout, x_dimid, x_NAME, dimlen2))       ! information about dimensionID

call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF2DArrayInfo
!==================2-D netcdf file reading ends  ==================================================================================
!==================2-D netcdf file reading starts and reads its diemension,variables etc.===========================================
Subroutine nCDF2DArrayInfo2(FILE_NAME,dimlen2,dimlen1,WatershedVARID,WsMissingValues,WsFillValues)
!Task: provides the length of x- and y-coordinate and time
!FILE_NAME (in) 2-D netccdf file
!Dimlen2 (out) length of y-coordinate
!dimlen1 (out) length of x-coordinate
!WatershedVARID (in) vaiable ID (in this case watershed variable ID)
!WsMissingValues (out) missing value attribute for a variable in a netCDF 
!WsFillValues (out) missing value attribute in a netCDF 
Implicit None
use netcdf
integer, parameter :: NDIMS = 2
integer:: dimlen1, dimlen2,WSVarId
character (20) :: dimname1,  dimname2
character (50) :: FILE_NAME
character (50) :: WatershedVARID
character (len = *), parameter :: missing_value = "missing_value",fillvalue="_FillValue"
integer::numAtts
character (len = 50):: AttName
real::WSMissingValues
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file
call check(nf90_inquire_dimension(ncidout, 1, dimname1, dimlen1))       ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout, 2, dimname2, dimlen2))       ! information about dimensionID 2 
call check(nf90_inq_varid(ncidout,WatershedVARID,WSVarId))
CALL check(nf90_inquire_variable(ncidout,WSVarId,natts=numAtts))

DO iii=1,numAtts
CALL check(nf90_inq_attname(ncidout,WSVarId,iii,AttName))
    if(AttName .eq. missing_value)THEN
        CALL check(nf90_get_att(ncidout,WSVarId,missing_value,WsMissingValues)) 
    ELSE
        WsFillValues=0
    end IF
    if(AttName .eq. fillvalue)THEN
        CALL check(nf90_get_att(ncidout,WSVarId,fillvalue,WsFillValues))  
    ELSE
        WsFillValues=0
    end IF
END DO

call check(nf90_close(ncidout))                                         ! Closing the netcdf file
end subroutine nCDF2DArrayInfo2
!==================2-D netcdf file reading ends  ==================================================================================

!==================Start hecking eavh of the netCDf command functions  ============================================================
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
!TIME.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT TIME VALUE
!CAN BE IN ANY UT-LIKE TIME SCALE (UTC, UT1, TT, ETC.) - OUTPUT
!JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
!VAN FLANDERN.
!SOURCE: http://aa.usno.navy.mil/software/novas/novas_f/novasf_intro.php
!I = YEAR (IN)
!M = MONTH NUMBER (IN)
!K = DAY OF MONTH (IN)
!H = UT HOURS (IN)
!TJD = JULIAN DATE (OUT)
Implicit None
DOUBLE PRECISION H,TJD
Integer:: I,M,K
!JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
JD = K-32075+1461*(I+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12-3*((I+4900+(M-14)/12)/100)/4
TJD = JD - 0.5D0 + H/24.D0
RETURN
END
!================================================================
!================================================================
SUBROUTINE CALDAT (TJD,I,M,K,H)
!THIS SUBROUTINE COMPUTES CALENDAR DATE AND TIME, GIVEN JULIAN
!DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE TIME SCALE
!(UTC, UT1, TT, ETC.) - OUTPUT TIME VALUE WILL HAVE SAME BASIS.
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
DOUBLE PRECISION TJD,H,DJD,DMOD
Integer:: I,M,K
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
!!==================3-D netcdf file reading starts and reads TIMESTEPS=============================================
Subroutine NetCDFTimeArray(file_name,time_out,timelength)
!Task: provides all the values in time diension and the legth of time dimension
!file_name (in) 2-D netccdf file
!time_out (out) array that holds all the time dimension values
!timelength (out) the length of time dimension
Implicit None
use netcdf
integer, parameter :: NDIMS = 3
integer :: dimlen3
character (50) :: FILE_NAME, Var_name="time"
integer :: start(NDIMS), count(NDIMS),VarId
integer:: i
character (len = *), parameter :: UNITS = "units"
integer:: dimlen2,dimlen1  
integer timelength
Double precision :: JDS
Double precision, dimension(:), allocatable :: time_in
double precision time_out(timelength)

!Double precision:: time_in(dimlen3),time_out(dimlen3)
integer, parameter :: MAX_ATT_LEN = 100
character*(MAX_ATT_LEN) :: time_units_in
integer :: att_len
character*50:: CharYear, CharMonth, CharDay, CharHour
real::  shour,as,ys,ms

count = (/ 1, 1, 1 /)
start = (/ 1, 1, 1 /)

!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file  
CALL nCDF3DArrayInfo (FILE_NAME,dimlen2,dimlen1,dimlen3)               
call check(nf90_inq_varid(ncidout, Var_name, VarId))                       ! information about variableID for a given VariableName
call check(nf90_get_att(ncidout, VarId, UNITS, time_units_in))
call check(nf90_inquire_attribute(ncidout, Varid, UNITS, len = att_len))
!days since 2010-03-10T00.00
allocate(time_in(dimlen3))
!allocate(time_out(dimlen3))
call check(nf90_get_var(ncidout,VarId,time_in))          ! Read the surface  Elevation Data from the file
call check(nf90_close(ncidout))                          ! Closing the netcdf file

CharYear=time_units_in(12:15)
CharMonth=time_units_in(17:18)
CharDay=time_units_in(20:21)
Charhour=time_units_in(23:26)
read (CharYear,*)syear
read (CharMonth,*)smonth
read (CharDay,*)sday
read (CharHour,*)shour
as=(14.0-Smonth)/12.0
ys=Syear+4800.-a
ms=Smonth+12.*as-3.
JDS=Sday+(153.*ms+2.)/5.+365.*ys+ys/4.-ys/100.+ys/400.-32045.+(Shour-12.)/24.

do i=1,dimlen3
    time_out(i)=time_in(i)+JDS
end do
Deallocate(time_in)
!Deallocate(time_out)
end subroutine NetCDFTimeArray
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

!!!==================3-D netcdf file reading starts and reads TIMESTEPS=============================================
Subroutine SpatialCoordinate(File_name,dimlen1,dimlen2,DimName1,DimName2,DimValue1,DimValue2,DimUnit1,DimUnit2)
! this will give us spatial coornate names, their unit and value
! file_name  (in)  Netcdf file to read from
! dimlen1 (out) length of dimension (lon) 1
! dimlen2(out) length of dimension (lat) 2
! DimName1(out) name of dimension (lon) 1
! DimName2(out) name of dimension (lat) 2
! DimValue1(out) values in dimension (lon) 1
! DimValue2(out) values in dimension (lat) 1
! DimUnit1(out) unit of dimension (lon) 1
! DimUnit2(out) unit of dimension (lat) 2
Implicit None
use netcdf
character (50) :: FILE_NAME
character (50):: DimName1,DimName2
integer, parameter :: MAX_ATT_LEN = 100
character(100) :: DimUnit1,DimUnit2
integer:: ncidout,dimlen1,dimlen2,len1,len2
character (len = *), parameter :: UNITS = "units",missing_value = "missing_value"
REAL:: DimValue1(dimlen1)
REAL:: DimValue2(dimlen2)
integer :: varid

!Open the file and see whats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file    
                 
call check(nf90_inquire_dimension(ncidout,1,dimname1,len1))       ! Information about dimensionID 1
call check(nf90_inquire_dimension(ncidout,2,dimname2,len2))       ! information about dimensionID 2

Varid=1  !  We require that time is the 1st dimension
call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit1)) 
call check(nf90_get_var(ncidout,Varid,DimValue1))   

Varid=2  !  We require that time is the 2nd dimension
call check(nf90_get_att(ncidout,Varid, UNITS,DimUnit2)) 
call check(nf90_get_var(ncidout,Varid,DimValue2)) 

call check(nf90_close(ncidout)) 

End Subroutine SpatialCoordinate

SubRoutine NCDFReadAllTS(file_name,Var_name,AllVal,iycoord,jxcoord,rec)
use netcdf
Implicit None
integer, parameter :: NDIMS = 3 
character (50) :: FILE_NAME, Var_name
integer :: start(NDIMS), count(NDIMS),VarId
integer:: iycoord,jxcoord, rec
real, dimension(rec):: AllVal  
count = (/ 1, 1, rec /)
start = (/ 1, 1, 1 /)
!Open the file and see hats inside
call check(nf90_open(File_name, nf90_nowrite, ncidout))                 ! open the netcdf file                      
call check(nf90_inq_varid(ncidout, Var_name, VarId))                    ! information about variableID for a given VariableName
start(1) = iycoord
start(2) = jxcoord
call check(nf90_get_var(ncidout,VarId, AllVal, start, count))           ! Read the surface  Elevation Data from the file
call check(nf90_close(ncidout))                                         ! Closing the netcdf file

End SubRoutine NCDFReadAllTS



