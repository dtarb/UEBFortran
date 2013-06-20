!  File  snowx.f  
!  Subroutines and function subprograms for the Utah Energy Balance
!  Snow Accumulation and Melt Model.
!  David G. Tarboton, Utah Water Research Laboratory, Utah State University
!  
!  Last Change 9/9/12 to accommodate glacier melt.
!
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

!**************READ Bristow Campbell Parameter values **************

!      SUBROUTINE bcparm(dtbar,bca,bcc,bcfile)
!      CHARACTER*200 bcfile
!      REAL dtbar(12)
!      OPEN(8,FILE=bcfile,STATUS='OLD')
!      READ(8,*)bca,bcc
! 1    read(8,*,end=20)month,val
!      dtbar(month)=val
!      go to 1
! 20   CLOSE(8)
!      RETURN
!      END

!**************READ PARAMETER or state variable VALUES**************

!      SUBROUTINE readvals(vals,vfile,n)
!      CHARACTER*200 vfile
!      REAL vals(1)
!      OPEN(8,FILE=vfile,STATUS='OLD')
!      READ(8,*) (vals(i),i=1,n)
!      CLOSE(8)
!      RETURN
!      END

!********************* READ THE INITIAL CONDITIONS *********************

!     SUBROUTINE FINITIAL(statev,icfile,nxv,YEAR,MONTH,DAY,HOUR,DT)
!     CHARACTER*20 icfile
!     INTEGER YEAR,DAY
!     REAL statev(1)
!     OPEN(8,FILE=icfile,STATUS='OLD')
!      READ(8,*)(statev(i),i=1,nxv)
!      READ(8,*)YEAR,MONTH,DAY,HOUR,DT
!      CLOSE(8)
!      RETURN
!      END

!********************* READ THE site variables  *********************

!      SUBROUTINE readsv(sitev,svfile,nsv,slope,azi,lat)
!      CHARACTER*200 svfile
!      REAL sitev(1),lat
!      OPEN(8,FILE=svfile,STATUS='OLD')
!      READ(8,*)(sitev(i),i=1,nsv)
!      READ(8,*)slope,azi,lat
!      CLOSE(8)
!      RETURN
!      END

!************************** UPDATEtime () ***************************
!                 Update time for each time step

      SUBROUTINE UPDATEtime(YEAR,MONTH,DAY,HOUR,DT)
 
      INTEGER*4 YEAR,DAY,DMON(12),DM, MONTH,I       ! 30/03/2004 ITB 
      INTEGER*4 LYEAR  ! 30/03/2004 ITB  
        real  hour, dt  ! DGT Dec 10, 2004.  Fixing ITB errors 
   
      DATA (DMON(I),I=1,12)/31,28,31,30,31,30,31,31,30,31,30,31/
      HOUR=HOUR+DT
        DM=DMON(MONTH)
!  check for leap years 
      if(month .eq. 2)dm=lyear(year)
 10   continue
      IF(HOUR.GE.24.0) THEN
        HOUR=HOUR-24.0
        DAY=DAY+1
        go to 10
      ENDIF
 20   continue
      IF(DAY.GT.DM) THEN
        DAY=day - dm
        MONTH=MONTH+1
        IF(MONTH.GT.12) THEN
          MONTH=1
          YEAR=YEAR+1
            DM=DMON(MONTH)
          if(month .eq. 2)dm=lyear(year)
          endif
          go to 20
      ENDIF
      RETURN
      END

! ************************** lyear () ***************************
!    function to return number of days in February checking for
!    leap years
      function lyear(year)   
      integer year,lyear
      IF(MOD(YEAR,4).GT.0 .or. &
          (mod(year,100) .eq.0 .and. mod(year,400) .ne. 0)) THEN
! Leap years are every 4 years 
! - except for years that are multiples of centuries (e.g. 1800, 1900)
! - except again that when the century is divisible by 4 (e.g. 1600, 2000)
!   then it is a leap year 
         lyear=28
      ELSE
         lyear=29
      ENDIF
        return
        end

!**************************** atf () ****************************
!    to get the atmospheric transmissivity using the Bristow and Campbell
!    (1984) approach

      Subroutine atf(atff,trange,month,dtbar,a,c)
      DIMENSION dtbar(12)
      b=0.036*exp(-0.154*dtbar(month))
      atff=a*(1-exp(-b*trange**c))
!     write(6,*)trange,month,a,c,dtbar(month),atf
      RETURN
      END


!************************** hourlyRI () **********************
!               To get hourly radiation index

      SUBROUTINE hyri(YEAR,MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT, &
                         HRI,COSZEN)
      INTEGER YEAR,DAY
      REAL LP,LAT1,LAT   
! lp= latitude of equivalent plane in radians
! lat1 = latitude in radians
! lat = latitude in degrees
! a number that speaks for itself - every kissable digit
      PI=3.141592653589793238462643383279502884197169399375105820974944592308
      CRAD=PI/180.0
!  crad = degree to radian conversion factor
!    CONVERT timeS TO RADIANS FROM NOON
      T=(HOUR-12.0)*PI/12.0
      DELT1=DT*PI/12.0
!    CONVERT angles TO RADIANS
      SLOPE1=SLOPE*CRAD
      AZI1=AZI*CRAD
      LAT1=LAT*CRAD
        FJULIAN=FLOAT(JULIAN(year,MONTH,DAY))
      D=CRAD*23.5*SIN((FJULIAN-82.0)*0.017214206321)  
! 0.017214206321 is 2 pi / 365  
! D is solar declination
      LP=ASIN(SIN(SLOPE1)*COS(AZI1)*COS(LAT1) &
        +COS(SLOPE1)*SIN(LAT1))
! LP is latitude of equivalent plane
!     TD=ACOS(-TAN(LAT1)*TAN(D))  This formula abandoned 1/8/04 
!     to make the code work for polar conditions
! TD is half day length, i.e. the time from noon to sunset.  Sunrise is at -TD
      tanprod=TAN(LAT1)*TAN(D)
        if(tanprod .gt. 1.)then
          td=pi  ! This is the condition for perpetual light
        else if(tanprod .lt. -1.)then
          td=0   ! The condition for perpetual night
        else
          td=acos(-tanprod)  ! The condition where there is a sunrise and set
        endif
!  Equivalent longitude offset.  Modified on 1/8/04
!  so that it correctly accounts for shift in longitude if equivalent 
!  plane slope goes over a pole.  Achieved using atan2.
!     DDT=ATAN(SIN(AZI1)*SIN(SLOPE1)/(COS(SLOPE1)*COS(LAT1)
!    *    -COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))
      ddt=atan2(SIN(AZI1)*SIN(SLOPE1),(COS(SLOPE1)*COS(LAT1)-COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))  

!  Now similar logic as before needs to be repeated for equivalent plane
!  but with times reflecting
      TPeqp=TAN(LP)*TAN(D)
! Keep track of beginning and end of exposure of equiv plane to sunlight
      IF(tpeqp .gt. 1.0) THEN
          TPbeg=-pi   ! perpetual light
            tpend=pi
      ELSEif(tpeqp .lt. -1.)then
          TPbeg=0.0  ! perpetual dark
          tpend=0.0 
        else
            tpbeg = -acos(-tpeqp)-ddt
            tpend = acos(-tpeqp)-ddt
      ENDIF 

!  Start and end times for integration of radiation exposure
!  need to account for both horizon, slope and time step
      T1=AMAX1(T,tpbeg,-TD)
      T2=AMIN1(T+DELT1,TD,tpend)
!     write(6,*)t1,t2
      IF(T2.LE.T1) THEN
        HRI=0.0
      ELSE
        HRI=(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT) &
            -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)
!  In the above the divide by cos slope normalizes illumination to per unit horizontal area
      ENDIF
!  There is a special case if tpbeg is less than -pi that occurs in polar regions
!  where a poleward facing slope may be illuminated at night more than the day.
!  Add this in
      if(tpbeg .lt. -pi)then
          T1=AMAX1(T,-tpbeg+2*pi,-TD)
          T2=AMIN1(T+DELT1,TD)
          if(t2 .gt. t1)then
            hri=hri+(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT) &
            -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)
          endif
        endif
! for the purposes of calculating albedo we need a cosine of the
! illumination angle.  This does not have slope correction so back
! this out again.  This is an average over the time step
      COSZEN = HRI*COS(SLOPE1)
!     write(6,*)hri,coszen

      RETURN
      END


!***************************** JULIAN () ****************************

!             To convert the real date to julian date
! YJS The Julian are change to a new version to take the Lean Yean
! into consideration
!    in the old version, there are 365 days each year.
!     FUNCTION JULIAN(MONTH,DAY)
      function julian(yy,mm,dd)
      integer yy,dd
      dimension mmstrt(1:12)
      data (mmstrt(i),i=1,12)/0,31,59,90,120,151,181,212,243,273, &
                             304,334/
      jday = mmstrt(mm) + dd
      ileap = yy - int(yy/4)*4 
      if(ileap.eq.0.and.mm.ge.3) jday = jday + 1
      julian = jday
      return 
      end

!******************** For cloudiness fraction cf *********************
!    Computes the incoming longwave radiation using satterlund Formula
!    Modified 10/13/94 to account for cloudiness.  Emissivity of cloud
!    cover fraction is assumed to be 1.
!
      subroutine cloud(param,atff,cf)
        REAL param(*),atff,cf

      as     = param(28)                        ! Fraction of extraterrestaial radiation on cloudy day,Shuttleworth (1993)  
        bs     = param(29)                      ! (as+bs):Fraction of extraterrestaial radiation on clear day, Shuttleworth (1993) 
 
        IF (atff.GE.(as+bs)) THEN 
                        cf=0                                                             ! Cloudiness fraction
                ELSE IF (atff.LE.as) THEN
                        cf=1
                ELSE 
                        cf    = 1.-(atff -as)/bs
        END IF
   
      RETURN
      END

!************************************ QLIF ()*********************************
!
      subroutine qlif(TA,RH,TK,SBC,Ema,Eacl,cf,qliff )

        REAL TA,RH,TK,SBC,Ema,Eacl,cf,qliff

        TAK  = TA + TK
          EA   = SVPW(TA)*RH
!******************************************************  old option
!    
      EAcl   =  1.08*(1.0-EXP(-(EA/100.0)**((TAK)/2016.0)))   ! Clear sky emissivity
      Ema   =  (cf+(1.-cf)*EAcl)                              ! Emissivity for cloudy sky
      QLIFf =  Ema*SBC*TAK**4                                 ! Incoming longwave 

      RETURN
      END
!************************************NEW SUBROUTINE and FUNCTION*********************************


