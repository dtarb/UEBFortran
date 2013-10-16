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

!
      SUBROUTINE SNOWUEB2(dt,nt,input,sitev,statev, &
        tsprevday, taveprevday, nstepday, param,iflag, &
        cump,cumes,cumEc,cummr,cumGM,outv,mtime &
       ,atff,cf,OutArr,CumEG)

!      PARAMETER(nsv=9,npar=32,nxv=6,niv=8,nov=14) 
      parameter(niv=8) 
      REAL:: Ms,R,MG
      REAL MR,k,lc,ks
	real lans,lang,rkinc,Tac
      integer pflag,ounit,iflag(*),nstepday
    real:: cumGM,CumEG
!	For canopy variables
	REAL int,ieff,Inmax,Mc
	DOUBLE PRECISION Uc
    REAL AvaiableWater
    REAL:: MRX,MSX
      !CHANGES TO ACCOMODATE GLACIER
      real:: WGT ! WGT=WATER EQUIVALENT GLACIER THICKNESS
      real:: WGM ! WGT=WATER EQUIVALENT GLACIER MELTED
!	data pi/3.141592654/ !define PI forangle conversion

! Definitions
!  dt  time step in hours
!  nt number of time steps
! input  -- input forcing
!	  input(1,*) air temperature (C)
!	  input(2,*) precipitation (m/hr)
!	  input(3,*) wind speed  (m/s)
!	  input(4,*) relative humidity (on a 0-1 scale)
!	  input(5,*) incoming short wave  (kJ/m^2/h)
!	  input(6,*) net radiation  (kJ/m^2/h)
!	  input(7,*) Cosine of Zenith angle  
! SITEV -- site variables
!        site variables (1-5)
!        sitev(1)  drift factor  (No detailed information give 1)
!        sitev(2)  air pressure (Pa)
!        sitev(3) ground heat flux  Kj/m^2/hr (3.6 = 1 W/m^2)
!        sitev(4) albedo extinction parameter (m)
! STATEV
!        statev(1)  Snow Energy Content  (KJ/m^2)
!        statev(2)  Snow Water Equivalent (m) relative to T = 0 C solid phase
!        statev(4)  Canopy Snow Water Equivalent (m) relative to T = 0 C solid phase
!        statev(3)  Dimensionless age of snow surface (or albedo - depending on flag 4)
!        statev(5)  Refreezing depth (m) used as refdepth
!        statev(6)  Refreezing depth (m) used as totalrefdepth
!  totalrefdepth is a misnomer.  These are the same quantity - the refreezing depth.  They are repeated because 
!  when refdepth exceeds the diurnal penetration depth it is set to 0 to tell the code to use the regular 
!  surface temperature functions while totalrefdepth is only reset to 0 when there is actual melt generated
!  or energy content becomes negative indicating freezing of all liquid phase.  This ensures that the regular 
!  surface temperature functions persist when there is liquid water present but the refreezing front has penetrated
!  to depth greater than diurnal temperature fluctuations.
!        TsPrevday(1:nstepday)   Surface temperature over the last 24 hours
! 	   TavePrevday(1:nstepday)   Depth average te

! Imperature over the last 24 hours
!
! PARAM  --  snowmelt model parameters (see below)
! iflag  -- flags 
!	   iflag(1) 0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7
!        iflag(2)        no 0 (/yes 1) printing
!        iflag(3)  unit to which to print
!        iflag(4)  how albedo calculations are done (a value 1 means albedo is calculated, otherwise statev(3) is albedo
! 	   iflag(5)  model option for surface temperature approximation 
!              1) the Simple Gradient, almost the same as Original UEB,
!              2) Simple Gradient approach with shallow snow correction. 
!              3) The classical force-restore approach.
!              4) Modified force-restore approach.
! cump,cume,cummr  -- cumulative precipitation (with df), cumulative evaporation, cumulative melt over time step in m
! outv  -- output variables 
!       outv(1)=prain   rain  m/hr
!       outv(2)=ps     snow  m/hr
!       outv(3)=a     albedo
!       outv(4)=qh    sensible heat (kJ/m2/hr) 
!       outv(5)=qe    latent heat (kJ/m2/hr) 
!       outv(6)=e     sublimation m/hr
!       outv(7)=mr    melt outflow m/hr
!       outv(8)=qm    heat advected with melt
!       outv(9)=q     Energy rate of change (kJ/m2/hr) 
!       outv(10)=fm   Mass rate of change (m/hr)
!       outv(11)=tave  Average temperature (C)
!       outv(12)=tsurf  Surface temperature C
!       outv(13)=qnet  Net Radiation (kJ/m2/hr) 
!       outv(14)=smelt   Surface melt  m/hr
!

!
!yjs Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow 
!yjs drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
!yjs snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes 
!yjs back. 


    real input(niv,nt)
    real sitev(*)
    real outv(*)
    real statev(*)
    real param(*)   ! inherit dimension from calling program
    real mtime(*)  
    real tsprevday(*)
    real taveprevday(*)
    integer iTsMethod       !yjs Add model time initialization 09/19/2000
    integer windfl
    REAL::SWIT,SWISM,SWIR,SWIGM
    Real OutArr(53)
	common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

!yjs  Constant data set 
      data to /0.0/        !  Temperature of freezing (0 C)
      data tk /273.15/     !  Temperature to convert C to K (273.15)
      data sbc /2.0747e-7/ !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data hf /333.5/      !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cw /4.18/       !  Water Heat Capacity (4.18 KJ/kg/C)
      data cs /2.09/       !  Ice heat capacity (2.09 KJ/kg/C)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data Rag /287.0/     !  Ideal Gas constant for dry air (287 J/kg/K)
      data k  /0.4/        !  Von Karmans constant (0.4)
      data hff /3600.0/    !  Factor to convert /s into /hr (3600)
      data rhoi /917.0/    !  Density of Ice (917 kg/m^3)
      data rhow /1000.0/   !  Density of Water (1000 kg/m^3)
      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)
      data w1day /0.261799/!  Daily frequency (2pi/24 hr 0.261799 radians/hr) 
      data pi /3.141592654/!  Pi
!yjs  End of constant declaration
    
!  Parameters
    tr=param(1)     !  Temperature above which all is rain (3 C)
    ts=param(2)     !  Temperature below which all is snow (-1 C)
    Ems=param(3)    !  emmissivity of snow (nominally 0.99)
    cg =param(4)    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
    z=param(5)      !  Nominal meas. height for air temp. and humidity (2m)
    zo=param(6)     !  Surface aerodynamic roughness (m)
    rho=param(7)    !  Snow Density (Nominally 450 kg/m^3)
    rhog=param(8)   !  Soil Density (nominally 1700 kg/m^3)
    lc=param(9)     !  Liquid holding capacity of snow (0.05)
    ks=param(10)    !  Snow Saturated hydraulic conductivity (20 m/hr)
    de=param(11)    !  Thermally active depth of soil (0.1 m)
    abg=param(12)   !  Bare ground albedo  (0.25)
    avo=param(13)   !  Visual new snow albedo (0.95)
    anir0=param(14) !  NIR new snow albedo (0.65)
	lans= param(15) !  the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr) (Vinod! 0.36 Ref: Snow and Climate :Richard L Armstrong and Eric Brun ) 
	lang= param(16) !  the thermal conductivity of soil (:9.68 kJ/m/k/hr) (TK of ice or wet soil(2.22~ 3.48W/m/k):Vinod)
	wlf= param(17)  !  Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr) 
	rd1= param(18)  !  Apmlitude correction coefficient of heat conduction (1)
    fstab=param(19) !  Stability correction control parameter 0 = no corrections, 1 = full corrections
	Tref=param(20)  !  Reference temperature of soil layer in ground heat calculation input
	dNewS=param(21) !  The threshold depth of for new snow (0.001 m)
 	gsurf=param(22) !  The fraction of surface melt that runs off (e.g. from a glacier)

! 7 Parameters added for canopy
	EmC      = param(23)			! Emissivity of canopy 
	alpha    = param(24)			! Scattering coefficient for solar radiation
	alphaL   = param(25)		    ! Scattering coefficient for long wave radiation
	G        = param(26)            ! leaf orientation with respect to zenith angle
    Uc       = param(27)		    ! Unloading rate coefficient (Per hour) (Hedstrom and pomeroy, 1998) 
	as       = param(28)			! Fraction of extraterrestaial radiation on cloudy day,Shuttleworth (1993)  
	bs       = param(29)		  	! (as+bs):Fraction of extraterrestaial radiation on clear day, Shuttleworth (1993) 
	Lambda   = param(30)		    ! Ratio of direct atm radiation to diffuse,worked out from Dingman (1993)
	Rimax    = param(31)            ! Maximum value of Richardsion number for stability corretion
	Wcoeff   = param(32)            ! Wind decay coefficient for the forest
  
  ! Some initializations
  rkinc=0.
  tac=0.

!  State Variables - These serve as initial conditions
    Us = statev(1)				    	! Snow Energy Content  (KJ/m^2)
    Ws = statev(2)						! Snow Water Equivalent (m) relative to T = 0 C solid phase
    Wc = statev(4)                      ! Added for canopy

	If(Us.le.0.0) THEN 
	  refDepth = 0.0
	  totalRefDepth = 0.0
	ELSE
	  refDepth      = statev(5)
	  totalRefDepth = statev(6)
    ENDIF

!  Save Old Value 07/23/01     
	Us_old        =  Us
	refDepth_old  =  refDepth

!  Site Variables
    df = sitev(1)     !  Drift factor
    APr= sitev(2)     !  Atmospheric Pressure (Pa)
    qg = sitev(3)     !  Ground heat flux (KJ/m^2/hr)  This is more logically an
                        !  input variable,but is put here because it is never known at each
                        !  time step. Usually it will be assigned a value 0.
    Aep= sitev(4)     !  Albedo extinction parameter to smooth
                        !  transition of albedo when snow is shallow. Depends on Veg. height (m)
! 7 Site Variables added for canopy
	Cc    = sitev(5)   ! Canopy Coverage
	Hcan  = sitev(6)   ! Canopy height  
	LAI   = sitev(7)   ! Leaf Area Index  
	Sbar  = sitev(8)  ! Maximum snow load held per unit branch area(Kg/m2 for Pine)
	Ycage = sitev(9)  ! Parameters for wind speed transformation
                                 ! Ycage=1 for young pine   Should be in parameters
                                 ! Ycage=2 for Leafed deciduous
                                 ! Ycage=3 for Old pine with logn stems (Paw U and Meyers, 1987) 
			        		   ! Requires for wind speed transformation		
!  Control flags
    iradfl= iflag(1)
    pflag = iflag(2)
    ounit = iflag(3)      !iflag(4) albedo caculation
	                 
	iTsMethod = iflag(5)  ! the method to approximate the surface temperature
						  ! 1 normal old snow melt model
						  ! 2 revised direct gradient method (New method for ke) and qe
						  ! 3 force restore approach
						  ! 4 modified force restore approach
    windfl = iflag(6)

!  Model time step information
	yy = mtime(1)
	mm = mtime(2)
	dd = mtime(3)
	hr = mtime(4)

!**************************************************
!   Different densities and their ratio   
      RRHOI= RHOI/RHOW
      RRHO = RHO/RHOW
      RID  = 1.0/RRHO-1.0/RRHOI    !used to compute melt water flux (Fmelt)
	rhom = lc*rho

	! for Glacier melting calculation
    IF(SITEV(10) .EQ. 0 .OR. SITEV(10) .EQ. 3)THEN
        WGT=0.0
    ELSE
        WGT=1.0
    END IF
	Ws=Ws+WGT
!   Loop for each time step
!
      DO 2 i = 1,nt				
!   Input variables
        Ta =input(1,i)			! Air temperature input (Degrees C)
        P  =input(2,i)			! Precipitation rate input (m/hr)
        V  =input(3,i)			! Wind Speed (m/s)
        RH =input(4,i)			! Relative humidity (fraction 0-1)
!   IF(iradfl.eq.0)THEN		! input is incoming short and longwave
        Qsi=input(5,i)			! Incoming shortwave radiation (KJ/m^2/hr)
        Qli=input(6,i)			! Incoming longwave radiation (KJ/m^2/hr)
!   ELSE
        Qnetob = input(7,i) ! Net allwave radiation (KJ/m^2/hr)
!   ENDIF
        coszen = input(8,i)   ! Cos(angle between direct sunlight and surface normal)

!         Representative value over time step used in albedo calculation.  
!         We neglect the difference between direct and diffuse albedo.
!DGT Daily average temperatures handled internally so that multiple time steps will work
	  Tssk_ave = daily_ave(tsprevday, nstepday, -100.)+tk ! (C)  Surface temperature average over last 24 hours
	  Tsavek_ave = daily_ave(taveprevday, nstepday, -100.)+tk ! (C)  Depth averaged temperature average over last 24 hours
      Tssk_old = tsprevday(nstepday)+tk ! (C) Surface temperature from previous time step
      Tsavek_old = taveprevday(nstepday)+tk ! (C) Average temperature from previous time step
!   If any of these variables are out of range due to any problem set them back to freezing
	  if(Tssk_old .lt. 0.)then
	    If (snowdgtvariteflag .EQ. 1)then
	        write(66,*)"Invalid previous time step surface temperature ", &
          	Tssk_old," set to 273 K"
        end if
	        Tssk_old =tk  
	  endif
	  if(Tsavek_old .lt. 0.)then
	    If (snowdgtvariteflag .EQ. 1)then
	        write(66,*)"Invalid previous time step average temperature ", &
            Tsavek_old," set to 273 K"
        end if
	    Tsavek_old=tk
	  endif
	  if(Tssk_ave .lt. 0.)then
	    If (snowdgtvariteflag .EQ. 1)then
	        write(66,*)"Invalid last 24 hr average surface temperature ", &
          	Tssk_ave," set to 273 K"
	    end if
	    Tssk_ave=tk
	  endif
	  if(Tsavek_ave .lt. 0.)then
	    If (snowdgtvariteflag .EQ. 1)then
	        write(66,*)"Invalid last 24 hr average temperature ", &
          	Tsavek_ave," set to 273 K"
	    end if
	    Tsavek_ave=tk   
	  endif

!  Separate rain and snow      
        Ps = PARTSNOW(P,TA,TR,TS)
        Pr = P - PS
!  Increase precipitation as snow by drift multiplication factor
        PS = PS*dF

!  Calculate Albedo	
        ! for Glacier melting calculation
        
        IF(iflag(4).eq.1)THEN
            IF(SITEV(10) .EQ. 1 .OR. SITEV(10) .EQ. 2)THEN  
            !  Here we are in glacier so need to remove WGT from Ws before computing depth and guard against negative
                snowdepth=max((Ws-WGT)/RRHO,0.0)
                A=albedo(statev(3),coszen,snowdepth,aep,abg,avo,anir0)      
	        ELSE    
	            snowdepth=Ws/RRHO    ! Snow depth (Ws/RRho)            
                A=albedo(statev(3),coszen,snowdepth,aep,abg,avo,anir0)
            ENDIF
            ! Use of this albedo throughout time step neglects the
            ! changes due to new snow within a time step.
        ELSE
            A=statev(3)                                                 
        END IF   
												

!************  Maximum Interception, Literature  ************************
!  This model assumes that rainfall interception and thoroughfall occurs in a similar way that 
!  snowfall interception and unloading does. However this model is not suitable to model the
!  rainfall interception and throughfall

	Rhofs = 67.92+ 51.25* exp(Ta/2.59)		! Density of fresh snow from air temp, Army Corps 1956.
	S     = Sbar*(0.27+46/Rhofs)			! canopy holding capacity depends on density of snow (Headstrom and Pomeroy (1998))
	InMax = S*LAI							! kg/m2 =lit/m2 = .001m3/m2=.001 m =1/1000 m ()
	If (Cc .GT. 0)then
	Inmax = Inmax/1000*Cc 					! convert to m for per sq m of canopy
	Else
	Inmax = 0
	End if
	

	
!*************************************************************************
!   Call predictor corrector subroutine to do all the work for Canopy

	 CALL PREDICORRc(Us,Ws,Wc,A,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,Qli,atff, &
     	    COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime,  &
          	                                   ! Following variables are output
        Tc,QHc,QEc,Ec,Qpc,Qmc,Mc,FMc,int,Inmax,ieff,Ur,&
        Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc,&
         Rkinsc,Vz,Tac,  &                             ! Just for testing
               QHs,QEs,Es,QPs, &
               MR,QMs,Q,FM,TSURFs,tave,qnet,refDepth, totalRefDepth, &
               smelt,gsurf,Qlis)
                 
!  DGT 4/2/05   Despite all checks in predicor It can (and does) occur that 
!   we still have Us so large that it results in tave greater than 0, which implies that all the 
!   snow is liquid.  In these cases - just force the snow to disappear and add the energy involved to Qm.
 
     Tave = tavg(Us,Ws,rhow,cs,to,rhog,de,cg,hf) !  this call 
!   necessary to keep track of average internal temperature used in some surface energy algorithms.
	 IF(Tave .gt. 0.)THEN   !  all is liquid so snow must disappear
		mr = mr+Ws/dt
		qms = qms+Ws/dt*rhow*hf
		q  = q-Ws/dt*rhow*hf
		Us = Us-Ws/dt*rhow*hf
		Ws = 0.
     ENDIF
!  since Us, Ws was changed need to re-evaluate Tave
    Tave  = tavg(Us,Ws,rhow,cs,to,rhog,de,cg,hf)   

! DGT 7/25/05  To guard against unreasonable Us when there is no snow do not allow bulk temperature to go above 10 C
    if(Tave .gt. 10.)then
         Us=rhog*de*cg*10.0
         Tave = 10.0
	endif
! Update snow surface age based on snowfall in time step
    if(iflag(4).eq.1) call agesn(statev(3),dt,ps,Tsurfs,tk,dNewS)  

!  2013 - Introduced glacier substrate into model
!  Partition melt outflow into the portion that is from rain, snow and glacier
!  Inputs are condensation (C=max(-Es,0)), precipitation as rain (pr) and precipitation as snow (Ps)
!  Outputs are melt (Mr) and sublimation (Sub=Max(Es,0))
!  DW is change in snow water equivalent (positive if decrease)
!  DG is change in glacier (positive if decrease)
!  Mass balance
!  DW+DG=Mr+Sub-Pr-Ps-C
!  In the case of snow left on the ground/glacier - any rain is assumed to have added to the snow
    SnowLeft = MAX(Ws-WGT,0.0)
    DW = (statev(2) - SnowLeft)/dt   ! change in seasonal snow rate
    ! Reduction in snow storage is snow at the beginning minus the snow at the end in excess of glacier.  Can be negative if snow increases
    DG = Max(WGT-Ws,0.0)/dt   !  Reduction in glacier storage rate
 
    if(SnowLeft > 0.)then
      R=0. ! No outflow due to rain
!  Inequalities that hold
!  0<= Ms <= Ps+Pr+C+DW
!  0<= MG <= DG
!  The slack in these inequalities is due to Sub
!  Compute Ms and MG using proportioning
      MSX = Ps+Pr+DW+Max(-Es,0.0)   !  Maximum melt contrib from snow.  Note that condensation is taken to add to snow as is rain
      MRX = MSX+DG     !  Maximum melt overall
    ! Mr is less than MRX due to sublimation losses.  Do the rain proportioning accordingly  
      IF(MRX <= 0.0)THEN
        Ms=0.
        MG=0.
      ELSE
        Ms=Mr*MSX/MRX
        MG=Mr*DG/MRX
      END If
    else  ! The case where all seasonal snow is gone and precipitation may contribute to total outflow
!  Inequalities that hold
!  0<= Ms <= Ps+C+DW
!  0<= R <= Pr
!  0<= MG <= DG
!  The slack in these inequalities is due to Sub.  Note that here Pr may contribute to R but that C is still put to snow
      MSX = Ps+DW+Max(-Es,0.0)   !  Maximum melt contrib from snow.  Note that condensation is taken to add to snow as is rain
      MRX = MSX+DG+Pr     !  Maximum melt overall 
      IF(MRX <= 0.0)THEN
        Ms=0.
        R=0.
        MG=0.
      ELSE
        Ms=Mr*MSX/MRX
        R=Mr*Pr/MRX
        MG=Mr*DG/MRX
      END If
    endif     
    Eglacier=DG-MG
    SWISM=Ms
    SWIGM=MG
    SWIR=R
    SWIT=Mr
   IF (Ws .LT. WGT)Then                  !  There was glacier melt
       Ws=0                              !  Reset      
   ELSE                                  !  Here ws is greater than the glacier thickness so no glacier melt
       Ws=Ws-WGT  
   END IF          
    
!  accumulate for mass balance
   cump  = cump+(Ps+Pr)*dt
   cumes  = cumes+Es*dt                 
   cumEc = cumEc+Ec*dt                ! Evaporation from canopy
   cumMr = cumMr+Mr*dt                ! canopy melt not added
   cumGM = cumGM+SWIGM*dt             ! Cumulative glacier melt
   CumEG= CumEG+Eglacier*dt
!  yjs update the total depth of the refreezing depth in the snow pack according the 
!  yjs the refreezing depth at time step and the positive energy input. 07/22/01
!  DGT's revised logic  1/13/05
      if(lc.gt.0.0) then
       if(refDepth.gt. 0.0) then
	    totalRefDepth=refDepth  ! if refreezing depth increases totalrefdepth keeps track
	 else  ! here refdepth has gone to 0 somehow
	  if(mr.gt.0.0 .or. (Us.gt.Us_old .and. Us .gt.0.0)) &
!   If there is melt or an increase in energy refdepth is reset 
     	    totalRefDepth = 0.0 
	 endif
	elseif(mr.gt.0.0 .or. (Us.gt.Us_old .and. Us .gt.0.0)) then
!   Here lc=0.  If there is melt or an increase in energy refdepth is reset
!   This is likely redundant because if lc=0 then there is no meltwater to refreeze
	 totalRefDepth =0.0	 
	endif

	if(totalRefDepth .lt. 0.0)  totalRefDepth=0.0
!
!yjs update tsbackup and tavebackup
	 do 50 ii = 1 , nstepday-1
		tsprevday(ii)= tsprevday(ii+1)
	    taveprevday(ii)= taveprevday(ii+1)
 50	 continue
	   tsprevday(nstepday)= Tsurfs
	   taveprevday(nstepday)= tave

       IF(PFLAG.eq.1) THEN 
            write(ounit,*) US,Ws,statev(3),&
                     Pr,ps,A,qhs,qes,es,mr,qms,q,fm,tave,Tsurfs&
                    ,cump,cumes,cummr,qnet,smelt,refdepth,totalrefdepth&
       ,Cf,Taufb,Taufd,Qsib,Qsid &                  ! Atm Trans.,Solar,Dir. & DIff. (Atm Transm & Radiation)
       ,Taub,Taud &                                 ! Solar Dir, and Diff.( canopy transmittance)
       ,Qsns,Qsnc &                                 ! Solar Radiation Canopy & Sub Canopy    
       ,Qlns,Qlnc &                                 ! Longwave Radiation Canopy & Sub Canopy 
       ,Vz &                                        ! Wind speed Sub-Canopy
       ,Rkinsc, Rkinc &                             ! heat and water vap Cond, Sub-Canopy, Canopy
       ,Inmax,int,ieff ,Ur &                        ! Interception 
       ,Wc,Tc,Tac,QHc,QEc,Ec,Qpc,Qmc,Mc,FMc,&       ! Fluxes and Energies
       SWIGM,SWISM,SWIR
       ENDIF
        
      OutArr(1)=US
      OutArr(2)=Ws
      OutArr(3)=statev(3)
      OutArr(4)=Pr
      OutArr(5)=ps
      OutArr(6)=A
      OutArr(7)=qhs
      OutArr(8)=qes
      OutArr(9)=es
      OutArr(10)=SWIT
      OutArr(11)=qms
      OutArr(12)=q
      OutArr(13)=fm
      OutArr(14)=tave
      OutArr(15)=Tsurfs
      OutArr(16)=cump
      OutArr(17)=cumes
      OutArr(18)=cummr
      OutArr(19)=qnet
      OutArr(20)=smelt
      OutArr(21)=refdepth
      OutArr(22)=totalrefdepth
      OutArr(23)=Cf
      OutArr(24)=Taufb
      OutArr(25)=Taufd
      OutArr(26)=Qsib
      OutArr(27)=Qsid
      OutArr(28)=Taub
      OutArr(29)=Taud
      OutArr(30)=Qsns
      OutArr(31)=Qsnc
      OutArr(32)=Qlns
      OutArr(33)=Qlnc  
      OutArr(34)=Vz
      OutArr(35)=Rkinsc
      OutArr(36)=Rkinc 
      OutArr(37)=Inmax
      OutArr(38)=int
      OutArr(39)=ieff
      OutArr(40)=Ur
      OutArr(41)=Wc
      OutArr(42)=Tc
      OutArr(43)=Tac
      OutArr(44)=QHc
      OutArr(45)=QEc
      OutArr(46)=Ec
      OutArr(47)=Qpc
      OutArr(48)=Qmc
      OutArr(49)=Mc
      OutArr(50)=FMc
      OutArr(51)=SWIGM
      OutArr(52)=SWISM
      OutArr(53)=SWIR
       
  2    continue

       statev(1)=Us
       statev(2)=Ws
	   statev(5)=refDepth
	   statev(6)=totalRefDepth
	   statev(4)= Wc                      !  Added vinod
       outv(1)=Pr
       outv(2)=ps
       outv(3)=a
       outv(4)=qhs
       outv(5)=qes
       outv(6)=es
       outv(7)=mr
       outv(8)=qms
       outv(9)=q
       outv(10)=fm
       outv(11)=tave
       outv(12)=Tsurfs
       outv(13)=qnet
	   outv(14)=smelt   ! dgt 5/4/04

       RETURN
       END

!************************  CORRECTOR ***********************************

      SUBROUTINE PREDICORRc(Us,Ws,Wc,A,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,Qli &
      ,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime,  & 
          	                                   ! Following variables are output
        Tc,QHc,QEc,Ec,Qpc,Qmc,Mc,FMc,int,Inmax,ieff,Ur, &
        Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, &
         Rkinsc,Vz,Tac, &                             ! Just for testing
               QHs,QEs,Es,QPs, &
               MR,QMs,Q,FM,TSURFs,tave,qnet,refDepth, totalRefDepth, &
               smelt,gsurf, Qlis)

	REAL  k,Ta,Mc,LAI,int,int1,Inmax,ieff,i1,Mc1

	REAL MR,mr1,Lans,LanG
      real sitev(*)
      real param(*)
	real mtime(*)        !yjs add model time 
	Integer iradfl
    !real:: WGT
    
!	common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

!yjs  Constant data set 
      data to /0.0/        !  Temperature of freezing (0 C)
      data tk /273.15/     !  Temperature to convert C to K (273.15)
      data sbc /2.0747e-7/ !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data hf /333.5/      !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cw /4.18/       !  Water Heat Capacity (4.18 KJ/kg/C)
      data cs /2.09/       !  Ice heat capacity (2.09 KJ/kg/C)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data Rag /287.0/     !  Ideal Gas constant for dry air (287 J/kg/K)
      data k  /0.4/        !  Von Karmans constant (0.4)
      data hff /3600.0/    !  Factor to convert /s into /hr (3600)
      data rhoi /917.0/    !  Density of Ice (917 kg/m^3)
      data rhow /1000.0/   !  Density of Water (1000 kg/m^3)
      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)
      data w1day /0.261799/!  Daily frequency (2pi/24 hr 0.261799 radians/hr) 
      data pi /3.141592654/!  Pi
!yjs  End of constant declaration
 	Apr  = sitev(2)
	Cc   = sitev(5)
      LAI  = sitev(7)	

	data wtol,utol / 0.025,2000./
      data ncall / 0/
      ncall=ncall+1

      CALL  QFMc(Us,Ws,Wc,A,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, &
     	    COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &  
          	                                   ! Following variables are output
        Tc,QHc,QEc,Ec,Qpc,Qps,Qmc,Mc,FMc,int,Inmax,ieff,Ur, &
        Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, &
         Rkinsc,Vz,TSURFs,Tac,  &                              ! Just for testing
      tave,qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, &
      Smelt,smeltc)


!     PREDICTOR FOR CANOPY SNOW
      Wc1 = Wc + DT*FMc

		IF(Wc1.LT.0.0) THEN
			Wc1 = 0.0
			FMc = (Wc1- Wc)/DT
			Ur  = max(0.,(int-Fmc-Ec-Mc))
			Ec  = max(0.,(int-Fmc-Mc-Ur))
			Mc  = (int-Fmc-Ec-Ur)
	        FM  = PR+PS-int+Ur+Mc-MR-Es         ! int,ur,Mc added here
		ENDIF
!   Save values so that they can be averaged for output
      FMc1 = FMc 
	int1 = int
	Ur1  = Ur
	Mc1  = Mc
	Ec1  = Ec
      Tc1  = Tc

!     PREDICTOR FOR SUB-CANOPY SNOW
       Ws1 = Ws + DT*FM
       
	 IF(Ws1.LT.0.0) THEN
         Ws1=0.0
         CALL PREHELP(Ws1,Ws,DT,FM,0.,1.,PS,PR,Es,RHOW,HF,Q,QMs,MR,qes, &
              hneu)
       ENDIF

       Us1 = Us + DT*Q
       Q1  = Q
       FM1 = FM 
!   Save values so that they can be averaged for output
       qhs1	=qhs
       qes1	=qes
       es1	=es
       mr1	=mr
       smelt1	=smelt                         !cdgt 5/4/04 surface melt smelt
       qms1	=qms
       Tsurfs1=Tsurfs
       qnet1	=qnet

      CALL  QFMc(Us1,Ws1,Wc1,A,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, &
     	    COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime,   &
          	                                   ! Following variables are output
        Tc,QHc,QEc,Ec,Qpc,Qps,Qmc,Mc,FMc,int,Inmax,ieff,Ur, &
        Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, &
         Rkinsc,Vz,TSURFs,Tac,  &                              ! Just for testing
      tave,qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, &
      Smelt,smeltc)


!     CORRECTOR FOR CANOPY SNOW

      Wc2 = Wc + DT/2.0*(FMc1 + FMc)


		IF(Wc2.LT.0.0) THEN
		    Wc2 = 0.0
			FMc = (Wc2- Wc)/DT*2-FMc1
			Ur  = max(0.,(int-Fmc-Ec-Mc))
			Ec  = max(0.,(int-Fmc-Mc-Ur))
			Mc  = (int-Fmc-Ec-Ur)
              FM  = PR+PS-int+Ur+Mc-MR-Es         ! int,ur,Mc added here  
		ENDIF
 	int = (int1 +int)/2
	Ur  = (Ur1+Ur)/2
	Mc  = (Mc1+Mc)/2
	Ec  = (Ec1+Ec)/2                            ! But Mc and Ec don't change
	Wc  =  Wc2
	Tc  = (Tc1+ Tc)/2
	FMc = (FMc1+FMc)/2


!     CORRECTOR FOR SUB_CANOPY SNOW

       Ws2 = Ws + DT/2.0*(FM1 + FM)
       IF(Ws2.LT.0.0) THEN
         Ws2=0.0
         CALL PREHELP(Ws2,Ws,DT,FM,FM1,2.,PS,Pr,Es,RHOW,HF,Q,QMs,MR,qes,&
        hneu)
       ENDIF
       Us2 = Us + DT/2.0*(Q1 + Q)   


!   iterate to convergence to enhance stability
       niter=-3
       imax=1
 1     if((abs(Ws2-Ws1).gt.wtol .or. abs(Us2-Us1).gt.utol) .and. &
            niter .lt. imax)then
!      write(6,10)ncall,"_ncall"
! 10   format(1x,i5,a6)
          Ws1 = Ws2
          Us1 = Us2
          WC1 = Wc2

      CALL  QFMc(Us1,Ws1,Wc1,A,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, &
     	    COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime,  & 
          	                                   ! Following variables are output
        Tc,QHc,QEc,Ec,Qpc,Qps,Qmc,Mc,FMc,int,Inmax,ieff,Ur, &
        Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, &
         Rkinsc,Vz,TSURFs,Tac,  &                              ! Just for testing
      tave,qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, &
      Smelt,smeltc)


!     CORRECTOR FOR CANOPY SNOW AGAIN

      Wc2 = Wc + DT/2.0*(FMc1 + FMc)

		IF(Wc2.LT.0.0) THEN
		    Wc2 = 0.0
			FMc = (Wc2- Wc)/DT*2-FMc1
			Ur  = max(0.,(int-Fmc-Ec-Mc))
			Ec  = max(0.,(int-Fmc-Mc-Ur))
			Mc  = (int-Fmc-Ec-Ur)
              FM  = PR+PS-int+Ur+Mc-MR-Es         ! int,ur,Mc added here
!			FM  = FM+ (FMc-Ec)  
		ENDIF
 	int = (int1 +int)/2
	Ur  = (Ur1+Ur)/2
	Mc  = (Mc1+Mc)/2
	Ec  = (Ec1+Ec)/2                         
	Wc  =  Wc2
	Tc  = (Tc1+ Tc)/2
	FMc = (FMc1+FMc)/2

!dgt 5/4/04 surface melt smelt
!    corrector again

         Ws2 = Ws + DT/2.0*(FM1 + FM)
         IF(Ws2.LT.0.0) THEN
           Ws2=0.0
           CALL PREHELP(Ws2,Ws,DT,FM,FM1,2.,PS,Pr,Es,RHOW,HF,Q,QMs,MR,&
                qes,hneu)
         ENDIF
         Us2 = Us + DT/2.0*(Q1 + Q)
         niter=niter+1
         if(niter .ge. 1)then  ! had * steps to converge now hit it.
!  What follows is a fix to numerical instability that results from
!  nonlinearity when the snowpack is shallow and melting rapidly.  If
!  convergence does not occur when the snowpack is not melting (a very
!  rare thing) I just accept the predictor corrector solution.
!
!  Modified by DGT 7/27/05 to add complete meltout condition
!  The quantities that this changes are w2, ub2, mr and qms
!
!  The fix first checks whether the added energy is sufficient to melt the snow
!  completely.  If this is the case then the snow disappears.
!  In other cases we assume that the liquid fraction of the snow remains constant.
!  This to some extent bypasses the melt outflow estimates. 
!  ae is added energy during the time step.
          ae=(q1+q+qms1+qms)*.5*dt
!   This fix is only physically sensible under melting conditions
!   and when ae is positive and there is snow
          if(Us .gt. 0. .and. ae .gt. 0. .and. Ws .gt. 0.)then
	        es2=(es+es1)*.5   !  This is the average sublimation
!   Check liquid fraction with added energy.  If this is greater than 1 then all snow melts
!   Otherwise implement a solution assuming that the liquid fraction remains constant
	        rlf=(Us+ae)/(rhow*Ws*hf)
			if(rlf .ge. 1.)then
				mr=Ws/dt+(Ps+Pr-es2)   ! Here snow disappears
				if(mr .lt. 0.)then
					mr=0.   !  Force this to not go negative
!      This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
!      There is a resulting energy discrepancy due to a limitation on sublimation and latent heat flux
!      This is ignored because once the snow is gone energy computations are not pertinent.
                  endif
				qms=mr*rhow*hf
				Ws2=0.
				Us2=Us+ae-qms*dt
			else
!   Determine the w/ub ratio at the beginning of the time step.  
!   Physically liquid fraction = ub/(rhow*w*hf) and since rhow and hf are constants
!   keeping r=w/ub constant is equivalent to keeping liquid fraction constant.
!   If liquid fraction is constant this is constant.
				r=Ws/Us
!   Solve for ub2 and w2 that satisfy the three equations
!            r=w2/ub2 
!            ub2=ub+ae-rhow*hf*mr*dt     Energy balance the last term being the energy advected by melt
!            w2=w+(ps+prain-e2-mr)*dt    Mass balance 
!   The unknowns in the above are ub2, w2 and m and the equations are linear 
!   once the first eqn is multiplied by ub2
!   The solution is          
				Us2=(rhow*hf*(Ws+(Ps+Pr-es2)*dt)-ae-Us)/(rhow*hf*r-1)
				Ws2=r*Us2
				if(Ws2 .lt. 0.)then  ! Avoid negative W 
					Ws2=0.
				endif
				mr=(Ws-Ws2)/dt-es2+ps+Pr
				if(mr .lt. 0.)then   ! Avoid negative mr
					mr=0.
					Ws2=Ws+(ps+Pr-es2)/dt
					if(Ws2 .lt. 0)then
						Ws2=0.
!      This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
!      There is a resulting energy discrepancy due to a limitation on sublimation and latent heat flux
!      This is ignored because once the snow is gone energy computations are not pertinent.
					endif
				endif
				qms=mr*rhow*hf
				Us2=Us+ae-qms*dt   ! redundant most of the time but recalc due to exceptions
			endif
!    Check that nothing went wrong
              if(mr .lt. 0.)then
	            write(6,*)'Error - negative melt rate in snow'
			endif
	        if(Ws2 .lt. 0.)then
	            write(6,*)'Error - negative w2 in snow'
			endif
			q=ae/dt-qms
!   Now set first pass values equal to final values to fake the averages below
			qms1=qms
			mr1=mr
			q1=q
          endif
         endif
         go to 1
       endif
       Ws = Ws2
       Us=Us2
!  average values from two time steps for output.  This is done for mr
!  and e to ensure mass balance and the others for better physical
!  comparisons
       qh=(qhs+qhs1)*.5
       qe=(qes+qes1)*.5
       es=(es+es1)*.5
       mr=(mr+mr1)*.5
       qms=(qms+qms1)*.5
       Tsurfs=(Tsurfs+Tsurfs1)*.5
       qnet=(qnet+qnet1)*.5
       q=(q+q1)*.5
	 smelt=(smelt+smelt1)*.5/(hf*rhow)  !cdgt 5/4/04 surface melt smelt
!  convert from energy KJ/m^2/hr to water depth m/hr of melt.
       RETURN
       END 

!*************************************************************************
!				CALCULATE CANOPY MASS FLUX AT ANY INSTANT


 	SUBROUTINE QFMc(Us,Ws,Wc,A,dt,rid,P,Pr,Ps,Ta,V,RH,Qsi,atff,Qli, &
     	    COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime,  & 
          	                                   ! Following variables are output
        Tc,QHc,QEc,Ec,Qpc,Qps,Qmc,Mc,FMc,int,Inmax,ieff,Ur, &
        Cf,Taufb,Taufd,Qsib,Qsid,Taub,Taud,Qsnc,Qsns,Qlnc,Qlns,Rkinc, &
         Rkinsc,Vz,TSURFs,Tac,  &                              ! Just for testing
     
      tave,qnet,QHs,QEs,Es,MR,QMs,Q,FM,refDepth,totalRefDepth, &
      Smelt,smeltc)


	REAL  Mc,LAI,int,Inmax,ieff,k
	REAL param(*)
      real sitev(*)
	real mtime(*)                   !yjs 09/19/2000 add the model time
      real mr,lc,ks,Lans,LanG,refdepth,refdep, smeltc
!	DOUBLE PRECISION Ur

!yjs  Constant data set 
      data to /0.0/        !  Temperature of freezing (0 C)
      data tk /273.15/     !  Temperature to convert C to K (273.15)
      data sbc /2.0747e-7/ !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data hf /333.5/      !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cw /4.18/       !  Water Heat Capacity (4.18 KJ/kg/C)
      data cs /2.09/       !  Ice heat capacity (2.09 KJ/kg/C)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data Rag /287.0/      !  Ideal Gas constant for dry air (287 J/kg/K)
      data k  /0.4/        !  Von Karmans constant (0.4)
      data hff /3600.0/    !  Factor to convert /s into /hr (3600)
      data rhoi /917.0/    !  Density of Ice (917 kg/m^3)
      data rhow /1000.0/   !  Density of Water (1000 kg/m^3)
      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)
	data w1day /0.261799/!  Daily frequency (2pi/24 hr 0.261799 radians/hr) 
	data pi /3.141592654/!  Pi

	common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

!yjs  End of constant declaration

!  Parameters
      Ems=param(3)    !  emmissivity of snow (nominally 0.99)
      cg =param(4)    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      z=param(5)      !  Nominal meas. height for air temp. and humidity (2m)
      zo=param(6)     !  Surface aerodynamic roughness (m)
      rho=param(7)    !  Snow Density (Nominally 450 kg/m^3)
      rhog=param(8)   !  Soil Density (nominally 1700 kg/m^3)
      lc=param(9)     !  Liquid holding capacity of snow (0.05)
      ks=param(10)    !  Snow Saturated hydraulic conductivity (160 m/hr)
      de=param(11)    !  Thermally active depth of soil (0.4 m)
      abg=param(12)   !  Bare ground albedo  (0.25)
      avo=param(13)   !  Visual new snow albedo (0.95)
      anir0=param(14) ! NIR new snow albedo (0.65)
	lans= param(15) ! the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	lang= param(16) ! the thermal conductivity of soil (9.68 kJ/m/k/hr)
	wlf= param(17)  ! Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr) 
	rd1= param(18)  ! Apmlitude correction coefficient of heat conduction (1)
      fstab=param(19) ! Stability correction control parameter 0 = no corrections, 1 = full corrections
	Tref=param(20)  !  Reference temperature of soil layer in ground heat calculation input
	dNewS=param(21) !  The threshold depth of for new snow (0.001 m)
	gsurf=param(22) !  The fraction of surface melt that runs off (e.g. from a glacier)
                      !  cdgt gsurf added for modeling surface runoff from a glacier

!   Site variables
      df    =  sitev(1)      !  Drift factor
      Apr   =  sitev(2)      !  Atmospheric Pressure (Pa)
      qg    =  sitev(3)      !  Ground heat flux (KJ/m^2/hr)  This is more logically an
	Cc    =  sitev(5)      !  Canopy Coverage
	LAI   =  sitev(7)	   !  Leaf Area Index


!  To ensure all temperatures in kelvin
      TAK    = TA+TK
      TAVEK  = TAVE+TK
      RHOA   = APr/(RAg*TAK)     ! Density of Air in kg/m3 if APr in Pa
	rhom   = lc*rho

!   Some computations for conductivity 			   
	fKappaS = lans/(rho*cs)                           ! lans= lamda= thermao conductivity
	ds		= sqrt(2*fKappaS/w1day)                   ! fkappas =k= diffusivity
	TherC	= Lans/(ds*rd1*rho*cs)					  ! In QFM, related to thermal conductivity , lambda, of snow (Rs in Ori)
	Zs      = Ws*rhow/rho							  ! snow depth

!yjs save the old values
	tave_old=tave	

       EA=SVPW(TA)*RH              ! The saturation vapour pressure over water is used 
!                                   because most instruments report relative humidity relative to water.

!  Fraction of snow on canopy [Dickinson et. al (1993) eqn (50a)]
	IF (Inmax .EQ.0.)THEN             !(When LAI is zero for open area model)
	 Fs= 0.0
	ELSE

	 Fs   = (Wc/Inmax)**(2./3)  
	 IF (Fs.GT.1) Fs=1
!	 Fs   = (Wc/Inmax)**(.05)  
      ENDIF 

	 CALL INTERCEPT(Ta,LAI,p,Wc,dt,Inmax,param,sitev, &
                ieff,Ur,int)          ! Output variables
	                
!  Total heat advected by precipitation is 
	 Qp   = QPF(PR,TA,TO,PS,RHOW,HF,CW,CS)
!  Heat advected by precipitaion on canopy snow and sub-canopy snow is 
       IF (p.GT.0) THEN
	    IF(Wc.GT.0.)THEN
          Qpc  = int/p*Qp          ! Energy from precipitaiton that hits the canopy snow
          Qps  = Qp-Qpc            ! Energy advected by precipitaion that passes through canopy gap and thoroughfall from canopy 
	    ELSE
	    Qpc = 0.
          Qps = ((p-int)/p)*Qp 
	    ENDIF
	 ELSE 
	    Qpc=0.
      	Qps=0.
	 ENDIF  
                  
       Tave = tavg(Us,Ws,rhow,cs,to,rhog,de,cg,hf)

!yjs as described as below, the original UEB does not consider the refreezing. in this
! change, the model will consider the refreezing effect starting from the top to the bottom.
! Here is the predictor.

!yjs the method is: refreezing can be modeled if the Us is greater than 0.0 without any precipiation, 
! meanwhile the total refreezing depth in the snowpack is less than the Refreezing depth times the daily damping of wet snow.
 
!       IF (Wc.EQ.0.and.p.EQ.0.) THEN   ! When there is no snow in the canopy
!		 Tc  =  Ta                   
!  	 ELSE
!		Tc  = 0.0            ! First assume snow is melting and then if result is negative - calculate canopy temperature having determined snow is not melting
!	 ENDIF


!**********   To Calculate Refr Depth	    
 	if(Us .gt.0 .and. Ps.le.0.0 .and. Pr.le.0.0 .and. &
       totalRefDepth .le. rd1*ds .and. Ws.gt.0.0) then


	Qc1=  QcEst(Ws,p,Tk,Tck,V,Zm,d,Z0c,Rimax,Rcastar,Cf,Fs, &     
     		Qli,Hcan,Vz,Ta,Rh,RKINsc,Qps,To,Ps,Qsi,atff,COSZEN, &
          APr,TAK, EA,A,Ac,Wc,Inmax, Qnetob,Iradfl,param,sitev )

	Qc2=  QcEst(Ws,p,Tk-.01,Tck,V,Zm,d,Z0c,Rimax,Rcastar,Cf,Fs, &
     		Qli,Hcan,Vz,Ta,Rh,RKINsc,Qps,To,Ps,Qsi,atff,COSZEN, &
           APr,TAK, EA,A,Ac,Wc,Inmax, Qnetob,Iradfl,param,sitev )

	  call Grad(qc1,qc2,0.0,-0.01, var_a, var_b)
          x1=refDep(lans,var_a, var_b, hf, rhom, dt, refDepth)       !refreezing depth
	  refDepth=x1

	else
        refDepth=0.0
      endif

!**********************

! call Temperature

      Tsurfs= SRFTMPsc(Tssk, Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Qsi,atff,Cf, &
     	Qli,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
       Qpc,Qps, Inmax, Rkinc,Rkinsc,Vz,Tc,Tk,Tak,EA,RHOA, &
         fkappaS,RHO,TherC,Fs,  &                                          
         tave,refDepth,Smelt, SmeltC)    ! This will give me Tsurf, Tc, and Smeltc
 

	IF(IRADFL.EQ.1) Go to 13   !  To avoid these steps if the net radiation is input
		CALL PSOLRAD(Qsi,atff,param,cf  &                ! Input Variables 
     	                    ,Taufb,Taufd,Qsib,Qsid)    ! Output variables:




		CALL TRANSRADCAN (COSZEN,sitev,param, &               ! Input variables , leaf parameters
     	         Betab,Betad,Taub,Taud)   ! Output variables:


	 CALL NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, &
                  Inmax,Qsib,Qsid,param,Fs, &
                  Qsns,Qsnc ) !  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 

	  CALL NETLONGRAD(RH,Ta,Tsurfs,Tc,Tk,Fs,EmC,EmS,SBC,cf,sitev, &
                  	Qli,param, &
                        Qlis,Qlns,Qlnc )                     !  Output: Qsns,Qsnc 


13	Ess    = SVPI(Tsurfs)
	Esc    = SVPI(Tc) 

	 CALL TURBFLUX(Ws,Wc,A,TK,Tc,Ta,Tsurfs,RH,V,EA,p,param,sitev, &
     	              d,Z0c,Vz,RKINc,RKINsc,Tac,Fs,Ess,Esc, &                     ! Output variables
                          QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E)          

    
  	CALL INTERCEPT(Ta,LAI,p,Wc,dt,Inmax,param,sitev, &
     	                            ieff,Ur,int)    ! Output variables

	
	! melt from the canopy			                                    	
       Qmc  =   smeltC
	 Mc   =   Qmc/(Rhow*HF)
 
!     Ridistribution of intercepted snow subtracting in mass balance
!     We are applying same amount as unloading is lost from the canopy from the wind redistribution
 !      Rid = Ur

!*****************  UPDATING THE MASS BALANCE AT CANOPY

12	FMc = int-Ec-Mc-Ur

!*****************  UPDATING THE MASS and ENERGY BALANCE AT SUB-CANOPY


	 MR=FMELT(Us,RHOW,Ws,HF,LC,RID,KS,Pr)   !yjs Add a fraction to reduce the evaporation after snow gone
       
!                                                 MR in m/hr
       QMs=MR*RHOW*(HF+(Tave-to)*Cw)           !  Include advection of
!                                               meltwater/rain that is at tave so that the model does not
!                                               crash when there is no snow and it rains.
!                                               QMs in kj/m2/hr 
!  dgt 5/4/04  Add surface melt that runs off to melt runoff and melt energy so that energy and mass balances
!   are still correct
       Mr=mr+smelt/(hf*rhow)*gsurf  
	 Qms=qms+smelt*gsurf

       IF(IRADFL.EQ.0) THEN
!          QNET = QSI*(1.0-A)+QLNET
          QNET  = Qsns+Qlns                     ! At sub canopy snow
       ELSE
          QNET = QNETOB
       ENDIF

!      IF (Ws.EQ.0.and. p.EQ.0.) THEN   ! vinod added to orgi UEB
!		 Es  =  0.
! 		 Mr =  0.	
!	ENDIF

       Q = QNET+QPs+QG+QHs+QEs-QMs       
       FM=1.0*(PR+PS)-int+Ur+Mc-MR-Es                                 ! int,ur,Mc added here

	RETURN
	END

!*********************************************************************


!************************* SRFTMP () *********************************
!                COMPUTE THE SURFACE TEMPERATURE OF SNOW
 
      FUNCTION SRFTMPsc(Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Qsi,atff,Cf, &
      Qli,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
       Qpcin,Qpsin, Inmax, Rkinc,Rkinsc,Vz,Tc,Tk,Tak,EA,RHOA, &
         fkappaS,RHO,TherC,Fs,tave,refDepth,Smelt,SmeltC) 

  	
   						        	          
	real k,Inmax, LAI
	real F1, F2, F1ts,F1tc,F2ts,F2tc,J11, J12, J21, J22, delTs, delTc
	REAL param(*)
	REAL sitev(*)
	real mtime(*)
	integer iTsMethod

!     common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)   
      data Rag /287.0/     !  Ideal Gas constant for dry air (287 J/kg/K)
      data k  /0.4/        !  Von Karmans constant (0.4)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
!  Parameters
      Ems  = param(3)      !  emmissivity of snow (nominally 0.99)
      cg   = param(4)	   !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      z    = param(5)	   !  Nominal meas. height for air temp. and humidity (2m)
      zo   = param(6)	   !  Surface aerodynamic roughness (m)
      rho  = param(7)	   !  Snow Density (Nominally 450 kg/m^3)
      rhog = param(8)	   !  Soil Density (nominally 1700 kg/m^3)
      lans = param(15)	   !  thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
      lang = param(16)	   !  the thermal conductivity of soil (9.68 kJ/m/k/hr)
      wlf  = param(17)	   !  Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr) 
      rd1  = param(18)	   !  Apmlitude correction coefficient of heat conduction (1)
      Apr  = sitev(2)      !  Atmospheric Pressure (Pa)
     
      Cc    =   sitev(5)   !  Canopy Coverage
      LAI   = sitev(7)     ! Leaf Area Index  

!   This version written on 4/23/97 by Charlie Luce solves directly the
!   energy balance equation using Newtons method - avoiding the linearizations
!   used previously.  The derivative is evaluated numerically over the range
!   ts to fff*ts  where fff = 0.999

      data tol,nitermax/0.005,20/   ! changed 10 to 20 and .05 to .0005
	fff    = (273-tol)/273.                 ! Scale the difference used for derivatives by the tolerance
      TAK    = TA+TK
      TAVEK  = TAVE+TK

      Qps  =  Qpsin                          ! store input variable locally without changing global value
      Qpc  =  Qpcin
       if(Ws.le.0 .and. Qps.gt.0.0)Qps=0.0


!******** FOR OPEN AREA MODLING ********************************
		IF ((LAI .EQ. 0.0) .OR. (cc .EQ. 0.0)) THEN
		!IF (LAI .EQ. 0.0) THEN
		 Tssk=Tak                 
		 go to 18
 		ENDIF

!************  SOLVING NON LINEAR EQUATION IN TWO DIMENSTION (newton's method using Jacobian matrix)*********
   
       Tssk1   = TaK                             ! first approximation
       Tck1    = TaK                             ! first approximation  

       niter = 0
 15     IF(niter.LT. nitermax)THEN 

         F1 = SURFEBsc(Tssk1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
       Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck1,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC, &                                        
       TSURFs,tave,refDepth)   
                                                        	!yjs add three value to reflect model control changes)) 


        F2 = SURFEBc(Tck1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk1,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC, &                                           
       TSURFs,tave,refDepth)   !yjs add three value to reflect model control changes))


!     assumed small increament to estimate the Jocobian matrix

      delTs = 0.01
	delTc = 0.01
                   
        F1ts = SURFEBsc(Tssk1+delTs,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf, &
      Qli,Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob, &
      iTsMethod,mtime,Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck1, &
      Tk,Tak,EA,RHOA, FkappaS,RHO,TherC, &                        
       TSURFs,tave,refDepth)


         F1tc = SURFEBsc(Tssk1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck1+delTc,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC, &                                    
       TSURFs,tave,refDepth)      
                                                        					 
	
	    F2ts = SURFEBc(Tck1,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk1+delTs,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,   &                                         
       TSURFs,tave,refDepth)
     	
				 
	F2tc = SURFEBc(Tck1+delTc,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk1,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)	
				                                      	 
!      Jacobian matrix
       J11 = (F1ts-F1)/delTs 
	 J12 = (F1tc-F1)/delTc
	 J21 = (F2ts-F2)/delTs
	 J22 = (F2tc-F2)/delTc  

       s1 = (F1*J22-F2*J12)/(J12*J21-J22*J11)
       s2 = (F2*J11-F1*J21)/(J12*J21-J22*J11)

      IF((J12*J21).EQ.0.0 .and. (J22*J11).EQ.0.0) Go To 14     ! It gives NAN or infinity value for s1 and s2,  occues when Ws is zero

      IF(S1.GT.5.) S1= 5.
	IF(S1.LT.-5.) S1= -5.
	IF(S2.GT.5.) S2= 5.
	IF(S2.LT.-5.) S2= -5.

	Tssk1 = Tssk1+s1
	Tck1  = Tck1+s2

      
		IF (abs(S1).GT.0.001 .or. abs(S2).GT.0.001)THEN	
			niter =niter+1
			Go TO 15       ! loop back and iterate

! When there is snow and temperature is positive,we need to iterate for again putting temp zero
 
	    ELSEIF (Ws .GT. 0 .and. Tssk1 .GT. Tk) THEN   ! Iteration doesnot fail, 
    	             Tssk1 = Tk                                ! Just for first assumption, doesnot estimate melt
	             Tssk  = Tssk1
	             Tck   = Tck1 
				 go to 17


	    ELSEIF (Wc.GT.0 .and. Tck1.GT.Tk) THEN
	             Tck1 = Tk                                  ! Just for first assumption, doesnot estimate melt
	             Tck  = Tck1
                   Tc   = Tck-Tk
	             Tssk  = Tssk1 
				 go to 17
          ELSE
			SRFTMPsc = Tssk1-Tk 
			smelt    = 0.0  
			Tc       = Tck1-Tk 
			smeltC   = 0.0
			go to 21

	    ENDIF

       
	 ELSE
!     We consider the iteration is not complete since the values of S1 and S2 are not satisfied

	  GO TO 14  ! Use another method :one dimension approach

      ENDIF

!************************  Doubt the iteration When the temperature difference between snow surface and air is more thar 20 C  *******************
	  	 IF (Tssk1.LT.(Tak-20.))THEN
		  GO TO 14  ! Use another method  :one dimension approach
            ENDIF

          IF (Ws.LE.0.and. Tssk1.GT.(Tak+20.))THEN        
		  GO TO 14  ! Use another method  :one dimension approach       
		ENDIF

!***************************************
	  	 IF (Tck1.LT.(Tak-10.))THEN
		  GO TO 14  ! Use another method  :one dimension approach
            ENDIF
        
	    IF (Wc.LE.0.and. Tck1.GT.(Tak+10.))THEN        
		  GO TO 14  ! Use another method  :one dimension approach       
		ENDIF
 

!************  SOLVING NON LINEAR EQUATION IN ONE DIMENSTION *********************************************

!      ignore the effect of precip advected
!      energy on the calculation of surface temperature when there is no snow.
!      Without this ridiculously high temperatures can result as the model
!      tries to balance outgoing radiation with precip advected energy.    
!


 14      Tssk = Tak                         ! first approximation
	   Tck  = Tak

 17       ERc    = tol*2                           ! so that it does not end on first loop
         iterC  =  0.
! Estimate Tc based on assumed snow surface temperature
 
13     Tclast = Tck 
       Tc    = CanTemp(Tck, Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Qsi,atff,Cf, &
     	 Qli,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
       Qpc, Inmax, Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
         fkappaS,RHO,TherC,Fs,  &                                          
         tave,refDepth,SmeltC)                            ! Reduced parametre later 


                         
 18      ER    = tol*2                           ! so that it does not end on first loop ! This is the start when LAI is zero
       niter = 0
 1     if(ER.gt.tol.and.niter.lt. nitermax)then  
           Tslast = Tssk

        F1 = SURFEBsc(Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC, &                                           
       TSURFs,tave,refDepth)   
                                                        	!yjs add three value to reflect model control changes)) 
           
		 
 	 F2 = SURFEBsc(fff*Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC, &                                           
       TSURFs,tave,refDepth)          	!yjs add three value to reflect model control changes)) 
          
		
	 Tssk = Tssk - ((1.-fff) * Tssk * F1) / (F1 - F2)
	   if(Tssk .lt. Tak - 50)go to 11                      !If it looks like it is getting unstable go straight to bisection
	   ER = abs(Tssk - Tslast)
         niter=niter+1
	   go to 1        ! loop back and iterate
	   endif

	if(er.le.tol) goto 10                              ! The solution has converged
  
!   If still not converged use bisection method
 11    Tlb = TaK - 20.                             ! First guess at a reasonable range                 
	 Tub = Tak + 10.

	 Flb = SURFEBsc(Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,      &                                      
       TSURFs,tave,refDepth)      
     	
	 Fub = SURFEBsc(Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,    &                                        
       TSURFs,tave,refDepth)  
     	
	ibtowrite=0
       if(Flb*fub .gt. 0.)then     ! these are of the same sign so the range needs to be enlarged
		Tlb= TaK - 150.          ! an almost ridiculously large range - solution should be in this if it exists
		Tub= Tak + 100.  

	 Flb = SURFEBsc(Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)  
                                              
     	
	 Fub = SURFEBsc(Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)  
     

	 ibtowrite=1

          if(Flb*fub .gt. 0.)then   ! these are of the same sign so no bisection solution
              If (snowdgtvariteflag .EQ. 1)then
			        write(66,*)  &
              'Bisection surface temperature solution failed with large range'
			        write(66,*)'Date: ',mtime(1),mtime(2),mtime(3)
			        write(66,*)'time: ',mtime(4)
			        write(66,*)'Model element: ',mtime(5)
			        write(66,*) &
              'A surface temperature of 273 K assumed'
	          end if
              Tssk=tk
	          go to 10 
	    else
	        If (snowdgtvariteflag .EQ. 1)then
		        write(66,*) &
             'Bisection surface temperature solution with large range'
		        write(66,*)'Date: ',mtime(1),mtime(2),mtime(3)
		        write(66,*)'time: ',mtime(4)
		        write(66,*)'Model element: ',mtime(5)
		        write(66,*) &
 	        'This is not a critical problem unless it happens frequently'
                 write(66,*) &
              'and solution below appears incorrect'
            end if
	     endif
        else
        endif

!     Here do the bisection
       niter = log((Tub-Tlb)/tol)/log(2.)   ! Number of iterations needed for temperature tolerance
	 do iter  =1,niter
           Tssk = 0.5*(tub+tlb)
           F1   = SURFEBsc(Tssk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,   &                                         
       TSURFs,tave,refDepth)  
                                              	!yjs add three value to reflect model control changes)) 
		 if(f1.gt.0.0) then  ! This relies on the assumption (fact) that this is a monotonically decreasing function
			tlb=tssk
	     else
			tub=tssk
	     endif
       enddo
	if(ibtowrite .eq. 1)then
	    If (snowdgtvariteflag .EQ. 1)then
		    write(66,*)'Surface temperature: ',Tssk,' K'
		    write(66,*)'Energy closure: ',f1
		    write(66,*)'Iterations: ',niter
	    end if
	endif


 10    Tss = Tssk - Tk


!**************   check if snow is melting  ********************************

       IF(Ws.GT.0..AND.Tss.GT.0.) THEN

!dgt 5/4/04 surface melt smelt
          SRFTMPsc = 0.0
	smelt= SURFEBsc(SRFTMPsc+Tk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf, &
      Qli,Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod, &
      mtime, Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,    &                                        
       TSURFs,tave,refDepth) 
     
!dgt 5/4/04 surface melt smelt is the energy not accommodated by conduction into the snow
!  so it results in surface melt which then infiltrates and adds energy to the snowpack
!  through refreezing
!  This change added to capture this quantity when the model is used for a glacier where
!  surface melt may run off rather than infiltrate.  
!  For modeling of glaciers in Antarctica by Mike Gooseff

       ELSE
          SRFTMPsc = Tss 
	    smelt=0.0           !  No surface melt this case
	ENDIF
	
     IF ((LAI .EQ. 0.0) .OR. (cc .EQ. 0.0)) THEN
       !IF(LAI.EQ.0.)THEN
	  GO TO 21
  	 ENDIF
!****************  end melt check    ***************************************

! ITERATION FOR THE NEW Tc, FROM ESTIMATED Tss
        Tssk= SRFTMPsc+Tk
        Tck = Tc+Tk

      IF(ERc.gt.tol.and.iterC.LT. 10.)THEN 
          ERc    = abs(Tck -  Tclast)
         iterC = iterC+1
	  go to 13                  ! To estimate the new TC for the estimated Tss. loop back
      ENDIF

	 
21	 RETURN
       END

!************  END OF ITERATION TO SOLVE LINEAR EQUATION IN ONE DIMENSTION ************************
!*************************************************************************************************
!                  COMPUTE THE CANOPY TEMPERATURE

!************************************************************************************************** 
      FUNCTION CanTemp(Tck, Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Qsi,atff,Cf, &
     	Qli, COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
       Qpcin, Inmax, Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
         fkappaS,RHO,TherC,Fs,   &                                         
         tave,refDepth,SmeltC) 
	
  						        	          
	real Tssk,k,Inmax, SmeltC
	REAL param(*)
	REAL sitev(*)
	real mtime(*)
	integer iTsMethod

!     common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)   
      data Rag /287.0/     !  Ideal Gas constant for dry air (287 J/kg/K)    
      data k  /0.4/        !  Von Karmans constant (0.4)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
!  Parameters
      Ems  = param(3)      !  emmissivity of snow (nominally 0.99)
      cg   = param(4)	   !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      z    = param(5)	   !  Nominal meas. height for air temp. and humidity (2m)
      zo   = param(6)	   !  Surface aerodynamic roughness (m)
      rho  = param(7)	   !  Snow Density (Nominally 450 kg/m^3)
      rhog = param(8)	   !  Soil Density (nominally 1700 kg/m^3)
      lans = param(15)	   !  thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
      lang = param(16)     !  the thermal conductivity of soil (9.68 kJ/m/k/hr)
      wlf  = param(17)     !  Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr) 
      rd1  = param(18)     !  Apmlitude correction coefficient of heat conduction (1)
      Apr  = sitev(2)      !  Atmospheric Pressure (Pa)



      data tol,nitermax/0.05,20/   ! changed 10 to 20 and .05 to .0005
      fff    = (273-tol)/273.                 ! Scale the difference used for derivatives by the tolerance
      TAK    = TA+TK
      TAVEK  = TAVE+TK

      Qpc  =  Qpcin                          ! store input variable locally without changing global value

      if(Wc.le.0 .and. Qpc.gt.0.0) Qpc=0.0 
!
!      ignore the effect of precip advected
!      energy on the calculation of surface temperature when there is no snow.
!      Without this ridiculously high temperatures can result as the model
!      tries to balance outgoing radiation with precip advected energy.    
!
!       Tck  = TaK                              ! first approximation

       ER    = tol*2                           ! so that it does not end on first loop
       niter = 0
 1     if(ER.gt.tol.and.niter.lt. nitermax)then  
           Tclast = Tck

        F1 = SURFEBc(Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)             	!yjs add three value to reflect model control changes)) 

		 
 	 F2 = SURFEBc(fff*Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC, &                                           
       TSURFs,tave,refDepth)          	!yjs add three value to reflect model control changes)) 
          
		
	 Tck = Tck - ((1.-fff) * Tck * F1) / (F1 - F2)
	   if(Tck .lt. Tak - 20)go to 11                      !If it looks like it is getting unstable go straight to bisection
	   ER = abs(Tck - Tclast)
         niter=niter+1
	   go to 1        ! loop back and iterate
	   endif

!	if(abs(F1-F2).GT.100000 ) goto 11 
	if(er.le.tol) goto 10                              ! The solution has converged
  
!   If still not converged use bisection method
 11    Tlb = TaK - 20.                             ! First guess at a reasonable range                 
	 Tub = Tak + 10.

	 Flb = SURFEBc(Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)      
     	
	 Fub = SURFEBc(Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)  
     	
	ibtowrite=0
       if(Flb*fub .gt. 0.)then     ! these are of the same sign so the range needs to be enlarged
		Tlb= TaK - 150.          ! an almost ridiculously large range - solution should be in this if it exists
		Tub= Tak + 100.  

	 Flb = SURFEBc(Tlb,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,     &                                       
       TSURFs,tave,refDepth)  
                                              
     	
	 Fub = SURFEBc(Tub,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,      &                                      
       TSURFs,tave,refDepth)  
     

	 ibtowrite=1

          if(Flb*fub .gt. 0.)then   ! these are of the same sign so no bisection solution
              If (snowdgtvariteflag .EQ. 1)then
			        write(66,*) &
              'Bisection canopy temperature solution failed with large range'
			        write(66,*)'Date: ',mtime(1),mtime(2),mtime(3)
			        write(66,*)'time: ',mtime(4)
			        write(66,*)'Model element: ',mtime(5)
			        write(66,*) &
              'A canopy temperature of 273 K assumed'
              end if
              Tck=tk
		      go to 10 
	     else
	        If (snowdgtvariteflag .EQ. 1)then
			        write(66,*) &
                 'Bisection canopy temperature solution with large range'
			        write(66,*)'Date: ',mtime(1),mtime(2),mtime(3)
			        write(66,*)'time: ',mtime(4)
			        write(66,*)'Model element: ',mtime(5)
			        write(66,*) &
     	        'This is not a critical problem unless it happens frequently'
	                 write(66,*) &
                  'and solution below appears incorrect'
            END IF 
	     endif
        else
        endif

!     Here do the bisection
       niter = log((Tub-Tlb)/tol)/log(2.)   ! Number of iterations needed for temperature tolerance
	 do iter  =1,niter
           Tck = 0.5*(tub+tlb)
           F1   = SURFEBc(Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,           &                                 
       TSURFs,tave,refDepth)  
                                              	!yjs add three value to reflect model control changes)) 
		 if(f1.gt.0.0) then  ! This relies on the assumption (fact) that this is a monotonically decreasing function
			tlb=tck
	     else
			tub=tck
	     endif
       enddo
	if(ibtowrite .eq. 1)then
	    If (snowdgtvariteflag .EQ. 1)then
		    write(66,*)'Surface temperature: ',Tck,' K'
		    write(66,*)'Energy closure: ',f1
		    write(66,*)'Iterations: ',niter
	    END IF
	endif


 10    Tc = Tck - Tk


       IF(Wc.GT.0..AND.Tc.GT.0.) THEN

!dgt 5/4/04 surface melt smelt
          CanTemp = 0.0
	smeltC= SURFEBc(CanTemp+Tk,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf, &
     	Qli,Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob, &
     	iTsMethod,mtime,Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz, &
     	Tssk,Tk,Tak,EA,RHOA, FkappaS,RHO,TherC,   &                  
       TSURFs,tave,refDepth) 


       CanTemp =  (Tck*(1-Fs)+(CanTemp+Tk)*Fs)-Tk      ! Equiv canopy temp during melt
  

!dgt 5/4/04 surface melt smelt is the energy not accommodated by conduction into the snow
!  so it results in surface melt which then infiltrates and adds energy to the snowpack
!  through refreezing
!  This change added to capture this quantity when the model is used for a glacier where
!  surface melt may run off rather than infiltrate.  
!  For modeling of glaciers in Antarctica by Mike Gooseff
       ELSE
          CanTemp = Tc
	    smeltC  = 0.0
!  No surface melt this case
       ENDIF
       RETURN
       END


!*************************************************************************************************************************
!*****      FUNCTION TO EVALUATE THE SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR SURFACE TEMPERATURE     ****** 
!*                 DGT and C Luce 4/23/97
!************************************************************************************************************************** 
      FUNCTION SURFEBsc(Tssk, Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tck,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,  &                                          
       TSURFs,tave,refDepth)                                                                 ! Heat and vapor conductance for neutral

	real param(*),sitev(*),mtime(*), Inmax
	integer iTsMethod,Iradfl
	real LanS, LanG,LanE,LanE_Ze, LanE_de, LanE_ze2,LanE_de2
	real cg,rho,rhog,w1f,rd1,fkappas,ds,TherC
      real tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave,TsAvek

!	common /ts_save/ ts_old, tave_old, Ts_ave, Tave_ave 
      	 common /tsk_save/ tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave

!     Constant data set 
      data cs /2.09/				    !  Ice heat capacity (2.09 KJ/kg/C)
      data rhow /1000.0/				!  Density of Water (1000 kg/m^3)
	data w1day /0.261799/	        !  Daily frequency (2pi/24 hr 0.261799 radians/hr) 
      data sbc /2.0747e-7/            !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data Rag /287.0/                !  Ideal Gas constant for dry air (287 J/kg/K)
	data hneu /2834.0/              !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cp /1.005/                 !  Air Heat Capacity (1.005 KJ/kg/K)
	data Rag /287.0/                !  Ideal Gas constant for dry air (287 J/kg/K) (name changed)

!     Parameters
      cg   = param(4)				    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
	z    = param(5)					!  Nominal meas. height for air temp. and humidity (2m)
      zo   = param(6)					!  Surface aerodynamic roughness (m)
      rho  = param(7)				    !  Snow Density (Nominally 450 kg/m^3)
      rhog = param(8)					!  Soil Density (nominally 1700 kg/m^3)
      lans = param(15)					! the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	lang = param(16)					! the thermal conductivity of soil (9.68 kJ/m/k/hr)
	wlf  = param(17)					! Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr) 
	rd1  = param(18)					! Apmlitude correction coefficient of heat conduction (1)
  
      APr   =   sitev(2)     !  Atmospheric Pressure (Pa)
	Cc    =   sitev(5)     !  Canopy Coverage
      Hcan  =   sitev(6)     !  
      LAI   =   sitev(7)     !  

	Zs     =  Ws*rhow/rho                   ! Snow Depth  
	Tsavek = Tave + Tk
      Qp     =  Qps


!     07/25/02   at Boise.  To make the UEB2 work under Unix/linux system the fancy stuff like "Select case" shall not be used
	    !select case(iTsMethod)				!Four choice of the surface temperature modeling
	  

	!case (1)									!1. Old snow, use calibrated snow surface heat conductance
         	    if(iTsMethod .eq. 1) then
			qcs = RHO*CS*TherC*(Tssk-Tsavek)	!2. Revised scheme LanE/Ze of the snow surface
												!3. Force restore approach
	
	!case (2)									!4. Modified force restore approach.
          elseif (iTsMethod .eq. 2) then
			fKappaS = LanS/(rho*cs)

			fKappaG = LanG/(rhog*cg)

			d1 = sqrt(2*fKappaS/w1day)
			if(zs .ge. rd1*d1) then

				LanE_Ze=LanS/(rd1*d1)
	        else
			
				LanE_Ze=LanE(LanS, LanG, Zs, rho, rhog, cs, cg, rd1,  &
     				ze, w1day)						!Cyjs   call the subroutine to update the heat conductivity. LanE()
				LanE_Ze=LanE_Ze/ze
			end if

      			qcs= LanE_Ze*(Tssk-Tsavek)

	!case (3)
		elseif (iTsMethod .eq. 3) then 
			fKappaS = LanS/(rho*cs)
			fKappaG = LanG/(rhog*cg)

			d1 = sqrt(2*fKappaS/w1day)

			if(zs .ge. rd1*d1) then
				LanE_Ze=LanS/(rd1*d1)
				Ze=rd1*d1
	        else
			
				LanE_Ze=LanE(LanS, LanG, Zs, rho, rhog, cs, cg, rd1,   &
     				ze, w1day) !Cyjs   call the subroutine to update the heat conductivity. LanE()
				LanE_Ze=LanE_Ze/ze
			end if
			
			de=ze/rd1
			LanE_de=LanE_ze/de*ze
		qcs= LanE_de*(Tssk-Tssk_old)/(w1day*dt)+LanE_Ze*(Tssk-Tsavek)
		
		else 


	!case (4)						!change to all default cases. If not for model comparison

			fKappaS = LanS/(rho*cs)
			fKappaG = LanG/(rhog*cg)

			d1  = sqrt(2*fKappaS/w1day)
			dlf = sqrt(2*fKappaG/wlf) 	
			 
			if(zs .ge. rd1*d1) then
				LanE_Ze=LanS/(rd1*d1)
				Ze=rd1*d1
              else
				LanE_Ze=LanE(LanS,LanG,Zs,rho,rhog,cs,cg,rd1,ze,w1day)   !Cyjs   call the subroutine to update the heat conductivity. 
				LanE_Ze=LanE_Ze/ze
			end if  

			if(zs .ge. rd1*dlf) then
				LanE_Ze2=LanS/(rd1*dlf)
				Ze2=rd1*dlf
			else
				LanE_Ze2=LanE(LanS,LanG,Zs,rho,rhog,cs,cg,rd1,ze2,wlf)  !Cyjs   call the subroutine to update the heat conductivity.                   
			    LanE_Ze2=LanE_Ze2/ze2  
			endif

			de=ze/rd1
			LanE_de=LanE_ze/de*ze
			de2=ze2/rd1
			LanE_de2=LanE_ze2/de2*ze2

	if(Us.le.0.0 .or. refDepth.le.0.0)then
	  qcs= LanE_de*(Tssk-Tssk_old)/(w1day*dt)+LanE_Ze*(Tssk-Tssk_Ave)+ &
     			LanE_de2*(Tssk_ave-Tsavek_Ave)
	       elseif(refDepth .gt. rd1*d1) then
	  qcs= LanE_de*(Tssk-Tssk_old)/(w1day*dt)+LanE_Ze*(Tssk-Tssk_Ave)+ &
     			LanE_de2*(Tssk_ave-Tsavek_Ave)
             else
                qcs=lanE_ze*ze*(tssk-tk)/refDepth
	      endif


		!End select
            endif

	Ess    = SVP(Tssk-Tk)
	Esc    = SVP(Tck-Tk)  

	 CALL TURBFLUX(Ws,Wc,A,TK,Tck-Tk,Ta,Tssk-Tk,RH,V,EA,p,param,sitev, &
     	              d,Z0c,Vz,RKINc,RKINsc,Tac, Fs,Ess,Esc,  &  ! Output variables
                          QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E)          

       SURFEBsc = QP+QHs+QEs- Qcs



	IF (iradfl.EQ.0.0) THEN

	 CALL PSOLRAD(Qsi,atff,param,cf,  &                  ! Input Variables 
     	                    Taufb,Taufd,Qsib,Qsid)    ! Output variables:


	 CALL TRANSRADCAN (COSZEN,sitev,param, &               ! Input variables , leaf parameters
     	         Betab,Betad,Taub,Taud)   ! Output variables:

	 CALL NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, &
             Inmax,Qsib,Qsid,param,Fs, &
                            Qsns,Qsnc ) !  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
	 CALL NETLONGRAD(RH,Ta,Tssk-Tk,Tck-Tk,Tk,Fs,EmC,EmS,SBC,cf,sitev, &
                  	Qli,param, &
     	                Qlis,Qlns,Qlnc )                     !  Output: Qsns,Qsnc


		SURFEBsc = SURFEBsc + QSNs + QLNs

      ELSE
 		SURFEBsc = SURFEBsc + qnetob

	ENDIF

	RETURN
	END
!*****************************************************************************************************************************************	

!*************  FUNCTION TO EVALUATE THE CANPPY SURFACE ENERGY BALANCE FOR USE IN SOLVING FOR CANOPY TEMPERATURE *********************** 
                 
!*****************************************************************************************************************************************
 
      FUNCTION SURFEBc(Tck,Us,Ws,Wc,A,dt,P,Pr,Ps,Ta,V,RH,Fs,Cf,Qli, &
     Qsi,atff,COSZEN,EmC,EmS,param,sitev,iradfl,qnetob,iTsMethod,mtime, &
        Qpc,Qps,Inmax, Rkinc,Rkinsc,Vz,Tssk,Tk,Tak,EA,RHOA, &
        FkappaS,RHO,TherC,    &                                        
       TSURFs,tave,refDepth)                             ! Reduced parametre later                                    ! Heat and vapor conductance for neutral

	real param(*),sitev(*),mtime(*), Inmax
	integer iTsMethod,Iradfl
	real LanS, LanG,LanE,LanE_Ze, LanE_de, LanE_ze2,LanE_de2
	real cg,rho,rhog,w1f,rd1,fkappas,ds,TherC
      real tssk_old, tsavek_old, Tsavek_Ave, Tssk_ave,TsAvek


!     Constant data set 
      data cs /2.09/				    !  Ice heat capacity (2.09 KJ/kg/C)
      data rhow /1000.0/				!  Density of Water (1000 kg/m^3)
	data w1day /0.261799/	        !  Daily frequency (2pi/24 hr 0.261799 radians/hr) 
      data sbc /2.0747e-7/            !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
      data Rag /287.0/                !  Ideal Gas constant for dry air (287 J/kg/K)
	data hneu /2834.0/              !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data cp /1.005/                 !  Air Heat Capacity (1.005 KJ/kg/K)
	data Rag /287.0/                !  Ideal Gas constant for dry air (287 J/kg/K) (name changed)

!     Parameters
      cg   = param(4)				    !  Ground heat capacity (nominally 2.09 KJ/kg/C)
	z    = param(5)					!  Nominal meas. height for air temp. and humidity (2m)
      zo   = param(6)					!  Surface aerodynamic roughness (m)
      rho  = param(7)				    !  Snow Density (Nominally 450 kg/m^3)
      rhog = param(8)					!  Soil Density (nominally 1700 kg/m^3)
      lans = param(15)					! the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	lang = param(16)					! the thermal conductivity of soil (9.68 kJ/m/k/hr)
	wlf  = param(17)					! Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr) 
	rd1  = param(18)					! Apmlitude correction coefficient of heat conduction (1)
  
      APr  = sitev(2)                 !  Atmospheric Pressure (Pa)


       IF(IRADFL.EQ.1) Go to 13 

	 CALL PSOLRAD(Qsi,atff,param,cf,  &                  ! Input Variables 
     	                    Taufb,Taufd,Qsib,Qsid)        ! Output variables:

	 CALL TRANSRADCAN (COSZEN,sitev,param, &               ! Input variables , leaf parameters
     	         Betab,Betad,Taub,Taud)                     ! Output variables:

	 CALL NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, &
     	 Inmax,Qsib,Qsid,param,Fs, &
                            Qsns,Qsnc )                       !  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
    
	 CALL NETLONGRAD(RH,Ta,Tssk-Tk,Tck-Tk,Tk,Fs,EmC,EmS,SBC,cf,sitev, &
                  	Qli,param, &
     	                Qlis,Qlns,Qlnc )                     !  Output: Qsns,Qsnc



13	Ess    = SVP(Tssk-Tk)
	Esc    = SVP(Tck-Tk)  

	 CALL TURBFLUX(Ws,Wc,A,TK,Tck-tk,Ta,Tssk-Tk,RH,V,EA,p,param,sitev, &
     	              d,Z0c,Vz,RKINc,RKINsc,Tac, Fs,Ess,Esc,   &                    ! Output variables
                          QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E)          

         	    IF(Iradfl.eq.0)THEN

			    SURFEBc = QSNc+QLNc+QPc+QHc+QEc
              
			ELSE
     			    surfeb = surfeb + qnetob
          	ENDIF

	RETURN
	END


!********** FUNCTION TO ESTIMATE SURFACE HEAT CONDUCTION FOR USE IN SOLVING SURF TEMP *********
!                            ********  QcEst() *******

       FUNCTION QcEst(Ws,p,Tssk,Tck,V,Zm,d,Z0c,Rimax,Rcastar,Cf,Fs, &
     		Qli,Hcan,Vz,Ta,Rh,RKINsc,Qps,To,Ps,Qsi,atff,COSZEN, &
          APr,TAK, EA,A,Ac,Wc,Inmax, Qnetob,Iradfl,param,sitev)

	REAL k,n,Inmax,RH, Qps,Qnetob, RHOA,param(*),sitev(*)
	INTEGER Iradfl

	data tk /273.15/     !  Temperature to convert C to K (273.15)
	data Rag /287.0/     !  Ideal Gas constant for dry air (287 J/kg/K) (name changed)
      data k  /0.4/        !  Von Karmans constant (0.4)
      data hff /3600.0/    !  Factor to convert /s into /hr (3600)
      data rhow /1000.0/   !  Density of Water (1000 kg/m^3)
      data g    /9.81/     !  Gravitational acceleration (9.81 m/s^2)
      data cw /4.18/       !  Water Heat Capacity (4.18 KJ/kg/C)
      data cs /2.09/       !  Ice heat capacity (2.09 KJ/kg/C)
      data cp /1.005/      !  Air Heat Capacity (1.005 KJ/kg/K)
      data hf /333.5/      !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0/   !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)
      data sbc /2.0747e-7/ !  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
     
      Ems   =   param(3)   !  emmissivity of snow (nominally 0.99)
	cg    =   param(4)   !  Ground heat capacity (nominally 2.09 KJ/kg/C)
      z     =   param(5)   !  Nominal meas. height for air temp. and humidity (2m)
      zo    =   param(6)   !  Surface aerodynamic roughness (m)
      fstab =   param(19)  !  Stability correction control parameter 0 = no corrections, 1 = full corrections
	EmC   =   param(23)	 !  Emissivity of canopy 
  	Cc    =   sitev(5)   !  Canopy Coverage
	LAI   =   sitev(7)   !  Leaf Area Index  


	If (Wc.EQ.0.)THEN                          
	 Tck=Tak                                     ! Tc = Ta assumed here when there is no snow in canopy
	ELSE
	Tck= Tk                                     ! When snow at gound is melting, canopy snow is melting (assumed)
	ENDIF
	Tss= Tssk-Tk


	Ess    = SVP(Tssk-Tk)
	Esc    = SVP(Tck-Tk) 

	 CALL TURBFLUX(Ws,Wc,A,TK,Tck-Tk,Ta,Tssk-Tk,RH,V,EA,p,param,sitev, &
     	              d,Z0c,Vz,RKINc,RKINsc,Tac, Fs,Ess,Esc,   &   ! Output variables
                      QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E)          

              QcEst = QPs+QHs+QEs     
 
                                                     
	IF(Iradfl.eq.0)THEN 

		CALL PSOLRAD(Qsi,atff,param,cf  &                  ! Input Variables 
     	                    ,Taufb,Taufd,Qsib,Qsid)    ! Output variables:


		CALL TRANSRADCAN (COSZEN,sitev,param,  &              ! Input variables , leaf parameters
     	         Betab,Betad,Taub,Taud)   ! Output variables:

	 CALL NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, &
     	     Inmax,Qsib,Qsid,param,Fs, &
                            Qsns,Qsnc ) !  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 

	CALL NETLONGRAD(RH,Ta,Tssk-Tk,Tck-Tk,Tk,Fs,EmC,EmS,SBC,cf,sitev, &
                  	Qli,param, &
     	                Qlis,Qlns,Qlnc )                     !  Output: Qsns,Qsnc 


	qcEst = qcEst + QSNs + QLNs

      ELSE
 	qcEst = qcEst + qnetob

     	ENDIF

      RETURN
	END

!*********************** FMELT () ****************************************************************
!     Calculates the melt rate and melt outflow

      FUNCTION FMELT(UB,RHOW,W,HF,LC,RID,KS,PRAIN)
      REAL LC,KS

!       write(*,*) 'I am in FMELT!!!!'

      UU=0.0
      IF(UB.LT.0.) THEN
         FMELT=0.0
      ELSEIF(W.LE.0.0) THEN
         FMELT=PRAIN
         IF(PRAIN.LE.0.) THEN
            FMELT=0.0
         ENDIF
      ELSE
         UU=UB/(RHOW*W*HF)
!                              liquid fraction
         IF(UU.GT.0.99) THEN
            UU=0.99
!                    TO ENSURE HIGH MELT RATE WHEN ALL LIQUID
      ENDIF  

      IF((UU/(1-UU)).LE.LC) THEN
         SS=0.0
      ELSE
         SS=(UU/((1- UU)*RID)-LC/RID)/(1-LC/RID)
      ENDIF
      FMELT=KS*SS**3
      ENDIF
      IF(FMELT.LT.0.0) THEN
          WRITE(*,*)'FMELT is NEGATIVE!'
          STOP
       ENDIF

!       write(*,*) ' I am leaving FMELT!!!'

       RETURN
       END


!******************************************** Grad () ********************************************
! Linear simple function to calculate the gradient of Q and T

      subroutine Grad(qc1,qc2,t1,t2, a, b)
	if((t2-t1) .ne. 0.0) then
	  b=(qc2-qc1)/(t1-t2)
	  a=qc1+b*t1
	endif

	return
	end


!*************************** refDep() ************************************************************
!     function to calculate the refreezing depth
	real function refDep(flans, a, b, hf, rhom, dt, x1 )

	temp=flans*flans+2*b*(-a*flans*dt/(hf*rhom)+flans*x1+0.5*b*x1*x1)

	if(temp .lt. 0.0 .or. a.gt.0.0 .or. b.eq.0.0) then
	  refDep=0
	else 
	  refDep=(-flans+sqrt(temp))/b
  	 
	endif        	
	
	return 
	end



!******************** CALCULATE THE AVE.TEMP.OF SNOW AND INTERACTING SOIL LAYER ********************

      FUNCTION  TAVG(UB,W,RHOW,CS,TO,RHOG,DE,CG,HF)

      SNHC = RHOW*W*CS							 ! SNHC = Snow heat capacity
      SHC  = RHOG*DE*CG							 ! SHC  = Soil heat capacity
      CHC  = SNHC+SHC								 ! CHC  = Combined heat capacity

      IF(UB.LE.0.) THEN
         TS=UB/CHC
      ELSE 
         AL=UB/(RHOW*HF)
         IF(W.GT.AL) THEN
            TS=TO                                   !DE   = Therm. active depth of soil (param(11)) 
         ELSE
            TS=(UB-W*RHOW*HF)/CHC
         ENDIF
      ENDIF
      TAVG=TS

      RETURN
      END


!**************************************************************************************************
!********************* Caculate the daily average  *********************
!yjs  Calculate the daily average value
!yjs  n number of records, a minimum value -100 or some

	REAL FUNCTION daily_ave(backup, n, a) 
	REAL backup(*)
	sum = 0.
	count=0.
	DO 1 i = 1, n
	IF(backup(i).gt.a) THEN
		sum=sum+backup(i)
		count=count+1.
	END IF

 1	CONTINUE
	IF(count.ne.0) THEN 
		daily_ave=sum/count
	ELSE 
		daily_ave=0.
	END IF
	RETURN
	END


!*************************************************************************************************
!***************************** LanE () *************************************
!	Function to get the LanE which is the thermal conductivity by ze

	real function LanE(LanS,LanG,Zs,rho,rhog,cs,cg,r,ze,w1day)

	real LanS, LanG,Zs,rho,rhog,cs,cg,r,ze,w1day
!	real fKappaS, fKappG, d1, d2
	

	 fKappaS = LanS/(rho*cs)
	 fKappaG = LanG/(rhog*cg)

	 d1 = sqrt(2*fKappaS/w1day)
	 d2 = sqrt(2*fKappaG/w1day)	 

	 LanE=amin1(r,zs/d1)*d1/LanS + amax1(0.,(r-zs/d1))*d2/LanG

	 ze=amin1(r,zs/d1)*d1+amax1(0.,(r-zs/d1))*d2
	 
	 if (LanE .ne. 0.) LanE = 1/LanE*ze
	 return
	end


!******************************** ALBEDO () *****************************************************
!     Function to calculate Albedo
!     BATS Albedo Model (Dickinson et. al P.21)

      FUNCTION ALBEDO(tausn,coszen,d,aep,abg,avo,airo)     ! snow age, zenith angle, snow depth
                                                           ! albedo extinpara,bare Gr albedo, visoble and NIR for fresh snow
      B  = 2.0                                             ! aep.input file
      CS = 0.2
      CN = 0.5


      FAGE = tausn/(1.0+tausn)

      IF(coszen.LT.0.5) THEN
         FZEN = 1.0/B*((B+1.0)/(1.0+2.0*B*coszen)-1.0)
      ELSE
         FZEN = 0.0
      ENDIF

      AVD = (1.0-CS*FAGE)*AVO
      AVIS = AVD+0.4*FZEN*(1.0-AVD)
      AIRD = (1.0-CN*FAGE)*AIRO
      ANIR = AIRD+0.4*FZEN*(1.0-AIRD)
      ALBEDO = (AVIS+ANIR)/2.0
      if(d.lt.aep)then                               ! need to transition albedo to a bare ground value
        rr=(1.-d/aep)*exp(-d*.5/aep)
        albedo=rr*abg+(1.-rr)*albedo
      endif

      RETURN
      END

!******************************** AGESN () *******************************************************
!     Function to calculate Dimensionless age of snow for use in
!     BATS Albedo Model (Dickinson et. al P.21)

      SUBROUTINE AGESN(tausn,dt,Ps,tsurf,tk,dNewS) 
	     
!     Input Variables  : dt, ps,tsurf,tk,dnewS (threshold depth of new snow, param 21)
!     Output Variable  : Agesn

      tsk= tsurf+tk                                         ! Express surface temperature in kelvin
      R1 = EXP(5000.0*(1.0/TK - 1.0/TSK))
      R2 = R1**10

      IF(R2.GT.1.0) THEN
         R2 = 1.0
      ENDIF

      R3 = 0.03
      tausn = tausn+0.0036*(R1+R2+R3)*DT
	if(ps.gt.0.0) then
	  if(dNewS .gt. 0.0) then
	    tausn = max(tausn*(1-ps*dt/dNewS),0.)
	  else
	    tausn = 0.0
	  endif
	endif

      RETURN
      END

!************************** PARTSNOW () **********************************************************
!     Partitioning of precipitation into rain and snow
      
      FUNCTION PARTSNOW(P,TA,TR,TS)

      IF(TA.LT.TS) THEN
        PARTSNOW=P
      ELSEIF(TA.GT.TR) THEN
        PARTSNOW=0.0
      ELSE
        PARTSNOW=P*(TR-TA)/(TR-TS)
      ENDIF

      RETURN
      END