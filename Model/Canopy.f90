!*******************************************************************************************************************
!*		               CANOPY RADIATION TRANSMISSION PROCESSES
!*******************************************************************************************************************
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
!********** PARTITION OF INCOMING SOLAR RADIATION INTO DIRECT AND DIFFUSE BEFORE IT HITS THE CANOPY *************       
!     Partition the incoming solar radiation into direct and diffuse components

   SUBROUTINE PSOLRAD(Qsi,atff,param, cf, &       ! Input Variables 
        Taufb,Taufd,Qsib,Qsid)     ! Output variables:


        REAL Qsi,cf,Taufb,Taufd,Qsib,Qsid,Taufc
        REAL param(*),as,bs,Lambda
		
        as     = param(28)			! Fraction of extraterrestaial radiation on cloudy day,Shuttleworth (1993)  
	bs     = param(29)			! (as+bs):Fraction of extraterrestaial radiation on clear day, Shuttleworth (1993) 
	Lambda = param(30)		    ! Ratio of diffuse atm radiation to direct for cler sky,from Dingman (1993)



       Taufc =max(atff,as+bs)     ! Mac atmosphereci transmittance for clear sky
   
	IF (cf.EQ.0) THEN 
		Taufb = Lambda*Taufc				 ! Direct solar radiation atm.transmittance
		Taufd = (1.-Lambda)*Taufc			 ! Diffuse solar radiation atm.transmittance
	ELSE IF (cf.EQ.1) THEN
		Taufb = 0.0						     
		Taufd = atff
	ELSE 
		Taufb = Lambda*(1.-cf)*(as+bs)
		Taufd = atff- Taufb
	END IF
    
!     Direct canopy solar radiation incident at the canopy 

		Qsib  = Taufb*Qsi/(Taufb+Taufd)
!	Diffuse canopy solar radiation incident at the canopy   
     
		Qsid  = Taufd*Qsi/(Taufb+Taufd)    
	
	RETURN
	END

!*********** TRANSMISSION OF DIRECT AND DIFFUSE SOLAR RADIATION THROUGH THE CANOPY ***************
!     Estimates the direct and diffuse solar radiation fractions transmitted through the canopy      
	SUBROUTINE  TRANSRADCAN (COSZEN,sitev,param, &      ! Input variables : zenith angle, leaf parameters etc.
     	         Betab,Betad,Taub,Taud)                     ! Output variables: Transmission and Reflection fractions

    
  REAL sitev(*),param(*),LAI,COSZEN,EXPI,G, kk
  REAL Taub,Taud,Taubh,Taudh, Beta, Betab, Betad


	Cc    = sitev(5)     ! Leaf Area Index 
	Hcan  = sitev(6)     ! Canopy height  
	LAI   = sitev(7)     ! Leaf Area Index  
	alpha = param(24)    ! 
      G     = param(26)    ! 

!      For open area where no vegetaion canopy is present                                       
	IF (LAI.EQ.0.)THEN 
	Taub = 1.
	Taud = 1.
      Betab = 0.
      Betad = 0.
 	ElSE 

!      Transmission without scattering 
      LAI   = Cc*LAI                      ! Effective leaf area index
      Rho   = LAI/Hcan                    ! leaf density per unit height of the canopy

!     Multiple scattering light penetration
	 kk    = sqrt(1-alpha)                             ! k-prime in the paper

	EXPI  = EXPINT(kk*G*LAI)                           ! Exponential integral function	

!      Deep Canopy Solution
!   Avoid divide by 0 error when COSZEN is 0
       if(COSZEN .LE. 0.)then
           Taubh=0.
       else
           Taubh =  EXP(- kk*G*Rho*Hcan/COSZEN)                 ! Transmission function for deep canopy : Direct
       endif
	 Taudh =  (1-kk*G*Rho*Hcan)*EXP(-kk*G*Rho*Hcan)+ &
     	      (kk*G*Rho*Hcan)**2*EXPI                     ! Transmission function for deep canopy : Diffuse
     	 Beta   = (1-kk)/(1+kk)                                ! BetaPrime : reflection coefficient for infinitely deep canopy

!     Finite Canopy Solution
       Taub   = Taubh*(1-Beta**2)/(1-Beta**2*Taubh**2)     ! Direct fraction of radiation transmitted down through the canopy 
       Taud   = Taudh*(1-Beta**2)/(1-Beta**2*Taudh**2)         ! Diffuse fraction of radiation transmitted down through the canopy 
       Betab  = Beta*(1-Taubh**2) /(1-Beta**2*Taubh**2)    ! Direct fraction of radiation scattered up from the canopy 
       Betad  = Beta*(1-Taudh**2) /(1-Beta**2*Taudh**2)       ! Diffuse fraction of radiation scattered up from the canopy 

	ENDIF
	RETURN
	END

!****************** NET CANOPY & SUB CANOPY SOLAR RADIATION *********************************
!      Computes the net canopy and sub-canopy solar radiation
!      considering the multiple scatterings of radiation by the canopy and multile reflections 
!      between the canopy and the snow surface
	SUBROUTINE NETSOLRAD(Ta,A,Betab,Betad,Wc,Taub,Taud, &
     	 Inmax,Qsib,Qsid,param,Fs, &
                            Qsns,Qsnc ) !  Output: Qsns,Qsnc (Net subcanopy, canopy solar radiation) 
      
	REAL A,Betab,Betad,Wc,Inmax,Taub,Taud,Qsib,Qsid
	REAL param(*),Qsns, Qsnc


!	Net radiation at snow surface below canopy

		  f1b = (1-A)*Taub/(1-Betad*(1-Taud)*A) 
            f1d = (1-A)*Taud/(1-Betad*(1-Taud)*A)              
            Qsns =  f1b*Qsib + f1d*Qsid

	IF (Taub.EQ.1.0)THEN              ! (For LAI=0, Qsnc gives the numerical errors, so)
	Qsnc =0.0
	ELSE	                              
      
!     Net canopy solar radiation	

		 f3b = ((1-A)*Betab + A*Taub*Taud)/(1-Betad*(1-Taud)*A)
		 f3d = ((1-A)*Betad + A*Taud*Taud)/(1-Betad*(1-Taud)*A) 
		 f2b = 1- f1b-f3b  
		 f2d = 1- f1d-f3d
		 Qsnc = f2b*Qsib + f2d*Qsid  

	ENDIF  

	RETURN
	END

!****************** NET CANOPY & SUB CANOPY LONGWAVE RADIATION ******************************
!     Computes the net canopy and beneath canopy longwave radiation
	SUBROUTINE NETLONGRAD(RH,Ta,Tss,Tc,Tk,Fs,EmC,EmS,SBC,cf,sitev,Qli, &
     	                       param, &
                                     Qlis,Qlns,Qlnc )                     !  Output: Qsns,Qsnc 
 
	REAL RH,Ta,Tss,Tc,EmC,EmS,SBC,Fs,TaudL,sitev(*),LAI,G
	REAL EmA,ea, Qlis, Qlns,Qlnc,EXPI,kk,Beta,Betad,param(*)


	Cc    = sitev(5)     ! Leaf Area Index 
	Hcan  = sitev(6)     ! Canopy height  
	LAI   = sitev(7)     ! Leaf Area Index  
	alphaL= param(25)    ! 
      G     = param(26)    ! 	                              

!	ea   = SVPW(TA)*RH          	! Air vapor pressure 
      TAK  = TA + TK
	Tssk = Tss +Tk  
      EmCS = EmC*(1-Fs)+EmS*Fs

	Tck  = Tc+Tk
      Qle  = EmS*SBC*Tssk**4

	IF (LAI.EQ.0.)THEN          !  For Open area where no vegetation caonpy is present
	  Qlcd = 0.
	  Qlcu = 0.
        Qlns = Qli*EmS -Qle                    ! To avoid numerical errors
        Qlnc = 0.
	ElSE 

!      Transmission without scattering      
	LAI      = Cc*LAI
      Rho      = LAI/Hcan  

!     Downward scattered 
	 kk    = sqrt(1-alphaL)
!     Upward scattered 
     	 Beta = (1-kk)/(1+kk)

	EXPI  = EXPINT(kk*G*LAI)                      	

!     Deep canopy solution
	 Taudh =  (1-kk*G*Rho*Hcan)*EXP(-kk*G*Rho*Hcan)+ &
     	      (kk*G*Rho*Hcan)**2*EXPI                               ! Transmission function for deep canopy : Diffuse

!	 Finite canopy solution
         Taud   = (Taudh-Beta**2*Taudh) /(1-Beta**2*Taudh**2)         ! Transmission function for finite canopy depth
         Betad  = (Beta-Taudh*Beta*Taudh) /(1-Beta**2*Taudh**2)       !  Scatterd up  from top of the canopy at 0 depth  
!      Canopy emitted longwave radiation 
	  Qlcd = (1-Taud)*EmC*SBC*Tck**4
	  Qlcu = (1-Taud)*EmCs*SBC*Tck**4

! Net sub canopy longwave radiation
       Qlns = Taud*Qli*EmS + Ems*Qlcd-Qle+ (1-Taud)*(1.-1.0)*Qle

! Net canopy longwave radiation
      Qlnc = (((1-Taud)*1)+((1-Ems)*Taud*1))*Qli + (1-Taud)*Qle*1   &      ! Putting 1 for Emc
      + (1-Taud)*(1-Ems)*1*Qlcu - Qlcu - Qlcd      

	ENDIF
      RETURN
      END

!********************************************************************************************************************

!*		               	TURBULENT TRANSFER PROCESSES

!********************************************************************************************************************

!****************  BULK AERODYNAMIC AND CANOPY BOUNDARY LAYER RESISTANCES  *******************
!      Calculates resistances for the above and beneath canopy turbulent heat fluxes
 
	SUBROUTINE AeroRes(p,Wc,V,Ta,Tss,Tc,Fs,param,sitev,Tk, &
                  d,Z0c,Vz,Rc,Ra,Rbc,Rl,RKINc,RKINa,Rkinbc,Rkinl)         ! Output variables

	REAL z,d,k,Zm,Rimax,YcAge,Hcan,LAI,Vz,Ta,Tss,RKINsc,Kh, &
     	 Lbmean,Rkinca,Ru,Ra,Rc,Rs,Rca,Ri ,Rcastar, Rustar,ndecay, &
        param(*),sitev(*), Tac, Dc,Rastar,Rkinaa,V

	data g   /9.81/        ! Gravitational acceleration (9.81 m/s^2)
      data k  /0.4/          ! Von Karmans constant (0.4)
      data hff /3600.0/      ! Factor to convert /s into /hr (3600)
      z     = param(5)       ! Nominal meas. height for air temp. and humidity (2m)
      zo    = param(6)       ! Surface aerodynamic roughness (m)
	Rimax = param(31)      ! Maximum value of Richardsion number for stability corretion
	fstab = param(19)      ! Stability correction control parameter 0 = no corrections, 1 = full corrections
      Wcoeff= param(32)      ! 

	Cc    = sitev(5)       ! Canopy Cover 
	Hcan  = sitev(6)       ! Canopy height  
	LAI   = sitev(7)       ! Leaf Area Index 

	Zm    = Hcan+z         ! observed wind speed above the canopy

      LAI = Cc*LAI
              
      IF(V.LE.0.)THEN
               Ra      = 0.
               Rc      = 0.
               Rs      = 0.
               Rbc     = 0.
               Rl      = 0.
        	   RKINc   = 0.                       
               RKINa   = 0.  
			 RKINbc  = 0.                       
               RKINl   = 0. 
			 Vz      = V 
			  !  If no wind, all resistances are zero; giving fluxes zero

      ELSE

        IF(LAI.EQ.0.)THEN      ! For open area
                       
          Vz      = V
          Rbc     = (1.0/(k**2*Vz)*(log(z/zo))**2.0 )*1/hff            ! sub layer snow resistance; hff: to convert resistance second to hour

	          Ri = g*(Ta-Tss)*Z /(Vz**2.0*(0.5*(Ta+Tss)+273.15))     
                                                    
	    		IF (Ri.GT.Rimax) THEN                                 
				Ri= Rimax
                  ENDIF                                               	    		
               
             		IF (Ri.GT.0.) THEN
					Rbcstar = Rbc/(1-5*Ri)**2      ! Bulk
				ELSE
      				Rbcstar = Rbc/((1-5*Ri)**0.75)
 		    	ENDIF

             RKINbc  = 1.0/Rbcstar                     ! Bulk conductance
	  ELSE


      CALL WINDTRANab(V,Zm,param,sitev, &
                             Vstar,Vh,d,z0c, Vz,Vc )               ! Output wind speed within canopy at height 2 m           
      
!      Aerodynamic Resistances for Canopy
          ndecay =   Wcoeff*LAI              !(Wcoeff=0.5 for TWDEF, 0.6 for Niwot) ! wind speed decay coefficient inside canopy 

			Kh = k**2.0*V*(Hcan-d)/(log((Zm-d)/Z0c))
           
!        Above canopy aerodynamic resistance 
			Ra = (1.0* 1.0/(k**2.0*V)*log((Zm-d)/(z0c))*log((Zm-d)/  &        ! Increased ten times (Calder 1990, Lundberg and Halldin 1994)
     		 (Hcan-d))+ Hcan/(Kh*ndecay)* &
                (exp(ndecay-ndecay*(d+z0c)/Hcan)-1) )*1/hff

!        Below canopy aerodynamic resistance 
			Rc =(Hcan*exp(ndecay)/(Kh*ndecay)*(exp(-ndecay*z/Hcan) &
     		      -exp(-ndecay*(d+z0c)/Hcan))) *1/hff

			Rs = (1.0/(k**2*Vz)*(log(z/zo))**2.0)*1/hff	               
              Rc = Rc+ Rs                                                     ! Total below canopy resistance

!        Bulk aerodynamic resistance (used when there is no snow in the canopy)
		   Rbc = Rc+Ra                                               ! Total aerodynamic resistance from ground to above canopy used when there is no snow on canopy
                 
!     Boundary layer resistance of canopy(Rl)
        Dc      = .04    	                                             ! Deckinson et.al (1986);Bonan(1991): Characterstics dimension (m) of leaves (average width) ! put in parameters later
        Lbmean  = 0.02/ndecay*sqrt(Vc/Dc)* (1-exp(-ndecay/2))               
	  Rl      = 1/(LAI*Lbmean)*1/hff 

!        Leaf boundary layer conductance: reciprocal of resistance
	     RKINl  = 1.0/Rl  
	    
!      Correction for Ra, Rc and Rac for stable and unstable condition

	          Ri = g*(Ta-Tss)*Z /(Vz**2.0*(0.5*(Ta+Tss)+273.15))    ! wind speed unit should not be changed.
                                                    
	    		IF (Ri.GT.Rimax) THEN                                 
				Ri= Rimax
	    		ENDIF
                 
             		IF (Ri.GT.0.) THEN
					Rastar  = Ra  !/((1-5*Ri)**2)       ! Above canopy
					Rcstar  = Rc /((1-5*Ri)**2)      ! Below canpy
					Rbcstar = Rbc !  /((1-5*Ri)**2)     ! Bulk
                                                                      ! Rlstar  = Rl/((1-5*Ri)**2)         ! leaf boundary layer
				ELSE
					Rastar  = Ra  !/(1-5*Ri)**0.75 
					Rcstar  = Rc /(1-5*Ri)**0.75              
      				Rbcstar = Rbc  ! /(1-5*Ri)**0.75
				                                                    !	Rlstar  = Rl/(1-5*Ri)**0.75 
 				ENDIF
  

!     Turbulent (heat and vapor) transfer coefficient (Conductance) ! 
               RKINa   = 1.0/Rastar              ! Above canopy
	 		   RKINc   = 1.0/Rcstar              ! Below canpy                  
               RKINbc  = 1.0/Rbcstar             ! Bulk
!	         RKINl   = 1.0/Rlstar              ! leaf boundary layer
   
      ENDIF 
	ENDIF           
	RETURN
	END

!********************* WIND PROFILES ****************************************************************
!     Calculates wind at the 2 m above the surface and at the sink using the input of measured wind at 2 m above the canopy

	SUBROUTINE WINDTRANab(V,Zm,param,sitev, &
                             Vstar,Vh,d,z0c,Vz,Vc )               ! Out put wind speed within canopy at height 2 m

	REAL k,V,Vz,Hcan,LAI,Beta,Tac,RHc,Z,Zm,Z0c,d,Vstar,Vh
	REAL sitev(*),param(*) 

      data k  /0.4/           !  Von Karmans constant (0.4)
      Z      = param(5)       !  Nominal meas. height for air temp. and humidity (2m)
      Wcoeff= param(32)       ! 

	Cc     = sitev(5)       ! Canopy Cover 
	Hcan   = sitev(6)       ! Canopy height  
 	LAI    = sitev(7)       ! Leaf Area Index 
	Ycage  = sitev(9)      ! Parameters for wind speed transformation

      LAI = Cc*LAI
      
!     Roughness length (z0c) and zero plane displacement(d) for canopy [Periera, 1982]

      IF(V.LE.0.)THEN
	   Vz= 0.
	   Vc= 0.
	ELSE

	 d = Hcan*(0.05 + (LAI**0.20)/2 + (YcAge-1)/20)
                            

	IF (LAI.LT.1) THEN
		z0c= 0.1*Hcan
	ELSE
		z0c = Hcan*(0.23 - (LAI**0.25)/10 - (YcAge-1)/67)
	END IF

	Vstar  = k*V/log((Zm-d)/z0c)       !Friction velocity
        
!	Wind speed at height,Z from the ground/snow surface below canopy (Exponential profile from canopy height to below canopy)

      Vh  = 1./k*Vstar*log((Hcan-d)/z0c)  !Wind speed at the height of canopy (Logarithmic profile above canopy to canopy)
	  
	   Vz   = 1.0* Vh*exp(-Wcoeff*LAI*(1-z/Hcan))
	   Vc   = 1.0* Vh*exp(-Wcoeff*LAI*(1-(d+Z0c)/Hcan))     !To estimate canopy boundary layer conductance
!	   Vc   = 1.0* Vh*exp(-Wcoeff*LAI*(1-6./Hcan)**Beta)       

	ENDIF
	                          
	RETURN
	END

!************************ TURBULENT FLUXES (ABOVE AND BELOW THE CANOPY *******************************

!     Calculates the turbulent heat fluxes (sensible and latent
!     heat fluxes) and condensation/sublimation.

      SUBROUTINE TURBFLUX(Ws,Wc,A,TK,Tc,Ta,Tss,RH,V,EA,p,param,sitev, &
     	              d,Z0c,Vz,RKINc,RKINbc,Tac,Fs,Ess,Esc, &   ! Output variables
                      QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E)          
                                                      
      REAL Wc,RHOA,ESc,k,Apr, param(*),sitev(*),A,EA, LAI, EHo, EEo
      REAL QHc,QEc,Ec,QHs,QEs,Es,QH,QE,E
    

      z     =   param(5)      !  Nominal meas. height for air temp. and humidity (2m)
      zo    =   param(6)      !  Surface aerodynamic roughness (m)
      fstab =   param(19)     !  Stability correction control parameter 0 = no corrections, 1 = full corrections

      APr   =   sitev(2)      !  Atmospheric Pressure (Pa)
      Cc    =   sitev(5)      ! Canopy Coverage
	  LAI   =   sitev(7)      ! Leaf Area Index 

	  data Rag /287.0/        !  Ideal Gas constant for dry air (287 J/kg/K) (name changed)
      data k  /0.4/           !  Von Karmans constant (0.4)
      data hff /3600.0/       !  Factor to convert /s into /hr (3600)
      data rhow /1000.0/      !  Density of Water (1000 kg/m^3)
      data g    /9.81/        !  Gravitational acceleration (9.81 m/s^2)
      data cp /1.005/         !  Air Heat Capacity (1.005 KJ/kg/K)
      data hf /333.5/         !  Heat of fusion (333.5 KJ/kg)
      data hneu /2834.0/      !  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)

      data tol,nitermax/0.001,20/   
 
      Tak   = Ta+Tk
      Tck   = Tc+Tk

	  RHOA   = APr/(RAg*(TAK))	  ! RHOA in kg/m3,  APr: Atm Press
	  Ea     = SVPW(Ta)*RH          ! Actual vapor pressure sorrounding canopy

!     Wind less coefficient:
     EHoA  = 0.0       ! for sensibleheat flux
	 EHoC  = 0.0
	 EHoS  = 0.0
	 EEoA  = 0.0       ! for latent heat flux
     EEoC  = 0.0
     EEoS  = 0.0

 	CALL AeroRes(p,Wc,V,Ta,Tss,Tc,Fs,param,sitev,Tk, &
               d,Z0c,Vz,Rc,Ra,Rbc,Rl,RKINc,RKINa,Rkinbc,Rkinl)         ! Output variables
       

	IF (V.LE.0.) THEN
		   QHc  = 0.0
		   QEc  = 0.
		   Ec   = 0.
		   QHs  = 0.0
		   QEs  = 0.
		   Es   = 0.
             Tac  = Ta        ! Though We don't need Tac when V=0.
		
	ELSE

		IF (LAI.EQ.0.)THEN
						        
!     Flux estimation from the open areas when no vanopy is present
			QHs = (RHOA*cp*RKINbc+EHos)* (TA-TSs)                        !	0.0 added for zero wind condition  
			QEs = (0.622*hneu/(RAg*(TAK))*RKINbc+EEos)*(EA-ESs) 
			Es  = -QEs/(RHOW*HNEU)                                       ! Es in  m/hr
			QHc = 0.
			QEc = 0.
			Ec  = 0.

		ELSE 

		   Tac = (Tc*Rkinl + Tss*Rkinc + Ta*Rkina) &
     		     	/(1*Rkina + 1*Rkinl + 1*Rkinc)
		   Eac = (ESc*Rkinl + ESs*Rkinc + Ea*Rkina) &
     		     	/(1*Rkina + 1*Rkinl + 1*Rkinc)
             Tack    = Tac+Tk   
    	       RHOAc   =  APr/(RAg*(TAcK))	              
                 
!	 Flux from canopy
		   QHc  = (RHOAc*cp*RKINl+EHoC)*(Tac-Tc) 
			       
             IF (Wc.EQ.0.and. p.EQ.0.)THEN 
                      QEc = 0.                          ! No latent heat flux when there is no snow in the caanopy
	                Ec  = 0.
		   ELSE	
               QEc  = (0.622*hneu/(RAg*(TAcK))*RKINl+ EEoC)*(EAc-ESc)*Fs
               Ec   = -QEc/(RHOW*HNEU)   
             ENDIF 			
!     Flux from snow at ground
			QHs  =(RHOAc*cp*RKINc +EHoS)* (Tac-TSs)                  !	EHoS added for zero wind condition [Jordan et al.(1999); Koivusalo (2002); Herrero et al.(2009) ] 
			QEs  =(0.622*hneu/(RAg*(TAcK))*RKINc+ EEoS)*(EAc-ESs) 
			Es   =- QEs/(RHOW*HNEU)
	ENDIF
!	Total flux
			QH  =  QHs + QHc 
			QE  =  QEs + QEc 
			 E  =  Es  + Ec  
	ENDIF             
      RETURN
      END

!************************** RATE OF INTERCEPTION AND MAX INTERCEPTION ***********************
!     Calculates amount of snow intercepted by the canopy

	SUBROUTINE INTERCEPT(Ta,LAI,p,Wc,dt,Inmax,param,sitev, &
     	                            ieff,Ur,int)    ! Output variables

	REAL p,LAI,CC,Inmax,ieff, Rhofs, Sbar,S,int,dt,Ur 
	REAL Intma 
	REAL param(*),sitev(*)
	DOUBLE PRECISION Uc

    	Uc    = param(27)		    ! Unloading rate coefficient (Per hour) (Hedstrom and pomeroy, 1998) 
	Cc    = sitev(5)            ! Canopy Coverage

	IF(LAI.EQ.0.) 	CC=0.

!	Interception rate (di/dt)
      IF(p.EQ.0.) THEN
	    int  = 0.0
          ieff = 0.0
      ELSE
          Intma=0.
       IF (Wc.LT.Inmax) THEN
	    int  = CC*(1.-Wc*CC/(Inmax))*p
	    ieff = CC*(1.-Wc*CC/(Inmax))          ! Interception efficiency (di/dp) 
	    Intma =Inmax-Wc
		IF(int.GE.Intma) int= Intma

	  ELSE
	    int  = 0.0
          ieff = 0.0
	  END IF	

      END IF	

!	Unloading rate
	
	Ur = Uc*(Wc+(int*dt)/2)       ! Melt and Evap. are subtracted 
                                          ! half of the current interception is also considered for unloading
	RETURN	
	END

!******************************   ENERGY ADVECTED BY RAIN   **********************************


!     Calculates the heat advected to the snowpack due to rain

      FUNCTION QPF(PR,TA,TO,PS,RHOW,HF,CW,CS)
      IF(TA.GT.TO) THEN
         TRAIN=TA
         TSNOW=TO													! TO : Freezing Temp (int.e.0C)
      ELSE
         TRAIN=TO
         TSNOW=TA
      ENDIF
      QPF=PR*RHOW*(HF+CW*(TRAIN-TO))+PS*RHOW*CS*(TSNOW-TO)
      RETURN
      END


!****************************** PREHELP () ***************************************************
!      Routine to correct energy and mass fluxes when 
!      numerical overshoots dictate that W was changed in 
!      the calling routine - either because W went negative
!      or due to the liquid fraction being held constant.

       SUBROUTINE PREHELP(W1,W,DT,FM,FM1,fac,PS,PRAIN,E,RHOW,HF,Q,QM,MR, &
            qe,hneu)
       REAL MR,DT

       FM = (W1-W)/DT*fac-FM1
 
       MR = MAX( 0.0 , (PS + PRAIN - FM - E))
       E = Ps + Prain - FM - MR
!    Because melt rate changes the advected energy also changes.  Here
!     advected and melt energy are separated,
       QOTHER = Q + QM - QE
!     then the part due to melt recalculated
       QM = MR*RHOW*HF
!     then the part due to evaporation recalculated
       QE = -E*RHOW*HNEU
!     Energy recombined
       Q = QOTHER - QM + QE
       RETURN
       END

!****************** EXPONENTIAL INTEGRAL FUNCTION *****************************************
!     Computes the exponential integral function for the given value      
	FUNCTION EXPINT (LAI)
	REAL LAI
	DOUBLE PRECISION a0,a1,a2,a3,a4,a5,b1,b2,b3,b4
	IF (LAI.EQ.0)THEN
	    EXPINT=1.                          

	ELSEIF (LAI.LE.1.0) THEN
		a0=-.57721566
		a1=.99999193
		a2=-.24991055
		a3=.05519968
		a4=-.00976004
		a5=.00107857

		EXPINT = a0+a1*LAI+a2*LAI**2+a3*LAI**3+a4*LAI**4+a5*LAI**5  &
     		 -log(LAI)

	ELSE
		a1=8.5733287401
		a2=18.0590169730
		a3=8.6347637343
		a4=.2677737343
		b1=9.5733223454
		b2=25.6329561486
		b3=21.0996530827
		b4=3.9584969228

		EXPINT=(LAI**4+a1*LAI**3+a2*LAI**2+a3*LAI+a4)/ &
          	((LAI**4+b1*LAI**3+b2*LAI**2+b3*LAI+b4)*LAI*exp(LAI))

	END IF
	RETURN
      END

!      Reference Book:Handbook of mathematical functions with Formulas, Graphs, and Mathematical Tables
!      Edited by Milton Abramowitz and Irene A. Stegun, Volume 55, Issue 1972, page no. 231
!      Dover Publications, Inc, Mineola, New York.

!*************************** SVP () ***********************************************************
!     Calculates the vapour pressure at a specified temperature over water or ice
!     depending upon temperature.  Temperature is celsius here.

      FUNCTION SVP(T)
      IF(t .ge. 0.)THEN
        svp=svpw(t)
      ELSE
        svp=svpi(t)

      ENDIF
	RETURN
	END



!*************************** SVPW () **********************************************************
!     Calculates the vapour pressure at a specified temperature over water
!     using polynomial from Lowe (1977).

      FUNCTION SVPW(T)

      SVPW=6.107799961 + t * (0.4436518521 + t * (0.01428945805 +  &
           t * (0.0002650648471 + t * (3.031240936e-06 + t *  &
           (2.034080948e-08 + t * 6.136820929e-11)))))
      SVPW=SVPW*100											! convert from mb to Pa
      RETURN
      END



!*************************** SVPI () *********************************************************
!     Calculates the vapour pressure at a specified temperature over ice.
!     using polynomial from Lowe (1977).

      FUNCTION SVPI(T)
      SVPI=6.109177956 + t * (0.503469897 + t * (0.01886013408 +  &
           t * (0.0004176223716 + t * (5.82472028e-06 + t *  &
           (4.838803174e-08 + t * 1.838826904e-10)))))
      SVPI=SVPI*100											! convert from mb to Pa
      RETURN
      END


!**********************************************************************************************
!     Estimates reflection and scattering coefficient

       FUNCTION Tau1(Rho,G,h,COSZEN,kk)
       REAL Rho,G,h,COSZEN,kk, Tau1

		IF (h.EQ.0.) THEN
		  Tau1=1.0         ! To avoid computational error when coszen =0
		ELSE
	      IF (COSZEN.EQ. 0.) THEN
	       Tau1 = 0
            ELSE
	      Tau1  = EXP(- kk*G*Rho*h/COSZEN) 
		  ENDIF                  ! is the penetration function for deep canopy: direct
          ENDIF

       RETURN
       END
!************************************************************************************************
       FUNCTION Tau2(Rho,G,h,kk,EXPI)
       REAL Rho,G,h,kk,EXPI,Tau2

 	     Tau2  = (1-kk*G*Rho*h)*EXP(-kk*G*Rho*h)+(kk*G*Rho*h)**2*EXPI     ! is the penetration function for deep canopy: diffuse

       RETURN
       END
!************************************************************************************************