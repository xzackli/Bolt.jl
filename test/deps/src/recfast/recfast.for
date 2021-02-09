cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C Integrator for Cosmic Recombination of Hydrogen and Helium,
C developed by Douglas Scott (dscott@astro.ubc.ca)
C based on calculations in the papers Seager, Sasselov & Scott
C (ApJ, 523, L1, 1999; ApJS, 128, 407, 2000)
C and "fudge" updates in Wong, Moss & Scott (2008).
C
C Permission to use, copy, modify and distribute without fee or royalty at
C any tier, this software and its documentation, for any purpose and without
C fee or royalty is hereby granted, provided that you agree to comply with
C the following copyright notice and statements, including the disclaimer,
C and that the same appear on ALL copies of the software and documentation,
C including modifications that you make for internal use or for distribution:
C
C Copyright 1999-2010 by University of British Columbia.  All rights reserved.
C
C THIS SOFTWARE IS PROVIDED "AS IS", AND U.B.C. MAKES NO 
C REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  
C BY WAY OF EXAMPLE, BUT NOT LIMITATION,
C U.B.C. MAKES NO REPRESENTATIONS OR WARRANTIES OF 
C MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT 
C THE USE OF THE LICENSED SOFTWARE OR DOCUMENTATION WILL NOT INFRINGE 
C ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.   
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CN	Name:	RECFAST
CV	Version: 1.5.2
C 
CP	Purpose:  Calculate ionised fraction as a function of redshift.
CP		  Solves for H and He simultaneously, and includes
CP		  H "fudge factor" for low z effect, as well as
CP	          HeI fudge factor.
C
CD	Description: Solves for ionisation history since recombination
CD	using the equations in Seager, Sasselov & Scott (ApJ, 1999).
CD	The Cosmological model can be flat or open.
CD	The matter temperature is also followed, with an update from
CD	Moss & Scott (2009).
CD	The values for \alpha_B for H are from Hummer (1994).
CD	The singlet HeI coefficient is a fit from the full code.
CD	Additional He "fudge factors" are as described in Wong, Moss
CD	and Scott (2008).
CD	Extra fitting function included (in optical depth) to account
CD	for extra H physics described in Rubino-Martin et al. (2010).
CD	Care is taken to use the most accurate constants.
CD	Note that some parameters are fixed (e.g. N_nu=3, nu's are
CD	massless, w=-1, etc.) - some users may want to explictly
CD	imput their own H(z) to account for extra physics.
CD	This is provided as a PROGRAM, which can be easily converted
CD	to a SUBROUTINE for use in CMB Boltzmann codes.
C		
CA	Arguments:
CA	Name, Description
CA	Double precision throughout
CA
CA	z is redshift - W is sqrt(1+z), like conformal time
CA	x is total ionised fraction, relative to H
CA	x_H is ionized fraction of H - y(1) in R-K routine
CA	x_He is ionized fraction of He - y(2) in R-K routine
CA	(note that x_He=n_He+/n_He here and not n_He+/n_H)
CA	Tmat is matter temperature - y(3) in R-K routine
CA	f's are the derivatives of the Y's
CA	alphaB is case B recombination rate
CA	alpHe is the singlet only HeII recombination rate
CA	a_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	b_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	c_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	d_PPB is Pequignot, Petitjean & Boisson fitting parameter for Hydrogen
CA	a_VF is Verner and Ferland type fitting parameter for Helium
CA	b_VF is Verner and Ferland type fitting parameter for Helium
CA	T_0 is Verner and Ferland type fitting parameter for Helium
CA	T_1 is Verner and Ferland type fitting parameter for Helium
CA	Tnow is the observed CMB temperature today
CA	OmegaT is the total Omega_0
CA      OmegaL is the Omega_0 contribution from a Cosmological constant
CA      OmegaK is the Omega_0 contribution in curvature (1-O_T-O_L)
CA      OmegaB is Omega in baryons today
CA	OmegaC is the Omega_0 in (cold) dark matter: OmegaT=OmegaC+OmegaB
CA	Yp is the primordial helium abundace
CA	fHe is He/H number ratio = Yp/4(1-Yp)
CA	Trad and Tmat are radiation and matter temperatures
CA	epsilon is the approximate difference (=Trad-Tmat) at high z
CA	OmegaB is Omega in baryons today
CA	H is Hubble constant in units of 100 km/s/Mpc
CA	HOinp is input value of Hubble constant in units of 100 km/s/Mpc
CA	HO is Hubble constant in SI units
CA	bigH is 100 km/s/Mpc in SI units
CA	Hz is the value of H at the specific z (in ION)
CA	G is grvitational constant
CA	n is number density of hydrogen
CA	Nnow is number density today
CA	x0 is initial ionized fraction
CA	x_H0 is initial ionized fraction of Hydrogen
CA	x_He0 is initial ionized fraction of Helium
CA	rhs is dummy for calculating x0
CA	zinitial and zfinal are starting and ending redshifts
CA	fnu is the contribution of neutrinos to the radn. energy density
CA	zeq is the redshift of matter-radiation equality
CA	zstart and zend are for each pass to the integrator
CA	w0 and w1 are conformal-time-like initial and final zi and zf's
CA	Lw0 and Lw1 are logs of w0 and w1
CA	hw is the interval in W
CA	C,k_B,h_P: speed of light, Boltzmann's and Planck's constants
CA	m_e,m_H: electron mass and H atomic mass in SI
CA	not4: ratio of 4He atomic mass to 1H atomic mass
CA	sigma: Thomson cross-section
CA	a: radiation constant for u=aT^4
CA	Pi: Pi
CA	Lambda: 2s-1s two photon rate for Hydrogen
CA	Lambda_He: 2s-1s two photon rate for Helium
CA	DeltaB: energy of first excited state from continuum = 3.4eV
CA	DeltaB_He: energy of first excited state from cont. for He = 3.4eV
CA	L_H_ion: level for H ionization in m^-1
CA	L_H_alpha: level for H Ly alpha in m^-1
CA	L_He1_ion: level for HeI ionization
CA	L_He2_ion: level for HeII ionization
CA	L_He_2s: level for HeI 2s
CA	L_He_2p: level for He 2p (21P1-11S0) in m^-1
CA	Lalpha: Ly alpha wavelength in SI
CA	Lalpha_He: Helium I 2p-1s wavelength in SI
CA	mu_H,mu_T: mass per H atom and mass per particle
CA	H_frac: follow Tmat when t_Compton / t_Hubble > H_frac
CA	dHdz is the derivative of H at the specific z (in ION)
CA	CDB=DeltaB/k_B			Constants derived from B1,B2,R
CA	CDB_He=DeltaB_He/k_B		n=2-infinity for He in Kelvin
CA	CB1=CDB*4.			Lalpha and sigma_Th, calculated
CA	CB1_He1: CB1 for HeI ionization potential
CA	CB1_He2: CB1 for HeII ionization potential
CA	CR=2*Pi*(m_e/h_P)*(k_B/h_P)	once and passed in a common block
CA	CK=Lalpha**3/(8.*Pi)
CA	CK_He=Lalpha_He**3/(8.*Pi)
CA	CL=C*h_P/(k_B*Lalpha)
CA	CL_He=C*h_P/(k_B*Lalpha_He)
CA	CT=(8./3.)*(sigma/(m_e*C))*a
CA	Bfact=exp((E_2p-E_2s)/kT)	Extra Boltzmann factor
CA	fu is a "fudge factor" for H, to approximate low z behaviour
CA	b_He is a "fudge factor" for HeI, to approximate higher z behaviour
CA	Heswitch is an integer for modifying HeI recombination
CA	Parameters and quantities to describe the extra triplet states
CA	 and also the continuum opacity of H, with a fitting function
CA	 suggested by KIV, astro-ph/0703438
CA	a_trip: used to fit HeI triplet recombination rate
CA	b_trip: used to fit HeI triplet recombination rate
CA	L_He_2Pt: level for 23P012-11S0 in m^-1
CA	L_He_2St: level for 23S1-11S0 in m^-1
CA	L_He2St_ion: level for 23S1-continuum in m^-1
CA	A2P_s: Einstein A coefficient for He 21P1-11S0
CA	A2P_t: Einstein A coefficient for He 23P1-11S0    
CA	sigma_He_2Ps: H ionization x-section at HeI 21P1-11S0 freq. in m^2
CA	sigma_He_2Pt: H ionization x-section at HeI 23P1-11S0 freq. in m^2
CA	CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
CA	CfHe_t: triplet statistical correction
CA	Hswitch is an integer for modifying the H recombination
CA	AGauss1 is the amplitude of the 1st Gaussian for the H fudging
CA	AGauss2 is the amplitude of the 2nd Gaussian for the H fudging
CA	zGauss1 is the ln(1+z) central value of the 1st Gaussian
CA	zGauss2 is the ln(1+z) central value of the 2nd Gaussian
CA	wGauss1 is the width of the 1st Gaussian
CA	wGauss2 is the width of the 2nd Gaussian
CA	tol: tolerance for the integrator
CA	cw(24),w(3,9): work space for DVERK
CA	Ndim: number of d.e.'s to solve (integer)
CA	Nz: number of output redshitf (integer)
CA	I: loop index (integer)
CA	ind,nw: work-space for DVERK (integer)
C
CG	Global data (common blocks) referenced:
CG	/zLIST/zinitial,zfinal,Nz
CG	/Cfund/C,k_B,h_P,m_e,m_H,not4,sigma,a,Pi
CG	/data/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
CG		fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He
CG      /Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
CG	/Hemod/b_He,A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt,
CG		L_He_2p,L_He_2Pt,L_He_2St,L_He2St_ion
CG	/Hmod/AGauss1,AGauss2,zGauss1,zGauss2,wGauss1,wGauss2
CG	/Switch/Heswitch,Hswitch
C
CF	File & device access:
CF	Unit	/I,IO,O	/Name (if known)
C
CM	Modules called:
CM	DVERK (numerical integrator)
CM	GET_INIT (initial values for ionization fractions)
CM	ION (ionization and Temp derivatives)
C
CC	Comments:
CC	none
C
CH	History:
CH	CREATED		(simplest version) 19th March 1989
CH	RECREATED	11th January 1995
CH			includes variable Cosmology
CH			uses DVERK integrator
CH			initial conditions are Saha
CH	TESTED		a bunch, well, OK, not really
CH	MODIFIED	January 1995 (include Hummer's 1994 alpha table)
CH			January 1995 (include new value for 2s-1s rate)
CH			January 1995 (expand comments)
CH			March 1995 (add Saha for Helium)
CH			August 1997 (add HeII alpha table)
CH			July 1998 (include OmegaT correction and H fudge factor)
CH			Nov 1998 (change Trad to Tmat in Rup)
CH			Jan 1999 (tidied up for public consumption)
CH			Sept 1999 (switch to formula for alpha's, fix glitch)
CH			Feb 2000 (fixed overflow problem in He_Boltz)
CH			Oct 2001 (fixed OmegaT in z_eq)
CH			June 2003 (fixed error in Rdown_He formula)
CH			June 2003 (fixed He recombination coefficients)
CH			June 2003 (comments to point out fixed N_nu etc.)
CH			Oct 2006 (included new value for G)
CH			Oct 2006 (improved m_He/m_H to be "not4")
CH			Oct 2006 (fixed error, x for x_H in part of f(1))
CH			Jan 2008 (improved HeI recombination effects,
CH			              including HeI rec. fudge factor)
CH			Feb 2008 (avoid calculating gamma_2Ps and
CH			               gamma_2Pt when x_H close to 1.0)
CH			Aug 2008 (correction for x_H when Heflag=2
CH				     and Helfag>=5 to be smoother)
CH			Sept 2008 (added extra term to make transition
CH				     smoother for Tmat evolution)
CH			Jan 2010 (added fitting function to modify K
CH				to match x_e(z) for new H physics)
CH			July 2012 (modified fudge factors for better
CH				match to codes with more detailed physics)
CH			Sept 2012 (fixed "fu" at low z to match modifications)
C-
C	===============================================================

	subroutine get_xe(in_OmegaB, in_OmegaC, in_OmegaL, in_HOinp, in_Tnow, in_Yp, in_Hswitch, in_Heswitch, in_Nz, in_zinitial, in_zfinal, out_xe)

	real*8 in_OmegaB, in_OmegaC, in_OmegaL, in_HOinp, in_Tnow, in_Yp, in_zinitial, in_zfinal
	integer in_Hswitch, in_Heswitch, in_Nz
	real*8 out_xe(in_Nz)
	
C	--- Arguments
	real*8 Trad,Tmat
	real*8 OmegaT,OmegaB,H,HO,HOinp,bigH,G,OmegaL,OmegaK,OmegaC
	real*8 z,n,x,x0,rhs,x_H,x_He,x_H0,x_He0
	real*8 Tnow,zinitial,zfinal,Nnow,z_eq,fnu
	real*8 zstart,zend,w0,w1,Lw0,Lw1,hw
	real*8 C,k_B,h_P,m_e,m_H,not4,sigma,a,Pi
	real*8 Lambda,DeltaB,DeltaB_He,Lalpha,mu_H,mu_T,H_frac
	real*8 Lambda_He,Lalpha_He,Bfact,CK_He,CL_He
	real*8 L_H_ion,L_H_alpha,L_He1_ion,L_He2_ion,L_He_2s,L_He_2p
	real*8 CB1,CDB,CR,CK,CL,CT,Yp,fHe,CB1_He1,CB1_He2,CDB_He,fu,b_He
	real*8 A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt
	real*8 L_He_2Pt,L_He_2St,L_He2St_ion
	real*8 AGauss1,AGauss2,zGauss1,zGauss2,wGauss1,wGauss2

	real*8 tol
	real*8 cw(24),w(3,9)
	real*8 y(3)

	integer Ndim,Nz,I
	integer ind,nw
	integer Heswitch,Hswitch

	character*80 fileout

C	--- Parameter statements
	parameter(bigH=100.0D3/(1.0D6*3.0856775807D16))	!Ho in s-1
	parameter(tol=1.D-5)				!Tolerance for R-K

	external ION

C	--- Commons
	common/zLIST/zinitial,zfinal,Nz
	common/Cfund/C,k_B,h_P,m_e,m_H,not4,sigma,a,Pi
	common/Cdata/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
	1	fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He,fu
	common/Hemod/b_He,A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt,
	1	L_He_2p,L_He_2Pt,L_He_2St,L_He2St_ion
	common/Hmod/AGauss1,AGauss2,zGauss1,zGauss2,wGauss1,wGauss2
	common/Switch/Heswitch,Hswitch

	common/Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
C	===============================================================

C	--- Data
	data	C,k_B,h_P	/2.99792458D8,1.380658D-23,6.6260755D-34/
	data	m_e,m_H		/9.1093897D-31,1.673575D-27/	!av. H atom
c	note: neglecting deuterium, making an O(e-5) effect
	data	not4		/3.9715D0/		!mass He/H atom
	data	sigma,a		/6.6524616D-29,7.565914D-16/
	data	Pi		/3.141592653589d0/
	data	G		/6.6742D-11/ 			!new value
C	Fundamental constants in SI units
C	("not4" pointed out by Gary Steigman)

	data	Lambda		/8.2245809d0/
	data	Lambda_He	/51.3d0/	!new value from Dalgarno
	data	L_H_ion		/1.096787737D7/	!level for H ion. (in m^-1)
	data	L_H_alpha	/8.225916453D6/ !averaged over 2 levels
	data	L_He1_ion	/1.98310772D7/	!from Drake (1993)
	data	L_He2_ion	/4.389088863D7/	!from JPhysChemRefData (1987)
	data	L_He_2s		/1.66277434D7/	!from Drake (1993)
	data	L_He_2p		/1.71134891D7/	!from Drake (1993)
C	2 photon rates and atomic levels in SI units

	data	A2P_s		/1.798287D9/    !Morton, Wu & Drake (2006)
	data	A2P_t		/177.58D0/      !Lach & Pachuski (2001)
	data	L_He_2Pt	/1.690871466D7/ !Drake & Morton (2007)
	data	L_He_2St	/1.5985597526D7/ !Drake & Morton (2007)
	data	L_He2St_ion	/3.8454693845D6/ !Drake & Morton (2007)
	data	sigma_He_2Ps	/1.436289D-22/  !Hummer & Storey (1998)
	data	sigma_He_2Pt	/1.484872D-22/  !Hummer & Storey (1998)
C	Atomic data for HeI 

	data	AGauss1		/-0.14D0/	!Amplitude of 1st Gaussian
	data	AGauss2		/0.079D0/	!Amplitude of 2nd Gaussian
	data	zGauss1		/7.28D0/	!ln(1+z) of 1st Gaussian
	data	zGauss2		/6.73D0/	!ln(1+z) of 2nd Gaussian
	data	wGauss1		/0.18D0/	!Width of 1st Gaussian
	data	wGauss2		/0.33D0/	!Width of 2nd Gaussian
C	Gaussian fits for extra H physics (fit by Adam Moss, modified by
C	Antony Lewis)

C   Copy function argument parameters to common block data
	OmegaB = in_OmegaB
	OmegaC = in_OmegaC
	OmegaL = in_OmegaL
	HOinp = in_HOinp
	Tnow = in_Tnow
	Yp = in_Yp
	Hswitch = in_Hswitch
	Heswitch = in_Heswitch
	Nz = in_Nz
	zinitial = in_zinitial
	zfinal = in_zfinal

C	dimensions for integrator
	Ndim = 3


	z = zinitial

	OmegaT=OmegaC+OmegaB            !total dark matter + baryons
	OmegaK=1.d0-OmegaT-OmegaL	!curvature

c	convert the Hubble constant units
	H = HOinp/100.d0
	HO = H*bigH

c	sort out the helium abundance parameters
	mu_H = 1.d0/(1.d0-Yp)			!Mass per H atom
	mu_T = not4/(not4-(not4-1.d0)*Yp)	!Mass per atom
	fHe = Yp/(not4*(1.d0-Yp))		!n_He_tot / n_H_tot

	Nnow = 3.d0*HO*HO*OmegaB/(8.d0*Pi*G*mu_H*m_H)
	n = Nnow * (1.d0+z)**3
	fnu = (21.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
c	(this is explictly for 3 massless neutrinos - change if N_nu.ne.3)
	z_eq = (3.d0*(HO*C)**2/(8.d0*Pi*G*a*(1.d0+fnu)*Tnow**4))*OmegaT
	z_eq = z_eq - 1.d0

C	Set up some constants so they don't have to be calculated later
	Lalpha = 1.d0/L_H_alpha
	Lalpha_He = 1.d0/L_He_2p
	DeltaB = h_P*C*(L_H_ion-L_H_alpha)
	CDB = DeltaB/k_B
	DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)	!2s, not 2p
	CDB_He = DeltaB_He/k_B
	CB1 = h_P*C*L_H_ion/k_B
	CB1_He1 = h_P*C*L_He1_ion/k_B	!ionization for HeI
	CB1_He2 = h_P*C*L_He2_ion/k_B	!ionization for HeII
	CR = 2.d0*Pi*(m_e/h_P)*(k_B/h_P)
	CK = Lalpha**3/(8.d0*Pi)
	CK_He = Lalpha_He**3/(8.d0*Pi)
	CL = C*h_P/(k_B*Lalpha)
	CL_He = C*h_P/(k_B/L_He_2s)	!comes from det.bal. of 2s-1s
	CT = (8.d0/3.d0)*(sigma/(m_e*C))*a
	Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

C	Matter departs from radiation when t(Th) > H_frac * t(H)
C	choose some safely small number
	H_frac = 1.D-3


C	Fudge factor to approximate the low z out of equilibrium effect
	if (Hswitch .eq. 0) then
	  fu=1.14d0
	else
	  fu=1.125d0
	end if

	 
c	Set the He fudge factor
cc	if((Heswitch.eq. 2).or.(Heswitch.eq. 5).or.(Heswitch.eq.6))then
cc	  write(*,*)'Enter the fudge factor b_He'
cc	  read(*,*)b_He
cc	endif
	b_He = 0.86

c	Set initial matter temperature
	y(3) = Tnow*(1.d0+z)            !Initial rad. & mat. temperature
	Tmat = y(3)

	call get_init(z,x_H0,x_He0,x0)

	y(1) = x_H0
	y(2) = x_He0


	w0=1.d0/ dsqrt(1.d0 + zinitial)	!like a conformal time
	w1=1.d0/ dsqrt(1.d0 + zfinal)
	Lw0 = dLog(w0)
	Lw1 = dLog(w1)
	hW=(Lw1-Lw0)/dfloat(Nz)		!interval in log of conf time

c	Set up work-space stuff for DVERK
	ind  = 1
	nw   = 3
	do i = 1,24
	  cw(i) = 0.d0
	end do

	do i = 1,Nz
C       calculate the start and end redshift for the interval at each z
C	or just at each z
	  zstart = zinitial + dfloat(i-1)*(zfinal-zinitial)/dfloat(Nz)
	  zend   = zinitial + dfloat(i)*(zfinal-zinitial)/dfloat(Nz)

C Use Saha to get x_e, using the equation for x_e for ionized helium
C and for neutral helium.
C Everyb_trip ionized above z=8000.  First ionization over by z=5000.
C Assume He all singly ionized down to z=3500, then use He Saha until
C He is 99% singly ionized, and *then* switch to joint H/He recombination.

	  z = zend

	  if (zend.gt.8000.d0) then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    x0 = 1.d0+2.d0*fHe
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

	  else if(z.gt.5000.d0)then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He2/(Tnow*(1.d0+z)) ) / Nnow
	    rhs = rhs*1.d0		!ratio of g's is 1 for He++ <-> He+
	    x0 = 0.5d0 * ( dsqrt( (rhs-1.d0-fHe)**2
	1	+ 4.d0*(1.d0+2.d0*fHe)*rhs) - (rhs-1.d0-fHe) )
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

	  else if(z.gt.3500.d0)then

	    x_H0 = 1.d0
            x_He0 = 1.d0
            x0 = x_H0 + fHe*x_He0
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

          else if(y(2).gt.0.99)then

	    x_H0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He1/(Tnow*(1.d0+z)) ) / Nnow
	    rhs = rhs*4.d0		!ratio of g's is 4 for He+ <-> He0
	    x_He0 = 0.5d0 * ( dsqrt( (rhs-1.d0)**2 + 4.d0*(1.d0+fHe)*rhs )
	1	- (rhs-1.d0))
	    x0 = x_He0
	    x_He0 = (x0 - 1.d0)/fHe
	    y(1) = x_H0
	    y(2) = x_He0
	    y(3) = Tnow*(1.d0+z)

	  else if (y(1).gt.0.99d0) then

	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1/(Tnow*(1.d0+z)) ) / Nnow
	    x_H0 = 0.5d0 * (dsqrt( rhs**2+4.d0*rhs ) - rhs )

	    call DVERK(nw,ION,zstart,y,zend,tol,ind,cw,nw,w)
	    y(1) = x_H0
	    x0 = y(1) + fHe*y(2)

	  else

	    call DVERK(nw,ION,zstart,y,zend,tol,ind,cw,nw,w)

	    x0 = y(1) + fHe*y(2)

	  end if

	  Trad = Tnow * (1.d0+zend)
	  Tmat = y(3)
	  x_H = y(1)
	  x_He = y(2)
	  x = x0

	  out_xe(i)=x

	end do

	end

C	===============================================================
	subroutine GET_INIT(z,x_H0,x_He0,x0)

C	Set up the initial conditions so it will work for general,
C	but not pathological choices of zstart
C	Initial ionization fraction using Saha for relevant species

	implicit none

	real*8 OmegaT,HO,OmegaL,OmegaK
	real*8 z,x0,rhs,x_H0,x_He0
	real*8 Tnow,Nnow,z_eq
	real*8 Lambda,H_frac
	real*8 Lambda_He,Bfact,CK_He,CL_He
	real*8 CB1,CDB,CR,CK,CL,CT,fHe,CB1_He1,CB1_He2,CDB_He,fu

	common/Cdata/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
	1	fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He,fu
	common/Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
C	===============================================================

	if(z.gt.8000.d0)then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    x0 = 1.d0+2.d0*fHe

	else if(z.gt.3500.d0)then

	    x_H0 = 1.d0
	    x_He0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He2/(Tnow*(1.d0+z)) ) / Nnow
	rhs = rhs*1.d0		!ratio of g's is 1 for He++ <-> He+
	x0 = 0.5d0 * ( dsqrt( (rhs-1.d0-fHe)**2
	1	+ 4.d0*(1.d0+2.d0*fHe)*rhs) - (rhs-1.d0-fHe) )

	else if(z.gt.2000.d0)then

	x_H0 = 1.d0
	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1_He1/(Tnow*(1.d0+z)) ) / Nnow
	rhs = rhs*4.d0		!ratio of g's is 4 for He+ <-> He0
	    x_He0 = 0.5d0 * ( dsqrt( (rhs-1.d0)**2 + 4.d0*(1.d0+fHe)*rhs )
	1	- (rhs-1.d0))
	    x0 = x_He0
	    x_He0 = (x0 - 1.d0)/fHe

	else

	    rhs = dexp( 1.5d0 * dLog(CR*Tnow/(1.d0+z))
	1	- CB1/(Tnow*(1.d0+z)) ) / Nnow
	    x_H0 = 0.5d0 * (dsqrt( rhs**2+4.d0*rhs ) - rhs )
	    x_He0 = 0.d0
	    x0 = x_H0

	end if

	return

	end

C	===============================================================
	subroutine ION(Ndim,z,Y,f)

	implicit none

	integer Ndim,Heflag,Heswitch,Hswitch

	real*8 z,x,n,n_He,Trad,Tmat,x_H,x_He
	real*8 y(Ndim),f(Ndim)
	real*8 C,k_B,h_P,m_e,m_H,not4,sigma,a,Pi
	real*8 Lambda,H_frac,Lambda_He
	real*8 Tnow,HO,Nnow,z_eq,Hz,OmegaT,OmegaL,OmegaK
	real*8 Rup,Rdown,K,K_He,Rup_He,Rdown_He,He_Boltz
	real*8 timeTh,timeH,factor
	real*8 CB1,CDB,CR,CK,CL,CT,fHe,CB1_He1,CB1_He2,CDB_He,fu,b_He
	real*8 Bfact,CK_He,CL_He
	real*8 a_VF,b_VF,T_0,T_1,sq_0,sq_1,a_PPB,b_PPB,c_PPB,d_PPB
	real*8 tauHe_s,pHe_s
	real*8 A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt
	real*8 Doppler,gamma_2Ps,pb,qb,AHcon
	real*8 L_He_2p,L_He_2Pt,L_He_2St,L_He2St_ion
	real*8 a_trip,b_trip,Rdown_trip,Rup_trip
	real*8 tauHe_t,pHe_t,CL_PSt,CfHe_t,gamma_2Pt
	real*8 AGauss1,AGauss2,zGauss1,zGauss2,wGauss1,wGauss2
	real*8 dHdz,epsilon

	common/Cfund/C,k_B,h_P,m_e,m_H,not4,sigma,a,Pi
	common/Cdata/Lambda,H_frac,CB1,CDB,CR,CK,CL,CT,
	1	fHe,CB1_He1,CB1_He2,CDB_He,Lambda_He,Bfact,CK_He,CL_He,fu
	common/Hemod/b_He,A2P_s,A2P_t,sigma_He_2Ps,sigma_He_2Pt,
	1	L_He_2p,L_He_2Pt,L_He_2St,L_He2St_ion
	common/Hmod/AGauss1,AGauss2,zGauss1,zGauss2,wGauss1,wGauss2
	common/Switch/Heswitch,Hswitch
	common/Cosmo/Tnow,HO,Nnow,z_eq,OmegaT,OmegaL,OmegaK
C	===============================================================

c	the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen	
	a_PPB = 4.309d0
	b_PPB = -0.6166d0
	c_PPB = 0.6703d0
	d_PPB = 0.5300d0
c	the Verner and Ferland type fitting parameters for Helium
c	fixed to match those in the SSS papers, and now correct
	a_VF = 10.d0**(-16.744d0)
	b_VF = 0.711d0
	T_0 = 10.d0**(0.477121d0)	!3K
	T_1 = 10.d0**(5.114d0)
c	fitting parameters for HeI triplets
c	(matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)
	a_trip = 10.d0**(-16.306d0)
	b_trip = 0.761D0

	x_H = y(1)
	x_He = y(2)
	x = x_H + fHe * x_He
	Tmat = y(3)

	n = Nnow * (1.d0+z)**3
	n_He = fHe * Nnow * (1.d0+z)**3
	Trad = Tnow * (1.d0+z)
	Hz = HO * dsqrt((1.d0+z)**4/(1.d0+z_eq)*OmegaT + OmegaT*(1.d0+z)**3
	1	+ OmegaK*(1.d0+z)**2 + OmegaL)

c	Also calculate derivative for use later
	dHdz = (HO**2/2.d0/Hz)*(4.d0*(1.d0+z)**3/(1.d0+z_eq)*OmegaT 
	1	+ 3.d0*OmegaT*(1.d0+z)**2 + 2.d0*OmegaK*(1.d0+z) )

c	Get the radiative rates using PPQ fit (identical to Hummer's table)
	Rdown=1.d-19*a_PPB*(Tmat/1.d4)**b_PPB
	1	/(1.d0+c_PPB*(Tmat/1.d4)**d_PPB)
	Rup = Rdown * (CR*Tmat)**(1.5d0)*dexp(-CDB/Tmat)

c	calculate He using a fit to a Verner & Ferland type formula
	sq_0 = dsqrt(Tmat/T_0)
	sq_1 = dsqrt(Tmat/T_1)
c	typo here corrected by Wayne Hu and Savita Gahlaut
	Rdown_He = a_VF/(sq_0*(1.d0+sq_0)**(1.d0-b_VF))
	Rdown_He = Rdown_He/(1.d0+sq_1)**(1.d0+b_VF)
	Rup_He = Rdown_He*(CR*Tmat)**(1.5d0)*dexp(-CDB_He/Tmat)
	Rup_He = 4.d0*Rup_He	!statistical weights factor for HeI
c	Avoid overflow (pointed out by Jacques Roland)
	if((Bfact/Tmat).gt.680.d0)then
	  He_Boltz = dexp(680.d0)
	else
	  He_Boltz = dexp(Bfact/Tmat)
	end if

c	now deal with H and its fudges
	if (Hswitch.eq.0) then 
	  K = CK/Hz		!Peebles coefficient K=lambda_a^3/8piH
	else
c	fit a double Gaussian correction function
	  K = CK/Hz*(1.0d0
	1	+AGauss1*dexp(-((log(1.0d0+z)-zGauss1)/wGauss1)**2.d0)
	2	+AGauss2*dexp(-((log(1.0d0+z)-zGauss2)/wGauss2)**2.d0))
	end if

c	add the HeI part, using same T_0 and T_1 values
	Rdown_trip = a_trip/(sq_0*(1.d0+sq_0)**(1.0-b_trip))
	Rdown_trip = Rdown_trip/((1.d0+sq_1)**(1.d0+b_trip))
	Rup_trip = Rdown_trip*dexp(-h_P*C*L_He2St_ion/(k_B*Tmat))
	Rup_trip = Rup_trip*((CR*Tmat)**(1.5d0))*(4.d0/3.d0)
c	last factor here is the statistical weight

c       try to avoid "NaN" when x_He gets too small
	if ((x_He.lt.5.d-9) .or. (x_He.gt.0.980)) then 
	  Heflag = 0
	else
	  Heflag = Heswitch
	end if
	if (Heflag.eq.0)then		!use Peebles coeff. for He
	  K_He = CK_He/Hz
	else	!for Heflag>0 		!use Sobolev escape probability
	  tauHe_s = A2P_s*CK_He*3.d0*n_He*(1.d0-x_He)/Hz
	  pHe_s = (1.d0 - dexp(-tauHe_s))/tauHe_s
	  K_He = 1.d0/(A2P_s*pHe_s*3.d0*n_He*(1.d0-x_He))
c	smoother criterion here from Antony Lewis & Chad Fendt
	  if (((Heflag.eq.2).or.(Heflag.ge.5)).and.(x_H.lt.0.9999999d0))then
c	use fitting formula for continuum opacity of H
c	first get the Doppler width parameter
	    Doppler = 2.D0*k_B*Tmat/(m_H*not4*C*C)
	    Doppler = C*L_He_2p*dsqrt(Doppler)
	    gamma_2Ps = 3.d0*A2P_s*fHe*(1.d0-x_He)*C*C
	1	/(dsqrt(Pi)*sigma_He_2Ps*8.d0*Pi*Doppler*(1.d0-x_H))
	2	/((C*L_He_2p)**2.d0)
	    pb = 0.36d0  !value from KIV (2007)
	    qb = b_He
c	calculate AHcon, the value of A*p_(con,H) for H continuum opacity
	    AHcon = A2P_s/(1.d0+pb*(gamma_2Ps**qb))
	    K_He=1.d0/((A2P_s*pHe_s+AHcon)*3.d0*n_He*(1.d0-x_He))
	  end if
	  if (Heflag.ge.3) then		!include triplet effects
	    tauHe_t = A2P_t*n_He*(1.d0-x_He)*3.d0
	    tauHe_t = tauHe_t /(8.d0*Pi*Hz*L_He_2Pt**(3.d0))
	    pHe_t = (1.d0 - dexp(-tauHe_t))/tauHe_t
	    CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
	  if ((Heflag.eq.3) .or. (Heflag.eq.5).or.(x_H.gt.0.99999d0)) then
c	no H cont. effect
	    CfHe_t = A2P_t*pHe_t*dexp(-CL_PSt/Tmat)
	    CfHe_t = CfHe_t/(Rup_trip+CfHe_t)	!"C" factor for triplets
	  else					!include H cont. effect
	    Doppler = 2.d0*k_B*Tmat/(m_H*not4*C*C)
	    Doppler = C*L_He_2Pt*dsqrt(Doppler)
	    gamma_2Pt = 3.d0*A2P_t*fHe*(1.d0-x_He)*C*C
	1	/(dsqrt(Pi)*sigma_He_2Pt*8.d0*Pi*Doppler*(1.d0-x_H))
	2	/((C*L_He_2Pt)**2.d0)
c	use the fitting parameters from KIV (2007) in this case
	    pb = 0.66d0
	    qb = 0.9d0
	    AHcon = A2P_t/(1.d0+pb*gamma_2Pt**qb)/3.d0
	    CfHe_t = (A2P_t*pHe_t+AHcon)*dexp(-CL_PSt/Tmat)
	    CfHe_t = CfHe_t/(Rup_trip+CfHe_t)	!"C" factor for triplets
	  end if
	  end if
	end if

c	Estimates of Thomson scattering time and Hubble time
	timeTh=(1.d0/(CT*Trad**4))*(1.d0+x+fHe)/x	!Thomson time
	timeH=2.d0/(3.d0*HO*(1.d0+z)**1.5)		!Hubble time

c	calculate the derivatives
c	turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
c	(clunky, but seems to work)
	if (x_H.gt.0.99d0) then			!don't change at all
		f(1) = 0.d0
cc	else if ((x_H.gt.0.98d0).and.(Heflag.eq.0)) then	!don't modify
	  else if (x_H.gt.0.985d0) then		!use Saha rate for Hydrogen
		f(1) = (x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat))
	1	/(Hz*(1.d0+z))
c	for interest, calculate the correction factor compared to Saha
c	(without the fudge)
		factor=(1.d0 + K*Lambda*n*(1.d0-x_H))
	1	/(Hz*(1.d0+z)*(1.d0+K*Lambda*n*(1.d0-x)
	2	+K*Rup*n*(1.d0-x)))
	  else					!use full rate for H
		f(1) = ((x*x_H*n*Rdown - Rup*(1.d0-x_H)*dexp(-CL/Tmat))
	1		*(1.d0 + K*Lambda*n*(1.d0-x_H)))
	2	/(Hz*(1.d0+z)*(1.d0/fu+K*Lambda*n*(1.d0-x_H)/fu
	3	+K*Rup*n*(1.d0-x_H)))
	end if
c	turn off the He once it is small
	if (x_He.lt.1.d-15) then
		f(2)=0.d0
	else
		f(2) = ((x*x_He*n*Rdown_He 
	1	- Rup_He*(1.d0-x_He)*dexp(-CL_He/Tmat))
	2		*(1.d0+ K_He*Lambda_He*n_He*(1.d0-x_He)*He_Boltz))
	3	/(Hz*(1.d0+z)
	4	 * (1.d0 + K_He*(Lambda_He+Rup_He)*n_He*(1.d0-x_He)*He_Boltz))
c	Modification to HeI recombination including channel via triplets
	  if (Heflag.ge.3) then
		f(2) = f(2)+ (x*x_He*n*Rdown_trip
	1	- (1.d0-x_He)*3.d0*Rup_trip*dexp(-h_P*C*L_He_2st/(k_B*Tmat)))
	2	*CfHe_t/(Hz*(1.d0+z))
	  end if
	end if

c	follow the matter temperature once it has a chance of diverging

	if (timeTh.lt.H_frac*timeH) then
c		f(3)=Tmat/(1.d0+z)	!Tmat follows Trad
c	additional term to smooth transition to Tmat evolution,
c	(suggested by Adam Moss)
		epsilon = Hz*(1.d0+x+fHe)/(CT*Trad**3*x)
		f(3) = Tnow
	1	+ epsilon*((1.d0+fHe)/(1.d0+fHe+x))*((f(1)+fHe*f(2))/x)
	2	- epsilon* dHdz/Hz + 3.0d0*epsilon/(1.d0+z) 
	else
		f(3)= CT * (Trad**4) * x / (1.d0+x+fHe)
	1		* (Tmat-Trad) / (Hz*(1.d0+z)) + 2.d0*Tmat/(1.d0+z)
	end if

	return

	end
	end

C===============================================================================
      subroutine DVERK (n, fcn, x, y, xend, tol, ind, c, nw, w)
      integer n, ind, nw, k
      real*8 x, y(n), xend, tol, c(24), w(nw,9), temp
c
      external fcn
c
c     ******************************************************************
c     * begin initialization, parameter checking, interrupt re-entries *
c     ******************************************************************
c
c  ......abort if ind out of range 1 to 6
         if (ind.lt.1 .or. ind.gt.6) go to 500
c
c        cases - initial entry, normal re-entry, interrupt re-entries
         go to (5, 5, 45, 1111, 2222, 2222), ind
c        case 1 - initial entry (ind .eq. 1 or 2)
c  .........abort if n.gt.nw or tol.le.0
    5       if (n.gt.nw .or. tol.le.0.d0) go to 500
            if (ind.eq. 2) go to 15
c              initial entry without options (ind .eq. 1)
c              set c(1) to c(9) equal to 0
               do 10 k = 1, 9
                  c(k) = 0.d0
   10          continue
               go to 35
   15       continue
c              initial entry with options (ind .eq. 2)
c              make c(1) to c(9) non-negative
               do 20 k = 1, 9
                  c(k) = dabs(c(k))
   20          continue
c              make floor values non-negative if they are to be used
               if (c(1).ne.4.d0 .and. c(1).ne.5.d0) go to 30
                  do 25 k = 1, n
                     c(k+30) = dabs(c(k+30))
   25             continue
   30          continue
   35       continue
c           initialize rreb, dwarf, prev xend, flag, counts
            c(10) = 2.d0**(-56)
            c(11) = 1.d-35
c           set previous xend initially to initial value of x
            c(20) = x
            do 40 k = 21, 24
               c(k) = 0.d0
   40       continue
            go to 50
c        case 2 - normal re-entry (ind .eq. 3)
c  .........abort if xend reached, and either x changed or xend not
   45       if (c(21).ne.0.d0 .and.
     +                        (x.ne.c(20) .or. xend.eq.c(20))) go to 500
c           re-initialize flag
            c(21) = 0.d0
            go to 50
c        case 3 - re-entry following an interrupt (ind .eq. 4 to 6)
c           transfer control to the appropriate re-entry point..........
c           this has already been handled by the computed go to        .
c        end cases                                                     v
   50    continue
c
c     end initialization, etc.
c
c     ******************************************************************
c     * loop through the following 4 stages, once for each trial  step *
c     * until the occurrence of one of the following                   *
c     *    (a) the normal return (with ind .eq. 3) on reaching xend in *
c     *        stage 4                                                 *
c     *    (b) an error return (with ind .lt. 0) in stage 1 or stage 4 *
c     *    (c) an interrupt return (with ind  .eq.  4,  5  or  6),  if *
c     *        requested, in stage 1 or stage 4                        *
c     ******************************************************************
c
99999 continue
c
c        ***************************************************************
c        * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
c        * and some parameter  checking,  and  end  up  with  suitable *
c        * values of hmag, xtrial and htrial in preparation for taking *
c        * an integration step.                                        *
c        ***************************************************************
c
c***********error return (with ind=-1) if no of fcn evals too great
            if (c(7).eq.0.d0 .or. c(24).lt.c(7)) go to 100
               ind = -1
               return
  100       continue
c
c           calculate slope (adding 1 to no of fcn evals) if ind .ne. 6
            if (ind .eq. 6) go to 105
               call fcn(n, x, y, w(1,1))
               c(24) = c(24) + 1.d0
  105       continue
c
c           calculate hmin - use default unless value prescribed
            c(13) = c(3)
            if (c(3) .ne. 0.d0) go to 165
c              calculate default value of hmin
c              first calculate weighted norm y - c(12) - as specified
c              by the error control indicator c(1)
               temp = 0.d0
               if (c(1) .ne. 1.d0) go to 115
c                 absolute error control - weights are 1
                  do 110 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  110             continue
                  c(12) = temp
                  go to 160
  115          if (c(1) .ne. 2.d0) go to 120
c                 relative error control - weights are 1/dabs(y(k)) so
c                 weighted norm y is 1
                  c(12) = 1.d0
                  go to 160
  120          if (c(1) .ne. 3.d0) go to 130
c                 weights are 1/max(c(2),abs(y(k)))
                  do 125 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(2))
  125             continue
                  c(12) = dmin1(temp, 1.d0)
                  go to 160
  130          if (c(1) .ne. 4.d0) go to 140
c                 weights are 1/max(c(k+30),abs(y(k)))
                  do 135 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  135             continue
                  c(12) = dmin1(temp, 1.d0)
                  go to 160
  140          if (c(1) .ne. 5.d0) go to 150
c                 weights are 1/c(k+30)
                  do 145 k = 1, n
                     temp = dmax1(temp, dabs(y(k))/c(k+30))
  145             continue
                  c(12) = temp
                  go to 160
  150          continue
c                 default case - weights are 1/max(1,abs(y(k)))
                  do 155 k = 1, n
                     temp = dmax1(temp, dabs(y(k)))
  155             continue
                  c(12) = dmin1(temp, 1.d0)
  160          continue
               c(13) = 10.d0*dmax1(c(11),c(10)*dmax1(c(12)/tol,dabs(x)))
  165       continue
c
c           calculate scale - use default unless value prescribed
            c(15) = c(5)
            if (c(5) .eq. 0.d0) c(15) = 1.d0
c
c           calculate hmax - consider 4 cases
c           case 1 both hmax and scale prescribed
               if (c(6).ne.0.d0 .and. c(5).ne.0.d0)
     +                                    c(16) = dmin1(c(6), 2.d0/c(5))
c           case 2 - hmax prescribed, but scale not
               if (c(6).ne.0.d0 .and. c(5).eq.0.d0) c(16) = c(6)
c           case 3 - hmax not prescribed, but scale is
               if (c(6).eq.0.d0 .and. c(5).ne.0.d0) c(16) = 2.d0/c(5)
c           case 4 - neither hmax nor scale is provided
               if (c(6).eq.0.d0 .and. c(5).eq.0.d0) c(16) = 2.d0
c
c***********error return (with ind=-2) if hmin .gt. hmax
            if (c(13) .le. c(16)) go to 170
               ind = -2
               return
  170       continue
c
c           calculate preliminary hmag - consider 3 cases
            if (ind .gt. 2) go to 175
c           case 1 - initial entry - use prescribed value of hstart, if
c              any, else default
               c(14) = c(4)
               if (c(4) .eq. 0.d0) c(14) = c(16)*tol**(1./6.)
               go to 185
  175       if (c(23) .gt. 1.d0) go to 180
c           case 2 - after a successful step, or at most  one  failure,
c              use min(2, .9*(tol/est)**(1/6))*hmag, but avoid possible
c              overflow. then avoid reduction by more than half.
               temp = 2.d0*c(14)
               if (tol .lt. (2.d0/.9d0)**6*c(19))
     +                            temp = .9d0*(tol/c(19))**(1./6.)*c(14)
               c(14) = dmax1(temp, .5d0*c(14))
               go to 185
  180       continue
c           case 3 - after two or more successive failures
               c(14) = .5d0*c(14)
  185       continue
c
c           check against hmax
            c(14) = dmin1(c(14), c(16))
c
c           check against hmin
            c(14) = dmax1(c(14), c(13))
c
c***********interrupt no 1 (with ind=4) if requested
            if (c(8) .eq. 0.d0) go to 1111
               ind = 4
               return
c           resume here on re-entry with ind .eq. 4   ........re-entry..
 1111       continue
c
c           calculate hmag, xtrial - depending on preliminary hmag, xend
            if (c(14) .ge. dabs(xend - x)) go to 190
c              do not step more than half way to xend
               c(14) = dmin1(c(14), .5d0*dabs(xend - x))
               c(17) = x + dsign(c(14), xend - x)
               go to 195
  190       continue
c              hit xend exactly
               c(14) = dabs(xend - x)
               c(17) = xend
  195       continue
c
c           calculate htrial
            c(18) = c(17) - x
c
c        end stage 1
c
c        ***************************************************************
c        * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
c        * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
c        * stage 3. w(*,9) is temporary storage until finally it holds *
c        * ytrial.                                                     *
c        ***************************************************************
c
            temp = c(18)/1398169080000.d0
c
            do 200 k = 1, n
               w(k,9) = y(k) + temp*w(k,1)*233028180000.d0
  200       continue
            call fcn(n, x + c(18)/6.d0, w(1,9), w(1,2))
c
            do 205 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*74569017600.d0
     +                                + w(k,2)*298276070400.d0  )
  205       continue
            call fcn(n, x + c(18)*(4.d0/15.d0), w(1,9), w(1,3))
c
            do 210 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*1165140900000.d0
     +                                - w(k,2)*3728450880000.d0
     +                                + w(k,3)*3495422700000.d0 )
  210       continue
            call fcn(n, x + c(18)*(2.d0/3.d0), w(1,9), w(1,4))
c
            do 215 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*3604654659375.d0
     +                                + w(k,2)*12816549900000.d0
     +                                - w(k,3)*9284716546875.d0
     +                                + w(k,4)*1237962206250.d0 )
  215       continue
            call fcn(n, x + c(18)*(5.d0/6.d0), w(1,9), w(1,5))
c
            do 220 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*3355605792000.d0
     +                                - w(k,2)*11185352640000.d0
     +                                + w(k,3)*9172628850000.d0
     +                                - w(k,4)*427218330000.d0
     +                                + w(k,5)*482505408000.d0  )
  220       continue
            call fcn(n, x + c(18), w(1,9), w(1,6))
c
            do 225 k = 1, n
               w(k,9) = y(k) + temp*( - w(k,1)*770204740536.d0
     +                                + w(k,2)*2311639545600.d0
     +                                - w(k,3)*1322092233000.d0
     +                                - w(k,4)*453006781920.d0
     +                                + w(k,5)*326875481856.d0  )
  225       continue
            call fcn(n, x + c(18)/15.d0, w(1,9), w(1,7))
c
            do 230 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*2845924389000.d0
     +                                - w(k,2)*9754668000000.d0
     +                                + w(k,3)*7897110375000.d0
     +                                - w(k,4)*192082660000.d0
     +                                + w(k,5)*400298976000.d0
     +                                + w(k,7)*201586000000.d0  )
  230       continue
            call fcn(n, x + c(18), w(1,9), w(1,8))
c
c           calculate ytrial, the extrapolated approximation and store
c              in w(*,9)
            do 235 k = 1, n
               w(k,9) = y(k) + temp*(   w(k,1)*104862681000.d0
     +                                + w(k,3)*545186250000.d0
     +                                + w(k,4)*446637345000.d0
     +                                + w(k,5)*188806464000.d0
     +                                + w(k,7)*15076875000.d0
     +                                + w(k,8)*97599465000.d0   )
  235       continue
c
c           add 7 to the no of fcn evals
            c(24) = c(24) + 7.d0
c
c        end stage 2
c
c        ***************************************************************
c        * stage 3 - calculate the error estimate est. first calculate *
c        * the  unweighted  absolute  error  estimate vector (per unit *
c        * step) for the unextrapolated approximation and store it  in *
c        * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
c        * specified by the error  control  indicator  c(1).  finally, *
c        * modify  this result to produce est, the error estimate (per *
c        * unit step) for the extrapolated approximation ytrial.       *
c        ***************************************************************
c
c           calculate the unweighted absolute error estimate vector
            do 300 k = 1, n
               w(k,2) = (   w(k,1)*8738556750.d0
     +                    + w(k,3)*9735468750.d0
     +                    - w(k,4)*9709507500.d0
     +                    + w(k,5)*8582112000.d0
     +                    + w(k,6)*95329710000.d0
     +                    - w(k,7)*15076875000.d0
     +                    - w(k,8)*97599465000.d0)/1398169080000.d0
  300       continue
c
c           calculate the weighted max norm of w(*,2) as specified by
c           the error control indicator c(1)
            temp = 0.d0
            if (c(1) .ne. 1.d0) go to 310
c              absolute error control
               do 305 k = 1, n
                  temp = dmax1(temp,dabs(w(k,2)))
  305          continue
               go to 360
  310       if (c(1) .ne. 2.d0) go to 320
c              relative error control
               do 315 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/y(k)))
  315          continue
               go to 360
  320       if (c(1) .ne. 3.d0) go to 330
c              weights are 1/max(c(2),abs(y(k)))
               do 325 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))
     +                             / dmax1(c(2), dabs(y(k))) )
  325          continue
               go to 360
  330       if (c(1) .ne. 4.d0) go to 340
c              weights are 1/max(c(k+30),abs(y(k)))
               do 335 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))
     +                             / dmax1(c(k+30), dabs(y(k))) )
  335          continue
               go to 360
  340       if (c(1) .ne. 5.d0) go to 350
c              weights are 1/c(k+30)
               do 345 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2)/c(k+30)))
  345          continue
               go to 360
  350       continue
c              default case - weights are 1/max(1,abs(y(k)))
               do 355 k = 1, n
                  temp = dmax1(temp, dabs(w(k,2))
     +                             / dmax1(1.d0, dabs(y(k))) )
  355          continue
  360       continue
c
c           calculate est - (the weighted max norm of w(*,2))*hmag*scale
c              - est is intended to be a measure of the error  per  unit
c              step in ytrial
            c(19) = temp*c(14)*c(15)
c
c        end stage 3
c
c        ***************************************************************
c        * stage 4 - make decisions.                                   *
c        ***************************************************************
c
c           set ind=5 if step acceptable, else set ind=6
            ind = 5
            if (c(19) .gt. tol) ind = 6
c
c***********interrupt no 2 if requested
            if (c(9) .eq. 0.d0) go to 2222
               return
c           resume here on re-entry with ind .eq. 5 or 6   ...re-entry..
 2222       continue
c
            if (ind .eq. 6) go to 410
c              step accepted (ind .eq. 5), so update x, y from xtrial,
c                 ytrial, add 1 to the no of successful steps, and set
c                 the no of successive failures to zero
               x = c(17)
               do 400 k = 1, n
                  y(k) = w(k,9)
  400          continue
               c(22) = c(22) + 1.d0
               c(23) = 0.d0
c**************return(with ind=3, xend saved, flag set) if x .eq. xend
               if (x .ne. xend) go to 405
                  ind = 3
                  c(20) = xend
                  c(21) = 1.d0
                  return
  405          continue
               go to 420
  410       continue
c              step not accepted (ind .eq. 6), so add 1 to the no of
c                 successive failures
               c(23) = c(23) + 1.d0
c**************error return (with ind=-3) if hmag .le. hmin
               if (c(14) .gt. c(13)) go to 415
                  ind = -3
                  return
  415          continue
  420       continue
c
c        end stage 4
c
      go to 99999
c     end loop
c
c  begin abort action
  500 continue
c
      write(*,505) ind, tol, x, n, c(13), xend, nw, c(16), c(20),
     +      c(22), c(23), c(24), (y(k), k = 1, n)
  505 format( /// 1h0, 58hcomputation stopped in DVERK with the followin
     +g values -
     +   / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =,
     +          1pd22.15
     +   / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =,
     +          1pd22.15
     +   / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =,
     +          1pd22.15
     +   / 1h0, 14x, 27hno of successful steps    =, 0pf8.0
     +   / 1h , 14x, 27hno of successive failures =, 0pf8.0
     +   / 1h , 14x, 27hno of function evals      =, 0pf8.0
     +   / 1h0, 23hthe components of y are
     +   // (1h , 1p5d24.15)                                           )
c
      stop
c
c  end abort action
c
      end
