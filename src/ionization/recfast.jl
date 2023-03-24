abstract type IonizationIntegrator end
abstract type AbstractIonizationHistory{T, IT<:AbstractInterpolation{T,1}} end

const H0_natural_unit_conversion = ustrip(unnatural(u"s", 1.0*unit(natural(1u"s"))))
const Kelvin_natural_unit_conversion = ustrip(unnatural(1u"K", 1.0*unit(natural(1u"K")) ))

struct IonizationHistory{T, IT} <: AbstractIonizationHistory{T, IT}
	τ₀::T
    Xₑ::IT
    τ::IT
    τ′::IT
    τ′′::IT
    g̃::IT
    g̃′::IT
    g̃′′::IT
    Tmat::IT
    csb²::IT
    # Trad::IT #This is never used
end


@with_kw struct RECFAST{T, AB<:AbstractBackground{T}} <: IonizationIntegrator @deftype T
    bg::AB  # a RECFAST has an associated background evolution
    C  = 2.99792458e8  # Fundamental constants in SI units
    k_B = 1.380658e-23
    h_P = 6.6260755e-34
    m_e = 9.1093897e-31
    m_H = 1.673575e-27  #	av. H atom
    # note: neglecting deuterium, making an O(e-5) effect
    not4 = 3.9715e0  # mass He/H atom  ("not4" pointed out by Gary Steigman)
    sigma = 6.6524616e-29
    a = 7.565914e-16
    G = 6.6742e-11 	# new value

    Lambda = 8.2245809e0
    Lambda_He = 51.3e0              # new value from Dalgarno
    L_H_ion = 1.096787737e7         # level for H ion. (in m^-1)
    L_H_alpha = 8.225916453e6       # averaged over 2 levels
    L_He1_ion = 1.98310772e7        # from Drake (1993)
    L_He2_ion = 4.389088863e7       # from JPhysChemRefData (1987)
    L_He_2s	= 1.66277434e7          # from Drake (1993)
    L_He_2p	= 1.71134891e7          # from Drake (1993)
    # C	2 photon rates and atomic levels in SI units

    A2P_s = 1.798287e9              # Morton, Wu & Drake (2006)
    A2P_t = 177.58e0                # Lach & Pachuski (2001)
    L_He_2Pt = 1.690871466e7        # Drake & Morton (2007)
    L_He_2St = 1.5985597526e7       # Drake & Morton (2007)
    L_He2St_ion = 3.8454693845e6    # Drake & Morton (2007)
    sigma_He_2Ps = 1.436289e-22     # Hummer & Storey (1998)
    sigma_He_2Pt = 1.484872e-22     # Hummer & Storey (1998)
    # C	Atomic data for HeI

    AGauss1	= -0.14e0               # Amplitude of 1st Gaussian
    AGauss2 = 0.079e0               # Amplitude of 2nd Gaussian
    zGauss1 = 7.28e0                # ln(1+z) of 1st Gaussian
    zGauss2 = 6.73e0                # ln(1+z) of 2nd Gaussian
    wGauss1 = 0.18e0                # Width of 1st Gaussian
    wGauss2 = 0.33e0                # Width of 2nd Gaussian
    # Gaussian fits for extra H physics (fit by Adam Moss, modified by Antony Lewis)

    # the Pequignot, Petitjean & Boisson fitting parameters for Hydrogen
	a_PPB = 4.309
	b_PPB = -0.6166
	c_PPB = 0.6703
	d_PPB = 0.5300
    # the Verner and Ferland type fitting parameters for Helium
    # fixed to match those in the SSS papers, and now correct
	a_VF = 10^(-16.744)
	b_VF = 0.711
	T_0 = 10^(0.477121)	#!3K
	T_1 = 10^(5.114)
    # fitting parameters for HeI triplets
    # (matches Hummer's table with <1% error for 10^2.8 < T/K < 10^4)
	a_trip = 10^(-16.306)
	b_trip = 0.761

    # Set up some constants so they don't have to be calculated later
    Lalpha = 1/L_H_alpha
    Lalpha_He = 1/L_He_2p
    DeltaB = h_P*C*(L_H_ion-L_H_alpha)
    CDB = DeltaB/k_B
    DeltaB_He = h_P*C*(L_He1_ion-L_He_2s)	# 2s, not 2p
    CDB_He = DeltaB_He/k_B
    CB1 = h_P*C*L_H_ion/k_B
    CB1_He1 = h_P*C*L_He1_ion/k_B	# ionization for HeI
    CB1_He2 = h_P*C*L_He2_ion/k_B	# ionization for HeII
    CR = 2π * (m_e/h_P)*(k_B/h_P)
    CK = Lalpha^3/(8π)
    CK_He = Lalpha_He^3/(8π)
    CL = C*h_P/(k_B*Lalpha)
    CL_He = C*h_P/(k_B/L_He_2s)	# comes from det.bal. of 2s-1s
    CT = (8/3)*(sigma/(m_e*C))*a
    Bfact = h_P*C*(L_He_2p-L_He_2s)/k_B

    # Matter departs from radiation when t(Th) > H_frac * t(H)
    H_frac = 1e-3  # choose some safely small number

    # switches
    Hswitch::Int64 = 1
    Heswitch::Int64 = 6

    # Cosmology
    Yp::T = 0.24
    OmegaB::T = 0.046  # TODO: should replace during GREAT GENERALIZATION
    HO =  bg.H₀ / H0_natural_unit_conversion
	OmegaG::T = 5.0469e-5 #not sure this is the best way to do this
	Tnow = (15/ π^2 *bg.ρ_crit * OmegaG)^(1/4) * Kelvin_natural_unit_conversion #last thing is natural to K
	# This was hardcoded originally as: Tnow = 2.725, fixes issue downstream with Duals
	# Had to change RECFAST test temperature though - should double check this

    # sort out the helium abundance parameters
    mu_H = 1 / (1 - Yp)			 # Mass per H atom
    mu_T = not4/(not4-(not4-1)*Yp)	 # Mass per atom
    fHe = Yp/(not4*(1 - Yp))		# n_He_tot / n_H_tot

    Nnow = 3 * HO * HO * OmegaB / (8π * G * mu_H * m_H)  # TODO: should replace during GREAT GENERALIZATION
    fu = (Hswitch == 0) ? 1.14 : 1.125
    b_He = 0.86  # Set the He fudge factor
    tol = 1e-6
end

# helper constructor which dispatches on the background
RECFAST(bg::AB; kws...) where {T, AB<:AbstractBackground{T}} = RECFAST{T,AB}(bg=bg; kws...)


function recfast_init(𝕣::RECFAST, z)
    if z > 8000.
        x_H0 = 1.
        x_He0 = 1.
        x0 = 1. + 2 * 𝕣.fHe
    elseif z > 3500.
        x_H0 = 1.
        x_He0 = 1.
        rhs = exp( 1.5 * log(𝕣.CR*𝕣.Tnow/(1 + z)) - 𝕣.CB1_He2/(𝕣.Tnow*(1 + z)) ) / 𝕣.Nnow
	    rhs = rhs * 1.  # ratio of g's is 1 for He++ <-> He+
	    x0 = 0.5 * ( sqrt( (rhs - 1 - 𝕣.fHe)^2 + 4 * (1 + 2 * 𝕣.fHe) * rhs) - (rhs - 1 - 𝕣.fHe) )
    elseif z > 2000.
	    x_H0 = 1.
	    rhs = exp( 1.5 * log(𝕣.CR * 𝕣.Tnow / (1 + z)) - 𝕣.CB1_He1/(𝕣.Tnow*(1 + z)) ) / 𝕣.Nnow
	    rhs = 4rhs    # ratio of g's is 4 for He+ <-> He0
	    x_He0 = 0.5 * ( sqrt( (rhs-1)^2 + 4*(1 + 𝕣.fHe)*rhs) - (rhs-1))
	    x0 = x_He0
	    x_He0 = (x0 - 1.)/𝕣.fHe
    else
	    rhs = exp( 1.5 * log(𝕣.CR*𝕣.Tnow/(1 + z)) - 𝕣.CB1/(𝕣.Tnow*(1 + z)) ) / 𝕣.Nnow
	    x_H0 = 0.5 * (sqrt( rhs^2 + 4 * rhs ) - rhs )
	    x_He0 = 0.
	    x0 = x_H0
    end

    return x_H0, x_He0, x0
end

#RECFAST f' at particular z
function ion_recfast(y, 𝕣::RECFAST{T}, z) where {T}

    f1, f2, f3 = zero(T), zero(T), zero(T)

	x_H = y[1]
	x_He = y[2]
	x = x_H + 𝕣.fHe * x_He
	Tmat = y[3]

	n = 𝕣.Nnow * (1+z)^3
	n_He = 𝕣.fHe * 𝕣.Nnow * (1+z)^3
	Trad = 𝕣.Tnow * (1+z)

    a = 1 / (1+z)  # scale factor
    x_a = a2x(a)
	Hz = 𝕣.bg.ℋ(x_a) / a / H0_natural_unit_conversion
	dHdz = (-𝕣.bg.ℋ′(x_a) + 𝕣.bg.ℋ(x_a)) / H0_natural_unit_conversion

    # Get the radiative rates using PPQ fit (identical to Hummer's table)
	Rdown=1e-19*𝕣.a_PPB*(Tmat/1e4)^𝕣.b_PPB/(1. + 𝕣.c_PPB*(Tmat/1e4)^𝕣.d_PPB)
	Rup = Rdown * (𝕣.CR*Tmat)^(1.5)*exp(-𝕣.CDB/Tmat)

    # calculate He using a fit to a Verner & Ferland type formula
	sq_0 = sqrt(Tmat/𝕣.T_0)
	sq_1 = sqrt(Tmat/𝕣.T_1)
    # typo here corrected by Wayne Hu and Savita Gahlaut
	Rdown_He = 𝕣.a_VF/(sq_0*(1+sq_0)^(1-𝕣.b_VF))
	Rdown_He = Rdown_He/(1+sq_1)^(1+𝕣.b_VF)
	Rup_He = Rdown_He*(𝕣.CR*Tmat)^(1.5)*exp(-𝕣.CDB_He/Tmat)
	Rup_He = 4. * Rup_He # statistical weights factor for HeI
    # Avoid overflow (pointed out by Jacques Roland)
	if((𝕣.Bfact/Tmat) > 680.)
	  He_Boltz = exp(680.)
	else
	  He_Boltz = exp(𝕣.Bfact/Tmat)
	end

    # now deal with H and its fudges
	if (𝕣.Hswitch == 0)
	    K = 𝕣.CK / Hz # !Peebles coefficient K=lambda_a^3/8piH
	else
        # fit a double Gaussian correction function
        K = 𝕣.CK / Hz*(1.0
            + 𝕣.AGauss1*exp(-((log(1+z)-𝕣.zGauss1)/𝕣.wGauss1)^2)
            + 𝕣.AGauss2*exp(-((log(1+z)-𝕣.zGauss2)/𝕣.wGauss2)^2))
	end

    # add the HeI part, using same T_0 and T_1 values
	Rdown_trip = 𝕣.a_trip/(sq_0*(1+sq_0)^(1-𝕣.b_trip))
	Rdown_trip = Rdown_trip/((1+sq_1)^(1+𝕣.b_trip))
	Rup_trip = Rdown_trip*exp(-𝕣.h_P*𝕣.C*𝕣.L_He2St_ion/(𝕣.k_B*Tmat))
	Rup_trip = Rup_trip*((𝕣.CR*Tmat)^1.5)*(4/3)
    # last factor here is the statistical weight

    # try to avoid "NaN" when x_He gets too small
	if ((x_He < 5.e-9) || (x_He > 0.980))
        Heflag = 0
	else
	    Heflag = 𝕣.Heswitch
	end
	if (Heflag == 0)  # use Peebles coeff. for He
	    K_He = 𝕣.CK_He/Hz
	else # for Heflag>0 		!use Sobolev escape probability
        tauHe_s = 𝕣.A2P_s*𝕣.CK_He*3*n_He*(1-x_He)/Hz
        pHe_s = (1 - exp(-tauHe_s))/tauHe_s
        K_He = 1 / (𝕣.A2P_s*pHe_s*3*n_He*(1-x_He))
        # smoother criterion here from Antony Lewis & Chad Fendt
	    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.9999999))
            # use fitting formula for continuum opacity of H
            # first get the Doppler width parameter
            Doppler = 2*𝕣.k_B*Tmat/(𝕣.m_H*𝕣.not4*𝕣.C*𝕣.C)
            Doppler = 𝕣.C*𝕣.L_He_2p*sqrt(Doppler)
            gamma_2Ps = 3*𝕣.A2P_s*𝕣.fHe*(1-x_He)*𝕣.C*𝕣.C /(
                sqrt(π)*𝕣.sigma_He_2Ps*8π*Doppler*(1-x_H)) /((𝕣.C*𝕣.L_He_2p)^2)
            pb = 0.36 # value from KIV (2007)
            qb = 𝕣.b_He
            # calculate AHcon, the value of A*p_(con,H) for H continuum opacity
            AHcon = 𝕣.A2P_s/(1+pb*(gamma_2Ps^qb))
            K_He = 1/((𝕣.A2P_s*pHe_s+AHcon)*3*n_He*(1-x_He))
	    end
	    if (Heflag >= 3) # include triplet effects
            tauHe_t = 𝕣.A2P_t*n_He*(1. - x_He)*3
            tauHe_t = tauHe_t /(8π*Hz*𝕣.L_He_2Pt^3)
            pHe_t = (1 - exp(-tauHe_t))/tauHe_t
            CL_PSt = 𝕣.h_P*𝕣.C*(𝕣.L_He_2Pt - 𝕣.L_He_2St)/𝕣.k_B
            if ((Heflag == 3) || (Heflag == 5) || (x_H > 0.99999))
                # no H cont. effect
                CfHe_t = 𝕣.A2P_t*pHe_t*exp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t) # "C" factor for triplets
            else # include H cont. effect
                Doppler = 2*𝕣.k_B*Tmat/(𝕣.m_H*𝕣.not4*𝕣.C*𝕣.C)
                Doppler = 𝕣.C*𝕣.L_He_2Pt*sqrt(Doppler)
                gamma_2Pt = (3*𝕣.A2P_t*𝕣.fHe*(1-x_He)*𝕣.C*𝕣.C
                    /(sqrt(π)*𝕣.sigma_He_2Pt*8π*Doppler*(1-x_H))
                    /((𝕣.C*𝕣.L_He_2Pt)^2))
                # use the fitting parameters from KIV (2007) in this case
                pb = 0.66
                qb = 0.9
                AHcon = 𝕣.A2P_t/(1+pb*gamma_2Pt^qb)/3
                CfHe_t = (𝕣.A2P_t*pHe_t+AHcon)*exp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t)  # "C" factor for triplets
            end
	    end
	end

    # Estimates of Thomson scattering time and Hubble time
	timeTh=(1/(𝕣.CT*Trad^4))*(1+x+𝕣.fHe)/x	#!Thomson time
	timeH=2/(3*𝕣.HO*(1+z)^1.5)		#!Hubble time

    # calculate the derivatives
    # turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    # (clunky, but seems to work)
	if (x_H > 0.99)  # don't change at all
		f1 = 0.
    # else if ((x_H.gt.0.98d0).and.(Heflag.eq.0)) then	!don't modify
	elseif (x_H > 0.985)  # !use Saha rate for Hydrogen
		f1 = (x*x_H*n*Rdown - Rup*(1-x_H)*exp(-𝕣.CL/Tmat))/(Hz*(1+z))
        # for interest, calculate the correction factor compared to Saha
        # (without the fudge)
		factor=(1 + K*𝕣.Lambda*n*(1-x_H))/(Hz*(1+z)*(1+K*𝕣.Lambda*n*(1-x)+K*Rup*n*(1-x)))
    else  #!use full rate for H
		f1 = (((x*x_H*n*Rdown - Rup*(1.0-x_H)*exp(-𝕣.CL/Tmat))
			*(1.0 + K*𝕣.Lambda*n*(1.0-x_H)))
		    /(Hz*(1.0+z)*(1.0/𝕣.fu+K*𝕣.Lambda*n*(1.0-x_H)/𝕣.fu
		    +K*Rup*n*(1.0-x_H))))
	end
    # turn off the He once it is small
	if (x_He < 1e-15)
		f2 = 0.
	else
		f2 = (((x*x_He*n*Rdown_He - Rup_He*(1-x_He)*exp(-𝕣.CL_He/Tmat))
            *(1+ K_He*𝕣.Lambda_He*n_He*(1-x_He)*He_Boltz))
            / (Hz*(1+z)
            * (1 + K_He*(𝕣.Lambda_He+Rup_He)*n_He*(1-x_He)*He_Boltz)))
        # Modification to HeI recombination including channel via triplets
	    if (Heflag >= 3)
		    f2 = f2 + (x*x_He*n*Rdown_trip
                - (1-x_He)*3*Rup_trip*exp(-𝕣.h_P*𝕣.C*𝕣.L_He_2St/(𝕣.k_B*Tmat))
                ) * CfHe_t/(Hz*(1+z))
	    end
	end

    # follow the matter temperature once it has a chance of diverging
	if (timeTh < 𝕣.H_frac*timeH)
    # additional term to smooth transition to Tmat evolution,
    # (suggested by Adam Moss)
		epsilon = Hz*(1+x+𝕣.fHe)/(𝕣.CT*Trad^3*x)
		f3 = 𝕣.Tnow + epsilon*((1+𝕣.fHe)/(1+𝕣.fHe+x))*(
            (f1+𝕣.fHe*f2)/x) - epsilon* dHdz/Hz + 3*epsilon/(1+z)
	else
		f3 = 𝕣.CT * (Trad^4) * x / (1+x+𝕣.fHe)* (Tmat-Trad) / (Hz*(1+z)) + 2*Tmat/(1+z)
	end

	return SA[f1, f2, f3]
end


function Xe_He_evolution(𝕣, z, sol)
    y = sol(z, idxs=1)
    rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1/(𝕣.Tnow*(1+z))) / 𝕣.Nnow
    x_H0 = 0.5 * (sqrt(rhs^2+4*rhs) - rhs)
    return x_H0 + 𝕣.fHe * y
end

# get ionisation fraction out of 
function Xe_H_He_evolution(𝕣, z, sol)
    return sol(z, idxs=1) + 𝕣.fHe * sol(z, idxs=2)
end

function ion_recfast!(f, y, 𝕣::RECFAST, z)
    f1, f2, f3 = ion_recfast(y, 𝕣, z)
    f[1] = f1
    f[2] = f2
    f[3] = f3
	return
end


#For reionization - use only the late-time Tmat ode
function late_Tmat(Tm, p, z)
	𝕣,x = p #probably bad to mix types...
	a = 1 / (1+z)
	x_a = a2x(a)
	Hz = 𝕣.bg.ℋ(x_a) / a / H0_natural_unit_conversion
	Trad = 𝕣.Tnow * (1+z)
	dTm = 𝕣.CT * (Trad^4) * x / (1+x+𝕣.fHe)* (Tm-Trad) / (Hz*(1+z)) + 2*Tm/(1+z)
	return dTm
end


"""Tmat for z > 3500"""
Tmat_early(𝕣,z) = 𝕣.Tnow*(1+z)


"""Xe until joint H/He recombination"""
function Xe_early(𝕣, z)
    x0 = 1.0
    if (z > 8000.)
        x0 = 1 + 2*𝕣.fHe
    elseif (z > 5000.)
        rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He2/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
        rhs = rhs*1.  # ratio of g's is 1 for He++ <-> He+
        x0 = 0.5 * (sqrt( (rhs-1-𝕣.fHe)^2 + 4*(1+2𝕣.fHe)*rhs) - (rhs-1-𝕣.fHe) )
    elseif (z > 3500.)
        x0 = 1 + 𝕣.fHe
    else
        # attempt Helium Saha
        rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He1/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
        rhs = rhs*4  # ratio of g's is 4 for He+ <-> He0
        x_He0 = 0.5 * ( sqrt( (rhs-1)^2 + 4*(1+𝕣.fHe)*rhs ) - (rhs-1))
        x0 = x_He0
    end
    return x0
end


# determine redshift at which we have to stop He Saha
function end_of_saha_condition(z, 𝕣)
    rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He1/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
    rhs = rhs*4  # ratio of g's is 4 for He+ <-> He0
    x_He0 = 0.5 * ( sqrt( (rhs-1)^2 + 4*(1+𝕣.fHe)*rhs ) - (rhs-1))
    x0 = x_He0
    x_He0 = (x0 - 1) / 𝕣.fHe
    return x_He0 - 0.99
end


function ion_recfast_H_Saha(u, 𝕣, z)
    rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1/(𝕣.Tnow*(1+z))) / 𝕣.Nnow
    x_H0 = 0.5 * (sqrt(rhs^2+4*rhs) - rhs)
    u = SA[x_H0, u[1], u[2]]
    du = ion_recfast(u, 𝕣, z)
    return SA[du[2], du[3]]
end

function init_He_evolution(𝕣, z0)
    x_H0 = 1.
    rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z0)) - 𝕣.CB1_He1/(𝕣.Tnow*(1+z0)) ) / 𝕣.Nnow
    rhs = rhs*4.  # ratio of g's is 4 for He+ <-> He0
    x_He0 = 0.5 * ( sqrt( (rhs-1)^2 + 4*(1+𝕣.fHe)*rhs ) - (rhs-1))  # He Saha
    x0 = x_He0
    x_He0 = (x0 - 1) / 𝕣.fHe
    return SA[x_He0, Bolt.Tmat_early(𝕣, z0)]
end


function x_H0_H_Saha(𝕣, z)
    rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1/(𝕣.Tnow*(1+z))) / 𝕣.Nnow
    x_H0 = 0.5 * (sqrt(rhs^2+4*rhs) - rhs)
    return x_H0
end

end_He_evo_condition(z, 𝕣) = (x_H0_H_Saha(𝕣, z) - 0.985)



struct RECFASTHistory{T, R, HES, HS}
    𝕣::R
    zinitial::T
    zfinal::T
    z_He_evo_start::T
    z_H_He_evo_start::T
    sol_He::HES
    sol_H_He::HS
end


function Xe(rhist::RECFASTHistory, z) 
    𝕣 = rhist.𝕣
    if (z > 8000.)
        return 1 + 2*𝕣.fHe
    elseif (z > 5000.)
        rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He2/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
        rhs = rhs*1.  # ratio of g's is 1 for He++ <-> He+
        return 0.5 * (sqrt( (rhs-1-𝕣.fHe)^2 + 4*(1+2𝕣.fHe)*rhs) - (rhs-1-𝕣.fHe) )
    elseif (z > 3500.)
        return 1 + 𝕣.fHe
    elseif (z > rhist.z_He_evo_start)
        rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He1/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
        rhs = rhs*4  # ratio of g's is 4 for He+ <-> He0
        x_He0 = 0.5 * ( sqrt( (rhs-1)^2 + 4*(1+𝕣.fHe)*rhs ) - (rhs-1))
        return x_He0
    elseif (z > rhist.z_H_He_evo_start)
        return Bolt.Xe_He_evolution(𝕣, z, rhist.sol_He)
    else
        return Bolt.Xe_H_He_evolution(𝕣, z, rhist.sol_H_He)
    end
end

function Tmat(rhist::RECFASTHistory, z) 
    𝕣 = rhist.𝕣
    if (z > rhist.z_He_evo_start)
        return Bolt.Tmat_early(𝕣, z)
    elseif (z > rhist.z_H_He_evo_start)
        return rhist.sol_He(z, idxs=2)
    else
        return rhist.sol_H_He(z, idxs=3)
    end
end

function recfastsolve(𝕣, alg=Tsit5(), zinitial=10000., zfinal=0.)

    z_epoch_He_Saha_begin = min(zinitial, 3500.)

    # figure out when to start Helium and Hydrogen+Helium evolutions
    z_He_evo_start = assume_nondual(solve(
        IntervalNonlinearProblem(Bolt.end_of_saha_condition, 
        (z_epoch_He_Saha_begin, zfinal), 𝕣), 
        Falsi(), reltol = 1e-4).u)

    z_H_He_evo_start = assume_nondual(solve(
        IntervalNonlinearProblem(Bolt.end_He_evo_condition, 
        (z_He_evo_start, zfinal), 𝕣), 
        Falsi(), reltol = 1e-4).u)

    # evolve Helium and Tmat
    y2 = Bolt.init_He_evolution(𝕣, z_He_evo_start)
    prob2 = ODEProblem{false}(Bolt.ion_recfast_H_Saha, y2, 
        (z_He_evo_start, z_H_He_evo_start), 𝕣)
    sol_He = solve(prob2, alg, reltol=𝕣.tol)

    # evolve Hydrogen, Helium, and Tmat
    z3 = z_H_He_evo_start
    y3 = SA[Bolt.x_H0_H_Saha(𝕣, z3), sol_He(z3, idxs=1), sol_He(z3, idxs=2)]
    prob3 = ODEProblem{false}(Bolt.ion_recfast, y3, (z3, zfinal), 𝕣)
    sol_H_He = solve(prob3, alg, reltol=𝕣.tol)

    return RECFASTHistory(𝕣, zinitial, zfinal, 
        z_He_evo_start, z_H_He_evo_start, sol_He, sol_H_He)
end


function reionization_Xe(rh::RECFASTHistory, z)
    𝕣 = rh.𝕣
    X_fin = 1 + 𝕣.Yp / ( 𝕣.not4*(1-𝕣.Yp) ) #ionization frac today
    zre,α,ΔH,zHe,ΔHe,fHe = 7.6711,1.5,0.5,3.5,0.5,X_fin-1 #reion params, TO REPLACE
    x_orig = Bolt.Xe(rh, z)
    x_reio_H =  (X_fin - x_orig) / 2 * (
        1 + tanh(( (1+zre)^α - (1+z)^α ) / ( α*(1+zre)^(α-1) ) / ΔH)) + x_orig
    x_reio_He = fHe / 2 * ( 1 + tanh( (zHe - z) / ΔHe) )
    x_reio = x_reio_H + x_reio_He
    return x_reio
end


function reionization_Tmat_ode(Tm, rh::RECFASTHistory, z)
    𝕣 = rh.𝕣
    x_reio = reionization_Xe(rh, z)

	a = 1 / (1+z)
	x_a = a2x(a)
	Hz = 𝕣.bg.ℋ(x_a) / a / H0_natural_unit_conversion
	Trad = 𝕣.Tnow * (1+z)
	dTm = 𝕣.CT * Trad^4 * x_reio/(1 + x_reio + 𝕣.fHe) *
        (Tm - Trad) / (Hz * (1 + z)) + 2 * Tm / (1 + z)
	return dTm
end


struct TanhReionizationHistory{T, IH, TS}
    zre_ini::T
    ionization_history::IH
    sol_reionization_Tmat::TS
end


function Xe(trhist::TanhReionizationHistory, z) 
    if (z > trhist.zre_ini)
        return Xe(trhist.ionization_history, z)
    else
        return reionization_Xe(trhist.ionization_history, z)
    end
end

function Tmat(trhist::TanhReionizationHistory, z) 
    if (z > trhist.zre_ini)
        return Tmat(trhist.ionization_history, z)
    else
        return trhist.sol_reionization_Tmat(z)
    end
end

function tanh_reio_solve(ion_hist, zre_ini=50.0)

    𝕣 = ion_hist.𝕣
    zre_ini = 50.0
    reio_prob = ODEProblem(Bolt.reionization_Tmat_ode, 
        Bolt.Tmat(ion_hist, zre_ini), (zre_ini, ion_hist.zfinal), ion_hist)
    sol_reio_Tmat = solve(reio_prob, Tsit5(), reltol=𝕣.tol)
    trh = Bolt.TanhReionizationHistory(zre_ini, ion_hist, sol_reio_Tmat);

end


function old_recfast_xe(𝕣::RECFAST{T};
        Hswitch::Int=1, Heswitch::Int=6, Nz::Int=1000, zinitial=10000., zfinal=0.,
        alg=Tsit5()) where T

    z = zinitial
    n = 𝕣.Nnow * (1 + z)^3
    y = zeros(T,3)  # array is x_H, x_He, Tmat (Hydrogen ionization, Helium ionization, matter temperature)
    y[3] = 𝕣.Tnow * (1 + z)
    Tmat = y[3]

    x_H0, x_He0, x0 = recfast_init(𝕣, z)
    y[1] = x_H0
    y[2] = x_He0

    out_xe = zeros(T, Nz)
    out_Tmat = zeros(T, Nz)

    for i in 1:Nz
        # calculate the start and end redshift for the interval at each z
        # or just at each z
	    zstart = zinitial + float(i-1)*(zfinal-zinitial)/float(Nz)
	    zend   = zinitial + float(i)*(zfinal-zinitial)/float(Nz)

        # Use Saha to get x_e, using the equation for x_e for ionized helium
        # and for neutral helium.
        # Everyb_trip ionized above z=8000.  First ionization over by z=5000.
        # Assume He all singly ionized down to z=3500, then use He Saha until
        # He is 99% singly ionized, and *then* switch to joint H/He recombination.
	    z = zend
        if (zend > 8000.)
            x_H0 = 1.
            x_He0 = 1.
            x0 = 1 + 2*𝕣.fHe
            y[1] = x_H0
            y[2] = x_He0
            y[3] = 𝕣.Tnow*(1+z)
        elseif (z > 5000.)
            x_H0 = 1.
            x_He0 = 1.
            rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He2/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
            rhs = rhs*1.  # ratio of g's is 1 for He++ <-> He+
            x0 = 0.5 * (sqrt( (rhs-1-𝕣.fHe)^2 + 4*(1+2𝕣.fHe)*rhs) - (rhs-1-𝕣.fHe) )
			y[1] = x_H0
            y[2] = x_He0
            y[3] = 𝕣.Tnow*(1+z)
	    elseif (z > 3500.)
            x_H0 = 1.
            x_He0 = 1.
            x0 = x_H0 + 𝕣.fHe*x_He0
            y[1] = x_H0
            y[2] = x_He0
            y[3] = 𝕣.Tnow*(1+z)
        elseif (y[2] > 0.99)
            x_H0 = 1.
            rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1_He1/(𝕣.Tnow*(1+z)) ) / 𝕣.Nnow
            rhs = rhs*4.  # ratio of g's is 4 for He+ <-> He0
            x_He0 = 0.5 * ( sqrt( (rhs-1)^2 + 4*(1+𝕣.fHe)*rhs ) - (rhs-1))  # He Saha
            x0 = x_He0
            x_He0 = (x0 - 1) / 𝕣.fHe
            y[1] = x_H0
            y[2] = x_He0
            y[3] = 𝕣.Tnow*(1+z)
	    elseif (y[1] > 0.99)
            rhs = exp(1.5 * log(𝕣.CR*𝕣.Tnow/(1+z)) - 𝕣.CB1/(𝕣.Tnow*(1+z))) / 𝕣.Nnow
            x_H0 = 0.5 * (sqrt(rhs^2+4*rhs) - rhs)

            prob = ODEProblem{true}(ion_recfast!, y, (zstart, zend), 𝕣,
                save_everystep=false, save_start=false)
            sol = solve(prob, alg, reltol=𝕣.tol)
            y .= sol(zend)

            y[1] = x_H0
            x0 = y[1] + 𝕣.fHe*y[2]
	    else
            prob = ODEProblem{true}(ion_recfast!, y, (zstart, zend), 𝕣,
                save_everystep=false, save_start=false)
            sol = solve(prob, alg, reltol=𝕣.tol)
            y .= sol(zend)
            x0 = y[1] + 𝕣.fHe*y[2]
        end

        Trad = 𝕣.Tnow * (1+zend)
        Tmat = y[3]
        x_H = y[1]
        x_He = y[2]
        x = x0

        out_xe[i] = x
        out_Tmat[i] = Tmat

    end

	#Reionization
	# We need to
	# 0. Solve for x as usual at all z (as above)
	# 1. Apply the by-hand CAMB-style reionization function to obtain x_re from x, overwrite x (then τ is corrected)
	#FIXME I did this fast so it is probably not done well...
	X_fin = 1 + 𝕣.Yp / ( 𝕣.not4*(1-𝕣.Yp) ) #ionization frac today
	zre_ini,zre,α,ΔH,zHe,ΔHe,fHe = 50.,7.6711,1.5,0.5,3.5,0.5,X_fin-1 #reion params
	# 2. Feed x_re(z) and Tm(zre_ini) into dTm/dz - essentially re-solving for y[3] at z<zre_ini
	idx_zre_start = argmin( abs.(collect(zinitial:(zfinal-zinitial)/float(Nz):zfinal) .- zre_ini) ) #find first z s.t. z>=zre_ini
	Tmat = out_Tmat[idx_zre_start] #set initial Tmat, update in place as we go after
	for i in idx_zre_start:Nz
		# calculate the start and end redshift for the interval at each z
		# or just at each z
		zstart = zinitial + float(i-1)*(zfinal-zinitial)/float(Nz)
		zend   = zinitial + float(i)*(zfinal-zinitial)/float(Nz)
		z = zend

		#load non-reion values
		x = out_xe[i]

		#apply reionization tanh functions
		x_reio_H = (X_fin - x) / 2 *( 1 +
						tanh(( (1+zre)^α - (1+z)^α ) / ( α*(1+zre)^(α-1) ) / ΔH)
												 ) + x
		x_reio_He = fHe / 2 * ( 1 + tanh( (zHe - z) / ΔHe) )
		x_reio = (z > zre_ini) ? x : x_reio_H + x_reio_He
		out_xe[i] = x_reio #overwrite non-reion piece

		# use only the above "else" case for Tmat part, since t_T >> t_H after z_re_ini
		prob = ODEProblem(late_Tmat, Tmat, (zstart, zend), [𝕣 x_reio],
			save_everystep=false, save_start=false)
		sol = solve(prob, alg, reltol=𝕣.tol)
		# 3. Overwrite the Tm array before passing to csb function
		Tmat = sol(zend)
		out_Tmat[i] = Tmat
	end

    return out_xe, out_Tmat
end

RECFASTredshifts(Nz, zinitial, zfinal) =
    range(zinitial, stop=zfinal, length=Nz+1)[2:end]

function IonizationHistory(𝕣::RECFAST{T}, par::AbstractCosmoParams{T}, bg::AbstractBackground{T}) where
                           T#{T, ACP<:AbstractCosmoParams, AB<:AbstractBackground}
    x_grid = bg.x_grid
    # GRAFT RECFAST ONTO BOLT. TODO: MEGA-REFACTOR ==============
    # Nz = 100000 #add extra two zero for reion otherwise too low res, get wiggles (was spacing of Δz=1)
	#FIXME treat these wiggles better!
	# Xe_RECFAST, Tmat_RECFAST = recfast_xe(𝕣; Nz=Nz, zinitial=10000., zfinal=0.)
	# z_RECFAST = RECFASTredshifts(Nz, 10000., 0.)
    # RECFAST_Xₑ_z = spline(reverse(Xe_RECFAST), reverse(z_RECFAST))
    # RECFAST_Tmat_z = spline(reverse(Tmat_RECFAST), reverse(z_RECFAST))

    rhist = recfastsolve(𝕣)
    trhist = tanh_reio_solve(rhist)
    

    xinitial_RECFAST = z2x(rhist.zinitial)
    Xe_initial = Xe(rhist, rhist.zinitial)

    Xₑ_function = x -> (x < xinitial_RECFAST) ?
        Xe_initial : Xe(trhist, x2z(x))
    Trad_function = x -> 𝕣.Tnow * (1 + x2z(x))
    Tmat_function = x -> (x < xinitial_RECFAST) ?
        Trad_function(x) : Tmat(trhist, x2z(x))

    # =====================================================
	#j - do we really need bg to be passed to IonizationHistory separately from 𝕣.bg?
	#is there a reason not to just put par and bg into 𝕣?
	ℋ_function = bg.ℋ
    τ, τ′ = τ_functions(x_grid, Xₑ_function, par, ℋ_function)
    g̃ = g̃_function(τ, τ′)

    Xₑ_ = spline(Xₑ_function.(x_grid), x_grid)
    τ_ = spline(τ.(x_grid), x_grid)
    g̃_ = spline(g̃.(x_grid), x_grid)
    Tmat_ = spline(Tmat_function.(x_grid), x_grid)
	#sound speed
	csb²_pre = @.( 𝕣.C^-2 * 𝕣.k_B/𝕣.m_H * ( 1/𝕣.mu_T + (1-𝕣.Yp)*Xₑ_(x_grid) ) ) #not the most readable...
	#FIXME probably this is a bad way to do this...
	csb²_ = spline(csb²_pre .* (Tmat_.(x_grid) .- 1/3 *spline_∂ₓ(Tmat_, x_grid).(x_grid)),x_grid)

    return IonizationHistory(
		T(τ(0.)),
        Xₑ_,
        τ_,
        spline_∂ₓ(τ_, x_grid),
        spline_∂ₓ²(τ_, x_grid),
        g̃_,
        spline_∂ₓ(g̃_, x_grid),
        spline_∂ₓ²(g̃_, x_grid),
        Tmat_,
		csb²_,
    )
end
