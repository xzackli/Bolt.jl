using Bolt
include("../test/deps/deps.jl")

##
"""
Wrapper of RECFAST Fortran code with parameters as defined in that code.
Returns tuple of (z's, xe's)
"""
function get_xe(OmegaB::Float64, OmegaC::Float64, OmegaL::Float64,
                HOinp::Float64, Tnow::Float64, Yp::Float64;
                Hswitch::Int64=1, Heswitch::Int64=6,
                Nz::Int64=1000, zstart::Float64=10000., zend::Float64=0.)

    xe = Array{Float64}(undef,Nz)
    ccall(
        (:get_xe_, librecfast), Nothing,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64},
         Ref{Int64}, Ref{Int64}, Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        OmegaB, OmegaC, OmegaL, HOinp, Tnow, Yp, Hswitch, Heswitch, Nz, zstart, zend, xe
    )
    range(zstart,stop=zend,length=Nz+1)[2:end], xe
end


p = CosmoParams(Y_p = 0.24)
bg = Background(p)
T_cmb = 2.7255
z, xedat = get_xe(p.Ω_b, p.Ω_m, Bolt.Ω_Λ(p), p.h * 100, T_cmb, p.Y_p)

# clf()
# plt.plot(z, xedat, "-")
# xscale("log")
# gcf()


##



##
# using Parameters
# @with_kw struct RECFASTParams @deftype Float64
bigH = 100.0e3 / (1e6 * 3.0856775807e16)	 # H₀ in s-1
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
# end

##
# Gaussian fits for extra H physics (fit by Adam Moss, modified by
# Antony Lewis)

OmegaB = p.Ω_b
OmegaC = p.Ω_m
OmegaL = Bolt.Ω_Λ(p)
HOinp = p.h * 100
Tnow = T_cmb
Yp = p.Y_p
Hswitch=1
Heswitch=6
Nz=1000
zinitial=10000.
zfinal=0.

z = zinitial
OmegaT = OmegaC + OmegaB            # total dark matter + baryons
OmegaK = 1. - OmegaT - OmegaL	    # curvature

# convert the Hubble constant units
H = HOinp/100
HO = H*bigH

# sort out the helium abundance parameters
mu_H = 1 / (1 - Yp)			 # Mass per H atom
mu_T = not4/(not4-(not4-1)*Yp)	 # Mass per atom
fHe = Yp/(not4*(1 - Yp))		# n_He_tot / n_H_tot

Nnow = 3 * HO * HO * OmegaB / (8π * G * mu_H * m_H)
n = Nnow * (1 + z)^3
fnu = (21/8)*(4/11)^(4/3)
# 	(this is explictly for 3 massless neutrinos - change if N_nu.ne.3)  # TODO: should look into this
z_eq = (3 * (HO*C)^2 / (8π * G * a * (1+fnu)*Tnow^4))*OmegaT
z_eq = z_eq - 1

# 	Set up some constants so they don't have to be calculated later
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

#	Matter departs from radiation when t(Th) > H_frac * t(H)
#	choose some safely small number
H_frac = 1e-3

#	Fudge factor to approximate the low z out of equilibrium effect
if (Hswitch == 0)
    fu=1.14e0
else
    fu=1.125e0
end

b_He = 0.86  # Set the He fudge factor
y = zeros(3)  # array is x_H, x_He, Tmat (Hydrogen ionization, Helium ionization, matter temperature)
# y[3] = Tnow * (1 + z)
# x_H0, x_He0 = get_init(z, x0)
# y[1] = x_H0
# y[2] = x_He0

##
function get_init(z)
    x_H0, x_He0, x0 = [0.0], [0.0], [0.0]
    ccall(
        (:get_init_, librecfast), Nothing,
        (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        z, x_H0, x_He0, x0
    )
    return x_H0[1], x_He0[1], x0[1]
end

function recfast_init(z)
    if z > 8000.
        x_H0 = 1.
        x_He0 = 1.
        x0 = 1. + 2fHe
    elseif z > 3500.
        x_H0 = 1.
        x_He0 = 1.
        rhs = exp( 1.5 * log(CR*Tnow/(1. + z)) - CB1_He2/(Tnow*(1. + z)) ) / Nnow
	    rhs = rhs * 1.  #ratio of g's is 1 for He++ <-> He+
	    x0 = 0.5 * ( sqrt( (rhs - 1. - fHe)^2 + 4. * (1. + 2fHe)*rhs) - (rhs - 1. - fHe) )
    elseif z > 2000
	    x_H0 = 1.
	    rhs = exp( 1.5 * log(CR*Tnow/(1. + z)) - CB1_He1/(Tnow*(1. + z)) ) / Nnow
	    rhs = rhs*4.    # ratio of g's is 4 for He+ <-> He0
	    x_He0 = 0.5 * ( sqrt( (rhs-1.)^2 + 4*(1. + fHe)*rhs ) - (rhs-1.))
	    x0 = x_He0
	    x_He0 = (x0 - 1.)/fHe
    else
	    rhs = exp( 1.5 * log(CR*Tnow/(1. + z)) - CB1/(Tnow*(1. + z)) ) / Nnow
	    x_H0 = 0.5 * (sqrt( rhs^2 + 4 * rhs ) - rhs )
	    x_He0 = 0.
	    x0 = x_H0
    end

    return x_H0, x_He0, x0
end
using Test
@test all(get_init(9000.0) .≈ recfast_init(9000.0))
@test all(get_init(4000.0) .≈ recfast_init(4000.0))
@test all(get_init(3000.0) .≈ recfast_init(3000.0))
@test all(get_init(1000.0) .≈ recfast_init(1000.0))
@test all(get_init(500.0) .≈ recfast_init(500.0))
@test all(get_init(100.0) .≈ recfast_init(100.0))

##

function get_ion(z, y)
    # x_H0, x_He0, x0 = [0.0], [0.0], [0.0]
    Ndim = 3
    f = zeros(Ndim)
    ccall(
        (:ion_, librecfast), Nothing,
        (Ref{Int64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        Ndim, z, y, f
    )
    return f
end

function ion_recfast(z, y, f)
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

	x_H = y[1]
	x_He = y[2]
	x = x_H + fHe * x_He
	Tmat = y[3]

	n = Nnow * (1+z)^3
	n_He = fHe * Nnow * (1+z)^3
	Trad = Tnow * (1+z)
	Hz = HO * sqrt((1+z)^4/(1+z_eq)*OmegaT + OmegaT*(1+z)^3 + OmegaK*(1+z)^2 + OmegaL)

    # Also calculate derivative for use later
	dHdz = (HO^2 /2/Hz)*(4*(1+z)^3/(1+z_eq)*OmegaT + 3*OmegaT*(1+z)^2 + 2*OmegaK*(1+z) )

    # Get the radiative rates using PPQ fit (identical to Hummer's table)
	Rdown=1e-19*a_PPB*(Tmat/1e4)^b_PPB/(1. + c_PPB*(Tmat/1e4)^d_PPB)
	Rup = Rdown * (CR*Tmat)^(1.5)*exp(-CDB/Tmat)

    # calculate He using a fit to a Verner & Ferland type formula
	sq_0 = sqrt(Tmat/T_0)
	sq_1 = sqrt(Tmat/T_1)
    # typo here corrected by Wayne Hu and Savita Gahlaut
	Rdown_He = a_VF/(sq_0*(1+sq_0)^(1-b_VF))
	Rdown_He = Rdown_He/(1+sq_1)^(1+b_VF)
	Rup_He = Rdown_He*(CR*Tmat)^(1.5)*exp(-CDB_He/Tmat)
	Rup_He = 4. * Rup_He # statistical weights factor for HeI
    # Avoid overflow (pointed out by Jacques Roland)
	if((Bfact/Tmat) > 680.)
	  He_Boltz = exp(680.)
	else
	  He_Boltz = exp(Bfact/Tmat)
	end

    # now deal with H and its fudges
	if (Hswitch == 0)
	    K = CK/Hz # !Peebles coefficient K=lambda_a^3/8piH
	else
        # fit a double Gaussian correction function
        K = CK/Hz*(1.0
            + AGauss1*exp(-((log(1+z)-zGauss1)/wGauss1)^2)
            + AGauss2*exp(-((log(1+z)-zGauss2)/wGauss2)^2))
	end

    # add the HeI part, using same T_0 and T_1 values
	Rdown_trip = a_trip/(sq_0*(1+sq_0)^(1-b_trip))
	Rdown_trip = Rdown_trip/((1+sq_1)^(1+b_trip))
	Rup_trip = Rdown_trip*exp(-h_P*C*L_He2St_ion/(k_B*Tmat))
	Rup_trip = Rup_trip*((CR*Tmat)^1.5)*(4/3)
    # last factor here is the statistical weight

    # try to avoid "NaN" when x_He gets too small
	if ((x_He < 5.e-9) || (x_He > 0.980))
        Heflag = 0
	else
	    Heflag = Heswitch
	end
	if (Heflag == 0)  # use Peebles coeff. for He
	    K_He = CK_He/Hz
	else # for Heflag>0 		!use Sobolev escape probability
        tauHe_s = A2P_s*CK_He*3*n_He*(1-x_He)/Hz
        pHe_s = (1 - exp(-tauHe_s))/tauHe_s
        K_He = 1 / (A2P_s*pHe_s*3*n_He*(1-x_He))
        # smoother criterion here from Antony Lewis & Chad Fendt
	    if (((Heflag == 2) || (Heflag >= 5)) && (x_H < 0.9999999))
            # use fitting formula for continuum opacity of H
            # first get the Doppler width parameter
            Doppler = 2*k_B*Tmat/(m_H*not4*C*C)
            Doppler = C*L_He_2p*sqrt(Doppler)
            gamma_2Ps = 3*A2P_s*fHe*(1-x_He)*C*C /(
                sqrt(π)*sigma_He_2Ps*8π*Doppler*(1-x_H)) /((C*L_He_2p)^2)
            pb = 0.36 # value from KIV (2007)
            qb = b_He
            # calculate AHcon, the value of A*p_(con,H) for H continuum opacity
            AHcon = A2P_s/(1+pb*(gamma_2Ps^qb))
            K_He = 1/((A2P_s*pHe_s+AHcon)*3*n_He*(1-x_He))
	    end
	    if (Heflag >= 3) # include triplet effects
            tauHe_t = A2P_t*n_He*(1. - x_He)*3
            tauHe_t = tauHe_t /(8π*Hz*L_He_2Pt^3)
            pHe_t = (1 - exp(-tauHe_t))/tauHe_t
            CL_PSt = h_P*C*(L_He_2Pt - L_He_2st)/k_B
            if ((Heflag == 3) || (Heflag == 5) || (x_H > 0.99999))
                # no H cont. effect
                CfHe_t = A2P_t*pHe_t*exp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t) # "C" factor for triplets
            else # include H cont. effect
                Doppler = 2*k_B*Tmat/(m_H*not4*C*C)
                Doppler = C*L_He_2Pt*sqrt(Doppler)
                gamma_2Pt = (3*A2P_t*fHe*(1-x_He)*C*C
                    /(sqrt(π)*sigma_He_2Pt*8π*Doppler*(1-x_H))
                    /((C*L_He_2Pt)^2))
                # use the fitting parameters from KIV (2007) in this case
                pb = 0.66
                qb = 0.9
                AHcon = A2P_t/(1+pb*gamma_2Pt^qb)/3
                CfHe_t = (A2P_t*pHe_t+AHcon)*exp(-CL_PSt/Tmat)
                CfHe_t = CfHe_t/(Rup_trip+CfHe_t)  # "C" factor for triplets
            end
	    end
	end

    # Estimates of Thomson scattering time and Hubble time
	timeTh=(1/(CT*Trad^4))*(1+x+fHe)/x	#!Thomson time
	timeH=2/(3*HO*(1+z)^1.5)		#!Hubble time

    # calculate the derivatives
    # turn on H only for x_H<0.99, and use Saha derivative for 0.98<x_H<0.99
    # (clunky, but seems to work)
	if (x_H > 0.99)  # don't change at all
		f[1] = 0.
    # else if ((x_H.gt.0.98d0).and.(Heflag.eq.0)) then	!don't modify
	elseif (x_H > 0.985)  # !use Saha rate for Hydrogen
		f[1] = (x*x_H*n*Rdown - Rup*(1-x_H)*exp(-CL/Tmat))/(Hz*(1+z))
        # for interest, calculate the correction factor compared to Saha
        # (without the fudge)
		factor=(1 + K*Lambda*n*(1-x_H))/(Hz*(1+z)*(1+K*Lambda*n*(1-x)+K*Rup*n*(1-x)))
    else  #!use full rate for H
		f[1] = (((x*x_H*n*Rdown - Rup*(1.0-x_H)*exp(-CL/Tmat))
			*(1.0 + K*Lambda*n*(1.0-x_H)))
		    /(Hz*(1.0+z)*(1.0/fu+K*Lambda*n*(1.0-x_H)/fu
		    +K*Rup*n*(1.0-x_H))))
	end
    # turn off the He once it is small
	if (x_He < 1e-15)
		f[2] = 0.
	else
		f[2] = (((x*x_He*n*Rdown_He - Rup_He*(1-x_He)*exp(-CL_He/Tmat))
            *(1+ K_He*Lambda_He*n_He*(1-x_He)*He_Boltz))
            / (Hz*(1+z)
            * (1 + K_He*(Lambda_He+Rup_He)*n_He*(1-x_He)*He_Boltz)))
        # Modification to HeI recombination including channel via triplets
	    if (Heflag >= 3)
		    f[2] = f[2] + (x*x_He*n*Rdown_trip
                - (1-x_He)*3*Rup_trip*exp(-h_P*C*L_He_2st/(k_B*Tmat))
                ) * CfHe_t/(Hz*(1+z))
	    end
	end

    # follow the matter temperature once it has a chance of diverging

	if (timeTh < H_frac*timeH)
    # f(3)=Tmat/(1.d0+z)	!Tmat follows Trad
    # additional term to smooth transition to Tmat evolution,
    # (suggested by Adam Moss)
		epsilon = Hz*(1+x+fHe)/(CT*Trad^3*x)
		f[3] = Tnow + epsilon*((1+fHe)/(1+fHe+x))*(
            (f[1]+fHe*f[2])/x) - epsilon* dHdz/Hz + 3*epsilon/(1+z)
	else
		f[3] = CT * (Trad^4) * x / (1+x+fHe)* (Tmat-Trad) / (Hz*(1+z)) + 2*Tmat/(1+z)
	end

	return
end


z_TEST = 1400.0
x_H0, x_He0, x0 = recfast_init(z_TEST)
# print(x_H0, "\n")
get_ion(z_TEST, [x_H0, x_He0, Tnow * (1+z_TEST)] )
##

function test_fort(z)
    f_TEST = get_ion(z, [x_H0, x_He0, Tnow * (1+z_TEST)] )
    return f_TEST[1]
end

##
function test(z)
    f_TEST = zeros(3)
    ion_recfast(z, [x_H0, x_He0, Tnow * (1+z_TEST)], f_TEST)
    return f_TEST[1]
end

clf()
plot([abs(test_fort(z)) for z in 10:100:2000])
plot([abs(test(z)) for z in 10:100:2000])
yscale("log")
gcf()

##
clf()
plot([test_fort(z) ./ test(z) for z in 10:100:2000])
# plot([-test(z) for z in 10:100:2000])
# yscale("log")
gcf()

##
