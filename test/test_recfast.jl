using Bolt
include("../test/deps/deps.jl")

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

clf()
plt.plot(z, xedat, "-")
# xscale("log")
gcf()


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

z_TEST = 1500.0
x_H0, x_He0, x0 = recfast_init(z_TEST)
print(x_H0, "\n")
get_ion(z_TEST, [x_H0, x_He0, Tnow * (1+z_TEST)] )


##

##
