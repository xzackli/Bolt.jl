using Bolt
include("../test/deps/deps.jl")

using OrdinaryDiffEq
using Parameters

𝕣 = RECFASTIonization()

# OmegaB = p.Ω_b
# OmegaC = p.Ω_m
# OmegaL = Bolt.Ω_Λ(p)
# HOinp = p.h * 100
# Tnow = T_cmb
# Yp = p.Y_p
# Hswitch=1
# Heswitch=6
# Nz=1000
# zinitial=10000.
# zfinal=0.

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

z, xedat = get_xe(𝕣.OmegaB, 𝕣.OmegaC, 𝕣.OmegaL, 𝕣.HOinp, 𝕣.Tnow, 𝕣.Yp)

# clf()
# plt.plot(z, xedat, "-")
# xscale("log")
# gcf()


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

using Test
@test all(get_init(9000.0) .≈ recfast_init(𝕣, 9000.0))
@test all(get_init(4000.0) .≈ recfast_init(𝕣, 4000.0))
@test all(get_init(3000.0) .≈ recfast_init(𝕣, 3000.0))
@test all(get_init(1000.0) .≈ recfast_init(𝕣, 1000.0))
@test all(get_init(500.0) .≈ recfast_init(𝕣, 500.0))
@test all(get_init(100.0) .≈ recfast_init(𝕣, 100.0))

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



z_TEST = 1400.0
x_H0, x_He0, x0 = recfast_init(𝕣, z_TEST)
f_TEST = zeros(3)
# ion_recfast!(𝕣, z_TEST, [x_H0, x_He0, 𝕣.Tnow * (1+z_TEST)], f_TEST)

# print(x_H0, "\n")
# get_ion(z_TEST, [x_H0, x_He0, 𝕣.Tnow * (1+z_TEST)] )
##

function test_fort(z)
    f_TEST = get_ion(z, [x_H0, x_He0, 𝕣.Tnow * (1+z_TEST)] )
    return f_TEST[1]
end

function test(z)
    f_TEST = zeros(3)
    ion_recfast!(f_TEST, [x_H0, x_He0, 𝕣.Tnow * (1+z_TEST)], 𝕣, z)
    return f_TEST[1]
end

clf()
plot([abs(test_fort(z)) for z in 10:40:2000])
plot([abs(test(z)) for z in 10:40:2000], "-")
yscale("log")
gcf()

##
clf()
plot([test_fort(z) ./ test(z) for z in 10:100:2000])
# ylim(1-1e-6, 1+1e-6)
# plot([-test(z) for z in 10:100:2000])
# yscale("log")
gcf()

##

##
@time xe_bespoke = recfast_xe(𝕣);

##
clf()
plot(z, xedat, "-")
plot(z, xe_bespoke, "--")
# ylim(0.95, 1.01)
# xlim(1500,1600)
gcf()

##
clf()
plot(z, xe_bespoke ./ xedat, "-")
# xlim(0, 2000)

ylim(1 - 0.01^2, 1 + 0.01^2)
# plot(z, xe_bespoke, "--")
gcf()


##
