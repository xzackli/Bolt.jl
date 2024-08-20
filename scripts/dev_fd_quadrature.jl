# Script to get a feeling for speed vs accuracy of quadrature schemes
using QuadGK,FastGaussQuadrature
using Plots

# Get the "truth" function
# bare FD distribution (for decoupled neutrinos the form of the bare FD doesn't change with time)
# The analytic value of the density integral at early times is 7/8 π²/15 * Tν^4
# FD function multiplied by a reasonable pert (use the massive nu ICs)
# For "truth" use high-resolution integration
# Choose qmin, qmax, n_q

# Compare tolerances

#One-liner to map finite interval to unit interval (and back)
to_ui(lq,lqmi,lqma) = -1 + (1- (-1)) / (lqma-lqmi) * (lq-lqmi)
from_ui(x,lqmi,lqma) = lqmi + (lqma- lqmi) / (1- (-1)) * (x- (-1))

#how to choose this seems a bit up in the air, but CAMB uses qmax of 30T
T=2e-4
logqmin,logqmax=log10(T/30),log10(T*30)#1e-6,1e-1
dxdq(q) = (1+to_ui(1+logqmin,logqmin,logqmax))/(q*log(10)) #chovar factor, dx/dq=dx/dlogq dlogq/dq

plot()
f0(q) = 4π *2(2π)^-3 *q^2 / (1 + exp(q/T))
ζ=1.202056903159594
answer = 3 *T^3 *ζ / (2π^2)
logqq=collect(logqmin:.01:logqmax)
plot(logqq,f0.(10.0 .^logqq))
plot(to_ui.(logqq,logqmin,logqmax),
     f0.(10.0 .^from_ui.(to_ui.(logqq,logqmin,logqmax),logqmin,logqmax)))
ft(x) = f0(10.0 ^ from_ui(x,logqmin,logqmax)) / dxdq(10.0 ^ from_ui(x,logqmin,logqmax))
plot!(to_ui.(logqq,logqmin,logqmax),ft.(to_ui.(logqq,logqmin,logqmax)))


#QuadGK
# The lazy thing I was doing before - just run quadgk at low tolerance
@time Iqgk = quadgk(q ->  f0(q) ,10.0 .^logqmin, 10.0 .^logqmax,rtol=1e-6)[1]
println("Relative error on quadGK: ", abs.(1-Iqgk/answer))
#this is already very slow - it takes .73 seconds at rtol=1e-2
#going to 1e-6 (which seems about what fast quadrature is getting), gives ~.26 s
#Should get average timings b/c seeing big difference from first call
#Comparison with FastGaussQuadrature should take into account saved quantities.

#FastGaussQuadrature.jl routines
quad_sum(x,w,f) = sum(f.(x).*w)

# Legendre quadrature
n=10
@time xGLe, wGLe = gausslegendre( n )
@time IGLe = quad_sum(xGLe,wGLe,ft)
println("Relative error on GLeg: ", abs.(1-IGLe/answer), "  with n= ", n," points.")
#just from changing n it looks like rel error bottoms out at 1e-11
#^I guess this is the error we induce by not integrating the tails - truly negligible for bg...
#with 10 pts using the T/30 30T range get 1e-3 error, will do this for testing

# CC quadrature
ncc = n+3
@time xGCC, wGCC = gausslobatto( ncc )
@time IGCC = quad_sum(xGCC,wGCC,ft)
println("Relative error on CC: ", abs.(1-IGCC/answer), "  with n= ", ncc," points.")
#30 pts is pretty good - relative tolerance of 1e-6 for Leg, 1e-4 for cheby zeros

#FIXME: try using asymptotic expansion instead of just truncating the integral (standard practice)

#FIXME: Need to take into account possible integration over splines instead of f
#If integrate over spline is faster but probably introduces interp error, especially at ends
#Also we won't spline the perts - but should check this though for bg

# More realistic test on perts
#Look at a the result of the un-normalized density integral at z=0 using dlnf0dlnq to test q-dep
#"truth" as a higher-order low-tolerance quadgk over larger q-range
dlnf0dlnq0(q) = -q / T /(1 + exp(-q/T))

#high-tolerance QuadGK
a=1#10^-4
@time Iqgk_truth = quadgk(q ->  √(q^2 + a^2 * 0.06^2)* q^2 * f0(q)*dlnf0dlnq0(q) ,10.0 .^logqmin, 10.0 .^logqmax,rtol=1e-12)[1]

#fast QuadGK at 1e-3
@time Iqgk_fast = quadgk(q ->  √(q^2 + a^2 * 0.06^2)* q^2 * f0(q)*dlnf0dlnq0(q) ,10.0 .^logqmin, 10.0 .^logqmax,rtol=1e-3)[1]
println("Relative error on realistic 1e-3 gk: ", abs.(1-Iqgk_fast/Iqgk_truth))


#vs FastGaussQuadrature
function ftp(x)
    return (√((10.0 ^ from_ui(x,logqmin,logqmax))^2 + a^2 * 0.06^2)
            * f0(10.0 ^ from_ui(x,logqmin,logqmax)) * (10.0 ^ from_ui(x,logqmin,logqmax))^2
            *dlnf0dlnq0(10.0 ^ from_ui(x,logqmin,logqmax)) )/ dxdq(10.0 ^ from_ui(x,logqmin,logqmax))
end

#Gauss-Legendre
m = 15
@time xGLe2, wGLe2 = gausslegendre( m )
@time IGLe2 = quad_sum(xGLe2,wGLe2,ftp)
println("Relative error on realistic GLeg: ", abs.(1-IGLe2/Iqgk_truth), "  with m= ", m," points.")

#CC
mcc = m#+1
@time xGCC2, wGCC2 = gausslobatto( mcc )
@time IGCC2 = quad_sum(xGCC2,wGCC2,ftp)
println("Relative error on CC: ", abs.(1-IGCC2/Iqgk_truth), "  with n= ", mcc," points.")

#comparing to the f0 test, where we got 1e-3 error for 10 and 13 pts for GL and CC respectively, here we have
#12 and 13 pts, respectively, and similarly, we have that it is faster
#Need to do a careful comparison at some point where we can average over a number of times,
#as it is hard to know when the quad points are reusing p/re-computed values...
#Either way so long as the computation of quadrature points (which occurs in the background once)
#plus the sum operations that happen frequently are faster than quadgk * number of times it is called,
#then we are fine - and I think this is an easily attainable threshold for a small number of quad pts (<15)

#Would be interesting to see if cosmology affects this accuracy, because could reduce cost even further by
#just passing the q-pts (computed externally) once instead of computing in bg

#hmm I forgot a factor of q before, but even at higher redshift the performance
#is about the same on 15 pts, only around 1e-3 rerr

#------------------------------
#Check the background integrals that used to be quadgk against fixed gauss quadrature
#the old quadgk code - this is very lazy copying...
const km_s_Mpc_100 = 2.1331196424576403e-33  # [eV]
const G_natural = 6.708830858490363e-57
const ζ = 1.2020569 #Riemann ζ(3) for phase space integrals
function ρPold(a,par::AbstractCosmoParams)
    #Background phase space energy density and pressure
    Tν =  (par.N_ν/3)^(1/4) * (4/11)^(1/3) * (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    m = par.Σm_ν
    qmin=1e-18
    qmax=1e1
    ρ = 4π * a^(-4) * quadgk(q ->  q^2 * √( q^2 + (a*m)^2 ) * f0(q) ,qmin, qmax,rtol=1e-6)[1]
    P = 4π/3 * a^(-4) * quadgk(q -> q^2 * q^2 /√( q^2 + (a*m)^2) * f0(q), qmin, qmax,rtol=1e-6)[1]
    return ρ,P#,norm
end
H₀(par::AbstractCosmoParams) = par.h * km_s_Mpc_100
ρ_crit(par::AbstractCosmoParams) = (3 / 8π) * H₀(par)^2 / G_natural  # [eV⁴]
function Ω_Λ(par::AbstractCosmoParams)
    #Below can definitely be more streamlined, I am just making it work for now
    Tγ = (15/ π^2 *ρ_crit(par) *par.Ω_r)^(1/4)
    νfac = (90 * ζ /(11 * π^4)) * (par.Ω_r * par.h^2 / Tγ) *((par.N_ν/3)^(3/4))
    #^the factor that goes into nr approx to neutrino energy density, plus equal sharing ΔN_eff factor for single massive neutrino
    Ω_ν = par.Σm_ν*νfac/par.h^2
    return 1 - (par.Ω_r*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3))  # dark energy density
                                         + par.Ω_b + par.Ω_m
                                         + Ω_ν) #assume massive nus are non-rel today
end
function H_aold(a, par::AbstractCosmoParams)
    ρ_ν,_ = ρPold(a,par) # we don't atually need pressure?
    return H₀(par) * √((par.Ω_m + par.Ω_b ) * a^(-3)
                        + ρ_ν/ρ_crit(par)
                        + par.Ω_r*(1+(2/3)*(7par.N_ν/8)*(4/11)^(4/3)) * a^(-4)
                        + Ω_Λ(par))
end
ℋ_aold(a, par::AbstractCosmoParams) = a * H_aold(a, par)
ℋold(x, par::AbstractCosmoParams) = ℋ_aold(x2a(x), par)

function ηold(x, par::AbstractCosmoParams)
    return quadgk(a -> 1.0 / (a * ℋ_aold(a, par)), 0.0, x2a(x),rtol=1e-6)[1]
end

using Bolt
using Interpolations
𝕡 = CosmoParams()
x_grid=-20.0:0.1:0.0
bg = Background(𝕡;x_grid=x_grid,nq=n_q)
spline(f, x_grid) = scale(interpolate(f, BSpline(Cubic(Line(OnGrid())))), x_grid)
oldℋ_  = spline([ℋold(x, par) for x in x_grid], x_grid)

xmin,xmax=-20,0
xx = collect(xmin:.01:xmax)
#aa = exp.(xx)
#first look at just conformal ℋ, this is just a test of ρ
spline_fixed_ℋ =bg.ℋ(xx)
spline_quadgk_ℋ =oldℋ_(xx)
plot(xx,log10.(spline_fixed_ℋ))
plot!(xx,log10.(spline_quadgk_ℋ))
plot(xx,log10.(abs.(spline_fixed_ℋ ./ spline_quadgk_ℋ .- 1)))
plot(xx,abs.(spline_fixed_ℋ ./ spline_quadgk_ℋ .- 1))
#difference is always less than 4.5e-6

#conformal time integrand
plot(xx,log10.(exp.(xx) ./ spline_fixed_ℋ))
plot!(xx,log10.(exp.(xx) ./ spline_quadgk_ℋ))

#okay lets look at conformal spline time now
@time spline_fixed_η = bg.η(xx)
oldη_  = spline([ηold(x, 𝕡) for x in x_grid], x_grid)
@time spline_quadgk_η =oldη_(xx)
plot(xx,log10.(spline_fixed_η))
plot!(xx,log10.(spline_quadgk_η))
plot(xx,abs.(spline_fixed_η ./ spline_quadgk_η .- 1))
plot(xx,log10.(abs.(spline_fixed_η ./ spline_quadgk_η .- 1)))
#hmm always below 5e-4 which is what I would've expected from earlier accuracies with 15 pts
#EXCEPT for x<-15 in which case error rockets up to 5e-2
#This goes away when I drop the lower limit in quad in η on x from -10 to -15, accuracy slightly degrades at zero
#Going to minimum value of -11.75 gives error always less than 1e-3 - and we should not have expected
#better than this given the earlier tests, plus going to smaller amin
#results in larger error at x=0, so I will leave it here for now
