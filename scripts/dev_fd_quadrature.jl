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
#Either way, fast quad is at least 3 orders of magnitude faster
#Seems cheby is faster than Legendre by ~50x - but is only accurate at ~1%
#(which makes sense, lobatto pts much faster than Newton iteration to find Legendre zeros)
#though with 13 pts get 1e-3 error same as Legendre - this may be way to go, but should
#first check on perturbed f integrals

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

#Use asymptotic expansion instead of just truncating the integral (standard practice)

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
#as it is hard to know when the quad points are reusing pre-computed values...
#Either way so long as the computation of quadrature points (which occurs in the background once)
#plus the sum operations that happen frequently are faster than quadgk * number of times it is called,
#then we are fine - and I think this is an easily attainable threshold for a small number of quad pts (<15)

#Would be interesting to see if cosmology affects this accuracy, because could reduce cost even further by
#just passing the q-pts (computed externally) once instead of computing in bg

#hmm I forgot a factor of q before, but even at higher redshift the performance
#is about the same on 15 pts, only around 1e-3 rerr
