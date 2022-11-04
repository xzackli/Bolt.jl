using Bolt
using AdvancedHMC,ForwardDiff,Distributions
using Plots,StatsPlots
using LinearAlgebra
using DelimitedFiles
using Random
Random.seed!(0); # Set seed for reproducibility

D = 1
initial_θ=[1.0]

function TTTEEE(𝕡::DT) where DT
    #boilerplate
    bg = Background(𝕡; x_grid=-20.0:0.01:0.0, nq=15)
    𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b, OmegaG=𝕡.Ω_r)
    ih = IonizationHistory(𝕣, 𝕡, bg)
    #hardcode k limits for the moment
    kmin,kmax,nk = 2.335e-5,2.335e-1,100 #0.1bg.H₀, 1000bg.H₀
    k_grid = quadratic_k(kmin, kmax, nk)
    sf_t = source_grid(𝕡, bg, ih, k_grid, BasicNewtonian())
    sf_p = source_grid_P(𝕡, bg, ih, k_grid, BasicNewtonian())
    #spectra
    ells = 10:20:1200
    ellfac = @.ells*(ells+1)
    TT = cltt(ells, 𝕡, bg, ih, sf_t)
    TE = clte(ells, 𝕡, bg, ih, sf_t, sf_p)
    EE = clee(ells, 𝕡, bg, ih, sf_p)
    return hcat(TT,TE,EE).*ellfac
end

function toy_TTTEEE(θA::DT) where DT
    return fid_spec .* θA #(As/Afid)
end

function real_TTTEEE(θA::DT) where DT
    Afid = 2.097e-9
    As = θA*Afid
    𝕡 = CosmoParams{DT}(A=As)
    return vcat(TTTEEE(𝕡)...)
end


#FIXME qthreads dies on perlmutter
Bolt.bessel_interpolator(2, k_grid[end] * bg.η₀)(k_grid[end]*(bg.η₀ - bg.η(bg.x_grid[end])))
#^this is the wrong behavior? will open an issue...

#save the fiducial TTTEEE for fast checks
ells = 10:20:1200
ellfac = @.ells*(ells+1)
TTTEEE_fid = TTTEEE(CosmoParams())
# writedlim("fid_spec.dat",vcat(TTTEEE_fid...))
#check this
# plot(ells,TTTEEE_fid[:,3].*ellfac)
# plot(TTTEEE_fid.*ellfac)
# ylims!(-1e-11,1e-11)

#simplified example with just Aₛ
Afid = 2.097e-9
fid_spec=readdlm("fid_spec.dat")[:,1]
fid_cov = Diagonal(fid_spec.^2 .*I(length(fid_spec))) #diagonal covariance

# Define the target distribution
N=length(fid_spec)
d = fid_spec
C  =fid_cov
function ℓπ(θ)
    #TODO to make this a real example just need to put TTTEEE
    μ = toy_TTTEEE(θ)
    # μ = real_TTTEEE(θ)
    return logpdf(MvNormal(d, C), μ)
end

ℓπ(1.0)

#HMC stuff
n_samples, n_adapts = 2, 10 #only using 2 steps here to test "real" TTTEEE
metric = DiagEuclideanMetric(D)
hamiltonian = Hamiltonian(metric, ℓπ, ForwardDiff)
initial_ϵ = find_good_stepsize(hamiltonian, initial_θ)
integrator = Leapfrog(initial_ϵ)
proposal = NUTS{MultinomialTS, GeneralisedNoUTurn}(integrator)
adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(0.8, integrator)) #what is this?

#run
samples, stats = sample(hamiltonian, proposal, initial_θ, n_samples, adaptor, n_adapts; progress=true)

histogram(reduce(vcat,samples))
vline!([1.0])
logpdf(MvNormal(d-toy_TTTEEE(0.5), C.*(0.5).^2), d-toy_TTTEEE(0.5))
logpdf(MvNormal(d-toy_TTTEEE(0.5), C), toy_TTTEEE(0.5))
