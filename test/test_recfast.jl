using DelimitedFiles
using PyPlot
cd("test")

recfastdata = readdlm("data/test_recfast_1.dat", ',', Float64, '\n', header=true)[1]
z⃗, Xe = recfastdata[:,1], recfastdata[:,2]

# clf()
# plot(z, Xe, "-", label=raw"$X_e$")
# legend()
# gcf()

##
using Bolt

𝕡 = CosmoParams(Σm_ν = 0.0, N_ν = 3.0)
bg = Background(𝕡)
𝕣 = Bolt.RECFAST(bg=bg)
xe_bespoke = Bolt.recfast_xe(𝕣; Nz=1000, zinitial=10000., zfinal=0.);

##

clf()
plot(z⃗, Xe ./ xe_bespoke , "-", label=raw"RECFAST / recfast.jl")
ylim(1 - 0.01, 1 + 0.01)

# plot(z, Xe , "-")
# plot(z, xe_bespoke, "--")
xlabel(raw"redshift")
legend()
gcf()

##
using UnitfulAstro, NaturallyUnitful
x₁ = let z = 100.0
    Hz = 𝕣.HO * sqrt((1+z)^4/(1+𝕣.z_eq)*𝕣.OmegaT + 𝕣.OmegaT*(1+z)^3 + 𝕣.OmegaK*(1+z)^2 + 𝕣.OmegaL)
    (𝕣.HO^2 /2/Hz)*(4*(1+z)^3/(1+𝕣.z_eq)*𝕣.OmegaT + 3*𝕣.OmegaT*(1+z)^2 + 2*𝕣.OmegaK*(1+z))
end

##
x₂ =let z = 100.0
    H0_natural_unit_conversion = ustrip(u"s", unnatural(u"s", 1u"eV^-1"))

    a = 1 / (1+z)  # scale factor
    x_a = a2x(a)
	Hz = 𝕣.bg.ℋ(x_a) / a / 𝕣.H0_natural_unit_conversion
	dHdz = (-𝕣.bg.ℋ′(x_a) + 𝕣.bg.ℋ(x_a)) / 𝕣.H0_natural_unit_conversion
end

##
let z = 0.0
    a = 1.0
    x_a = a2x(a)
    H0_natural_unit_conversion = ustrip(u"s", unnatural(u"s", 1u"eV^-1"))
	Hz = 𝕣.bg.H₀ / 𝕣.H0_natural_unit_conversion

end
##
𝕣.HO
