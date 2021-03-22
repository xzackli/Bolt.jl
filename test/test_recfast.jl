using DelimitedFiles
using PyPlot
cd("test")

recfastdata = readdlm("data/test_recfast_1.dat", ',', Float64, '\n', header=true)[1]
z, Xe = recfastdata[:,1], recfastdata[:,2]

##
clf()
plot(z, Xe, "-", label=raw"$X_e$")
legend()
gcf()

##
using Bolt

𝕣 = Bolt.RECFASTIonization()
@time xe_bespoke = Bolt.recfast_xe(𝕣; Nz=1000, zinitial=10000., zfinal=0.);

##

clf()
plot(z, Xe ./ xe_bespoke , "-", label=raw"ratio")
ylim(1 - 0.001, 1 + 0.001)
# plot(z, xe_bespoke, "--")
legend()
gcf()

#
