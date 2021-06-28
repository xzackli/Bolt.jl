#Zack's sparsity pattern script
using Bolt
using BenchmarkTools
using OrdinaryDiffEq
using Plots
using SparsityDetection, SparseArrays
using Printf

par = CosmoParams()
ℓᵧ,ℓ_mν,nq = 10,10,15#2,2,2
ℓ_ν=ℓᵧ
bg = Background(par,nq=nq)
ih = IonizationHistory(Peebles(), par, bg)

k_grid = quadratic_k(0.1bg.H₀, 1000bg.H₀, 100)
k_i = 10
hierarchy = Hierarchy(BasicNewtonian(), par, bg, ih, k_grid[k_i],ℓᵧ,ℓ_mν,nq)
xᵢ = (hierarchy.bg.x_grid)[1]
u₀ = Bolt.initial_conditions(xᵢ, hierarchy)
du_dummy = deepcopy(u₀)

hf(du, u) = Bolt.hierarchy!(du, u, hierarchy, xᵢ)
sparsity_pattern = jacobian_sparsity(hf, du_dummy, u₀)
spy(sparsity_pattern, marker = (:square, 1))
lines = transpose([ℓᵧ+1,2(ℓᵧ+1),2(ℓᵧ+1)+ℓ_ν+1,2(ℓᵧ+1)+(ℓ_ν+1)+nq*(ℓ_mν+1)]).+0.5
labels = ["←γ" "←γP" "←ν" "←ℳ" ]
#regular lines
vline!(lines,label=false,color="green")
hline!(lines,label = false,color="green")
#nq lines
nq_lines = transpose([2(ℓᵧ+1)+ℓ_ν+1+nq*i for i in 1:ℓ_mν]).+0.5
hline!(nq_lines,label = false,color="green",ls=:dash)
vline!(nq_lines,label = false,color="green",ls=:dash)
plot!(xticks = (lines, labels))
ylabel!("X` depends on...")
xlabel!("X affects...")
title!("(ℓᵧ=ℓ_ν, ℓ_mν, nq) = ($(@sprintf("%d", ℓᵧ)),$(@sprintf("%d",ℓ_mν)),$(@sprintf("%d",nq)))")
xlims!(1,size(sparsity_pattern)[1])
ylims!(1,size(sparsity_pattern)[1])
savefig("../compare/sparse_nu_big.png")
#println(sparsity_pattern)

#back to the big one

#from trivial example
#so if there is a dot at row r and col c it means that the variable at column
#c impacts the output at row r


#github example
fcalls = 0
function f(dx,x)
  global fcalls += 1
  for i in 2:length(x)-1
    dx[i] = x[1] - 2x[i] + x[i+1]
  end
  dx[1] = -2x[1] + x[2]
  dx[end] = x[end-1] - 2x[end]
  nothing
end
nput = rand(10)
output = similar(nput)
sparsity_pattern = jacobian_sparsity(f,output,nput)
jac = Float64.(sparse(sparsity_pattern))
spy(sparsity_pattern, marker = (:square, 5))
