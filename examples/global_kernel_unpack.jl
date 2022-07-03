using CUDA, StaticArrays, OffsetArrays
using Pkg
Pkg.activate("/pscratch/sd/j/jsull/julia/Bolt.jl")
using Bolt
using BenchmarkTools
using Setfield


# bg/ion setup
𝕡 = CosmoParams()
bg = Background(𝕡; x_grid=-20.0:0.1:0.0, nq=6)
𝕣 = Bolt.RECFAST(bg=bg, Yp=𝕡.Y_p, OmegaB=𝕡.Ω_b)
ih = IonizationHistory(𝕣, 𝕡, bg)
k = 500bg.H₀
reltol=1e-5
#ℓᵧ = 20
#ℓ_ν = 20 
#ℓ_mν  = 4
#nq = 6
hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, 20,20,4,6)#ℓᵧ, ℓ_ν, ℓ_mν,nq)
let
	global const ℓᵧ = 20#4
	global const ℓ_ν = 20#4
	global const ℓ_mν  = 4#4
	global const nq = 6#4
    #hierarchy = Hierarchy(BasicNewtonian(), 𝕡, bg, ih, k, ℓᵧ, ℓ_ν, ℓ_mν,nq)
    xᵢ = first(hierarchy.bg.x_grid)
    u₀ = cu(Bolt.initial_conditions(xᵢ, hierarchy))
    #u₀ = cu(u₀)
    #u₀ = cu(zeros((2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5))
    println(typeof( u₀))
    du = cu([NaN])

	function gpu_unpack(u) 
        Θ = OffsetVector(SVector((u[i] for i=1:(ℓᵧ+1))...), 0:ℓᵧ)
        Θᵖ = OffsetVector(SVector((u[i] for i=(ℓᵧ+1)+1:(ℓᵧ+1)+(ℓᵧ+1) )...), 0:ℓᵧ)
        𝒩  = OffsetVector(SVector((u[i] for i=(2(ℓᵧ+1)+1):(2(ℓᵧ+1)+(ℓ_ν+1)) )...), 0:ℓ_ν)
        ℳ  = OffsetVector(SVector((u[i] for i=(2(ℓᵧ+1)+(ℓ_ν+1))+1:(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq) )...), 0:(ℓ_mν+1)*nq -1)
        Φ, δ, v, δ_b, v_b = SVector((u[i] for i=(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+1:(2(ℓᵧ+1)+(ℓ_ν+1)+(ℓ_mν+1)*nq)+5  )...)
        return Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b
    end

    function f_kernel!(du)
       Θ, Θᵖ, 𝒩, ℳ, Φ, δ, v, δ_b, v_b = gpu_unpack(u₀)
       du[1]=Θ[1]
       return nothing
   end
   @cuda f_kernel!(du)
   CUDA.@allowscalar du[1]
end
