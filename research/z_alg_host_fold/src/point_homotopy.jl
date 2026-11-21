# Standard multihomogeneous continuation on the same two point fibers.
# Julia 1.11 + HomotopyContinuation.jl. The separate Python certificate uses
# analytic host enclosures even though this candidate generator uses midpoints.
using HomotopyContinuation, DelimitedFiles, LinearAlgebra
BLAS.set_num_threads(1)
out = joinpath(@__DIR__, "..", "results")
@var c[1:3] d[1:3]
C = [1; c]
D = [1; d]
system_seconds = Dict()
for side in ["minus", "plus"]
    coeff = readdlm(joinpath(out, "hc_coeff_" * side * ".txt"))
    val(i,a,b) = complex(coeff[(i-1)*16+4*(a-1)+b,1],coeff[(i-1)*16+4*(a-1)+b,2])
    equation(i) = sum(val(i,a,b)*C[a]*D[b] for a=1:4,b=1:4)
    build = @elapsed sys = System([equation(i) for i=1:6]; variable_groups=[c,d])
    function track(sys)
        return solve(sys; start_system=:total_degree, seed=UInt32(17029),
                     threading=false, show_progress=false)
    end
    first = @elapsed r = track(sys)
    warm = @elapsed r = track(sys)
    sols = solutions(r)
    open(joinpath(out,"hc_roots_" * side * ".txt"),"w") do io
        for v in sols
            println(io,join([real.(v); imag.(v)]," "))
        end
    end
    system_seconds[side] = (build,first,warm,ntracked(r),length(sols),length(failed(r)))
end
open(joinpath(out,"hc_times.tsv"),"w") do io
    println(io,"side\tbuild_s\tfirst_s\twarm_s\tpaths\tsolutions\tfailed")
    for side in ["minus", "plus"]
        println(io,side,"\t",join(system_seconds[side],"\t"))
    end
end
println("Julia ", VERSION,"; HomotopyContinuation ",pkgversion(HomotopyContinuation))
