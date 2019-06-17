using Distributed, SharedArrays



if nprocs() < Sys.CPU_THREADS - 1
    addprocs(Sys.CPU_THREADS - 1)
end


@everywhere import LinearAlgebra: dot
@everywhere using CEC17
@everywhere using Metaheuristics

@everywhere const desired_accu = 1e-6
@everywhere const D = 10
@everywhere const max_NFEs = 10000*D
@everywhere const NRUNS = 31


@everywhere mutable struct approxProblem
    F::Function
    bounds_ul::Matrix{Float64}
end

@everywhere macro getBenchmark(D)
    fs = approxProblem[]

    bounds = Array([ -100.0ones($D) 100.0ones($D) ]')

    for i = 1:10
        fun = approxProblem(x -> abs(cec17_test_func(x, i) - 100i), bounds)
        push!(fs, fun)
    end

    fs
end

@everywhere const benchmark = @getBenchmark(D)


@everywhere function accuracy_termination(P::Array)
    if typeof(P[1]) <: Float64
        f = minimum(P)
    elseif typeof(P[1]) <: Metaheuristics.Bee
        f = minimum( map(x->x.sol.f, P) )
    else
        f = minimum( map(x->x.f, P) )
    end

    return abs(f) < desired_accu
end
