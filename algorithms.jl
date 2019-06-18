#                    N    Ne   limit
const ABC_bounds  = [10   0.1     0;
                     500  1.0   1000]

#                    N    K   η_max
const ECA_bounds  = [ 10  2    0.0;
                      500 10    4]

#                     N  F  CR
const DE_bounds   = [10  0  0;
                     500 2  1.0]

#                    N   C1   C2   omega
const PSO_bounds  = [10    0   0    0;
                     500   4   4    1.0]

function ABC_errors(Φ, benchmark)
    Errors_shared = SharedArray{Float64}(length(benchmark), NRUNS)
    Cost_shared = SharedArray{Int}(length(benchmark), NRUNS)

    Φ_ = SharedArray{Float64}([round(Φ[:N]), Φ[:Ne], round(Φ[:limit])])

    @sync @distributed for fn in 1:length(benchmark)
        N = Int(Φ_[1])
        Ne = round(Int, Φ_[2]*N)
        nruns = 0
        ff(x) = begin
            nruns += 1
            benchmark[fn].F(x)
        end
        for r = 1:NRUNS
            x, fx = ABC(ff, benchmark[fn].bounds_ul;
                                N = N,
                            limit = Int(Φ_[3]),
                               Ne = Ne,
                               No = N - Ne,
                                max_evals=max_NFEs,
                                termination=accuracy_termination,
                                )

            Errors_shared[fn, r] = fx < desired_accu ? 0.0 : fx
            Cost_shared[fn, r] = nruns 
            nruns *= 0
        end
    end

    Errors = Matrix(Errors_shared)
    Cost = Matrix(Cost_shared)

    return  [Errors, Cost]
    
end


function ECA_errors(Φ, benchmark)
    Errors_shared = SharedArray{Float64}(length(benchmark), NRUNS)
    Cost_shared = SharedArray{Int}(length(benchmark), NRUNS)
    
    Φ_ = SharedArray{Float64}([ round(Φ[:N]), round(Φ[:K]), Φ[:eta]  ])

    @sync @distributed for fn in 1:length(benchmark)
        nruns = 0
        ff(x) = begin
            nruns += 1
            benchmark[fn].F(x)
        end
        for r = 1:NRUNS
            x, fx = eca(ff, D;
                                N = Int(Φ_[1]),
                                K = Int(Φ_[2]),
                                η_max = Φ_[3],
                                max_evals=max_NFEs,
                                termination=accuracy_termination,
                                limits  =  benchmark[fn].bounds_ul,
                                showResults=false)

            Errors_shared[fn, r] = fx < desired_accu ? 0.0 : fx
            Cost_shared[fn, r] = nruns 
            nruns *= 0
        end
    end

    Errors = Matrix(Errors_shared)
    Cost = Matrix(Cost_shared)

    return  [Errors, Cost]
    
end

function DE_errors(Φ, benchmark)
    Errors_shared = SharedArray{Float64}(length(benchmark), NRUNS)
    Cost_shared = SharedArray{Int}(length(benchmark), NRUNS)
    Φ_ = SharedArray{Float64}([ round(Φ[:N]), Φ[:F], Φ[:CR] ])

    @sync @distributed for fn in 1:length(benchmark)
        nruns = 0
        ff(x) = begin
            nruns += 1
            benchmark[fn].F(x)
        end
        for r = 1:NRUNS
            x, fx = DE(ff, D;
                                N  = Int(Φ_[1]),
                                F  = Φ_[2],
                                CR = Φ_[3],
                                max_evals=max_NFEs,
                                termination=accuracy_termination,
                                limits  =  benchmark[fn].bounds_ul,
                                showResults=false)

            Errors_shared[fn, r] = fx < desired_accu ? 0.0 : fx
            Cost_shared[fn, r] = nruns 
            nruns *= 0
        end
    end

    Errors = Matrix(Errors_shared)
    Cost = Matrix(Cost_shared)

    return  [Errors, Cost]
    
end

function PSO_errors(Φ, benchmark)
    Errors_shared = SharedArray{Float64}(length(benchmark), NRUNS)
    Cost_shared = SharedArray{Int}(length(benchmark), NRUNS)
    Φ_ = SharedArray{Float64}( [ round(Φ[:N]), Φ[:C1], Φ[:C2], Φ[:omega] ] )

    @sync @distributed for fn in 1:length(benchmark)
        nruns = 0
        ff(x) = begin
            nruns += 1
            benchmark[fn].F(x)
        end
        for r = 1:NRUNS
            x, fx = pso(ff, D;
                                N  = Int(Φ_[1]),
                                C1 = Φ_[2],
                                C2 = Φ_[3],
                                ω  = Φ_[4],
                                max_evals=max_NFEs,
                                termination=accuracy_termination,
                                limits  =  benchmark[fn].bounds_ul,
                                showResults=false)

            Errors_shared[fn, r] = fx < desired_accu ? 0.0 : fx
            Cost_shared[fn, r] = nruns 
            nruns *= 0
        end
    end

    Errors = Matrix(Errors_shared)
    Cost = Matrix(Cost_shared)
    # Cost = Matrix(Cost_shared)

    return  [Errors, Cost]
    
end


function getErrors(Φ, benchmark, name)
    if name == "abc"
        return ABC_errors(Φ, benchmark)
    elseif name == "eca"
        return ECA_errors(Φ, benchmark)
    elseif name == "de"
        return DE_errors(Φ, benchmark)
    elseif name == "pso"
        return PSO_errors(Φ, benchmark)
    end
end
