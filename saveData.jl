using Metaheuristics, RCall, JLD
using Bilevel
using DataFrames
using HypothesisTests
using Statistics

include("distributed.jl")
include("algorithms.jl")

R""" 
library(irace)
getBestParms <- function(alg){
   dt <- data.frame()
   for(i in 1:10){
       fname = paste("data/irace/irace-", alg, "-", i, ".Rdata", sep="")
       load(fname)
       best.config <- getFinalElites(iraceResults = iraceResults, n = 1)
       dt = rbind(dt, best.config)
   }

   return(dt[,c(-1, -ncol(dt))])
}
"""

import Base.+
+(a::String, b::String) = string(a, b)

getBestParms = @rget getBestParms 



function getBCAPData(name, i)
    fname = "data/bcap/" + uppercase(name) + "_result_CEC17_D10_seed1560229400_run$(i).jld"
    load(fname, "result")
end

function getBestBCAP(pop)
    matrix = zeros(Int, length(pop), length(pop))
    for i =1:length(pop)
        for j=1:length(pop)
            i == j && continue
            for fn = 1:10
                s1 = pop[i].y[fn,:]
                s2 = pop[j].y[fn,:]
                matrix[i, j] += pvalue(SignedRankTest(s1, s2)) < 0.05 && median(s1) < median(s2)
            end
        end
    end

    score = sum(matrix, dims = 2)
   
    return pop[argmax(score)]
end

function getBestParmsBCAP(name)
    matrix = []
    for i =1:10
        res = getBCAPData(name, i)
        #best_sol = res.best_sol
        best_sol = getBestBCAP(res.population)
        
        s = mean(best_sol.y .== 0.0)
        x = best_sol.x
        if name == "abc"
            push!(matrix, [x[1:end-2]...,x[end], x[end-1], s])
        else
            push!(matrix, [x..., s])
        end
        
    end

    matrix = [ matrix[i][j] for i=1:length(matrix), j=1:length(matrix[1]) ]
    
    
    df = DataFrame(matrix)
    n = names(df)
    n[end] = :SRatio
    
    names!(df, n)
    sort!(df, :SRatio, rev=true)
    df
end


function mergeResults(name)
    df_irace =  getBestParms(name)

    df_bcap= getBestParmsBCAP(name)[:,1:end-1] # SRatio is not considered
    


    nms = String.(names(df_irace))
    names1 = "bcap_" .+ nms
    names2 = "irace_" .+ nms
    df = hcat(df_bcap, df_irace)
    names!(df, Symbol.([ names1..., names2...  ]))

    names!(df_bcap, names(df_irace))


    return Dict("merged" => df, "bcap" => df_bcap, "irace" => df_irace)
end

function saveBestParms()
    tables = Dict{String, Dict}()
    for name in ["abc", "de", "eca", "pso"]
        tables[name] = mergeResults(name)
    end

    save("dataFrame_parameters_bcap_irace.jld", "tables", tables)

end

function saveBestParmsResults()
    fname = "dataFrame_parameters_bcap_irace.jld"

    if !isfile(fname)
        saveBestParms()
    end

    bestParms = load(fname, "tables")

    tuners = ["bcap", "irace"]
    allRuns = Dict()
    for algname in ["eca", "de", "pso", "abc"]
        res = Dict()
        for tuner = tuners
            tab = bestParms[algname][tuner]
            
            errors = []
            costs = []
            for i = 1:10
                Φ = tab[i, :]
                r = getErrors(Φ, benchmark, algname)
                push!(errors, r[1])
                push!(costs,  r[2])
                println("Alg: $algname Tuner: $tuner Parms_id: $i ")
            end

            res[tuner] = Dict( "error" => errors, "cost" => costs )
        end

        allRuns[algname] = res

    end

    save("accuracy_of_parameters_bcap_irace.jld", "parameter_runs", allRuns)
end
