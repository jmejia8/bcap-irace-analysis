using JLD
import Statistics: median, mean
using HypothesisTests, DataFrames

function updateTable(runs, parms, algname, tuner1, tuner2  = tuner1)
 
    L = length(runs[algname][tuner1]["error"])
    matrix = zeros(Int, L, L)
    matrix2 = zeros(Int, L, L)

    for i =1:L
        for j=1:L
            i == j && tuner1==tuner2  && continue

            s1 = runs[algname][tuner1]["error"][i]
            s2 = runs[algname][tuner2]["error"][j]
            for fn = 1:10
                ξ1 = s1[fn,:]
                ξ2 = s2[fn,:]
                
                if pvalue(SignedRankTest(ξ1, ξ2)) < 0.05 
                    matrix[i, j] += median(ξ1) < median(ξ2)
                else
                    matrix2[i, j] += 1
                end
            end
        end
    end

    score = sum(matrix, dims = 2)[:,1]
    score2 = sum(matrix2, dims = 2)[:,1]
    score3 = ( [ sum(y .== 0.0) for y in runs[algname][tuner1]["error"]] )

    df = parms[algname][tuner1]
    insertcols!(df, ncol(df)+1, :win => score)
    insertcols!(df, ncol(df)+1, :tie => score2)
    insertcols!(df, ncol(df)+1, :lose => 100 .- (score2 + score))
    insertcols!(df, ncol(df)+1, :SR => score3)

    cost = [ mean(y) for y in runs[algname][tuner1]["cost"] ] ./ 100000
    insertcols!(df, ncol(df)+1, :NFEs => cost)
    
    sort!(df, [:lose])
   
    return df
end

function test(algname)
    bestParms = load("dataFrame_parameters_bcap_irace.jld", "tables")

    allRuns = load("accuracy_of_parameters_bcap_irace.jld", "parameter_runs")

    # algname = "eca"

    df = updateTable(allRuns, bestParms, algname, "irace", "bcap")
    df2 = updateTable(allRuns, bestParms, algname, "bcap", "irace")


    println("irace")
    display(df)
    println("bcap")
    display(df2)

    # return median(allRuns["pso"]["irace"]["error"][1], dims= 2)
end

test()
