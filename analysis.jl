using JLD
import Statistics: median, mean
using HypothesisTests, DataFrames, Bilevel

function getBCAPData(name, i)
    fname = "data/bcap/$(uppercase(name))_result_CEC17_D10_seed1560229400_run$(i).jld"
    load(fname, "result")
end

function table2latex(df)
    
    header = join(["Run", String.(names(df))...], " & ")

    println(header, " \\\\ \\hline")
    for i = 1:nrow(df)
        str = ["$i", sprint.(print, round.(df[i,:], digits=6))...]
        println(join(str, " & "), "\\\\ \\hline ")
    end

end


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
    score3 = ( [ sum(y .== 0.0)/length(y) for y in runs[algname][tuner1]["error"]] )

    df = parms[algname][tuner1]
    insertcols!(df, ncol(df)+1, :SR => round.(Int, 100score3))

    cost = [ mean(y) for y in runs[algname][tuner1]["cost"] ] ./ 100000
    cost = min.(1, cost)
    insertcols!(df, ncol(df)+1, :NFEs => (round.(Int, 100cost)))
    
    insertcols!(df, ncol(df)+1, :win => score)
    insertcols!(df, ncol(df)+1, :tie => score2)
    insertcols!(df, ncol(df)+1, :lose => 100 .- (score2 + score))
    sort!(df, [:lose])
   
    return df
end

function test(algname)
    bestParms = load("dataFrame_parameters_bcap_irace.jld", "tables")

    allRuns = load("accuracy_of_parameters_bcap_irace2.jld", "parameter_runs")

    # algname = "eca"

    df = updateTable(allRuns, bestParms, algname, "irace", "bcap")
    df2 = updateTable(allRuns, bestParms, algname, "bcap", "irace")


    println("irace")
    # display(df)
    display(table2latex(df))
    println("bcap")
    display(table2latex(df2))
    # display(df2)

    # return median(allRuns["pso"]["irace"]["error"][1], dims= 2)
end


using StatsPlots
import Random: shuffle
pyplot()

function convergence()


    name  = "eca"
 
    markers  = shuffle(Plots._shape_keys)
    l_styles = rand([:solid, :dash, :dot, :dashdot], 10)

    calls = Int[]

    l = @layout [a b; c d]

    p = plot(layout= l, size = (800, 500))

    algs = ["abc", "eca", "de", "pso"]

    for j = 1:4
        for i = 1:10
            println(algs[j], " ", i)
            res  = getBCAPData(algs[j], i)

            conv = map(s->s.best_sol.F, res.convergence )
            conv_calls = map(s->s.F_calls, res.convergence )

            plot!(p[j], conv_calls, (conv),
                        marker=markers[i],
                        linestyle=l_styles[i],
                        title="Convergence: $(uppercase(algs[j]))",
                        legend=false,
                        xlabel="UL calls")

        end
    end

    p

end

convergence()
