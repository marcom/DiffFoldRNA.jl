import Pkg
Pkg.activate(@__DIR__)

using Dates: now, Second
using BenchmarkTools, CSV, DataFrames
using DiffFoldRNA
using DiffFoldRNA: NTS

function main(args)
    @show outpath = joinpath(@__DIR__, "results-seqstruct-partition-$(trunc(now(), Second)).csv")
    df = DataFrame(n = Int[], model=String[],
                   time = Float64[],
                   val = Float64[])
    for em in [ViennaModel(), All1Model(), RandomModel(1)]
        println("seqstruct_partition timings, $(nameof(typeof(em))")
        @show nrange = [20, 30, 40, 50]
        n = 15
        p_seq = ones(n,NTS)/NTS
        seqstruct_partition(p_seq, em)
        for n in nrange
            p_seq = ones(n,NTS)/NTS
            t = @timed @time val = seqstruct_partition(p_seq, em)
            push!(df, (; n, model="$em", time=t.time, val))
            display(df)
        end
        println()
        println()
    end
    CSV.write(outpath, df)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
