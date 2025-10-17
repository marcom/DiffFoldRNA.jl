import Pkg
Pkg.activate(@__DIR__)

using Test

include("common-validation.jl")

function testcmd(cmd::Cmd)
    r = run(cmd)
    if r.exitcode != 0
        @warn "cmd failed: $cmd"
    end
    @test r.exitcode == 0
end

function main(args)
    exit_code, n_start, n_end, atol = parse_cmdline(args; n_start=2, n_end=2, atol=1e-1)
    exit_code == EXIT_SUCCESS || return exit_code
    scripts_to_check = [
        #"compare-to-python.jl",
        "compare-energy-to-ViennaRNA.jl",
        "check-seqstruct-partition-All1Model-vienna-vs-nussinov.jl",
        "check-seqstruct-partition-nussinov-vs-bruteforce.jl",
        "check-seqstruct-partition-vienna-vs-bruteforce.jl",
    ]
    println("n_start = $n_start, n_end = $n_end, atol = $atol")
    println("scripts_to_check =")
    display(scripts_to_check)
    println()
    @testset "quickcheck all scripts working" verbose=true begin
        jl = `julia --project=@.`
        for script in scripts_to_check
            @testset "$script" begin
                testcmd(`$jl $(joinpath(@__DIR__, script)) $n_start $n_end $atol`)
            end
        end
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
