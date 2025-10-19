import Pkg
Pkg.activate(@__DIR__)

using Test
using DiffFoldRNA: RandomModel, ViennaModel, test_seq_brute

include("common-validation.jl")

# TODO
# - collect table with results (absolute and relative errors)

function main(args)
    exit_code, n_start, n_end, atol = parse_cmdline(args; n_start=2, n_end=3, atol=1e-7)
    exit_code == EXIT_SUCCESS || return exit_code

    niter_per_n = 10
    @testset verbose=true "seq_partition vs bruteforce: random dbn, random pseq ($niter_per_n times per n)" begin
        @show n_start:n_end
        @testset verbose=true "RandomModel" begin
            test_seq_brute(; nrange=n_start:n_end, model=RandomModel(), hpmin=0, niter_per_n)
        end
        @testset verbose=true "ViennaModel" begin
            test_seq_brute(; nrange=n_start:n_end, model=ViennaModel(), hpmin=3, niter_per_n)
        end
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
