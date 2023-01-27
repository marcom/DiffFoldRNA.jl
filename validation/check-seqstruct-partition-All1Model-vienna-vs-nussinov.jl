using Test
using DiffFoldRNA: test_all_1_nussinov_vienna

include("common-validation.jl")

# TODO
# - collect table with results (absolute and relative errors)

function main(args)
    exit_code, n_start, n_end, atol = parse_cmdline(args; n_start=2, n_end=3, atol=1e-7)
    exit_code == EXIT_SUCCESS || return exit_code
    @testset "seqstruct All1Model, vienna vs nussinov, atol=$atol" verbose=true begin
        for n = n_start:n_end
            @testset "$n" begin
                test_all_1_nussinov_vienna(n; atol)
            end
        end
    end
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
