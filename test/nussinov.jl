using Test
using LogExpFunctions: softmax
#using DiffFoldRNA: NTS, nussinov_boltz, nussinov_boltz_pair,
#    nussinov_seqstruct_partition, seqstruct_partition_brute_force
using DiffFoldRNA: test_nussinov_brute

@testset "nussinov" begin
    showtestset()

    for hpmin = 0:3
        for n = 1:7
            # p_seq = rand(n, NTS)  # TODO: change to (NTS,n) when moved to column-major
            # p_seq = softmax(p_seq; dims=2)
            # pf = nussinov_seqstruct_partition(p_seq, nussinov_boltz_pair)
            # pf_brute = seqstruct_partition_brute_force(Float64, p_seq, nussinov_boltz)
            # @test pf ≈ pf_brute atol=1e-7
            test_nussinov_brute(n; hpmin)
        end
        #@test nussinov_test(2) isa Float64
    end
end
