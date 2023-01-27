using Test
using DiffFoldRNA
using DiffFoldRNA: HAIRPIN, allstructs_boltz_sum, boltz, nussinov_boltz, one_hot_seq,
    seqstruct_partition_brute_force, seq_partition_brute_force

@testset "brute_force" begin
    showtestset()

    @testset "allstructs_boltz_sum" begin
        T = Float64
        seq = "GAC"
        boltz_fn(seq, match) = one(T)
        s = allstructs_boltz_sum(T, seq, boltz_fn; hpmin=HAIRPIN)
        @test s isa T
        @test s == allstructs_boltz_sum(seq, boltz_fn)
    end

    @testset "seqstruct_partition_brute_force" begin
        n = 5
        p_seq = ones(n,NTS) / NTS
        T = Float64
        pf = seqstruct_partition_brute_force(T, p_seq, nussinov_boltz; hpmin=HAIRPIN)
        @test pf isa T
        @test pf == seqstruct_partition_brute_force(p_seq, nussinov_boltz)
    end

    @testset "seq_partition_brute_force" begin
        T = Float64
        model = ViennaModel()
        boltz(seq, dbn) = DiffFoldRNA.boltz(seq, dbn, model)

        n = 7
        p_seq = random_p_seq(n)
        dbn = "((...))"
        seq = "GCAUAGC"
        p = seq_partition_brute_force(T, p_seq, dbn, boltz)
        @test p isa T
        @test p == seq_partition_brute_force(p_seq, dbn, boltz)
        @test seq_partition_brute_force(one_hot_seq(seq), dbn, boltz) == boltz(seq, dbn)
    end
end
