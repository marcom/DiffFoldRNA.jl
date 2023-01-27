using Test, ProgressMeter, DataFrames
using DiffFoldRNA
using DiffFoldRNA: test_all_1_nussinov_vienna, test_model_brute_vienna,
    test_seq_one_hot, test_seq_brute

@testset "vienna" begin
    showtestset()
    T = Float64
    nts = DiffFoldRNA.NTS

    random_p_seq(n) = random_p_seq(Float64, n)
    function random_p_seq(T, n)
        p = rand(T, n, nts)
        p ./= sum(p; dims=2)
    end

    @testset "seqstruct_partition" begin
        # basic functionality test
        @showprogress "seqstruct_partition" for n in [1:12...,]
            @test seqstruct_partition(ones(T, n, nts)/nts, RandomModel{T}()) isa T
        end
        @test_throws ErrorException seqstruct_partition(ones(T, 0, nts)/nts, RandomModel{T}())
    end

    @testset "test_all_1_nussinov_vienna" begin
        @showprogress "test_all_1_nussinov_vienna" for n = 1:20
            # n = 12: first multiloop for HAIRPIN=3
            test_all_1_nussinov_vienna(n)
        end
    end

    @testset "test_model_brute_vienna(RandomModel)" begin
        @showprogress "test_model_brute_vienna(RandomModel)" for n = 1:8
            em = RandomModel()
            test_model_brute_vienna(em, n)
        end
    end

    @testset "test_model_brute_vienna(ViennaModel)" begin
        @showprogress "test_model_brute_vienna(ViennaModel)" for n = 1:8
            em = ViennaModel()
            test_model_brute_vienna(em, n)
        end
    end

    # test_seq_one_hot: RandomModel
    @testset "test_seq_one_hot" begin
        @showprogress "test_seq_one_hot, RandomModel/ViennaModel, rand seq/dbn, n = 1:20, 40, 50, 60" for n in 1:20
            for (em, hpmin) in [
                RandomModel()   => 0,
                RandomModel(5)  => 1,
                RandomModel(42) => 2,
                ViennaModel()   => 3,
                All1Model()     => 0,
                ]
                for _ = 1:5
                    test_seq_one_hot(n; model=em, hpmin)
                end
            end
        end
        em = RandomModel()
        @showprogress "test_seq_one_hot(RandomModel), fixed seq, fixed dbn" for (seq,dbn) in [
            "GGCGU" => "(.)..",
            "GAAGU" => "(...)",
            ]
            test_seq_one_hot(; seq, dbn, model=em)
        end
    end

    @testset "test_seq_brute" begin
        test_seq_brute(nrange=1:12, model=RandomModel(), hpmin=0)
        test_seq_brute(nrange=1:12, model=ViennaModel(), hpmin=3)
    end
end
