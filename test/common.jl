using Test
using DiffFoldRNA: RNA_ALPHA, ALL_PAIRS, NTS, one_hot_seq, valid_pair,
    is_valid_dbn, matching_to_dbn, dbn_to_matching, random_primary,
    random_seq_for_dbn, seq_prob, structure_tree, structure_tree_full,
    structure_list_postorder, random_p_seq, normalize_to_p_seq,
    normalize_to_p_seq!

const MATCH_DBN = [
    [2, 1, 3] => "().",
    [1, 3, 2] => ".()",
    [3, 2, 1] => "(.)",
    [1, 3, 2, 4] => ".().",
    [2, 1, 4, 3] => "()()",
    [6, 3, 2, 5, 4, 1] => "(()())",
    [1, 7, 4, 3, 6, 5, 2] => ".(()())",
    [2, 1, 8, 5, 4, 7, 6, 3] => "()(()())",
]

@testset "common" begin
    showtestset()

    @testset "one_hot_seq" begin
        for T in (Float32, Float64)
            for (seq, want_pseq) in [
                "A" => T.([1 0 0 0;]),
                "ACGU" => T.([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]),
                ]
                pseq = one_hot_seq(T, seq)
                @test pseq isa Matrix{T}
                @test pseq == want_pseq
                @test pseq == T.(one_hot_seq(seq))
            end
        end
    end

    @testset "valid_pair(::Char, ::Char)" begin
        for (ci, cj) in ALL_PAIRS
            @test valid_pair(ci, cj) == true
        end
        all_combinations = vec(join.(Iterators.product(collect(RNA_ALPHA), collect(RNA_ALPHA)) |> collect))
        NOT_PAIRS = setdiff(all_combinations, ALL_PAIRS)
        for (ci, cj) in NOT_PAIRS
            @test valid_pair(ci, cj) == false
        end
    end

    @testset "is_valid_dbn" begin
        @test is_valid_dbn("")
        for n = 1:10
            @test is_valid_dbn("."^n)
            @test is_valid_dbn("(" * "."^n * ")")
            @test is_valid_dbn("("^n * ")"^n)
            @test is_valid_dbn("("^(div(n,2)) * "."^n * ")"^(div(n,2)))
            @test ! is_valid_dbn("."^n * ")")
            @test ! is_valid_dbn("."^n * ")"^n)
            @test ! is_valid_dbn("("^(n-1) * "."^(n-1) * ")"^(n-1) * ")")
        end
        @test ! is_valid_dbn("(.")
    end
    
    @testset "matching_to_dbn" begin
        for n = 0:10
            @test matching_to_dbn(Int[1:n...]) == "."^n
        end
        for (match, dbn) in MATCH_DBN
            @test matching_to_dbn(match) == dbn
        end
        # These should throw Exceptions as the input is an illegal
        # match
        for match in [
            [2],
            [1,3],
            [3,1],
            ]
            @test (try matching_to_dbn(match) catch e; e; end) isa Exception
        end
    end

    @testset "dbn_to_matching" begin
        for n = 0:10
            @test dbn_to_matching("."^n) == Int[1:n...]
        end
        for (match, dbn) in MATCH_DBN
            @test dbn_to_matching(dbn) == match
        end
    end

    @testset "random_primary" begin
        for n = 0:10
            seq = random_primary(n)
            @test seq isa String
            @test length(seq) == n
            @test all(c -> c in RNA_ALPHA, seq)
        end
    end

    @testset "structure_tree" begin
        for (dbn, res) in [
            ""                => (Dict(-1 => Int[]),
                                  Int[]),
            "()"              => (Dict(-1 => [1]),
                                  [2, -1]),
            ".()"             => (Dict(-1 => [2]),
                                  [-1, 3, -1]),
            "(.)"             => (Dict(-1 => [1]),
                                  [3, -1, -1]),
            "()()"            => (Dict(-1 => [1, 3]),
                                  [2, -1, 4, -1]),
            "(()())"          => (Dict(-1 => [1], 1 => [2, 4]),
                                  [6, 3, -1, 5, -1, -1]),
            "...(.(.)(.)..)." => (Dict(-1 => [4], 4 => [6, 9]),
                                  [-1, -1, -1, 14, -1, 8, -1, -1, 11, -1, -1, -1, -1, -1, -1]),
            ]
            @test structure_tree(dbn) == res
        end
    end

    @testset "structure_tree_full" begin
        for (dbn, res) in [
            ""                => (Dict(-1 => Int[]),
                                  Int[]),
            "()"              => (Dict(-1 => [1], 1 => Int[]),
                                  [2, -1]),
            ".()"             => (Dict(-1 => [2], 2 => Int[]),
                                  [-1, 3, -1]),
            "(.)"             => (Dict(-1 => [1], 1 => Int[]),
                                  [3, -1, -1]),
            "()()"            => (Dict(-1 => [1, 3], 1 => [], 3 => Int[]),
                                  [2, -1, 4, -1]),
            "(()())"          => (Dict(-1 => [1], 1 => [2, 4], 2 => Int[], 4 => Int[]),
                                  [6, 3, -1, 5, -1, -1]),
            "...(.(.)(.)..)." => (Dict(-1 => [4], 4 => [6, 9], 6 => Int[], 9 => Int[]),
                                  [-1, -1, -1, 14, -1, 8, -1, -1, 11, -1, -1, -1, -1, -1, -1]),
            ]
            @test structure_tree_full(dbn) == res
        end
    end

    @testset "structure_list_postorder" begin
        for (dbn, res) in [
            ""                => ([-1 => Int[]],
                                  Int[]),
            "()"              => ([1 => Int[], -1 => [1]],
                                  [2, -1]),
            ".()"             => ([2 => Int[], -1 => [2]],
                                  [-1, 3, -1]),
            "(.)"             => ([1 => Int[], -1 => [1]],
                                  [3, -1, -1]),
            "()()"            => ([3 => Int[], 1 => Int[], -1 => [1, 3]],
                                  [2, -1, 4, -1]),
            "(()())"          => ([4 => Int[], 2 => Int[], 1 => [2,4], -1 => [1]],
                                  [6, 3, -1, 5, -1, -1]),
            "...(.(.)(.)..)." => ([9 => Int[], 6 => Int[], 4 => [6, 9], -1 => [4]],
                                  [-1, -1, -1, 14, -1, 8, -1, -1, 11, -1, -1, -1, -1, -1, -1]),
            ]
            @test structure_list_postorder(dbn) == res
        end
    end

    @testset "seq_prob" begin
        for n = 1:5
            @test seq_prob(ones(n,NTS)/NTS, "A"^n) == (1/NTS)^n
        end
    end

    @testset "random_p_seq" begin
        for T in (Float32, Float64)
            for n = 1:5
                p = random_p_seq(T, n)
                @test p isa Matrix{T}
                @test axes(p) == (1:n, 1:NTS)
                @test all(x -> 0.0 <= x <= 1.0, p)
            end
        end
    end

    @testset "normalize_to_p_seq" begin
        function test_p_seq(p1, p2)
            @test axes(p1) == axes(p2)
            @test all(p1 .== p2)
            @test all(x -> 0.0 <= x <= 1.0, p1)
            for psum in sum(p1; dims=2)
                @test psum ≈ 1.0 atol=1e-7
            end
        end
        for n = 5:10
            x = rand(n, NTS) .- 0.5
            x2 = copy(x)

            p_seq = normalize_to_p_seq(x)
            x2 .= x
            normalize_to_p_seq!(x2)
            test_p_seq(p_seq, x2)
            for (i, fn) in enumerate([x -> x^2,
                                      x -> exp(x),
                                      x -> exp(-x),
                                      x -> sin(x)^2,
                                      x -> cos(x)^2, ])
                x = rand(n, NTS) .- 0.5
                x2 .= x
                p_seq = normalize_to_p_seq(fn, x)
                normalize_to_p_seq!(fn, x2)
                test_p_seq(p_seq, x2)
            end
        end
    end

    @testset "random_seq_for_dbn" begin
        dbns = [
            ["."^n for n = 1:10]...,
            ["(" * "."^n * ")" for n = 0:5]...,
            ["..((" * "."^n * ")).." for n = 0:5]...,
            "..(.(..(((...)).).((...))).)...",
        ]
        for dbn in dbns, _ = 1:10
            seq = random_seq_for_dbn(dbn)
            @test length(seq) == length(dbn)
            @test all(c -> c in RNA_ALPHA, seq)
            match = dbn_to_matching(dbn)
            @test all(enumerate(match)) do (i,j)
                i >= j || "$(seq[i])$(seq[j])" in ALL_PAIRS
            end
        end
    end
end
