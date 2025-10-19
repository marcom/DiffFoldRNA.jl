using Test

function test_nussinov_brute(n::Int; atol::Float64=1e-7, verbose::Bool=true, hpmin::Int=HAIRPIN)
    p_seq = random_p_seq(n)
    nuss = nussinov_seqstruct_partition(p_seq, nussinov_boltz_pair; hpmin)
    brute = seqstruct_partition_brute_force(Float64, p_seq, nussinov_boltz; hpmin)
    if verbose
        @show n, nuss, brute, nuss-brute
    end
    @test nuss ≈ brute atol=atol
end

function test_all_1_nussinov_vienna(n::Int; atol::Float64=1e-7, verbose::Bool=true)
    em = All1Model()
    p_seq = random_p_seq(n)
    nuss = nussinov_seqstruct_partition(p_seq, nussinov_boltz_pair_1)
    vien = seqstruct_partition(p_seq, em)
    if verbose
        @show n, nuss, vien, nuss-vien
    end
    @test nuss ≈ vien atol=atol
end

function test_model_brute_vienna(em::AbstractModel, n::Integer; atol=1e-7, verbose=true)
    p_seq = random_p_seq(n)
    vien = seqstruct_partition(p_seq, em)
    boltzfn(seq, match) = boltz(seq, matching_to_dbn(match), em)
    brute = seqstruct_partition_brute_force(p_seq, boltzfn)
    if verbose
        @show n, brute, vien, brute-vien
    end
    @test brute ≈ vien atol=atol
end

function test_seq_one_hot(n::Int; model::AbstractModel=RandomModel(), atol::Float64=1e-7, hpmin::Int=HAIRPIN)
    seq = random_primary(n)
    return test_seq_one_hot_for_seq(seq; model, atol, hpmin)
end
function test_seq_one_hot_for_seq(seq::AbstractString; model::AbstractModel=RandomModel(), atol::Float64=1e-7, hpmin::Int=HAIRPIN)
    dbn = random_dbn(AllStructs(seq, hpmin))
    return test_seq_one_hot(; seq, dbn, model, atol)
end
function test_seq_one_hot_for_dbn(dbn::AbstractString; model::AbstractModel=RandomModel(), atol::Float64=1e-7)
    dbn = random_seq_for_dbn(dbn::AbstractString)
    return test_seq_one_hot(; seq, dbn, model, atol)
end
function test_seq_one_hot(; seq::AbstractString, dbn::AbstractString,
                          model::AbstractModel=RandomModel(), atol::Float64=1e-7)
    #@show seq, dbn
    n = length(seq)
    n == length(dbn) || throw(ArgumentError("seq and dbn must have same length, $n != $(length(dbn))"))
    p_seq = one_hot_seq(seq)
    b_calc = boltz(seq, dbn, model)
    b_spart = seq_partition(p_seq, dbn, model)
    #@show n, seq, dbn, b_calc-b_spart, b_calc, b_spart
    @test b_calc ≈ b_spart atol=atol
end

using ProgressMeter: @showprogress
using DataFrames
function test_quick_seq_one_hot()
    model = RandomModel()
    df = DataFrame(n=Int[], seq=String[], dbn=String[], b_calc=Float64[], b_spart=Float64[], delta=Float64[])
    @showprogress "test_seq_one_hot(RandomModel), fixed seq, fixed dbn" for (seq,dbn) in [
        "GGCGU" => "(.)..",
        "GAAGU" => "(...)",
        "UACGGCGCUGGCUGU" => "()().((.(.).)).",
        ]
        #test_seq_one_hot(; seq, dbn, model)
        atol = 1e-7
        n = length(seq)
        p_seq = one_hot_seq(seq)
        b_calc = boltz(seq, dbn, model)
        b_spart = seq_partition(p_seq, dbn, model)
        #@show n, seq, dbn, b_calc-b_spart, b_calc, b_spart
        if ! isapprox(b_calc, b_spart; atol)
            push!(df, (; n, seq, dbn, b_calc, b_spart, delta=b_calc-b_spart))
        end
        @test b_calc ≈ b_spart atol=atol
    end
    display(df)
end

function test_seq_brute(; nrange=1:12, niter_per_n::Int=5, model=RandomModel(), hpmin::Int=HAIRPIN, atol::Float64=1e-7)
    @showprogress "test_seq_brute" for it = 1:niter_per_n, n in nrange
        p_seq = random_p_seq(n)
        dbn = random_dbn(AllStructs(n, hpmin))
        boltz_fn(seq, dbn) = boltz(seq, dbn, model)
        b_brute = seq_partition_brute_force(p_seq, dbn, boltz_fn)
        b_spart = seq_partition(p_seq, dbn, model)
        @test b_spart ≈ b_brute atol=atol
    end
end
