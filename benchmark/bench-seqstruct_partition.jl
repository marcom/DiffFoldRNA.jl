#!/usr/bin/env julia

"""
Benchmark driver for `DiffFoldRNA.seqstruct_partition` and related helper kernels.

This script gathers runtime and GC statistics across a range of sequence lengths,
number types, and energy models.

The collected results are printed as a `DataFrame` and written to a CSV file
(`benchmark/seqstruct_partition_benchmarks.csv` by default).
"""

import Pkg
Pkg.activate(@__DIR__)

using BenchmarkTools
using CSV
using DataFrames
using DiffFoldRNA
using ForwardDiff
using OffsetArrays
using Random
using Statistics
using Dates

const RNG_SEED = 0x53524E41   # "SRNA" in hex
const RNG = MersenneTwister(RNG_SEED)
const NTS = DiffFoldRNA.NTS
const DEFAULT_SEQ_LENGTHS = (5, 10, 20, 30, 40, 50)
const DEFAULT_GRADIENT_LENGTHS = (5, 10, 20)
function default_output_path()
    now = Dates.now()
    stamp = Dates.format(now, dateformat"yyyy-mm-dd-HH:MM:SS")
    return joinpath(@__DIR__, "results-seqstruct_partition-$stamp.csv")
end
const PARTIAL_INDEX = (1, 1)  # entry receiving a dual seed for partial derivatives

const MODEL_SPEC = Dict{DataType, Vector{Pair{String, Function}}}(
    Float64 => [
        "All1Model"   => () -> DiffFoldRNA.All1Model{Float64}(),
        "ViennaModel" => () -> DiffFoldRNA.ViennaModel(),
        "RandomModel" => () -> DiffFoldRNA.RandomModel{Float64}(1),
    ],
    Float32 => [
        "All1Model"   => () -> DiffFoldRNA.All1Model{Float32}(),
        # ViennaModel stores Float64 parameters internally; include it for comparison.
        "ViennaModel" => () -> DiffFoldRNA.ViennaModel(),
        "RandomModel" => () -> DiffFoldRNA.RandomModel{Float32}(1),
    ],
)

# Cache reusable base probability matrices so different benchmarks see the same inputs.
const BASE_PSEQ_CACHE = Dict{Tuple{DataType, Int}, Matrix}()

function base_p_seq(::Type{T}, n::Int)::Matrix{T} where {T}
    key = (T, n)
    return get!(BASE_PSEQ_CACHE, key) do
        p = rand(RNG, T, n, NTS)
        p ./= sum(p; dims=2)
        copy(p)
    end
end

partial_dual_matrix(p_base::AbstractMatrix{T}, idx::Tuple{Int, Int}) where {T} =
    dualize_single_entry(p_base, idx)

function dualize_single_entry(p_base::AbstractMatrix{T}, idx::Tuple{Int, Int}) where {T}
    dual_mat = Array{ForwardDiff.Dual{Nothing, T, 1}}(undef, size(p_base))
    for j in axes(p_base, 2)
        for i in axes(p_base, 1)
            partial = (i == idx[1] && j == idx[2]) ? one(T) : zero(T)
            dual_mat[i, j] = ForwardDiff.Dual{Nothing, T, 1}(p_base[i, j],
                ForwardDiff.Partials{1, T}((partial,)))
        end
    end
    return dual_mat
end

function prepare_dp_state(p_core::AbstractMatrix{S}, ::Type{E}) where {S, E}
    n = size(p_core, 1)
    p_seq = OffsetMatrix(zeros(E, n + 2, NTS), 0:n+1, 1:NTS)
    p_seq[1:n, 1:NTS] .= E.(p_core)
    p_seq[begin, begin] = one(E)
    p_seq[end, begin] = one(E)

    E_tab  = OffsetArray(zeros(E, NTS, NTS, n + 2), 1:NTS, 1:NTS, 0:n+1)
    P_tab  = OffsetArray(zeros(E, NTS, NTS, n + 2, n + 2), 1:NTS, 1:NTS, 0:n+1, 0:n+1)
    ML_tab = OffsetArray(zeros(E, NTS, NTS, NTS, NTS, 3, n + 2, n + 2),
                         1:NTS, 1:NTS, 1:NTS, 1:NTS, 0:2, 0:n+1, 0:n+1)
    OMM    = OffsetArray(zeros(E, NTS, NTS, n + 2, n + 2), 1:NTS, 1:NTS, 0:n+1, 0:n+1)

    return (p_seq = p_seq, E = E_tab, P = P_tab, ML = ML_tab, OMM = OMM, hairpin = DiffFoldRNA.HAIRPIN)
end

function seqstruct_gradient(p_core::AbstractMatrix{T}, model; hpmin=DiffFoldRNA.HAIRPIN) where {T}
    n = size(p_core, 1)
    function f(v::AbstractVector)
        mat = reshape(v, n, NTS)
        DiffFoldRNA.seqstruct_partition(mat, model; hpmin=hpmin)
    end
    return ForwardDiff.gradient(f, vec(copy(p_core)))
end

mean_or_zero(x) = isempty(x) ? zero(eltype(x)) : mean(x)
std_or_zero(x) = (length(x) <= 1) ? zero(eltype(x)) : std(x)
minimum_or_zero(x) = isempty(x) ? zero(eltype(x)) : minimum(x)
maximum_or_zero(x) = isempty(x) ? zero(eltype(x)) : maximum(x)

function summarise_trial(trial::BenchmarkTools.Trial)
    times = trial.times
    gctimes = trial.gctimes
    return (
        samples = length(times),
        evals = trial.params.evals,
        time_min_ns = minimum(times),
        time_max_ns = maximum(times),
        time_mean_ns = mean(times),
        time_median_ns = median(times),
        time_std_ns = std_or_zero(times),
        gctime_total_ns = sum(gctimes),
        gctime_mean_ns = mean_or_zero(gctimes),
        gctime_min_ns = minimum_or_zero(gctimes),
        gctime_max_ns = maximum_or_zero(gctimes),
        allocs = trial.allocs,
        memory_bytes = trial.memory,
    )
end

function benchmark_helper_case!(rows, funcname::String, mode::String,
                                n::Int, elty::Type, model_name::String, bench_expr)
    trial = BenchmarkTools.run(bench_expr)
    stats = summarise_trial(trial)
    push!(rows, merge(stats, (
        function_name = funcname,
        mode = mode,
        sequence_length = n,
        element_type = string(elty),
        model = model_name,
    )))
end

function seqstruct_partition_cases!(rows, ::Type{T}, n::Int,
                                    model_name::String, model) where {T}
    p_base = base_p_seq(T, n)
    hairpin = DiffFoldRNA.HAIRPIN

    # Plain evaluation
    bench_plain = @benchmarkable DiffFoldRNA.seqstruct_partition($p_base, $model; hpmin=$hairpin)
    benchmark_helper_case!(rows, "seqstruct_partition", "value", n, T, model_name, bench_plain)

    # Single partial derivative via Dual numbers
    p_dual = partial_dual_matrix(p_base, PARTIAL_INDEX)
    bench_dual = @benchmarkable DiffFoldRNA.seqstruct_partition($p_dual, $model; hpmin=$hairpin)
    benchmark_helper_case!(rows, "seqstruct_partition", "partial_dual", n, eltype(p_dual), model_name, bench_dual)

    # Full gradient for shorter sequences
    if n in DEFAULT_GRADIENT_LENGTHS
        bench_grad = @benchmarkable seqstruct_gradient($p_base, $model; hpmin=$hairpin)
        benchmark_helper_case!(rows, "seqstruct_partition", "gradient", n, T, model_name, bench_grad)
    end
end

function helper_cases!(rows, ::Type{T}, n::Int,
                       model_name::String, model) where {T}
    p_base = base_p_seq(T, n)
    p_dual = partial_dual_matrix(p_base, PARTIAL_INDEX)
    modes = [
        ("value", T, p_base),
        ("partial_dual", eltype(p_dual), p_dual),
    ]
    bi, bj = 1, 2
    i = 1
    j = n

    for (mode_label, elty, source_matrix) in modes
        # fill_outer_mismatch!
        bench_outer = @benchmarkable DiffFoldRNA.fill_outer_mismatch!(OMM, p_seq, $model, $i, $n) setup = begin
            st = prepare_dp_state($source_matrix, $elty)
            (p_seq, E, P, ML, OMM) = (st.p_seq, st.E, st.P, st.ML, st.OMM)
        end evals = 1
        benchmark_helper_case!(rows, "fill_outer_mismatch!", mode_label, n, elty, model_name, bench_outer)

        # psum_bulges
        bench_bulges = @benchmarkable DiffFoldRNA.psum_bulges(P, p_seq, $model, $bi, $bj, $i, $j) setup = begin
            st = prepare_dp_state($source_matrix, $elty)
            (p_seq, E, P, ML, OMM) = (st.p_seq, st.E, st.P, st.ML, st.OMM)
        end evals = 1
        benchmark_helper_case!(rows, "psum_bulges", mode_label, n, elty, model_name, bench_bulges)

        # psum_internal_loops
        bench_internal = @benchmarkable DiffFoldRNA.psum_internal_loops(P, OMM, p_seq, $model, $bi, $bj, $i, $j) setup = begin
            st = prepare_dp_state($source_matrix, $elty)
            (p_seq, E, P, ML, OMM) = (st.p_seq, st.E, st.P, st.ML, st.OMM)
        end evals = 1
        benchmark_helper_case!(rows, "psum_internal_loops", mode_label, n, elty, model_name, bench_internal)

        # fill_paired!
        bench_paired = @benchmarkable DiffFoldRNA.fill_paired!(P, ML, OMM, p_seq, $model, $i, $n, DiffFoldRNA.HAIRPIN) setup = begin
            st = prepare_dp_state($source_matrix, $elty)
            (p_seq, E, P, ML, OMM) = (st.p_seq, st.E, st.P, st.ML, st.OMM)
        end evals = 1
        benchmark_helper_case!(rows, "fill_paired!", mode_label, n, elty, model_name, bench_paired)

        # fill_multi!
        bench_multi = @benchmarkable DiffFoldRNA.fill_multi!(ML, P, p_seq, $model, $i, $n) setup = begin
            st = prepare_dp_state($source_matrix, $elty)
            (p_seq, E, P, ML, OMM) = (st.p_seq, st.E, st.P, st.ML, st.OMM)
        end evals = 1
        benchmark_helper_case!(rows, "fill_multi!", mode_label, n, elty, model_name, bench_multi)

        # fill_external!
        bench_external = @benchmarkable DiffFoldRNA.fill_external!(E, P, p_seq, $model, $n, $i) setup = begin
            st = prepare_dp_state($source_matrix, $elty)
            (p_seq, E, P, ML, OMM) = (st.p_seq, st.E, st.P, st.ML, st.OMM)
        end evals = 1
        benchmark_helper_case!(rows, "fill_external!", mode_label, n, elty, model_name, bench_external)
    end
end

function run_all_benchmarks(; seq_lengths = DEFAULT_SEQ_LENGTHS,
                            output_csv::AbstractString = default_output_path())
    results = NamedTuple[]

    for T in (Float32, Float64)
        model_specs = MODEL_SPEC[T]
        for (model_name, constructor) in model_specs
            model = constructor()
            for n in seq_lengths
                seqstruct_partition_cases!(results, T, n, model_name, model)
                helper_cases!(results, T, n, model_name, model)
            end
        end
    end

    df = DataFrame(results)
    sort!(df, [:function_name, :mode, :element_type, :model, :sequence_length])

    mkpath(dirname(output_csv))
    CSV.write(output_csv, df; writeheader=true)

    show(stdout, MIME"text/plain", df; allrows=true, allcols=true, truncation=:none)
    println()

    return df
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_all_benchmarks()
end
