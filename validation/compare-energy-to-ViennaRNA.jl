import Pkg
Pkg.activate(@__DIR__)

using Test
using DiffFoldRNA: AllStructs, RNA_ALPHA, ViennaModel, boltz, count_structures, get_nth, matching_to_dbn
using ArgParse, CSV, DataFrames, Unitful, ProgressMeter, Random
import FoldRNA: randseq
import ViennaRNA

const DBN_HAIRPINS   = ["(" * "."^n * ")" for n = 3:32]
const DBN_INTLOOPS   = vec(["(" * "."^n1 * "(...)" * "."^n2 * ")"
                            for n1 in [0:5..., 29], n2 in [0:5..., 29]])
const DBN_EXTLOOPS   = ["(...)(...)", "(...).(...)", "(...)..(...)", ".(...).(...)(...).(...)"]
const DBN_MULTILOOPS = ["(.(...)(...))", "((...).(...)((...)(...).))"]
const DBN_MISC = [
    "(((...)))",
    "(((....)))",
    "((.(...)).)",
]
const DBN_ALL = vcat(DBN_HAIRPINS, DBN_INTLOOPS, DBN_EXTLOOPS, DBN_MULTILOOPS, DBN_MISC)

function test_energy_vrna(seq::AbstractString, dbn::AbstractString, model; atol=1e-7u"kcal/mol", df)
    RT = model.boltz.RT
    en_jl   = -RT * log(boltz(seq, dbn, model))
    en_vrna = ViennaRNA.energy(seq, dbn)
    @test en_jl ≈ en_vrna atol=atol
    if ! isapprox(en_jl, en_vrna; atol)
        push!(df, (; len=length(seq), seq, dbn, en_jl, en_vrna, delta=en_jl-en_vrna))
    end
end

function main(args)
    nseq = 50
    nseq_per_dbn = 20
    ndbn_per_seq = 20
    seqlen = 30
    atol = 1e-2u"kcal/mol"
    hpmin = 3

    Ten = typeof(1.0u"kcal/mol")
    @testset "compare energy(seq, dbn) vs ViennaRNA, atol=$atol" begin
        df = DataFrame(len=Int[], seq=String[], dbn=String[], en_jl=Ten[], en_vrna=Ten[], delta=Ten[])

        @testset "fixed dbn, random seq" verbose=true begin
            model = ViennaModel()
            @showprogress for dbn in DBN_ALL
                for _ in 1:nseq_per_dbn
                    seq = randseq(dbn)
                    test_energy_vrna(seq, dbn, model; atol, df)
                end
            end
        end

        @testset "random seq, random dbn" begin
            model = ViennaModel()
            @showprogress for i = 1:nseq
                seq = randstring(RNA_ALPHA, seqlen)
                as = AllStructs(seq, hpmin)
                ndbn = count_structures(as)
                for j = 1:ndbn_per_seq
                    match = get_nth(as, rand(0:ndbn-1))
                    dbn = matching_to_dbn(match)
                    test_energy_vrna(seq, dbn, model; atol, df)
                end
            end
        end
        display(df)
    end

    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
