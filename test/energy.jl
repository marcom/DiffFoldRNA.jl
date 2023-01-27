using Test
using DiffFoldRNA

SEQ_DBN = [
    ("GAGAGGAAACCACC",
     "(.(.((...)).))"),
    ("GAAGAAACAAC",
     "(..(...)..)"),
]

@testset "energy" begin
    showtestset()

    # All1Model, ViennaModel, RandomModel, RandomExtloopModel,
    # RandomMultiloopModel, RandomHairpinModel, RandomBulgeModel,
    # RandomStackModel, RandomILModel
    T = Float64
    for model in [All1Model(), All1Model{T}(), ViennaModel(),
                  [m for M in (RandomModel, RandomExtloopModel, RandomMultiloopModel,
                               RandomHairpinModel, RandomBulgeModel, RandomStackModel,
                               RandomILModel)
                       for m in (M(), M{T}(), M(42), M{T}(42)) ]...,
                  ]
        @testset "$model" begin
            for (seq, dbn) in SEQ_DBN
                @test boltz(seq, dbn, model) isa eltype(model)
            end
        end
    end

    # DebugModel, supress output
    model = DebugModel(RandomModel())
    @testset "$model" begin
        redirect_stdio(stdout=devnull, stderr=devnull) do
            for (seq, dbn) in SEQ_DBN
                @test boltz("GAAAC", "(...)", model) isa eltype(model)
            end
        end
    end

    # test other types
    T = Int64
    seq = "GAAGAAACAAC"
    dbn = "(..(...)..)"
    @test boltz(seq, dbn, All1Model{T}()) isa T

    # TODO: RandomModel broken for Int64
    @test_broken boltz(seq, dbn, RandomModel{T}()) isa T
end
