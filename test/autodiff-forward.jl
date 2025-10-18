using Test
using Random
using ForwardDiff
using DiffFoldRNA: All1Model, ViennaModel, NTS, get_RT_for_model, p_seq_to_str,
    make_loss_ptarget, make_grad!_fwd

@testset "autodiff-forward" begin
    showtestset()

    @testset "get_RT_for_model" begin
        all1_model = All1Model()
        @test get_RT_for_model(all1_model) == 1.0

        vienna_model = ViennaModel()
        @test get_RT_for_model(vienna_model) == vienna_model.boltz.RT_ustrip
    end

    @testset "p_seq_to_str" begin
        p_seq = [
            0.9 0.05 0.03 0.02;
            0.1 0.5 0.3 0.1;
            0.2 0.1 0.6 0.1;
        ]
        @test p_seq_to_str(p_seq) == "ACG"
    end

    @testset "make_loss_ptarget" begin
        model = All1Model()
        target_dbn = ".."
        n = length(target_dbn)
        x = fill(0.5, n * NTS)

        zero_loss = make_loss_ptarget(target_dbn, model, 0.0, 0.0)
        @test zero_loss(x) == 0.0

        loss = make_loss_ptarget(target_dbn, model)
        value = loss(x)
        @test isfinite(value)
    end

    @testset "make_grad!_fwd" begin
        Random.seed!(1)
        model = All1Model()
        target_dbn = ".."
        n = length(target_dbn)
        example_x = fill(0.1, n * NTS)
        loss = make_loss_ptarget(target_dbn, model)
        grad! = make_grad!_fwd(loss, example_x)
        x = rand(n * NTS)
        got = similar(x)
        grad!(got, x)
        expected = ForwardDiff.gradient(loss, x)
        @test got ≈ expected
    end
end
