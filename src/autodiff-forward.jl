import ForwardDiff
using Optim: Optim, optimize, GradientDescent, ConjugateGradient
using UnicodePlots: heatmap, lineplot
using Statistics: mean
using LogExpFunctions: xlogx
using LinearAlgebra: norm
import ViennaRNA

function get_RT_for_model(model::AbstractModel)
    RT = if model isa ViennaModel
        model.boltz.RT_ustrip::Float64
    else
        1.0
    end
    return RT
end

function p_seq_to_str(p_seq::AbstractMatrix)
    n = size(p_seq, 1)
    seq = repeat(['.'], n)
    for (i, I) in enumerate(last(findmax(p_seq; dims=2)))
        seq[i] = RNA_ALPHA[last(Tuple(I))]
    end
    return join(seq)
end

# TODO: maxintloop
function make_loss_ptarget(target_dbn::AbstractString, model::AbstractModel, αₚ=1.0, αₙ=1.0;
                           p_norm_fn = x -> x^2)
    RT = get_RT_for_model(model)
    n = length(target_dbn)
    nb = NTS
    function loss(x::Vector)
        p_seq = reshape(x, n, nb)
        normalize_to_p_seq!(p_norm_fn, p_seq)
        return (αₚ * (-RT) * log(seq_partition(p_seq, target_dbn, model))
                - αₙ * (-RT) * log(seqstruct_partition(p_seq, model; hpmin=3)))
    end
    return loss
end

function make_grad!_fwd(f, example_x)
    cfg = ForwardDiff.GradientConfig(f, example_x, ForwardDiff.Chunk(8))
    grad! = (out, x) -> ForwardDiff.gradient!(out, f, x, cfg)
    return grad!
end

function design_ptarget(target_dbn::AbstractString, model::AbstractModel; verbose=true, optim_options=(;),
                        p_norm_fn = x -> x^2)
    entropy_base = 2.0 # unit is bits
    RT = get_RT_for_model(model)
    n = length(target_dbn)
    loss = make_loss_ptarget(target_dbn, model; p_norm_fn)
    grad! = make_grad!_fwd(loss, zeros(n * NTS))
    initial_x = ones(n * NTS) / NTS
    @show loss(initial_x)

    nstep = 0
    function opt_callback(opt)
        @show nstep
        x = reshape(convert(Array{Float64}, opt.metadata["x"]), n, NTS)
        p_seq = normalize_to_p_seq(p_norm_fn, x)
        entropy = reshape(sum(xlogx, p_seq; dims=2), n) ./ -log(entropy_base)
        if verbose
            g_x = reshape(convert(Array{Float64}, opt.metadata["g(x)"]), n, NTS)
            @show opt.metadata["Current step size"], opt.metadata["time"]
            @show norm(g_x)
            @show p_target_p_seq = exp(opt.value/-RT)
            @show p_target_p_seq_argmax = ViennaRNA.prob_of_structure(p_seq_to_str(p_seq), target_dbn)
            println("entropy")
            @show mean(entropy)
            display(lineplot(entropy; ylim=(0.0,entropy_base), height=4, width=n))
            println("x")
            display(heatmap(x'; height=4, width=n))
            println("pseq")
            display(heatmap(p_seq'; height=4, width=n, xlabel=p_seq_to_str(p_seq)))
            println("g(x)")
            display(heatmap(g_x'; height=4, width=n, xlabel=target_dbn))


            #@show mean_dist_to_pure_seq(fold.pseq)
            #@show max_dist_to_pure_seq(fold.pseq)
            #return opt.value <= loss_stop ? true : false
        end
        nstep += 1
        return mean(entropy) < 0.4
    end

    optim_method = GradientDescent()
    #optim_method = ConjugateGradient()

    opt = optimize(loss, grad!, initial_x, optim_method,
                   Optim.Options(; callback=opt_callback, show_trace=false, extended_trace=true, optim_options...))
    display(opt)
    sol = Optim.minimizer(opt)
    @show sol_p_seq = normalize_to_p_seq(p_norm_fn, reshape(sol, n, NTS))
    @show exp(loss(sol)/-RT)
    @show p_seq_to_str(sol_p_seq)
    @show p_target_p_seq_argmax = ViennaRNA.prob_of_structure(p_seq_to_str(sol_p_seq), target_dbn)

    return sol_p_seq
end

