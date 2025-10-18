using Random: randstring

const RNA_ALPHA = "ACGU"

# For the RandomModel of julia and python to match, it must be +2
# because of 0-based (python) vs 1-based indexing (julia) as this is
# passed as a base type index to float_hash, and we subtract 1 from
# all base type indices in julia to account for the change from
# 0-based to 1-based indexing
#
# TODO: make this nicer
#
# const INVALID_BASE = length(RNA_ALPHA) + 1
const INVALID_BASE = length(RNA_ALPHA) + 2

const NTS = length(RNA_ALPHA)
const ALL_PAIRS = ["AU", "UA", "GC", "CG", "GU", "UG"]
const ALL_PAIRS_TO_BASETYPES = [
    (findfirst(b1, RNA_ALPHA), findfirst(b2, RNA_ALPHA))::Tuple{Int,Int} for (b1, b2) in ALL_PAIRS
]
const NBPS = length(ALL_PAIRS)
const HAIRPIN = 0
const SPECIAL_HAIRPINS = ["CAACG", "GUUAC"]  # Not complete

# TODO: some functions missing

one_hot_seq(seq::AbstractString) = one_hot_seq(Float64, seq)
function one_hot_seq(T::Type, seq::AbstractString)
    n = length(seq)
    one_hot = zeros(T, n, NTS)
    for i in 1:n
        one_hot[i, findfirst(seq[i], RNA_ALPHA)] = one(T)
    end
    return one_hot
end

function valid_pair(a::Char, b::Char)
    # TODO: no 'N' handling, but see also valid_pair(as::AllStructs, i, j)
    a == 'A' && return b == 'U'
    a == 'C' && return b == 'G'
    a == 'G' && return (b == 'U' || b == 'C')
    a == 'U' && return (b == 'A' || b == 'G')
    return false
end

function is_valid_dbn(dbn::AbstractString)
    stk = Int[]
    for i in 1:length(dbn)
        if dbn[i] == '('
            push!(stk, i)
        elseif dbn[i] == ')'
            if length(stk) == 0
                return false
            end
            pop!(stk)
        elseif dbn[i] != '.'
            return false
        end
    end
    return true
end

function matching_to_dbn(match::Vector{Int})
    dbn = Char[]
    for i in 1:length(match)
        if match[i] < i
            push!(dbn, ')')
        elseif match[i] > i
            push!(dbn, '(')
        else
            push!(dbn, '.')
        end
    end
    return join(dbn)
end

function dbn_to_matching(dbn::AbstractString)
    n = length(dbn)
    stk = Int[]
    match = Int[1:n...]
    for i in 1:n
        if dbn[i] == '('
            push!(stk, i)
        elseif dbn[i] == ')'
            j = pop!(stk)
            match[j] = i
            match[i] = j
        end
    end
    return match
end

random_primary(sz::Integer) = randstring(RNA_ALPHA, sz)

function structure_tree(dbn::AbstractString)
    n = length(dbn)
    ch = Dict{Int,Vector{Int}}(-1 => Int[])
    stk = Int[-1]
    right = repeat([-1], n)
    for i in 1:n
        if dbn[i] == '('
            push!(stk, i)
        elseif dbn[i] == ')'
            left = pop!(stk)
            if stk[end] ∉ keys(ch)
                ch[stk[end]] = Int[]
            end
            push!(ch[stk[end]], left)
            right[left] = i
        end
    end
    return ch, right
end

# like structure_tree, but also include leaves
function structure_tree_full(dbn::AbstractString)
    n = length(dbn)
    ch = Dict{Int,Vector{Int}}(-1 => Int[])
    stk = Int[-1]
    right = repeat([-1], n)
    for i in 1:n
        if dbn[i] == '('
            push!(stk, i)
            ch[i] = Int[]
        elseif dbn[i] == ')'
            length(stk) > 0 || error("Too many closing parens at position $i")
            left = pop!(stk)
            push!(ch[stk[end]], left)
            right[left] = i
        end
    end
    stk == [-1] || error("Missing closing parens, opening parens are at $stk")
    return ch, right
end

function structure_list_postorder(dbn::AbstractString)
    ch, rightmatch = structure_tree_full(dbn)
    loops_postorder = sort!(collect(ch); rev=true)
    return loops_postorder, rightmatch
end

function seq_prob(p_seq::AbstractMatrix, seq::AbstractString)
    # TODO: assert length(p_seq) == length(seq)
    n = length(seq)
    if axes(p_seq) != (1:n, 1:NTS)
        throw(ArgumentError("must have axes(p_seq) == (1:$n, 1:$NTS), axes(p_seq) = $(axes(p_seq))"))
    end
    p = one(eltype(p_seq))
    for i in 1:length(seq)
        p *= p_seq[i, findfirst(seq[i], RNA_ALPHA)]
    end
    return p
end

# call with fn=x->exp(-x) for softmax normalisation
random_p_seq(n::Integer; fn=identity) = random_p_seq(Float64, n; fn)
function random_p_seq(T::Type, n::Integer; fn=identity)
    p = rand(T, n, NTS)
    p .= fn.(p) ./ sum(fn, p; dims=2)
    return p
end

normalize_to_p_seq(p_seq::AbstractMatrix) = normalize_to_p_seq(x -> x^2, p_seq)
normalize_to_p_seq!(p_seq::AbstractMatrix) = normalize_to_p_seq!(x -> x^2, p_seq)
normalize_to_p_seq(f::Function, p_seq::AbstractMatrix) = f.(p_seq) ./ sum(f, p_seq; dims=2)
normalize_to_p_seq!(f::Function, p_seq::AbstractMatrix) = p_seq .= f.(p_seq) ./ sum(f, p_seq; dims=2)

function random_seq_for_dbn(dbn::AbstractString)
    n = length(dbn)
    match = dbn_to_matching(dbn)
    seq = repeat(['.'], n)
    axes(match) == (1:n,) || error("axes(match) != 1:n, axes(match) = $(axes(match))")
    for i in eachindex(match)
        j = match[i]
        if i == j
            seq[i] = rand(RNA_ALPHA)
        elseif i < j
            bi, bj = rand(ALL_PAIRS)
            seq[i] = bi
            seq[j] = bj
        end
    end
    if any(==('.'), seq)
        error("unassigned bases, seq = $seq, match = $match")
    end
    return join(seq)
end
