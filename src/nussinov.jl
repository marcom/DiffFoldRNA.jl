using OffsetArrays
using LogExpFunctions: softmax

function make_valid_pairs_table()
    valid_pairs = zeros(Bool, NTS, NTS)
    for i in 1:NTS, j in 1:NTS
        bp = join((RNA_ALPHA[i], RNA_ALPHA[j]))
        valid_pairs[i, j] = bp in ALL_PAIRS
    end
    return valid_pairs
end

nussinov_boltz_pair_1(bi, bj) = 1
nussinov_boltz_pair(bi, bj) = bi*37+bj

function nussinov_seqstruct_partition(p_seq, boltz_pair; hpmin::Int=HAIRPIN)
    HAIRPIN = hpmin
    n = size(p_seq,1)
    axes(p_seq,2) == 1:NTS || error("axes(p_seq,2) == $(axes(p_seq,2)), expected 1:$NTS")
    D = OffsetMatrix(ones(n+1, n), 0:n, 0:n-1)
    valid_pairs = make_valid_pairs_table()
    for i in n-1:-1:0, j in i:n-1  # py: i in range(n-1, -1, -1), j in range(i, n)
        sm = D[i+1, j]
        for k in i+HAIRPIN+1:j, bi in 1:NTS, bk in 1:NTS  # py: range(i+HAIRPIN+1, j+1)
            valid_pairs[bi, bk] || continue
            sm += D[i+1, k-1] *
                D[k+1, j] *
                p_seq[i+1, bi] *  # TODO: i+1 because p_seq is 1-based, D 0-based
                p_seq[k+1, bk] *  # TODO: i+1 because p_seq is 1-based, D 0-based
                boltz_pair(bi, bk)
        end
        D[i, j] = sm
    end
    return D[0, n-1]
end

function nussinov_boltz(seq::AbstractString, match::AbstractVector{Int})
    length(seq) == length(match) ||
        error("length and match must have same length, $(length(seq)) != $(length(match))")
    e = 1
    for i in 1:length(match)  # py: range(len(match)):
        if match[i] > i
            e *= nussinov_boltz_pair(findfirst(seq[i], RNA_ALPHA),
                                     findfirst(seq[match[i]], RNA_ALPHA))
        end
    end
    return e
end

function nussinov_test(n::Integer=9)
    p_seq = rand(n, NTS)  # TODO: change to (NTS,n)
    p_seq = softmax(p_seq; dims=2)
    @show nussinov_seqstruct_partition(p_seq, nussinov_boltz_pair)
    @show seqstruct_partition_brute_force(Float64, p_seq, nussinov_boltz)
end

