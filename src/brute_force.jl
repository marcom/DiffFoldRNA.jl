# brute-force calculation of single-sequence partition function
#
# Note: boltz_fn(seq, matching) is the Boltzmann factor,
# with boltz_fn(seq,matching) == exp(-beta * energy(seq, matching))
allstructs_boltz_sum(seq::AbstractString, boltz_fn::Function; hpmin::Int=HAIRPIN) =
    allstructs_boltz_sum(Float64, seq, boltz_fn; hpmin)
function allstructs_boltz_sum(T::Type, seq::AbstractString, boltz_fn::Function; hpmin::Int=HAIRPIN)
    as = AllStructs(seq, hpmin)
    sm = zero(T)
    # TODO: this assumes get_nth() counts structures starting from 0
    for i in 0:count_structures(as)-1
        sm += boltz_fn(seq, get_nth(as, i))
    end
    return sm
end

seqstruct_partition_brute_force(p_seq::AbstractMatrix, boltz_fn::Function; hpmin::Int=HAIRPIN) =
    seqstruct_partition_brute_force(Float64, p_seq, boltz_fn; hpmin)
function seqstruct_partition_brute_force(T::Type, p_seq::AbstractMatrix, boltz_fn::Function; hpmin::Int=HAIRPIN)
    n = size(p_seq, 1)

    function f(seq_list)
        if length(seq_list) == n
            seq = join(seq_list)
            return seq_prob(p_seq, seq) * allstructs_boltz_sum(T, seq, boltz_fn; hpmin)
        end
        sm = zero(T)
        for b in RNA_ALPHA
            # TODO: a lot of copying?
            new_seq_list = push!(copy(seq_list), b)
            sm += f(new_seq_list)
        end
        return sm
    end

    return f(Char[])
end

function seq_partition_brute_force(p_seq::AbstractMatrix, struc::AbstractString,
                                   boltz_fn::Function)
    return seq_partition_brute_force(Float64, p_seq, struc, boltz_fn)
end
function seq_partition_brute_force(Ten::Type, p_seq::AbstractMatrix{Tp},
                                   struc::AbstractString, boltz_fn::Function) where Tp
    n = size(p_seq, 1)
    axes(p_seq) == (1:n, 1:NTS) || throw(ArgumentError("axes(p_seq) != (1:$n, 1:$NTS)"))
    length(struc) == n || throw(ArgumentError("p_seq and structure must have same length ($n != $(length(struc)))"))
    match = dbn_to_matching(struc)

    function f(seq_list::Vector)
        if length(seq_list) == n
            seq = join(seq_list)
            return seq_prob(p_seq, seq) * boltz_fn(seq, struc)
        end
        sm = zero(Ten)
        idx = length(seq_list) + 1  # +1 because of 1-based indexing
        for b in RNA_ALPHA
            if match[idx] < idx &&  !valid_pair(b, seq_list[match[idx]])
                continue
            end
            # TODO: a lot of copying?
            new_seq_list = push!(copy(seq_list), b)
            sm += f(new_seq_list)
        end
        return sm
    end

    return f([])
end
