# TODO
# - faster
#   - get_bp_bases: lookups
#   - pr_special_hairpin: lookups
#   - make pseq more efficient (not vector of vectors), matrix?
#     i think SciML has some kind of VectorOfVectors
#   - arrays: better to do column-major access in julia, numpy is row-major by default
# - lots of hardcoded 1:NTS, 1:NBPS index ranges

using OffsetArrays

# bptype -> (basetype, basetype)
@inline get_bp_bases(bp::Integer) = ALL_PAIRS_TO_BASETYPES[bp]

# NOTE: python indexing starts at 0, so we use OffsetArrays to use the
# same indexing

function psum_hairpin_not_special(p_seq::AbstractMatrix{Tp}, em::M,
                                  bi::Integer, bj::Integer, i::Integer, j::Integer) where {Tp,Te,M<:AbstractModel{Te}}
    T = promote_type(Tp,Te)
    # Special case for HAIRPIN<=1
    # Necessary to respect conditional probability of the mismatch
    # Can be removed or made branchless/jax.where
    if i + 1 == j
        return en_hairpin_not_special(em, bi, bj, bj, bi, 0)::T
    end
    sm = zero(T)::T
    if i+1 == j-1
        @inbounds for bip1 in 1:NTS
            sm += p_seq[i+1, bip1] *
                en_hairpin_not_special(em, bi, bj, bip1, bip1, 1)
        end
        return sm::T
    end
    @inbounds for bip1 in 1:NTS, bjm1 in 1:NTS
        sm += p_seq[i+1, bip1] * p_seq[j-1, bjm1] *
            en_hairpin_not_special(em, bi, bj, bip1, bjm1, j-i-1)
    end
    return sm::T
end

function pr_special_hairpin(p_seq::AbstractMatrix{Tp}, em::M,
                            id::Integer, i::Integer, j::Integer) where {Tp,Te,M<:AbstractModel{Te}}
    T = promote_type(Tp, Te)
    SPECIAL_HAIRPINS = specialhairpins(em)
    pr = one(T)
    @inbounds for k in i+1:j-1
        # TODO: can be replaced with a table
        # TODO: changed k-i => k-i+1 (julia 1-based indexing)
        # Note: SPECIAL_HAIRPINS second index ∈ (2:j-i), length j-i-1
        pr *= p_seq[k, findfirst(SPECIAL_HAIRPINS[id][k-i+1], RNA_ALPHA)::Int]
    end
    return pr::T
end

function psum_hairpin_special(p_seq::AbstractMatrix{Tp}, em::M,
                              bi::Integer, bj::Integer, i::Integer, j::Integer) where {Tp,Te,M<:AbstractModel{Te}}
    T = promote_type(Tp, Te)
    SPECIAL_HAIRPINS = specialhairpins(em)
    sm = zero(T)::T
    @inbounds for id in LinearIndices(SPECIAL_HAIRPINS)
        # Can be made branchless
        if SPECIAL_HAIRPINS[id][begin] != RNA_ALPHA[bi]
            continue
        end
        if SPECIAL_HAIRPINS[id][end] != RNA_ALPHA[bj]
            continue
        end
        if length(SPECIAL_HAIRPINS[id]) != j - i + 1
            continue
        end
        sm += pr_special_hairpin(p_seq, em, id, i, j) * en_hairpin_special(em, id)
    end
    return sm::T
end

function psum_hairpin_special_correction(p_seq::AbstractMatrix{Tp}, em::M,
                                         bi::Integer, bj::Integer, i::Integer, j::Integer) where {Tp,Te,M<:AbstractModel{Te}}
    T = promote_type(Tp, Te)
    SPECIAL_HAIRPINS = specialhairpins(em)
    sm = zero(T)::T
    @inbounds for id in LinearIndices(SPECIAL_HAIRPINS)
        # Can be made branchless
        if SPECIAL_HAIRPINS[id][begin] != RNA_ALPHA[bi]
            continue
        end
        if SPECIAL_HAIRPINS[id][end] != RNA_ALPHA[bj]
            continue
        end
        if length(SPECIAL_HAIRPINS[id]) != j-i+1
            continue
        end
        # TODO: can be replace with a table
        # TODO: assumes length(SPECIAL_HAIRPINS[id]) >= 2, which should always be the case
        #       move this test to beginning of ss_partition
        bip1 = findfirst(SPECIAL_HAIRPINS[id][begin+1], RNA_ALPHA)::Int
        bjm1 = findfirst(SPECIAL_HAIRPINS[id][end-1], RNA_ALPHA)::Int
        sm += pr_special_hairpin(p_seq, em, id, i, j) *
            en_hairpin_not_special(em, bi, bj, bip1, bjm1, length(SPECIAL_HAIRPINS[id])-2)
    end
    return sm::T
end

# This can be precomputed
function psum_hairpin(p_seq::AbstractMatrix{Tp}, em::M,
                      bi::Integer, bj::Integer, i::Integer, j::Integer) where {Tp,Te,M<:AbstractModel{Te}}
    T = promote_type(Tp, Te)
    return psum_hairpin_not_special(p_seq, em, bi, bj, i, j)::T +
        psum_hairpin_special(p_seq, em, bi, bj, i, j)::T -
        psum_hairpin_special_correction(p_seq, em, bi, bj, i, j)::T
end

@inline bool_to_T(::Type{T}, cond::Bool) where T = cond ? one(T) : zero(T)

@inline function fill_external!(
    E::AbstractArray{T},
    P::AbstractArray{T},
    p_seq::AbstractMatrix{Tp},
    em::AbstractModel,
    n::Integer,
    i::Integer,
) where {T, Tp}
    @inbounds for bim1 in 1:NTS, bi in 1:NTS
        sm = zero(T)
        for bip1 in 1:NTS
            sm += (E[bi, bip1, i+1] + bool_to_T(T, i == n)) * p_seq[i+1, bip1]
        end
        for j in i+1:n, bj in 1:NTS, bjp1 in 1:NTS
            dangle5 = (i == 1) ? INVALID_BASE : bim1
            dangle3 = (j == n) ? INVALID_BASE : bjp1
            sm += P[bi, bj, i, j] *
                (E[bj, bjp1, j+1] + bool_to_T(T, j == n)) *
                en_ext_branch(em, dangle5, bi, bj, dangle3) *
                p_seq[j, bj] *
                p_seq[j+1, bjp1]
        end
        E[bim1, bi, i] = sm
    end
    return nothing
end

@inline function fill_multi!(
    ML::AbstractArray{T},
    P::AbstractArray{T},
    p_seq::AbstractMatrix{Tp},
    em::AbstractModel,
    i::Integer,
    n::Integer,
) where {T,Tp}
    @inbounds for bim1 in 1:NTS, bi in 1:NTS, bj in 1:NTS, bjp1 in 1:NTS, nb in 0:2, j in i:n
        nbm1_min0 = max(0, nb-1)
        sm = zero(T)
        for bip1 in 1:NTS
            sm += (ML[bi, bip1, bj, bjp1, nb, i+1, j] + bool_to_T(T, (i+1 > j) && (nb == 0))) *
                p_seq[i+1, bip1]
        end
        if nb <= 1
            sm += P[bi, bj, i, j] *
                en_multi_branch(em, bim1, bi, bj, bjp1)
        end
        for bk in 1:NTS
            if nb <= 1
                sm += P[bi, bk, i, j-1] *
                    en_multi_branch(em, bim1, bi, bk, bj) *
                    p_seq[j-1, bk]
            end
        end
        for k in i:j-2
            for bk in 1:NTS, bkp1 in 1:NTS
                sm += P[bi, bk, i, k] *
                    ML[bk, bkp1, bj, bjp1, nbm1_min0, k+1, j] *
                    en_multi_branch(em, bim1, bi, bk, bkp1) *
                    p_seq[k, bk] *
                    p_seq[k+1, bkp1]
            end
        end
        ML[bim1, bi, bj, bjp1, nb, i, j] = sm
    end
    return nothing
end

@inline function fill_outer_mismatch!(
    OMM::AbstractArray{T},
    p_seq::AbstractMatrix{Tp},
    em::AbstractModel,
    k::Integer,
    n::Integer,
) where {T,Tp}
    @inbounds for l in k+1:n, bpkl in 1:NBPS
        bk, bl = get_bp_bases(bpkl)
        for bkm1 in 1:NTS, blp1 in 1:NTS
            OMM[bk, bl, k, l] += en_il_outer_mismatch(em, bk, bl, bkm1, blp1) *
                p_seq[k-1, bkm1] *
                p_seq[l+1, blp1]
        end
    end
    return nothing
end

@inline function psum_internal_loops(
    P::AbstractArray{T},
    OMM::AbstractArray{T},
    p_seq::AbstractMatrix{Tp},
    em::AbstractModel,
    bi::Integer,
    bj::Integer,
    i::Integer,
    j::Integer,
) where {T,Tp}
    sm = zero(T)
    mmij = zero(T)
    @inbounds for bip1 in 1:NTS, bjm1 in 1:NTS
        mmij += p_seq[i+1, bip1] *
            p_seq[j-1, bjm1] *
            en_il_inner_mismatch(em, bi, bj, bip1, bjm1)
    end
    @inbounds for bpkl in 1:NBPS
        bk, bl = get_bp_bases(bpkl)
        for bip1 in 1:NTS, bjm1 in 1:NTS
            pr_ij_mm = p_seq[i+1, bip1] * p_seq[j-1, bjm1]
            sm += P[bk, bl, i+2, j-2] *
                p_seq[i+2, bk] *
                p_seq[j-2, bl] *
                pr_ij_mm *
                en_internal(em, bi, bj, bk, bl, bip1, bjm1, bip1, bjm1, 1, 1)
            for z in i+3:j-3
                for b in 1:NTS
                    il_en = en_internal(em, bi, bj, bk, bl, bip1, bjm1, bip1, b, 1, j-z-1)
                    sm += P[bk, bl, i+2, z] *
                        p_seq[i+2, bk] *
                        p_seq[z, bl] *
                        p_seq[z+1, b] *
                        pr_ij_mm *
                        il_en
                    il_en = en_internal(em, bi, bj, bk, bl, bip1, bjm1, b, bjm1, z-i-1, 1)
                    sm += P[bk, bl, z, j-2] *
                        p_seq[z, bk] *
                        p_seq[j-2, bl] *
                        p_seq[z-1, b] *
                        pr_ij_mm *
                        il_en
                end
            end
        end
        for k in i+2:j-3, l in k+1:j-2
            lup, rup = k-i-1, j-l-1
            if (lup <= 1) || (rup <= 1)
                continue
            end
            if (lup == 2 && rup == 2) || (lup == 2 && rup == 3) || (lup == 3 && rup == 2)
                for bip1 in 1:NTS, bjm1 in 1:NTS, bkm1 in 1:NTS, blp1 in 1:NTS
                    sm += P[bk, bl, k, l] *
                        p_seq[k, bk] *
                        p_seq[l, bl] *
                        en_internal(em, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup) *
                        p_seq[k-1, bkm1] *
                        p_seq[l+1, blp1] *
                        p_seq[i+1, bip1] *
                        p_seq[j-1, bjm1]
                end
            else
                init_and_pair = en_internal_init(em, lup+rup) *
                    en_internal_asym(em, abs(lup-rup)) *
                    P[bk, bl, k, l] *
                    p_seq[k, bk] *
                    p_seq[l, bl]
                sm += OMM[bk, bl, k, l] * mmij * init_and_pair
            end
        end
    end
    return sm
end

@inline function psum_bulges(
    P::AbstractArray{T},
    p_seq::AbstractMatrix{Tp},
    em::AbstractModel,
    bi::Integer,
    bj::Integer,
    i::Integer,
    j::Integer,
) where {T,Tp}
    sm = zero(T)
    @inbounds for bpkl in 1:NBPS
        bk, bl = get_bp_bases(bpkl)
        for kl in i+2:j-2
            sm += P[bk, bl, i+1, kl] *
                p_seq[i+1, bk] *
                p_seq[kl, bl] *
                en_bulge(em, bi, bj, bk, bl, j-kl-1)
            sm += P[bk, bl, kl, j-1] *
                p_seq[kl, bk] *
                p_seq[j-1, bl] *
                en_bulge(em, bi, bj, bk, bl, kl-i-1)
        end
    end
    return sm
end

@inline function fill_paired!(
    P::AbstractArray{T},
    ML::AbstractArray{T},
    OMM::AbstractArray{T},
    p_seq::AbstractMatrix{Tp},
    em::AbstractModel,
    i::Integer,
    n::Integer,
    hairpin_min::Integer,
) where {T,Tp}
    @inbounds for bpij in 1:NBPS, j in i+hairpin_min+1:n
        bi, bj = get_bp_bases(bpij)
        sm = psum_hairpin(p_seq, em, bi, bj, i, j)::T
        sm += psum_bulges(P, p_seq, em, bi, bj, i, j)
        sm += psum_internal_loops(P, OMM, p_seq, em, bi, bj, i, j)
        for bpkl in 1:NBPS
            bk, bl = get_bp_bases(bpkl)
            sm += P[bk, bl, i+1, j-1] *
                p_seq[i+1, bk] *
                p_seq[j-1, bl] *
                en_stack(em, bi, bj, bk, bl)
        end
        for bip1 in 1:NTS, bjm1 in 1:NTS
            sm += ML[bi, bip1, bjm1, bj, 2, i+1, j-1] *
                p_seq[i+1, bip1] *
                p_seq[j-1, bjm1] *
                en_multi_closing(em, bi, bip1, bjm1, bj)
        end
        P[bi, bj, i, j] = sm
    end
    return nothing
end


function seq_partition(p_seq::AbstractMatrix{Tp}, dbn::AbstractString,
                       em::M) where {Tp,Te,M<:AbstractModel{Te}}
    # Very cheesy hack to get a fast(er) O(N^2) sequence partition function
    # Linear time is possible
    n = length(dbn)
    axes(p_seq, 1) == 1:n || error("p_seq has wrong axes (axes(p_seq) = $(axes(p_seq))")
    dbn = OffsetArray(collect(dbn), 0:n-1)

    # This is all precomp that happens before compilation
    match = OffsetArray(collect(0:n), 0:n)
    stk = Int[]
    for i in LinearIndices(dbn)
        if dbn[i] == '('
            push!(stk, i)
        elseif dbn[i] == ')'
            j = pop!(stk)
            match[j+1] = i+1
        end
    end
    # Precomp unpaired ranges
    up = OffsetArray(ones(Bool, n+2, n+2), 0:n+1, 0:n+1)  # np.ones((n+2, n+2), dtype=bool)
    for i in 0:n  # range(n+1):
        for j in i+2:n  #  range(i+2, n+1):
            up[i, j] = (dbn[j-2] == '.') * up[i, j-1]
        end
    end

    # pad p_seq so indexing one out of bounds will be ok (left and right)
    #
    # NOTE: the python code uses 1-based indexing in a 0-based array
    # for this (presumably to avoid index -1, which in python means
    # "last array index"). We recreate this indexing with a
    # OffsetMatrix, but ONLY FOR THE FIRST INDEX, i.e. index `i` in
    # p_seq[i, j].
    p_padded = OffsetMatrix(zeros(Tp, n+2, NTS), 0:n+1, 1:NTS) # TODO: indices 0:n+1, 1:NTS
    p_padded[1:n, 1:NTS] .= p_seq
    # Base at 0 and n+1 is always the 0-th base
    p_padded[begin, begin] = one(Tp)
    p_padded[end, begin] = one(Tp)

    p_seq = p_padded

    # NOTE: NTS indices are 1:NTS, sequence indices are 0:n+1
    # external loop
    E = OffsetArray(zeros(Tp, NTS, NTS, n+2), 1:NTS, 1:NTS, 0:n+1)
    # P[i,j]: (i,j) are basepaired, i.e. close a loop
    P = OffsetArray(zeros(Tp, NTS, NTS, n+2, n+2), 1:NTS, 1:NTS, 0:n+1, 0:n+1)
    # ML: multiloop tables (3)
    # three tables needed: 0 stems, 1 stem, >= 2 stems (TODO: check)
    ML = OffsetArray(zeros(Tp, NTS, NTS, NTS, NTS, 3, n+2, n+2), 1:NTS, 1:NTS, 1:NTS, 1:NTS, 0:2, 0:n+1, 0:n+1)
    OMM = OffsetArray(zeros(Tp, NTS, NTS, n+2, n+2), 1:NTS, 1:NTS, 0:n+1, 0:n+1)

    # TODO: check again
    @inline bool_to_Te(cond::Bool) = cond ? one(Te) : zero(Te)

    function fill_external(i::Integer)
        for bim1 in 1:NTS, bi in 1:NTS
            sm = zero(Te)
            j = match[i]
            if j == i
                for bip1 in 1:NTS
                    # TODO: this assumes +(::Te, ::Bool) exists
                    sm += (E[bi, bip1, i+1] + bool_to_Te(i == n)) * p_seq[i+1, bip1]
                end
            else
                for bj in 1:NTS, bjp1 in 1:NTS  # j in i+1:n
                    # These can be where instead. Currently branchless arithmetic.
                    dangle5 = (i == 1)*INVALID_BASE + (i != 1)*bim1
                    dangle3 = (j == n)*INVALID_BASE + (j != n)*bjp1
                    sm += P[bi, bj, i, j] * (E[bj, bjp1, j+1] + bool_to_Te(j == n)) *
                        en_ext_branch(em, dangle5, bi, bj, dangle3) *
                        p_seq[j, bj] * p_seq[j+1, bjp1]
                end
            end
            E[bim1, bi, i] = sm
        end
        return
    end

    function fill_multi(i::Integer)
        # TODO: 0:2 hardcoded, 3 tables
        for bim1 in 1:NTS, bi in 1:NTS, bj in 1:NTS, bjp1 in 1:NTS, nb in 0:2, j in i:n
            # The max function can be made branchless
            nbm1_min0 = max(0, nb-1)
            sm = zero(Te)
            k = match[i]
            # Replace if/elif/else with where
            if k == i
                for bip1 in 1:NTS
                    sm += (ML[bi, bip1, bj, bjp1, nb, i+1, j] + bool_to_Te(i+1 > j && nb == 0)) *
                        p_seq[i+1, bip1]
                end
            elseif k == j
                # Special case for k==j
                if nb <= 1
                    sm += P[bi, bj, i, j] *
                        en_multi_branch(em, bim1, bi, bj, bjp1)
                end
            elseif k == j-1
                for bk in 1:NTS
                    # Special case for k==j-1
                    if nb <= 1
                        sm += P[bi, bk, i, j-1] *
                            en_multi_branch(em, bim1, bi, bk, bj) * p_seq[j-1, bk]
                    end
                end
            else
                for bk in 1:NTS, bkp1 in 1:NTS
                    sm += P[bi, bk, i, k] * ML[bk, bkp1, bj, bjp1, nbm1_min0, k+1, j] *
                        en_multi_branch(em, bim1, bi, bk, bkp1) *
                        p_seq[k, bk] *
                        p_seq[k+1, bkp1]
                end
            end
            ML[bim1, bi, bj, bjp1, nb, i, j] = sm
        end
        return
    end

    function fill_outer_mismatch(k::Integer)
        # TODO was: range(k+1, n+1):
        for l in k+1:n, bpkl in 1:NBPS
            bk, bl = get_bp_bases(bpkl)
            for bkm1 in 1:NTS, blp1 in 1:NTS
                OMM[bk, bl, k, l] += en_il_outer_mismatch(em, bk, bl, bkm1, blp1) *
                    p_seq[k-1, bkm1] *
                    p_seq[l+1, blp1]
            end
        end
        return
    end

    function psum_internal_loops(Te::Type, bi::Integer, bj::Integer, i::Integer, j::Integer)
        sm = zero(Te)
        mmij = zero(Te)
        for bip1 in 1:NTS, bjm1 in 1:NTS
            mmij += p_seq[i+1, bip1] *
                p_seq[j-1, bjm1] *
                en_il_inner_mismatch(em, bi, bj, bip1, bjm1)
        end
        for bpkl in 1:NBPS
            bk, bl = get_bp_bases(bpkl)
            # Both of these loops can be factored out
            for bip1 in 1:NTS, bjm1 in 1:NTS
                pr_ij_mm = p_seq[i+1, bip1] * p_seq[j-1, bjm1]
                # 1x1. TODO meaning: Don't need safe_P since we pad on both sides
                sm += P[bk, bl, i+2, j-2] *
                    p_seq[i+2, bk] *
                    p_seq[j-2, bl] *
                    pr_ij_mm *
                    en_internal(em, bi, bj, bk, bl, bip1, bjm1, bip1, bjm1, 1, 1) *
                    up[i, i+2] * up[j-2, j]
                # 1xn (n>=2)
                # TODO was: range(i+3, j-2):
                for z in i+3:j-3
                    # TODO: This loop could be optimised with the
                    # mismatch trick, though probably not worth it
                    for b in 1:NTS
                        il_en = en_internal(em, bi, bj, bk, bl, bip1, bjm1, bip1, b, 1, j-z-1)
                        sm += P[bk, bl, i+2, z] *
                            p_seq[i+2, bk] *
                            p_seq[z, bl] *
                            p_seq[z+1, b] *
                            pr_ij_mm *
                            il_en * up[i, i+2] * up[z, j]
                        il_en = en_internal(em, bi, bj, bk, bl, bip1, bjm1, b, bjm1, z-i-1, 1)
                        sm += P[bk, bl, z, j-2] *
                            p_seq[z, bk] *
                            p_seq[j-2, bl] *
                            p_seq[z-1, b] *
                            pr_ij_mm *
                            il_en * up[i, z] * up[j-2, j]
                    end
                end
            end
            # other internal loops
            for k in i+2:j-3 # , l in k+1:j-2
                l = match[k]
                actually_paired = l > k
                # Replace with where
                # Avoids out of bounds indexing
                if ! actually_paired
                    l = j-1
                end
                res = zero(Te)::Te
                lup, rup = k-i-1, j-l-1

                # TODO: Special cases. Can be replaced with wheres
                if lup <= 1 || rup <= 1
                    # 1xn already done
                    continue
                end
                if (lup == 2 && rup == 2) || (lup == 2 && rup == 3) || (lup == 3 && rup == 2)
                    # 2x2, 2x3, 3x2
                    # Could be optimised using the mismatch trick.
                    # Probably not worth it.
                    for bip1 in 1:NTS, bjm1 in 1:NTS, bkm1 in 1:NTS, blp1 in 1:NTS
                        res += P[bk, bl, k, l] *
                            p_seq[k, bk] *
                            p_seq[l, bl] *
                            en_internal(em, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup) *
                            p_seq[k-1, bkm1] *
                            p_seq[l+1, blp1] *
                            p_seq[i+1, bip1] *
                            p_seq[j-1, bjm1] * up[i, k] * up[l, j]
                    end
                else
                    # Optimised using the mismatch trick (mmij * OMM)
                    init_and_pair = en_internal_init(em, lup+rup) *
                        en_internal_asym(em, abs(lup-rup)) *
                        P[bk, bl, k, l] *
                        p_seq[k, bk] *
                        p_seq[l, bl]
                    res += OMM[bk, bl, k, l] * mmij * init_and_pair * up[i, k] * up[l, j]
                end
                # Replace with where
                if actually_paired
                    sm += res
                end
            end
        end
        return sm
    end

    function psum_bulges(Te::Type, bi::Integer, bj::Integer, i::Integer, j::Integer)
        sm = zero(Te)
        for bpkl in 1:NBPS
            bk, bl = get_bp_bases(bpkl)
            for kl in i+2:j-2  # TODO was: range(i+2, j-1):
                sm += P[bk, bl, i+1, kl] *
                    p_seq[i+1, bk] *
                    p_seq[kl, bl] *
                    en_bulge(em, bi, bj, bk, bl, j-kl-1) * up[kl, j]
                sm += P[bk, bl, kl, j-1] *
                    p_seq[kl, bk] *
                    p_seq[j-1, bl] *
                    en_bulge(em, bi, bj, bk, bl, kl-i-1) * up[i, kl]
            end
        end
        return sm
    end

    function fill_paired(i::Integer)
        j = match[i]
        # Replace with cond (not where)
        # I think you can actually do this with a cond because it's in the outer scan
        if j == i
            return
        end
        for bpij in 1:NBPS
            bi, bj = get_bp_bases(bpij)
            # TODO: shouldn't all other terms added to sm also be multiplied with up[i,j] ???
            sm = psum_hairpin(p_seq, em, bi, bj, i, j) * up[i, j]
            sm += psum_bulges(Te, bi, bj, i, j)
            sm += psum_internal_loops(Te, bi, bj, i, j)
            # Stacks
            for bpkl in 1:NBPS
                bk, bl = get_bp_bases(bpkl)
                sm += P[bk, bl, i+1, j-1] *
                    p_seq[i+1, bk] *
                    p_seq[j-1, bl] *
                    en_stack(em, bi, bj, bk, bl)
            end
            # Multi-loops
            for bip1 in 1:NTS, bjm1 in 1:NTS
                sm += ML[bi, bip1, bjm1, bj, 2, i+1, j-1] *
                    p_seq[i+1, bip1] *
                    p_seq[j-1, bjm1] *
                    en_multi_closing(em, bi, bip1, bjm1, bj)
            end
            P[bi, bj, i, j] = sm
        end
        return
    end

    for i in n:-1:1  # TODO was: range(n, 0, -1):
        fill_outer_mismatch(i)
        fill_paired(i)
        fill_multi(i)
        fill_external(i)
    end

    sm = zero(Te)
    for bim1 in 1:NTS, bi in 1:NTS
        sm += E[bim1, bi, 1] * p_seq[1, bi] * p_seq[0, bim1]
    end
    return sm
end

function seqstruct_partition(p_seq::AbstractMatrix{Tp}, em::M; hpmin::Int=HAIRPIN) where {Tp,Te,M<:AbstractModel{Te}}
    T = promote_type(Tp, Te)
    HAIRPIN = hpmin
    n = size(p_seq, 1)
    n > 0 || error("size(p_seq,1) == 0")
    axes(p_seq) == (1:n, 1:NTS) || error("expected axes(p_seq) == (1:n, 1:NTS)")

    # pad p_seq so indexing one out of bounds will be ok (left and right)
    #
    # NOTE: the python code uses 1-based indexing in a 0-based array
    # for this (presumably to avoid index -1, which in python means
    # "last array index"). We recreate this indexing with a
    # OffsetMatrix, but ONLY FOR THE FIRST INDEX, i.e. index `i` in
    # p_seq[i, j].
    p_padded = OffsetMatrix(zeros(Tp, n+2, NTS), 0:n+1, 1:NTS) # TODO: indices 0:n+1, 1:NTS
    p_padded[1:n, 1:NTS] .= p_seq
    p_padded[begin, begin] = one(Tp)
    p_padded[end, begin] = one(Tp)
    p_seq = p_padded

    # NOTE: NTS indices are 1:NTS, sequence indices are 0:n+1
    # external loop
    E = OffsetArray(zeros(T, NTS, NTS, n+2), 1:NTS, 1:NTS, 0:n+1)
    # P[i,j]: (i,j) are basepaired, i.e. close a loop
    P = OffsetArray(zeros(T, NTS, NTS, n+2, n+2), 1:NTS, 1:NTS, 0:n+1, 0:n+1)
    # ML: multiloop tables (3)
    # three tables needed: 0 stems, 1 stem, >= 2 stems (TODO: check)
    ML = OffsetArray(zeros(T, NTS, NTS, NTS, NTS, 3, n+2, n+2), 1:NTS, 1:NTS, 1:NTS, 1:NTS, 0:2, 0:n+1, 0:n+1)
    OMM = OffsetArray(zeros(T, NTS, NTS, n+2, n+2), 1:NTS, 1:NTS, 0:n+1, 0:n+1)

    for i in n:-1:1
        fill_outer_mismatch!(OMM, p_seq, em, i, n)
        fill_paired!(P, ML, OMM, p_seq, em, i, n, HAIRPIN)
        fill_multi!(ML, P, p_seq, em, i, n)
        fill_external!(E, P, p_seq, em, n, i)
    end

    sm = zero(T)::T
    for bim1 in 1:NTS, bi in 1:NTS
        sm += E[bim1, bi, 1] * p_seq[1, bi] * p_seq[0, bim1]
    end
    return sm
end
