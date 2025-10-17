using FoldRNA: FoldRNA, RNA_TURNER2004

abstract type AbstractModel{T} end
Base.eltype(::AbstractModel{T}) where T = T
specialhairpins(::AbstractModel) = SPECIAL_HAIRPINS

# energy contributions for the nearest-neighbour model
# em::AbstractModel
#
# en_bulge(em, bi, bj, bk, bl, j-kl-1)
# en_ext_branch(em, ee, bi, bj, ff)
# en_hairpin_not_special(bi, bj, bip1, bjm1, len)
# en_hairpin_special(em, id)
# en_il_inner_mismatch(em, bi, bj, bip1, bjm1)
# en_il_outer_mismatch(em, bk, bl, bkm1, blp1)
# en_internal_asym(em, abs(lup-rup)
# en_internal(em, bi, bj, bk, bl, bip1, bjm1, bip1, bjm1, 1, 1)
# en_internal_init(em, lup+rup)
# en_multi_branch(em, bim1, bi, bj, bjp1)
# en_multi_closing(em, bi, bip1, bjm1, bj)
# en_stack(em, bi, bj, bk, bl)

# comments on abstract base class methods
#
# en_hairpin_not_special(self, bi, bj, bip1, bjm1, nunpaired):
#     Note that vienna ignores the all-C case
# en_hairpin_special(self, id):
#     id is the index into specialhairpins(::AbstractModel)

function en_internal(em::AbstractModel, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup)
    en = en_internal_init(em, lup + rup) *
        en_internal_asym(em, abs(lup - rup)) *
        en_il_inner_mismatch(em, bi, bj, bip1, bjm1) *
        en_il_outer_mismatch(em, bk, bl, bkm1, blp1)
    return en
end

# ViennaModel: use ViennaRNA's nearest-neighbour parameters
# TODO: assumes all base type indices (bi, bj, bim1, etc.) have
#       1-based indexing
struct ViennaModel{T} <: AbstractModel{T}
    boltz :: FoldRNA.LoopModel{Float64, Int64, 5, 7, 30}
    special_hairpins :: Vector{String}
    special_hairpins_boltz :: Vector{Float64}
    function ViennaModel(dG::FoldRNA.LoopModel{Int64, Int64, 5, 7, 30})
        # TODO: hardcoded LoopModel types, ViennaModel{T}
        # TODO: filter wildcard chars?
        boltz = FoldRNA.LoopModel{Float64, Int64, 5, 7, 30}(; alphabet=dG.alphabet)
        m = new{Float64}(boltz, String[], Float64[])
        m.boltz.maxloop  = dG.maxloop
        m.boltz.bptype  .= dG.bptype
        make_exp_params!(m, dG)
        return m
    end
end
ViennaModel() = ViennaModel(RNA_TURNER2004)
@inline specialhairpins(m::ViennaModel) = m.special_hairpins
function make_exp_params!(m::ViennaModel, dG::FoldRNA.LoopModel{Int64, Int64, 5, 7, 30})
    # TODO: make temperature a function argument, but need dH and dS (not just dG)
    beta = dG.unit / dG.RT  # 1/RT in units that the ViennaRNA tables use
    b = m.boltz

    # arrays
    b.mismatch_extloop   .= exp.(-beta .* dG.mismatch_extloop)    # en_ext_branch
    b.mismatch_multiloop .= exp.(-beta .* dG.mismatch_multiloop)  # en_multi_branch, en_multi_closing
    b.mismatch_hairpin   .= exp.(-beta .* dG.mismatch_hairpin)    # en_hairpin_not_special
    b.mismatch_intloop   .= exp.(-beta .* dG.mismatch_intloop)    # en_il_inner_mismatch, en_il_outer_mismatch
    b.mismatch_intloop1n .= exp.(-beta .* dG.mismatch_intloop1n)  # en_internal
    b.mismatch_intloop23 .= exp.(-beta .* dG.mismatch_intloop23)  # en_internal
    b.stack              .= exp.(-beta .* dG.stack)               # en_stack, en_bulge
    b.hairpin_init       .= exp.(-beta .* dG.hairpin_init)        # en_hairpin_not_special
    b.bulge_init         .= exp.(-beta .* dG.bulge_init)          # en_bulge
    b.terminal_nonGC_bp  .= exp.(-beta .* dG.terminal_nonGC_bp)   # en_ext_branch, en_multi_branch, en_multi_closing, en_bulge
    b.dangle5            .= exp.(-beta .* dG.dangle5)             # en_ext_branch, en_multi_branch, en_multi_closing
    b.dangle3            .= exp.(-beta .* dG.dangle3)             # en_ext_branch, en_multi_branch, en_multi_closing
    b.intloop_init       .= exp.(-beta .* dG.intloop_init)        # en_internal_init, en_internal
    b.intloop11          .= exp.(-beta .* dG.intloop11)           # en_internal
    b.intloop12          .= exp.(-beta .* dG.intloop12)           # en_internal
    # intloop22 is the largest array and is size-reduced (wildcard
    # chars filtered)
    b.intloop22           = similar(dG.intloop22, eltype(m))
    b.intloop22          .= exp.(-beta .* dG.intloop22)           # en_internal

    # scalars
    b.ninio_max        = exp(-beta * dG.ninio_max)                # en_interal_asym, en_internal
    b.ninio_m          = exp(-beta * dG.ninio_m)                  # en_internal_asym, en_internal
    b.multiloop_branch = exp(-beta * dG.multiloop_branch)         # en_multi_branch, en_multi_closing
    b.multiloop_init   = exp(-beta * dG.multiloop_init)

    # no exponentiation, just multiply with -beta
    # because: hairpin init energies for len > maxloop are extrapolated as:
    #     en_max_init + lxc * log(len/maxloop)
    # for the Boltzmann factors x -> exp(-beta*x) we have
    # exp(-beta*en_max_init) * (len/maxloop)^(-beta*lxc)
    m.boltz.lxc = -beta * dG.lxc                                  # en_bulge, en_internal_init, en_hairpin_not_special

    # special hairpins
    for (i, (sh_seq, sh_en)) in enumerate(dG.specialhairpins)
        sh_str = join(map(k -> RNA_ALPHA[k], sh_seq))
        sh_boltz = exp(-beta * sh_en)
        push!(m.special_hairpins, sh_str)
        push!(m.special_hairpins_boltz, sh_boltz)
    end
    return m
end
@inline function _en_helper_branch(m::ViennaModel{T}, mismatch, bim1, bi, bj, bjp1) where T
    # TODO: make branchless
    b = m.boltz
    s = one(T)
    bptype_ij = b.bptype[bi, bj]
    if bptype_ij == 0
        return zero(T)
    end
    s *= b.terminal_nonGC_bp[bptype_ij]
    if (bim1 != INVALID_BASE) && (bjp1 != INVALID_BASE)
        #s += score_mismatch(fold, i, j, dangle5, dangle3, mismatch_extloop)
        s *= mismatch[bptype_ij, bim1, bjp1]
    else
        if bim1 != INVALID_BASE
            # s += fold.model.dangle5[bptype(fold, i, j), fold.seq[dangle5]]
            s *= b.dangle5[bptype_ij, bim1]
        end
        if bjp1 != INVALID_BASE
            # s += fold.model.dangle3[bptype(fold, i, j), fold.seq[dangle3]]
            s *= b.dangle3[bptype_ij, bjp1]
        end
    end
    return s
end
@inline function en_ext_branch(m::ViennaModel{T}, bim1, bi, bj, bjp1) where T
    return _en_helper_branch(m, m.boltz.mismatch_extloop, bim1, bi, bj, bjp1)
end
@inline function en_multi_branch(m::ViennaModel{T}, bim1, bi, bk, bkp1) where T
    return m.boltz.multiloop_branch *
        _en_helper_branch(m, m.boltz.mismatch_multiloop, bim1, bi, bk, bkp1)
end
@inline function en_multi_closing(m::ViennaModel{T}, bi, bip1, bjm1, bj) where T
    # Note: reversed order of (j,i), (bjm1,bip1) compared to en_multi_branch
    return m.boltz.multiloop_init * m.boltz.multiloop_branch *
        _en_helper_branch(m, m.boltz.mismatch_multiloop, bjm1, bj, bi, bip1)
end
@inline function en_hairpin_not_special(m::ViennaModel{T}, bi, bj, bip1, bjm1, nunpaired) where T
    b = m.boltz
    maxloop = b.maxloop
    bptype_ij = b.bptype[bi, bj]
    if bptype_ij == 0
        return zero(T)
    end
    s = one(T)
    # hairpin init
    if nunpaired <= b.maxloop
        # small loops have precomputed values
        s *= b.hairpin_init[nunpaired]
    else
        # large loops get extrapolated
        # TODO old energy: m.hairpin_init[m.maxloop] + trunc(Int, m.lxc * log(len / m.maxloop))
        s *= b.hairpin_init[maxloop] * (nunpaired / maxloop)^b.lxc
    end
    # hairpin mismatch
    if nunpaired == 3
        s *= b.terminal_nonGC_bp[bptype_ij]
    else
        s *= b.mismatch_hairpin[bptype_ij, bip1, bjm1]
    end
    return s
end
@inline en_hairpin_special(m::ViennaModel{T}, id) where T = m.special_hairpins_boltz[id]
@inline function en_stack(m::ViennaModel{T}, bi, bj, bk, bl) where T
    # TODO: make branchless
    bptype_ij = m.boltz.bptype[bi, bj]
    bptype_lk = m.boltz.bptype[bl, bk]  # note reversed order of (bl, bk)
    if (bptype_ij == 0) || (bptype_lk == 0)
        return zero(T)  # impossible base pairs have infinite energy / Boltzmann factors of 0
    end
    return m.boltz.stack[bptype_ij, bptype_lk]
end
@inline function en_bulge(m::ViennaModel{T}, bi, bj, bk, bl, nunpaired) where T
    # TODO: make branchless
    b = m.boltz
    maxloop = m.boltz.maxloop
    s = one(T)
    # initiation penalty
    if nunpaired < maxloop
        s *= b.bulge_init[nunpaired]
    else
        # large loops get extrapolated
        # TODO old: m.bulge_init[maxloop] + trunc(Int, m.lxc * log(nunpaired / m.maxloop)))
        s *= b.bulge_init[maxloop] * (nunpaired / maxloop)^b.lxc
    end
    if nunpaired == 1
        # short 0×1 or 1x0 bulge loop, so it gets stacking
        s *= en_stack(m, bi, bj, bk, bl)
    else
        bptype_ij = b.bptype[bi, bj]
        bptype_kl = b.bptype[bk, bl]
        if (bptype_ij == 0) || (bptype_kl == 0)
            return zero(T)  # impossible base pairs have infinite energy / Boltzmann factors of 0
        end
        s *= b.terminal_nonGC_bp[bptype_ij]
        s *= b.terminal_nonGC_bp[bptype_kl]  # TODO: why not l,k ?
    end
    return s
end
@inline function en_il_inner_mismatch(m::ViennaModel{T}, bi, bj, bip1, bjm1) where T
    b = m.boltz
    bptype_ij = b.bptype[bi, bj]
    if bptype_ij == 0
        return zero(T)
    end
    return b.mismatch_intloop[bptype_ij, bip1, bjm1]
end
@inline function en_il_outer_mismatch(m::ViennaModel{T}, bi, bj, bim1, bjp1) where T
    b = m.boltz
    bptype_ji = b.bptype[bj, bi]  # Note reversed order of bj, bi
    if bptype_ji == 0
        return zero(T)
    end
    return b.mismatch_intloop[bptype_ji, bjp1, bim1]  # Note reversed order of bjp1, bim1
end
@inline function en_internal_init(m::ViennaModel{T}, sz::Int) where T
    # TODO: make branchless
    b = m.boltz
    maxloop = b.maxloop
    return if sz <= maxloop
        b.intloop_init[sz]
    else
        # large loops get extrapolated
        # TODO energy form: m.intloop_init[m.maxloop] + trunc(Int, m.lxc * log(sz / m.maxloop))
        b.intloop_init[maxloop] * (sz / b.maxloop)^b.lxc
    end
end
@inline function en_internal_asym(m::ViennaModel{T}, asym::Int) where T
    # TODO: make branchless
    b = m.boltz
    # Note: min changed to max because exp(-x) reverses ordering
    # exp(-min(a, b)) == max(exp(-a), exp(-b))
    # TODO energy form: min(m.ninio_max, abs(asym) * m.ninio_m)
    return max(b.ninio_max, b.ninio_m^asym)
end
@inline function en_internal(m::ViennaModel{T}, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup) where T
    b = m.boltz
    bptype_ij = b.bptype[bi, bj]
    bptype_lk = b.bptype[bl, bk]  # Note reversed order of bl, bk
    if (bptype_ij == 0) || (bptype_lk == 0)
        return zero(T)  # impossible base pairs have infinite energy / Boltzmann factors of 0
    end
    maxloop = b.maxloop
    s = one(T)
    if lup == 0 && rup == 0
        # stacked base pairs
        s *= en_stack(m, bi, bj, bk, bl)
    elseif lup == 0 || rup == 0
        # bulge loop
        nunpaired = lup + rup
        s *= en_bulge(m, bi, bj, bk, bl, nunpaired)
    elseif lup == 1 && rup == 1
        # 1×1 interior loop
        s *= b.intloop11[bptype_ij, bptype_lk, bip1, bjm1]
    elseif lup == 1 && rup == 2
        # 1×2 interior loop
        s *= b.intloop12[bptype_ij, bptype_lk, bip1, blp1, bjm1]
    elseif lup == 2 && rup == 1
        # 2x1 interior loop
        # Note intloop12, and reversed order of bptype_lk, bptype_ij and blp1, bip1
        s *= b.intloop12[bptype_lk, bptype_ij, blp1, bip1, bkm1]
    elseif lup == 1 || rup == 1
        # 1xn, nx1 interior loop
        s *= en_internal_init(m, lup + rup)
        s *= en_internal_asym(m, abs(lup-rup))
        s *= b.mismatch_intloop1n[bptype_ij, bip1, bjm1]
        s *= b.mismatch_intloop1n[bptype_lk, blp1, bkm1]
    elseif lup == 2 && rup == 2
        # 2x2 interior loop
        s *= b.intloop22[bptype_ij, bptype_lk, bip1, bkm1, blp1, bjm1]
    elseif (lup, rup) == (2,3) || (lup, rup) == (3,2)
        # 2x3, 3x2 interior loop
        s *= en_internal_init(m, 5)
        s *= en_internal_asym(m, 1)  # abs(lup-rup) == 1
        s *= b.mismatch_intloop23[bptype_ij, bip1, bjm1]
        s *= b.mismatch_intloop23[bptype_lk, blp1, bkm1]
    else
        # generic interior loop
        s *= en_internal_init(m, lup + rup)
        s *= en_internal_asym(m, abs(lup-rup))
        s *= b.mismatch_intloop[bptype_ij, bip1, bjm1]
        s *= b.mismatch_intloop[bptype_lk, blp1, bkm1]
    end
    return s
end

# All1Model: all Boltzmann factors are 1, so the partition function is
# the number of a structures
struct All1Model{T} <: AbstractModel{T} end
All1Model() = All1Model{Float64}()
en_ext_branch(em::All1Model{T}, bim1, bi, bj, bjp1) where T = one(T)
en_multi_branch(em::All1Model{T}, bim1, bi, bk, bkp1) where T = one(T)
en_multi_closing(em::All1Model{T}, bi, bip1, bjm1, bj) where T = one(T)
en_hairpin_not_special(em::All1Model{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
en_hairpin_special(em::All1Model{T}, id) where T = one(T)
en_stack(em::All1Model{T}, bi, bj, bk, bl) where T = one(T)
en_bulge(em::All1Model{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
en_il_inner_mismatch(em::All1Model{T}, bi, bj, bip1, bjm1) where T = one(T)
en_il_outer_mismatch(em::All1Model{T}, bi, bj, bim1, bjp1) where T = one(T)
en_internal_init(em::All1Model{T}, sz) where T = one(T)
en_internal_asym(em::All1Model{T}, asym) where T = one(T)

struct DebugModel{T,M} <: AbstractModel{T}
    m :: M
end
DebugModel(m::AbstractModel{T}) where T = DebugModel{T,typeof(m)}(m)
# TODO: generate these functions automatically
# TODO: use @debug instead of @info ? how to show @debug messages when calling
function en_ext_branch(em::DebugModel, bim1, bi, bj, bjp1)
    val = en_ext_branch(em.m, bim1, bi, bj, bjp1)
    @info "en_ext_branch(em, bim1=$bim1, bi=$bi, bj=$bj, bjp1=$bjp1) = $val"
    return val
end
function en_multi_branch(em::DebugModel, bim1, bi, bk, bkp1)
    val = en_multi_branch(em.m, bim1, bi, bk, bkp1)
    @info "en_multi_branch(em, bim1=$bim1, bi=$bi, bk=$bk, bkp1=$bkp1) = $val"
    return val
end
function en_multi_closing(em::DebugModel, bi, bip1, bjm1, bj)
    @info "en_multi_closing(em, bi=$bi, bip1=$bip1, bjm1=$bjm1, bj=$bj) = $val"
    return val
end
function en_hairpin_not_special(em::DebugModel, bi, bj, bip1, bjm1, nunpaired)
    val = en_hairpin_not_special(em.m, bi, bj, bip1, bjm1, nunpaired)
    @info "en_hairpin_not_special(em, bi=$bi, bj=$bj, bip1=$bip1, bjm1=$bjm1, nunpaired=$nunpaired) = $val"
    return val
end
function en_hairpin_special(em::DebugModel, id)
    val = en_hairpin_special(em.m, id)
    @info "en_hairpin_special(em, id=$id) = $val"
    return val
end
function en_stack(em::DebugModel, bi, bj, bk, bl)
    val = en_stack(em.m, bi, bj, bk, bl)
    @info "en_stack(em, bi=$bi, bj=$bj, bk=$bk, bl=$bl) = $val"
    return val
end
function en_bulge(em::DebugModel, bi, bj, bk, bl, nunpaired)
    val = en_bulge(em.m, bi, bj, bk, bl, nunpaired)
    @info "en_bulge(em, bi=$bi, bj=$bj, bk=$bk, bl=$bl, nunpaired=$nunpaired) = $val"
    return val
end
function en_il_inner_mismatch(em::DebugModel, bi, bj, bip1, bjm1)
    val = en_il_inner_mismatch(em.m, bi, bj, bip1, bjm1)
    @info "en_il_inner_mismatch(em, bi=$bi, bj=$bj, bip1=$bip1, bjm1=$bjm1) = $val"
    return val
end
function en_il_outer_mismatch(em::DebugModel, bi, bj, bim1, bjp1)
    val = en_il_outer_mismatch(em.m, bi, bj, bim1, bjp1)
    @info "en_il_outer_mismatch(em, bi=$bi, bj=$bj, bim1=$bim1, bjp1=$bjp1) = $val"
    return val
end
function en_internal_init(em::DebugModel, sz)
    val = en_internal_init(em.m, sz)
    @info "en_internal_init(em, sz=$sz) = $val"
    return val
end
function en_internal_asym(em::DebugModel, asym)
    val = en_internal_asym(em.m, asym)
    @info "en_internal_asym(em, asym=$asym) = $val"
    return val
end
function en_internal(em::DebugModel, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup)
    val = en_internal(em.m, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup)
    @info "en_internal(em, bi=$bi, bj=$bj, bk=$bk, bl=$bl, bip1=$bip1, bjm1=$bjm1, bkm1=$bkm1, blp1=$blp1, lup=$lup, rup=$rup) = $val"
    return val
end

abstract type AbstractRandomModel{T} <: AbstractModel{T} end

struct RandomModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomModel{T}() where T = RandomModel{T}(1)
RandomModel() = RandomModel{Float64}(1)
RandomModel(seed::Integer) = RandomModel{Float64}(seed)
# Note: the numbers 1, 2, 3, ... act as a salt for the hash, to avoid
# different energy terms for the same indices hashing to the same
# value
# Note: the bim-1, etc is necessary to get the same hash value as
# python, where the indices bim, etc. come from 0-based arrays
en_ext_branch(em::AbstractRandomModel{T}, bim1, bi, bj, bjp1) where T =
    T(float_hash(em.seed, bim1-1, bi-1, bj-1, bjp1-1, 1))
en_multi_branch(em::AbstractRandomModel{T}, bim1, bi, bk, bkp1) where T =
    T(float_hash(em.seed, bim1-1, bi-1, bk-1, bkp1-1, 2))
en_multi_closing(em::AbstractRandomModel{T}, bi, bip1, bjm1, bj) where T =
    T(float_hash(em.seed, bi-1, bip1-1, bjm1-1, bj-1, 3))
en_hairpin_not_special(em::AbstractRandomModel{T}, bi, bj, bip1, bjm1, nunpaired) where T =
    T(float_hash(em.seed, bi-1, bj-1, bip1-1, bjm1-1, nunpaired, 4))
# TODO: id-firstindex(specialhairpins(em)) so we get the same values as python,
# where id is an index into an array starting from 0
en_hairpin_special(em::AbstractRandomModel{T}, id) where T =
    T(float_hash(em.seed, id-firstindex(specialhairpins(em)), 5))
en_stack(em::AbstractRandomModel{T}, bi, bj, bk, bl) where T =
    T(float_hash(em.seed, bi-1, bj-1, bk-1, bl-1, 6))
en_bulge(em::AbstractRandomModel{T}, bi, bj, bk, bl, nunpaired) where T =
    T(float_hash(em.seed, bi-1, bj-1, bk-1, bl-1, nunpaired, 7))
en_il_inner_mismatch(em::AbstractRandomModel{T}, bi, bj, bip1, bjm1) where T =
    T(float_hash(em.seed, bi-1, bj-1, bip1-1, bjm1-1, 8))
en_il_outer_mismatch(em::AbstractRandomModel{T}, bi, bj, bim1, bjp1) where T =
    T(float_hash(em.seed, bi-1, bj-1, bim1-1, bjp1-1, 9))
en_internal_init(em::AbstractRandomModel{T}, sz) where T =
    T(float_hash(em.seed, sz, 10))
en_internal_asym(em::AbstractRandomModel{T}, asym) where T =
    T(float_hash(em.seed, asym, 11))

@inline function float_hash(seed::Int, args::Int...)
    base = 97
    m = 10007 # modulus
    p = 1
    v = mod((seed * 7777777), m)
    for a in args
        #if isnothing(a) == 0
        #    continue
        #end
        v += p * a
        v = mod(v, m)
        p *= base
        p = mod(p, m)
    end
    v = mod((v * 9999999), m)
    return v / (m-1)
end

struct RandomExtloopModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomExtloopModel{T}() where T = RandomExtloopModel{T}(1)
RandomExtloopModel() = RandomExtloopModel{Float64}(1)
RandomExtloopModel(seed::Integer) = RandomExtloopModel{Float64}(seed)
# en_ext_branch(em::RandomExtloopModel{T}, bim1, bi, bj, bjp1) where T = one(T)
en_multi_branch(em::RandomExtloopModel{T}, bim1, bi, bk, bkp1) where T = one(T)
en_multi_closing(em::RandomExtloopModel{T}, bi, bip1, bjm1, bj) where T = one(T)
en_hairpin_not_special(em::RandomExtloopModel{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
en_hairpin_special(em::RandomExtloopModel{T}, id) where T = one(T)
en_stack(em::RandomExtloopModel{T}, bi, bj, bk, bl) where T = one(T)
en_bulge(em::RandomExtloopModel{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
en_il_inner_mismatch(em::RandomExtloopModel{T}, bi, bj, bip1, bjm1) where T = one(T)
en_il_outer_mismatch(em::RandomExtloopModel{T}, bi, bj, bim1, bjp1) where T = one(T)
en_internal_init(em::RandomExtloopModel{T}, sz) where T = one(T)
en_internal_asym(em::RandomExtloopModel{T}, asym) where T = one(T)

struct RandomMultiloopModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomMultiloopModel{T}() where T = RandomMultiloopModel{T}(1)
RandomMultiloopModel() = RandomMultiloopModel{Float64}(1)
RandomMultiloopModel(seed::Integer) = RandomMultiloopModel{Float64}(seed)
en_ext_branch(em::RandomMultiloopModel{T}, bim1, bi, bj, bjp1) where T = one(T)
# en_multi_branch(em::RandomMultiloopModel{T}, bim1, bi, bk, bkp1) where T = one(T)
# en_multi_closing(em::RandomMultiloopModel{T}, bi, bip1, bjm1, bj) where T = one(T)
en_hairpin_not_special(em::RandomMultiloopModel{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
en_hairpin_special(em::RandomMultiloopModel{T}, id) where T = one(T)
en_stack(em::RandomMultiloopModel{T}, bi, bj, bk, bl) where T = one(T)
en_bulge(em::RandomMultiloopModel{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
en_il_inner_mismatch(em::RandomMultiloopModel{T}, bi, bj, bip1, bjm1) where T = one(T)
en_il_outer_mismatch(em::RandomMultiloopModel{T}, bi, bj, bim1, bjp1) where T = one(T)
en_internal_init(em::RandomMultiloopModel{T}, sz) where T = one(T)
en_internal_asym(em::RandomMultiloopModel{T}, asym) where T = one(T)

struct RandomHairpinModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomHairpinModel{T}() where T = RandomHairpinModel{T}(1)
RandomHairpinModel() = RandomHairpinModel{Float64}(1)
RandomHairpinModel(seed::Integer) = RandomHairpinModel{Float64}(seed)
en_ext_branch(em::RandomHairpinModel{T}, bim1, bi, bj, bjp1) where T = one(T)
en_multi_branch(em::RandomHairpinModel{T}, bim1, bi, bk, bkp1) where T = one(T)
en_multi_closing(em::RandomHairpinModel{T}, bi, bip1, bjm1, bj) where T = one(T)
# en_hairpin_not_special(em::RandomHairpinModel{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
# en_hairpin_special(em::RandomHairpinModel{T}, id) where T = one(T)
en_stack(em::RandomHairpinModel{T}, bi, bj, bk, bl) where T = one(T)
en_bulge(em::RandomHairpinModel{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
en_il_inner_mismatch(em::RandomHairpinModel{T}, bi, bj, bip1, bjm1) where T = one(T)
en_il_outer_mismatch(em::RandomHairpinModel{T}, bi, bj, bim1, bjp1) where T = one(T)
en_internal_init(em::RandomHairpinModel{T}, sz) where T = one(T)
en_internal_asym(em::RandomHairpinModel{T}, asym) where T = one(T)

struct RandomStackModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomStackModel{T}() where T = RandomStackModel{T}(1)
RandomStackModel() = RandomStackModel{Float64}(1)
RandomStackModel(seed::Integer) = RandomStackModel{Float64}(seed)
en_ext_branch(em::RandomStackModel{T}, bim1, bi, bj, bjp1) where T = one(T)
en_multi_branch(em::RandomStackModel{T}, bim1, bi, bk, bkp1) where T = one(T)
en_multi_closing(em::RandomStackModel{T}, bi, bip1, bjm1, bj) where T = one(T)
en_hairpin_not_special(em::RandomStackModel{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
en_hairpin_special(em::RandomStackModel{T}, id) where T = one(T)
# en_stack(em::RandomStackModel{T}, bi, bj, bk, bl) where T = one(T)
en_bulge(em::RandomStackModel{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
en_il_inner_mismatch(em::RandomStackModel{T}, bi, bj, bip1, bjm1) where T = one(T)
en_il_outer_mismatch(em::RandomStackModel{T}, bi, bj, bim1, bjp1) where T = one(T)
en_internal_init(em::RandomStackModel{T}, sz) where T = one(T)
en_internal_asym(em::RandomStackModel{T}, asym) where T = one(T)

struct RandomBulgeModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomBulgeModel{T}() where T = RandomBulgeModel{T}(1)
RandomBulgeModel() = RandomBulgeModel{Float64}(1)
RandomBulgeModel(seed::Integer) = RandomBulgeModel{Float64}(seed)
en_ext_branch(em::RandomBulgeModel{T}, bim1, bi, bj, bjp1) where T = one(T)
en_multi_branch(em::RandomBulgeModel{T}, bim1, bi, bk, bkp1) where T = one(T)
en_multi_closing(em::RandomBulgeModel{T}, bi, bip1, bjm1, bj) where T = one(T)
en_hairpin_not_special(em::RandomBulgeModel{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
en_hairpin_special(em::RandomBulgeModel{T}, id) where T = one(T)
en_stack(em::RandomBulgeModel{T}, bi, bj, bk, bl) where T = one(T)
# en_bulge(em::RandomBulgeModel{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
en_il_inner_mismatch(em::RandomBulgeModel{T}, bi, bj, bip1, bjm1) where T = one(T)
en_il_outer_mismatch(em::RandomBulgeModel{T}, bi, bj, bim1, bjp1) where T = one(T)
en_internal_init(em::RandomBulgeModel{T}, sz) where T = one(T)
en_internal_asym(em::RandomBulgeModel{T}, asym) where T = one(T)

struct RandomILModel{T} <: AbstractRandomModel{T}
    seed :: Int
end
RandomILModel{T}() where T = RandomILModel{T}(1)
RandomILModel() = RandomILModel{Float64}(1)
RandomILModel(seed::Integer) = RandomILModel{Float64}(seed)
en_ext_branch(em::RandomILModel{T}, bim1, bi, bj, bjp1) where T = one(T)
en_multi_branch(em::RandomILModel{T}, bim1, bi, bk, bkp1) where T = one(T)
en_multi_closing(em::RandomILModel{T}, bi, bip1, bjm1, bj) where T = one(T)
en_hairpin_not_special(em::RandomILModel{T}, bi, bj, bip1, bjm1, nunpaired) where T = one(T)
en_hairpin_special(em::RandomILModel{T}, id) where T = one(T)
en_stack(em::RandomILModel{T}, bi, bj, bk, bl) where T = one(T)
en_bulge(em::RandomILModel{T}, bi, bj, bk, bl, nunpaired) where T = one(T)
# en_il_inner_mismatch(em::RandomILModel{T}, bi, bj, bip1, bjm1) where T = one(T)
# en_il_outer_mismatch(em::RandomILModel{T}, bi, bj, bim1, bjp1) where T = one(T)
# en_internal_init(em::RandomILModel{T}, sz) where T = one(T)
# en_internal_asym(em::RandomILModel{T}, asym) where T = one(T)

function boltz(seq::AbstractString, dbn::AbstractString, em::AbstractModel{T}) where {T}
    # Note: formerly called `calculate`
    n = length(dbn)
    length(seq) == n || throw(ArgumentError("seq and dbn must have same length ($(length(seq)) != $n)"))
    b_idx(i) = (i <= 0 || i >= n+1) ? INVALID_BASE : findfirst(seq[i], RNA_ALPHA)::Int
    ch, right = structure_tree(dbn)
    SPECIAL_HAIRPINS = specialhairpins(em)

    function calc_rec(atl)
        if atl == -1
            sm = one(T)
            for cl in ch[atl]
                sm *= calc_rec(cl) *
                    en_ext_branch(em, b_idx(cl-1), b_idx(cl), b_idx(right[cl]), b_idx(right[cl]+1))
            end
            return sm
        end
        if ! (atl in keys(ch))
            s = seq[atl:right[atl]]  # TODO was: seq[atl:right[atl]+1]
            return if s in SPECIAL_HAIRPINS
                # TODO: inefficient lookup (findfirst)
                idx = findfirst(x -> x == s, SPECIAL_HAIRPINS)::Int
                en_hairpin_special(em, idx)
            else
                en_hairpin_not_special(em, b_idx(atl), b_idx(right[atl]), b_idx(atl+1),
                                       b_idx(right[atl]-1), right[atl]-atl-1)
            end
        elseif length(ch[atl]) == 1
            cl, cr = ch[atl][begin], right[ch[atl][begin]]
            if cl == atl+1 && cr == right[atl]-1
                return en_stack(em, b_idx(atl), b_idx(right[atl]), b_idx(cl), b_idx(cr)) * calc_rec(cl)
            elseif cl == atl+1 || cr == right[atl]-1
                nunpaired = max(cl-atl-1, right[atl]-cr-1)
                return en_bulge(em, b_idx(atl), b_idx(right[atl]), b_idx(cl), b_idx(cr), nunpaired) * calc_rec(cl)
            else
                bi = b_idx(atl)
                bj = b_idx(right[atl])
                bip1 = b_idx(atl+1)
                bjm1 = b_idx(right[atl]-1)
                bk = b_idx(cl)
                bl = b_idx(cr)
                bkm1 = b_idx(cl-1)
                blp1 = b_idx(cr+1)
                lup = cl-atl-1
                rup = right[atl]-cr-1
                return en_internal(em, bi, bj, bk, bl, bip1, bjm1, bkm1, blp1, lup, rup) * calc_rec(cl)
            end
        else
            sm = en_multi_closing(em, b_idx(atl), b_idx(atl+1), b_idx(right[atl]-1), b_idx(right[atl]))
            for cl in ch[atl]
                sm *= calc_rec(cl) * en_multi_branch(em, b_idx(cl-1), b_idx(cl), b_idx(right[cl]), b_idx(right[cl]+1))
            end
            return sm
        end
    end
    return calc_rec(-1)
end
