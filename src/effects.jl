# This file contains the functions that how system state effected by reactions.
# There are only four possible effects in this model:
# 1. Cancer cell division
# 2. Cancer cell death
# 3. Effector cell recruitment
# 4. Effector cell death

using Random
using Distributions

# This file contains the functions that define the effects of reactions

"""
    mutations!(rng::AbstractRNG, counter::Base.RefValue, λ::Real, p_a::Real)

Generate mutations for a cell division event with a given poisson rate `λ`,
and a probability `p_a` that a mutation is antigenic.
`counter` is a reference to a counter that tracks the number of mutations,
All mutations will be assigned a unique ID from `counter[]+1` to `counter[]+n`.
The `counter` will be updated to `counter[]+n` after the function call.
The function returns two ranges, the first one for antigenic mutations,
and the second one for passenger mutations, respectively.
"""
function mutation!(rng::AbstractRNG, counter::Base.RefValue, λ::Real, p_a::Real)
    n = rand(rng, Poisson{typeof(λ)}(λ)) # number of new mutations
    id_start = counter[]
    id_end = id_start + n
    counter[] = id_end
    n_antigenic = rand(rng, Binomial{typeof(p_a)}(n, p_a))
    (id_start+1):(id_start+n_antigenic), (id_start+n_antigenic+1):(id_end)
end

growend!(v, n) = Base._growend!(v, n)

function cancer_div!(
    ref::ReactionRef{<:Real,<:CancerCellDivision},
    rs::NamedTuple, # reaction lists
    # A example of effector (used as a constructor when needed)
    effector_example::CLTPopulation,
    # parameters used in cell division
    ## random number generator
    rng::Random.AbstractRNG,
    counter::Base.RefValue{<:Integer},
    ## mutation parameters
    λ::Real,
    p_e::Real,
    p_a::Real,
    dist_i,
    dist_a,
)
    cbl = rs.cancer_cell_division
    cdl = rs.cancer_cell_death
    iil = rs.immune_inhibition

    cancer_parent::CancerCell, effectors_parent = children(ref)

    # Parent cell mutations and reactions
    antigenic_parent = antigenic(cancer_parent)
    passenger_parent = passenger(cancer_parent)
    deriveds_parent = refall(rs, deriveds(cancer_parent))

    CT = typeof(cancer_parent)

    # daughter cell 1, a new cancer cell
    antigenic_d1, passenger_d1 = mutation!(rng, counter, λ, p_a)
    if rand(rng) < p_e
        # mutation causes immune escape,
        # so the daughter cell 1 will not have any antigen and related effector cell
        cb1 = duplicate(deriveds_parent.cancer_cell_division)
        cb1_slot = popslot!(cbl)
        cd1 = duplicate(deriveds_parent.cancer_cell_death)
        cd1_slot = popslot!(cdl)
        # there is none neoantigen, so the inhibition rate is zero
        ii1 = duplicate_zero(deriveds_parent.immune_inhibition)
        ii1_slot = popslot!(iil)

        daughter1 = CT(
            [passenger_parent; passenger_d1; antigenic_parent; antigenic_d1],
            (
                cancer_cell_division = cb1_slot,
                cancer_cell_death = cd1_slot,
                immune_inhibition = ii1_slot,
            ),
        )

        # No effector cell related to this cancer cell
        effectors_d1 = similar(effectors_parent, 0)

        # push new reaction nodes to reaction lists
        unsafe_insert!(cbl, cb1_slot, cb1, (daughter1, effectors_d1))
        unsafe_insert!(cdl, cd1_slot, cd1, (daughter1, effectors_d1))
        unsafe_insert!(iil, ii1_slot, ii1, (daughter1, effectors_d1))
    else
        # if there is no immune escape, the daughter cell will
        # inherit all antigens from the parent cell and may have new antigens
        cb1 = duplicate(deriveds_parent.cancer_cell_division)
        cb1_slot = popslot!(cbl)
        cd1 = duplicate(deriveds_parent.cancer_cell_death)
        cd1_slot = popslot!(cdl)
        # rate should be same as the old one, because the daughter cell
        # will have same antigens as the parent cell
        ii1 = duplicate(deriveds_parent.immune_inhibition)
        ii1_slot = popslot!(iil)

        daughter1 = CT(
            [passenger_parent; passenger_d1],
            [antigenic_parent; antigenic_d1],
            (
                cancer_cell_division = cb1_slot,
                cancer_cell_death = cd1_slot,
                immune_inhibition = ii1_slot,
            ),
        )

        n_oldantigen = length(antigenic_parent)
        n_neoantigen = length(antigenic_d1)

        # effector cells related to this cancer cell
        effectors_d1 = similar(effectors_parent, n_oldantigen + n_neoantigen)

        # push new reaction nodes related to this cancer cell to reaction lists
        # it's ok to use effectors_d1 before it's filled
        unsafe_insert!(cbl, cb1_slot, cb1, (daughter1, effectors_d1))
        unsafe_insert!(cdl, cd1_slot, cd1, (daughter1, effectors_d1))
        unsafe_insert!(iil, ii1_slot, ii1, (daughter1, effectors_d1))

        # copy of effector cells of the parent cell
        # register inhibition reaction of daughter cell to all copied effector cells
        # and update the e_num of the effector cells
        for (i, effector::CLTPopulation) in enumerate(effectors_parent)
            effectors_d1[i] = effector
            push!(deriveds(effector).immune_inhibition, ii1_slot)
            setval_diff!(effector, (c = 1,), rs)
        end

        # create new effector cells and related reactions for new antigens
        # and push related reactions to reaction lists
        if n_neoantigen != 0
            effector_new!(
                effectors_d1,
                rs,
                ii1_slot,
                effector_example,
                n_oldantigen,
                rng,
                dist_i,
                dist_a,
            )
        end
    end

    # daughter cell 2, reuse the old cancer cell
    daughter2 = cancer_parent
    effectors_d2 = effectors_parent

    antigenic_d2, passenger_d2 = mutation!(rng, counter, λ, p_a)
    ii2::ReactionRef = deriveds_parent.immune_inhibition
    if rand(rng) < p_e # mutation causes immune escape
        # modify mutations of the daughter cell
        ## all mutations become passenger
        append!(daughter2.passenger, daughter2.antigenic, antigenic_d2, passenger_d2)
        empty!(daughter2.antigenic)
        # update effector cells
        # all related effector cells will lose a target cancer cell
        for effector in effectors_d2
            _delete!(
                deriveds(effector).immune_inhibition,
                (ii2.chunk_index, ii2.slot_index),
            )
            setval_diff!(effector, (c = -1,), rs)
        end
        # clear dependent of reactions
        # this will not only clear from the cbl, but also from the cdl, ccl, iil
        # because they share the same reference
        empty!(effectors_d2)
        ## update parent reactions of the daughter cell
        # we need to free this cancer cell from immune inhibition reactions
        # and create a new reaction node with zero inhibition rate
        replace!(ii2, duplicate_zero(ii2))
    else
        # modify mutations of the daughter cell
        ## append new antigenic and passenger mutations to the daughter cell respectively
        append!(daughter2.passenger, passenger_d2)
        append!(daughter2.antigenic, antigenic_d2)
        ## no need to update parent reactions of the daughter cell, because we reuse the old cell

        # create new effector cells and related reactions for new antigens
        n_neoantigen = length(antigenic_d2)
        if n_neoantigen != 0
            n_old_antigen = length(effectors_d2)
            growend!(effectors_d2, n_neoantigen)
            effector_new!(
                effectors_d2,
                rs,
                (ii2.chunk_index, ii2.slot_index),
                effector_example,
                n_old_antigen,
                rng,
                dist_i,
                dist_a,
            )
        end
    end
    return
end

function cancer_death!(ref::ReactionRef, rs::NamedTuple)
    cancer::CancerCell, effectors = children(ref)
    die!(cancer, rs)
    ii_index = deriveds(cancer).immune_inhibition
    for effector in effectors
        _delete!(deriveds(effector).immune_inhibition, ii_index)
        setval_diff!(effector, (c = -1,), rs)
    end
    return
end

function effector_birth!(ref::ReactionRef, rs::NamedTuple)
    effector::CLTPopulation = children(ref)
    setval_diff!(effector, (e = 1,), rs)
    return
end

function effector_death!(ref::ReactionRef, rs::NamedTuple)
    effector::CLTPopulation = children(ref)
    setval_diff!(effector, (e = -1,), rs)
    return
end

"""
    effector_new!(
        effectors::AbstractVector{<:EffectorCellPopulation},
        eml::ReactionCons{EM}, edl::ReactionCons{ED},
        ial::ReactionCons{IA}, iel::ReactionCons{IE},
        ii::ImmumeInhition,
        effector_example::EffectorCellPopulation,
        offset::Integer,
        rng::Random.AbstractRNG,
        fi, fa
    )

Create new effector cells and related reactions for new antigens.
Fill the `effectors` starting from `offset` with new effector cells,
and push related reactions to reaction lists `eml`, `edl`, `ial`, `iel`.
"""
function effector_new!(
    effectors::AbstractVector{<:CLTPopulation},
    rs,
    ii_slot,
    effector_example::CLTPopulation,
    offset::Integer, # the offset of the first new effector cell
    rng::Random.AbstractRNG,
    dist_i,
    dist_a,
)
    ET = typeof(effector_example)
    ds = refall(rs, deriveds(effector_example))
    eml = rs.clt_migration
    edl = rs.clt_death
    ial = rs.immune_activation
    iel = rs.immune_escape
    for i = offset+1:length(effectors)
        # A cell may not always have a related effector cell
        # So we need a example of effector as a constructor

        # all reactions are zero because the effector population is zero now
        em = duplicate_zero(ds.clt_migration)
        em_slot = popslot!(eml)
        ed = duplicate_zero(ds.clt_death)
        ed_slot = popslot!(edl)
        ia = duplicate_zero(ds.immune_activation)
        ia_slot = popslot!(ial)
        ie = duplicate_zero(ds.immune_escape)
        ie_slot = popslot!(iel)

        effector = ET(
            rand(rng, dist_i),
            rand(rng, dist_a),
            (
                clt_migration = em_slot,
                clt_death = ed_slot,
                immune_activation = ia_slot,
                immune_inhibition = [ii_slot],
                immune_escape = ie_slot,
            ),
        )
        effectors[i] = effector

        # push new reaction nodes related to this effector cell to reaction lists
        unsafe_insert!(eml, em_slot, em, effector)
        unsafe_insert!(edl, ed_slot, ed, effector)
        unsafe_insert!(ial, ia_slot, ia, effector)
        unsafe_insert!(iel, ie_slot, ie, effector)
    end
    return effectors
end

"""
    refall(rs::NamedTuple, inds::NamedTuple{names})

Create a reference to reaction nodes in `rs` with indices in `inds`.

This function is a generated function for better performance.
"""
@generated function refall(rs::NamedTuple, inds::NamedTuple{names}) where {names}
    ex = Expr(:tuple)
    for name in names
        push!(ex.args, :($name = unsafe_ref(rs.$name, inds.$name)))
    end
    return ex
end

function _delete!(indices, index)
    # there is one and only one index in the indices
    deleteat!(indices, findfirst(==(index), indices))
end
