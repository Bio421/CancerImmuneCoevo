# This file defines the reaction nodes in the model.
# For the details about the methods, please refer see `reactionlist.jl`

"""
A node represents cell division of a cancer cell, whose rate is a constant.
"""
struct CancerCellDivision{R} <: Reaction{R}
    coeff::R
end

Base.getindex(node::CancerCellDivision) = node.coeff

duplicate(node::N) where {N<:CancerCellDivision} = N(node.coeff)
duplicate_zero(node::CancerCellDivision) = duplicate(node)

remove_child!(::CancerCellDivision, ::CancerCell) = true

effect!(ref::ReactionRef{<:Real,<:CancerCellDivision}, args...) = cancer_div!(ref, args...)

"""
A node represents intrinsic cell death of a cancer cell, whose rate is a constant.
"""
struct CancerCellDeath{R} <: Reaction{R}
    coeff::R
end

Base.getindex(node::CancerCellDeath) = node.coeff

duplicate(node::N) where {N<:CancerCellDeath} = N(node.coeff)
duplicate_zero(node::CancerCellDeath) = duplicate(node)

remove_child!(::CancerCellDeath, ::CancerCell) = true

effect!(ref::ReactionRef{<:Real,<:CancerCellDeath}, args...) = cancer_death!(ref, args...)

"""
A node represents migration of a type of effector cells, whose rate is a constant.
"""
struct CLTMigration{R} <: Reaction{R}
    coeff::R
end

Base.getindex(node::CLTMigration) = node.coeff

duplicate(node::N) where {N<:CLTMigration} = N(node.coeff)
duplicate_zero(node::CLTMigration) = duplicate(node)

update!(::CLTMigration{R}, ::CLTPopulation, ::Any, ::Any) where {R} = zero(R)
remove_child!(::CLTMigration, ::CLTPopulation) = true

effect!(ref::ReactionRef{<:Real,<:CLTMigration}, args...) = effector_birth!(ref, args...)

"""
A node represents intrinsic cell death of a type of effector cells.
The rate is linear to the number of effector cells of this type.
"""
mutable struct CLTDeath{R} <: Reaction{R}
    rate::R
    const coeff::R
end

Base.getindex(node::CLTDeath) = node.rate

duplicate(node::N) where {N<:CLTDeath} = N(node.rate, node.coeff)
duplicate_zero(node::N) where {N<:CLTDeath} = N(zero(eltype(node)), node.coeff)

update!(::CLTDeath{R}, ::CLTPopulation, ::Any, ::Any) where {R} = zero(R)
function update!(node::CLTDeath, ::CLTPopulation, ::Any, diff::@NamedTuple{e::I}) where {I}
    diff = diff.e * node.coeff
    node.rate += diff
    return diff
end
remove_child!(::CLTDeath, ::CLTPopulation) = true

effect!(ref::ReactionRef{<:Real,<:CLTDeath}, args...) = effector_death!(ref, args...)

"""
A node representing activation effect of cancer cells to effector cells.
Note, the cancer cell number is stored in effector cell population node,
so this node is not dependent on cancer cell node.
"""
mutable struct Activation{R,F} <: Reaction{R}
    rate::R
    const func::F
end

Base.getindex(node::Activation) = node.rate

duplicate(node::N) where {N<:Activation} = N(node.rate, node.func)
duplicate_zero(node::N) where {N<:Activation} = N(zero(eltype(node)), node.func)

function update!(node::Activation, reason::CLTPopulation, ::Any, ::Any)
    orig = node.rate
    new = node.func(c_num(reason) * immunogenicity(reason), e_num(reason))
    node.rate = new
    return new - orig
end
remove_child!(::Activation, ::CLTPopulation) = true

effect!(ref::ReactionRef{<:Real,<:Activation}, args...) = effector_birth!(ref, args...)

"""
A node representing immune inhibition from effector cells to a cancer cell.
"""
mutable struct ImmumeInhition{R,F} <: Reaction{R}
    rate::R
    e_sum::R # sum(antigenicity * e_num)
    const func::F
end

Base.getindex(node::ImmumeInhition) = node.rate

duplicate(node::N) where {N<:ImmumeInhition} = N(node.rate, node.e_sum, node.func)
duplicate_zero(node::N) where {N<:ImmumeInhition} =
    (R = eltype(node); N(zero(R), zero(R), node.func))

update!(::ImmumeInhition{R}, ::CLTPopulation, ::Any, ::@NamedTuple{c::I}) where {R,I} =
    zero(R)
function update!(
    node::ImmumeInhition,
    reason::CLTPopulation,
    ::Any,
    diff::@NamedTuple{e::I}
) where {I}
    node.e_sum += diff.e * antigenicity(reason)
    orig = node.rate
    node.rate = node.func(1, node.e_sum)
    return node.rate - orig
end
remove_child!(::ImmumeInhition, ::CancerCell) = true
# This should never happen, because CLTPopulation dies only when all related cancer cells die.
# In this case, the immune inhibition node should already be removed
# remove_child!(::ImmumeInhition, ::CLTPopulation) = false

effect!(ref::ReactionRef{<:Real,<:ImmumeInhition}, args...) = cancer_death!(ref, args...)

# the same as Activation
"""
A node representing inhibition from cancer cells to effector cells.
"""
mutable struct ImmuneEscape{R<:Real,F} <: Reaction{R}
    rate::R
    const func::F
end

Base.getindex(node::ImmuneEscape) = node.rate

duplicate(node::N) where {N<:ImmuneEscape} = N(node.rate, node.func)
duplicate_zero(node::N) where {N<:ImmuneEscape} = N(zero(eltype(node)), node.func)

function update!(node::ImmuneEscape, reason::CLTPopulation, ::Any, ::Any)
    orig = node.rate
    new = node.func(e_num(reason), c_num(reason))
    node.rate = new
    return new - orig
end
remove_child!(::ImmuneEscape, ::CLTPopulation) = true

effect!(ref::ReactionRef{<:Real,<:ImmuneEscape}, args...) = effector_death!(ref, args...)
