# This file contains the definition of the substance nodes in the model.
# There are only two types of substance nodes in the model:
# - CancerCell: a node represents a cancer cell with mutations;
# - CLTPopulation: a node represents a type of effector cells.

"""
    Substance <: AbstractNode

An abstract type represents substance nodes.
Substance nodes are the basic nodes in the SSA graph,
which do not have a simple rate and children.
A substance can be a population, a type of cell, or a type of molecule.
The value of a substance node is the number of given substance.

## Required Methods:

- `Base.getindex(node)`: returns the state of the given `node`;
- `setval_impl!(node, num)`: raw `setval!` implementation, return a boolean value
    indicating whether the given `node` is dead;
- `deriveds(node)`: returns indices of derived nodes in the reaction system;

## Optional Methods:

- `setval!(node, num)`: sets the number of given substance `node` to `num`,
    and updates the parent node(s);
- `setval_diff!(node, diff)`: incremental version of `setval!`;
"""
abstract type Substance end

"""
    setval!(node::Substance, val, rs::NamedTuple)

Set the value of the given `node` to `val`, and update the parent nodes.
"""
function setval!(node::Substance, val, rs::NamedTuple)
    orig = node[]
    if setval_impl!(node, val)
        die!(node, rs)
    else
        diff = diff_impl(orig, val)
        foreach_ds(rs, deriveds(node)) do d
            update!(d, node, orig, diff)
        end
    end
    return node
end

"""
    setval_diff!(node::Substance, diff, rs::NamedTuple)

Incrementally set the value of the given `node` by `diff`, and update the parent nodes.
"""
function setval_diff!(node::Substance, diff, rs::NamedTuple)
    orig = node[]
    val = apply_impl(orig, diff)
    if setval_impl!(node, val)
        die!(node, rs)
    else
        foreach_ds(rs, deriveds(node)) do d
            update!(d, node, orig, diff)
        end
    end
    return node
end

update!(ds::Tuple, node, orig, diff) = foreach(d -> update!(d, node, orig, diff), ds)

"""
    die!(node::Substance, rs::NamedTuple)

Clean up the given `node` from the reaction system `rs`.
This function called in `setval!` when `setval_impl!` returns `true`.
"""
function die!(node::Substance, rs::NamedTuple)
    foreach_ds(rs, deriveds(node)) do d
        remove_child!(d, node)
    end
end

"""
    foreach_ds(f, rs::NamedTuple, ds_ind::NamedTuple{names})

Apply the function `f` to each derived node in the `rs` with the corresponding index in `ds_ind`.
"""
@generated function foreach_ds(f, rs::NamedTuple, ds_ind::NamedTuple{names}) where {names}
    ex = Expr(:block)
    for name in names
        push!(ex.args, :(f(unsafe_ref(rs.$name, ds_ind.$name))))
    end
    return ex
end

"""
    diff_impl(orig, val)

Calculate the difference between `val` and `orig`, default to `val - orig`.
For custom type of `orig` and `val`, this method should be overloaded.
"""
diff_impl(orig, val) = val - orig

"""
    apply_impl(orig, diff)

Apply the difference `diff` to the original value `orig`, default to `orig + diff`.
For custom type of `orig` and `diff`, this method should be overloaded.
"""
apply_impl(val, diff) = val + diff

"""
    CancerCell{M,D} <: Substance

A node that represents a cancer cell with mutations.
The mutation type (neoantigen or passenger) is not stored in the node,
but inferred from related parent reaction nodes.

`antigenic`, `passenger`, and `mutations` is used to get the antigenic mutations.
"""
struct CancerCell{M,D<:NamedTuple} <: Substance
    passenger::M # passenger mutations
    antigenic::M # antigenic mutations
    deriveds::D
    CancerCell{M,D}(passenger::M, antigenic::M, deriveds::D) where {M,D} =
        new{M,D}(passenger, antigenic, deriveds)
    CancerCell{M,D}(passenger::M, deriveds::D) where {M,D} =
        new{M,D}(passenger, M(), deriveds)
    CancerCell{M,D}(deriveds) where {M,D} = new{M,D}(M(), M(), deriveds)
end

"Get antigenic mutations of the given `cell`."
antigenic(cell::CancerCell) = cell.antigenic
"Get passenger mutations of the given `cell`."
passenger(cell::CancerCell) = cell.passenger
"Get all mutations of the given `cell`."
mutations(cell::CancerCell) = LazyVCat(antigenic(cell), passenger(cell))

deriveds(cell::CancerCell) = cell.deriveds
Base.eltype(::Type{<:CancerCell{T}}) where {T} = T

# Each cancer cell ships with a unique set of mutations,
# so a instance of CancerCell only has one element.
Base.getindex(::CancerCell) = 1

"""
    CLTPopulation{I,R,FA,FI,FE} <: Substance

A node that represents a population of CLTs, carrying the same TCR and can recognize the same antigen.

`e_num`, `c_num` is used to get the number of effector cells and cancer cells recognized by them.
`antigenicity` and `immunogenicity` is used to get the antigenicity and immunogenicity of the effector cells.
"""
mutable struct CLTPopulation{I<:Real,R<:Real,D<:NamedTuple} <: Substance
    e_num::I # number of this type of effector cells
    c_num::I # number of cancer cells recognized by this type of effector cells
    const immunogenicity::R
    const antigenicity::R
    const deriveds::D
    CLTPopulation{I,R,D}(e::I, c::I, i::R, a::R, parent::D) where {I,R,D} =
        new{I,R,D}(e, c, i, a, parent)
    CLTPopulation{I,R,D}(i::R, a::R, parent::D) where {I,R,D} =
        new{I,R,D}(zero(I), one(I), i, a, parent)
end

"Get the number of effector cells of the given `node`."
e_num(node::CLTPopulation) = node.e_num
"Get the number of cancer cells recognized by the effector cells of the given `node`."
c_num(node::CLTPopulation) = node.c_num

"Set the number of effector cells of the given `node`."
set_e_num!(node::CLTPopulation, e::I) where {I} = node.e_num = e
"Set the number of cancer cells recognized by the effector cells of the given `node`."
set_c_num!(node::CLTPopulation, c::I) where {I} = node.c_num = c

immunogenicity(node::CLTPopulation) = node.immunogenicity
antigenicity(node::CLTPopulation) = node.antigenicity

deriveds(node::CLTPopulation) = node.deriveds
Base.eltype(::Type{<:CLTPopulation{I}}) where {I} = @NamedTuple{e::I, c::I}
Base.getindex(node::CLTPopulation) = (e = e_num(node), c = c_num(node))

setval_impl!(node::CLTPopulation{I}, val::@NamedTuple{e::I}) where {I} =
    (set_e_num!(node, val.e); false)
setval_impl!(node::CLTPopulation{I}, val::@NamedTuple{c::I}) where {I} =
    (set_c_num!(node, val.c); iszero(val.c))

diff_impl(orig::@NamedTuple{e::I, c::I}, val::@NamedTuple{e::I}) where {I} =
    (e = val.e - orig.e,)
diff_impl(orig::@NamedTuple{e::I, c::I}, val::@NamedTuple{c::I}) where {I} =
    (c = val.c - orig.c,)

apply_impl(orig::@NamedTuple{e::I, c::I}, diff::@NamedTuple{e::I}) where {I} =
    (e = orig.e + diff.e,)
apply_impl(orig::@NamedTuple{e::I, c::I}, diff::@NamedTuple{c::I}) where {I} =
    (c = orig.c + diff.c,)
