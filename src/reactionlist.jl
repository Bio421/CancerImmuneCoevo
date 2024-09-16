# This file implements the ReactionList data structure,
# which is used to store a list of reactions and their children.
# Compared to normal lists, this data structure is have similar performance
# for iteration and insertion, but much faster for deletion (we don't need to shift elements).

# EDIT: -----------------------------------------------------------------------
# For whom want to implement SSA for a similar model:
# This data structure is suboptimal, a better implementation is just use a simple vector with
# `swap_deleteat!` function, which deletes an element by swapping it with the last element.
# We can use swap_deleteat! in our SSA, because the order of reactions does not matter.
# -----------------------------------------------------------------------------

using DataStructures: SortedSet

"""
    Reaction{R,C}

A abstract type represents a reaction.

## Required Methods:

- `Base.getindex(node::Reaction)`: returns the rate of the `node`;
- `duplicate(node, children)`: creates a clone of the `node` with given `children`;
- `duplicate_zero(node, children)`: creates a clone of the `node` with given `children` and zero rate;
- `remove_child!(node, child)`: removes the `child` node from the `node`'s children;
  the return value indicates whether this node should be removed from the reaction list.

## Other related methods:

- `effect(ref::ReactionRef{T}, args...) where {T <: Reaction}`: this function is called
    when this reaction is sampled, used to update the state of the system.
"""
abstract type Reaction{R<:Real} end

Base.eltype(::Type{<:Reaction{R}}) where {R} = R

const CHUNK_SIZE = 4096

"""
    Slots

Used to store available slots in a chunk, and available chunks in a list.
We use a sorted set to ensure that always pop the first available slot for
better iteration performance (because we make sure the chunk is more densely packed).
"""
const Slots = SortedSet{Int,Base.Order.ForwardOrdering}
slots(iter) = Slots(Base.Order.ForwardOrdering(), iter)

"""
    Chunk{T,R,C}

A chunk of reactions. It contains a fixed number of reactions and their children.
The `mask` indicates whether a slot is occupied and the `empty_slots` stores the available slots.

## Initialization

When a chunk is created, reactions and children are initialized with undefined values.
All slots are marked as empty and stored in the `empty_slots`.

## Insertion

To insert a new reaction, use `popslot!` to pop an available slot from the `empty_slots`.
Then use `unsafe_insert!` to insert the reaction and its children into the chunk at the given slot.
Note the `popslot!` and `unsafe_insert!` must be called in pairs to ensure the mark is correct.

## Deletion

When a reaction is removed, use `deleteat!` to delete the reaction from the chunk.
The slot is marked as empty and pushed back to the `empty_slots`.

"""
struct Chunk{T,R<:Reaction{T},C}
    reactions::Vector{R}
    children::Vector{C}
    mask::BitVector
    empty_slots::Slots
    function Chunk{T,R,C}(capacity::Int = CHUNK_SIZE) where {T,R,C}
        reactions = Vector{R}(undef, capacity)
        children = Vector{C}(undef, capacity)
        mask = falses(capacity)
        empty_slots = slots(1:capacity)
        return new{T,R,C}(reactions, children, mask, empty_slots)
    end
end

Base.length(chunk::Chunk) = length(chunk.reactions) - length(chunk.empty_slots)
isfull(chunk::Chunk) = isempty(chunk.empty_slots)

popslot!(chunk::Chunk) = pop!(chunk.empty_slots)

Base.@propagate_inbounds function unsafe_insert!(chunk::Chunk, i::Int, r, c)
    chunk.reactions[i] = r
    chunk.children[i] = c
    chunk.mask[i] = true
    return chunk
end

Base.@propagate_inbounds function Base.deleteat!(chunk::Chunk, i::Int)
    chunk.mask[i] = false
    push!(chunk.empty_slots, i)
    return chunk
end

"""
    ReactionList{T,R,C}

A list of reactions stored in a list of chunks.
The `available_chunks` stores the indices of chunks with empty slots.
The `sum` stores the sum of the rates of all reactions in the list,
which will be updated automatically when a reaction is updated, inserted, or deleted.

## Initialization

When a list is created, a chunk with the given capacity is created.

## Insertion

To insert a new reaction, use `popslot!` to pop an available slot from the `available_chunks`.
If there is no available chunk, create a new chunk to store the reactions.
Then use `unsafe_insert!` to insert the reaction and its children into the list at the given slot.

## Deletion

When a reaction is removed, use `deleteat!` to delete the reaction from the corresponding chunk.
Unless the chunk is available before, the chunk will be available after the deletion.

## Iteration

We always insert a new reaction into the first available slot of the first available chunk.
So the iteration is simple and efficient. We iterate over all chunks and all reactions in each chunk.
"""
struct ReactionList{T,R<:Reaction{T},C}
    chunks::Vector{Chunk{T,R,C}}
    available_chunks::Slots
    sum::Base.RefValue{T}
    chunk_size::Int
    function ReactionList{T,R,C}(capacity::Int = CHUNK_SIZE) where {T,R,C}
        chunks = [Chunk{T,R,C}(capacity)]
        available_chunks = slots(1)
        sum = Ref{T}(zero(T))
        return new{T,R,C}(chunks, available_chunks, sum, capacity)
    end
end

ratetype(x) = ratetype(typeof(x))
ratetype(::Type{<:ReactionList{T}}) where {T} = T

Base.eltype(::Type{<:ReactionList{T,R}}) where {T,R} = R
Base.length(l::ReactionList) = sum(length, l.chunks)

isfull(l::ReactionList) = isempty(l.available_chunks)

"""
    popslot!(l::ReactionList)

Pop an available slot to insert a new element.
If the list is full, double the size of the data array.
"""
function popslot!(l::ReactionList)
    if isfull(l)
        n_chunks = length(l.chunks)
        chunk = newchunk(l)
        push!(l.chunks, chunk)
        push!(l.available_chunks, n_chunks + 1)
        n_chunks + 1, popslot!(chunk)
    else
        i = first(l.available_chunks)
        chunk = @inbounds l.chunks[i]
        slot = popslot!(chunk)
        isfull(chunk) && pop!(l.available_chunks)
        i, slot
    end
end
newchunk(l::ReactionList{T,R,C}) where {T,R,C} = Chunk{T,R,C}(l.chunk_size)

"""
    unsafe_insert!(l::ReactionList, (i, slot), r, c)

Insert the reaction `r` and its children `c` into the list `l` at the `slot` of the chunk `i`.

# Safety

- The `i` must be an available slot.
- The `i` must be inbounds.

It's safe to call this function with the result of `popslot!(l::ReactionList)`.
"""
Base.@propagate_inbounds function unsafe_insert!(
    l::ReactionList,
    (i, j)::Tuple{Int,Int},
    r,
    c,
)
    chunk = l.chunks[i]
    unsafe_insert!(chunk, j, r, c)
    l.sum[] += getindex(r)
    return l
end

Base.@propagate_inbounds function Base.deleteat!(l::ReactionList, (i, j)::Tuple{Int,Int})
    chunk = l.chunks[i]
    deleteat!(chunk, j)
    push!(l.available_chunks, i) # no need to check the duplicate
    l.sum[] -= getindex(chunk.reactions[j])
    return l
end
Base.@propagate_inbounds Base.deleteat!(l::ReactionList, (i, slot)) =
    deleteat!(l, (convert(Int, i), convert(Int, slot)))

function Base.iterate(
    l::ReactionList,
    state::Tuple{Int,Int,Int} = (1, 1, length(l.chunks[1])),
)
    # i: current chunk index, j: current slot index, n: number of rest non-empty slots
    (i, j, n) = state
    n_chunks = length(l.chunks)
    chunk = @inbounds l.chunks[i]
    # iterate until no more chunks
    while i <= n_chunks
        # iterare until no more non-empty slots
        while n > 0
            # if the slot is not empty return a ref and next state
            @inbounds if chunk.mask[j]
                return unsafe_ref(l, (i, j)), (i, j + 1, n - 1)
            end
            j += 1
        end
        # if no more non-empty slots, move to the next chunk
        i += 1
        i > n_chunks && return nothing
        j = 1
        chunk = @inbounds l.chunks[i]
        # length(chunk) is the number of non-empty slots instead of CHUNK_SIZE
        n = length(chunk)
    end
    return nothing
end

# we can use the iterate function to implement the following functions
# but we implement them separately for better performance (it may be not obvious)
function sample(l::ReactionList, r::Real)
    subsum = zero(r)
    n_chunks = length(l.chunks)
    i = 1
    while i <= n_chunks
        chunk = @inbounds l.chunks[i]
        j = 1
        n = length(chunk)
        while n > 0
            @inbounds if chunk.mask[j]
                subsum += getindex(@inbounds chunk.reactions[j])
                r <= subsum && return unsafe_ref(l, (i, j))
                n -= 1
            end
            j += 1
        end
        i += 1
    end
    throw(ArgumentError("The random number $(r) out of sum
        ($(subsum), $(l.sum[]), $(sum(r -> getindex(reaction(r)), l)))"))
end

"""
    ReactionRef{T,R,C}

A reference to a reaction in a list of reactions.

The reaction and its children can be accessed by the `reaction` and `children` functions.
"""
struct ReactionRef{T,R<:Reaction{T},C}
    list::ReactionList{T,R,C}
    chunk_index::Int
    slot_index::Int
    ReactionRef{T,R,C}(
        list::ReactionList{T,R,C},
        chunk_index::Int,
        slot_index::Int,
    ) where {T,R,C} = new{T,R,C}(list, chunk_index, slot_index)
end
ReactionRef(l::ReactionList{T,R,C}, i::Int, j::Int) where {T,R,C} =
    ReactionRef{T,R,C}(l, i, j)

unsafe_ref(l::ReactionList, (i, j)::Tuple{Int,Int}) = ReactionRef(l, i, j)

chunk(ref::ReactionRef) = @inbounds ref.list.chunks[ref.chunk_index]

exists(ref::ReactionRef) = @inbounds chunk(ref).mask[ref.slot_index]
reaction(ref::ReactionRef) = @inbounds chunk(ref).reactions[ref.slot_index]
children(ref::ReactionRef) = @inbounds chunk(ref).children[ref.slot_index]

remove!(ref::ReactionRef) = @inbounds deleteat!(ref.list, (ref.chunk_index, ref.slot_index))

function replace!(ref::ReactionRef, r)
    reactions = chunk(ref).reactions
    r_old = @inbounds reactions[ref.slot_index]
    @inbounds reactions[ref.slot_index] = r
    ref.list.sum[] += (getindex(r) - getindex(r_old))
    return ref
end

function remove_child!(ref::ReactionRef, child)
    remove_self = remove_child!(reaction(ref), child)
    remove_self && remove!(ref)
    return remove_self
end
function update!(ref::ReactionRef, reason::Substance, orig, diff)
    rate_diff = update!(reaction(ref), reason, orig, diff)
    iszero(rate_diff) || (ref.list.sum[] += rate_diff)
    return rate_diff
end

duplicate(ref::ReactionRef) = duplicate(reaction(ref))
duplicate_zero(ref::ReactionRef) = duplicate_zero(reaction(ref))

"""
    ReactionView{T,R,C}

A view of a list of reactions, similar to a list of `ReactionRef`.
"""
struct ReactionView{T,R<:Reaction{T},C}
    list::ReactionList{T,R,C}
    indices::Vector{Tuple{Int,Int}}
    ReactionView{T,R,C}(
        list::ReactionList{T,R,C},
        indices::Vector{Tuple{Int,Int}},
    ) where {T,R,C} = new{T,R,C}(list, indices)
end
ReactionView(l::ReactionList{T,R,C}, indices::Vector{Tuple{Int,Int}}) where {T,R,C} =
    ReactionView{T,R,C}(l, indices)

unsafe_ref(l::ReactionList, indices::Vector{Tuple{Int,Int}}) = ReactionView(l, indices)

function Base.push!(view::ReactionView, ref::ReactionRef)
    push!(view.indices, (ref.chunk_index, ref.slot_index))
    return view
end

function Base.delete!(view::ReactionView, ref::ReactionRef)
    deleteat!(view.indices, findfirst(==((ref.chunk_index, ref.slot_index)), view.indices))
    return view
end

function Base.empty!(view::ReactionView)
    empty!(view.indices)
    return view
end

function remove_child!(view::ReactionView, child)
    for (i, j) in view.indices
        ref = unsafe_ref(view.list, (i, j))
        remove_child!(ref, child)
    end
    return length(view.indices) == 0
end
function update!(view::ReactionView, reason::Substance, orig, diff)
    diff_sum = zero(ratetype(view.list))
    filter!(view.indices) do (i, j)
        ref = unsafe_ref(view.list, (i, j))
        exists(ref) || return false
        diff_sum += update!(ref, reason, orig, diff)
        return true
    end
    return diff_sum
end
