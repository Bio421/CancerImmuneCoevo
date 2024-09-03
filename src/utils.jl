"""
    LazyVCat

A lazy vertical concatenation of AbstractVector{T}.

For the case of high mutation rate, each cells may contains a large number of mutations.
Use this to avoid too many memory allocations and RAM usage.
"""
struct LazyVCat{T,I} <: AbstractVector{T}
    inners::I # any iterable and indexable collection of AbstractVector{T}
    function LazyVCat{T,I}(inners::I) where {T,I}
        # these checks are necessary and no overhead
        # because for most cases, these checks are done at compile time
        Base.require_one_based_indexing(inners)
        foreach(inners) do inner
            Base.require_one_based_indexing(inner)
            @assert(
                isa(inner, AbstractVector{T}),
                "inner is not an abstract vector with elements of type $T"
            )
        end
        new{T,I}(inners)
    end
    LazyVCat{T,I}() where {T,I<:Vector} = new{T,I}(I())
end

# homogeneous inners constructor
LazyVCat(inners::V...) where {T,V<:AbstractVector{T}} = LazyVCat{T,Vector{V}}([inners...])

Base.copy(l::L) where {L<:LazyVCat} = L(copy(l.inners))

# If the `inners` is a homogeneous vector, we can implement the `append!` method
# and `vcat` method for `LazyVCat` to make it more efficient.
function Base.append!(l::LazyVCat{T,I}, v::V) where {T,V<:AbstractVector{T},I<:Vector{V}}
    isempty(v) && return l
    push!(l.inners, v)
    l
end
Base.append!(l::LazyVCat{T,I}, v::LazyVCat{T,I}) where {T,I<:Vector} =
    (append!(l.inners, v.inners); l)
Base.vcat(
    l::L,
    rest::Union{L,V}...,
) where {T,V<:AbstractVector{T},I<:Vector{V},L<:LazyVCat{T,I}} =
    reduce(append!, rest; init = L(copy(l.inners)))

Base.empty!(l::LazyVCat{T,I}) where {T,I<:Vector} = empty!(l.inners)

# interface methods for AbstractArray
Base.size(l::LazyVCat) = (sum(length, l.inners; init = 0),)
Base.IndexStyle(::Type{<:LazyVCat}) = IndexLinear()
function Base.getindex(l::LazyVCat, i::Int)
    i < 1 && throw(BoundsError(l, i))
    j = i
    for inner in l.inners
        len = length(inner)
        if j > len
            j -= len
        else
            # it is safe to use `@inbounds` here
            # because the inner is not an offset array
            return @inbounds inner[j]
        end
    end
    throw(BoundsError(l, i))
end
function Base.iterate(l::LazyVCat, state::Tuple{Int,Int} = (1, 1))
    nth, j = state
    n = length(l.inners)
    nth > n && return
    while j > length(@inbounds l.inners[nth])
        nth += 1
        nth > n && return
        j = 1
    end
    @inbounds(l.inners[nth][j]), (nth, j + 1)
end
