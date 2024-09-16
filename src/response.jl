"""
A module defines the three types of functional response.
"""
module FunctionalResponse

abstract type ResponseFunction <: Function end

"""
    Linear(a)

Type I functional response coefficient. The function is

```math
    f(x, y) = a * x
```
"""
struct Linear{R<:Real} <: ResponseFunction
    a::R
end
(f::Linear)(x::Real, ::Real) = f.a * x

struct Bilinear{R<:Real} <: ResponseFunction
    a::R
end
(f::Bilinear)(x::Real, y::Real) = f.a * x * y
const TypeI = Bilinear

"""
    TypeII(a, h)

Type II functional response coefficient. The function is

```math
    f(x, y) = a * x * y / (1 + k * x)
```

where `k = a * h`, and `1 / h` is the maximum of the function.
"""
struct TypeII{R<:Real} <: ResponseFunction
    a::R
    k::R # a * h, where 1/h is the maximum of the function
    function TypeII(a::Real, h::Real)
        a, h = promote(a, h)
        k = a * h
        return new{typeof(a)}(a, k)
    end
end
(f::TypeII)(x::Real, y::Real) = f.a * x * y / (1 + f.k * x)

"""
    TypeIII(a, h, n)

Type III functional response coefficient. The function is

```math
    f(x, y) = a * x^n * y / (1 + k * x^n)
```

where `k = a * h`, and `1 / h` is the maximum of the function.
"""
struct TypeIII{R<:Real,I<:Integer} <: ResponseFunction
    a::R
    k::R # a * h, where 1/h is the maximum of the function
    n::I
end

function TypeIII(a::Real, h::Real, n::Integer)
    a, h = promote(a, h)
    k = a * h
    return TypeIII(a, k, n)
end

function (f::TypeIII)(x::Real, y::Real)
    x_n = x^f.n
    return f.a * x_n * y / (1 + f.k * x_n)
end

"""
    Reversed(f)

Swap the `x` and `y` arguments of a `ResponseFunction`.
"""
struct Reversed{F<:ResponseFunction} <: ResponseFunction
    f::F
end

reversed(f::Reversed) = f.f
reversed(f::Bilinear) = f # Bilinear is symmetric
reversed(f::ResponseFunction) = Reversed(f)

(f::Reversed)(x::Real, y::Real) = f.f(y, x)

end
