module PowerNumbers

using Base, DualNumbers, LinearAlgebra

import Base: convert, *, +, -, ==, <, <=, >, !, !=, >=, /, ^, \, in, isapprox, promote_rule
import Base: exp, atanh, log1p, abs, log, inv, real, imag, conj, sqrt,
                sin, cos, cbrt, abs2, log10, log2, exp, exp2, expm1,
                sin, cos, tan, sec, csc, cot, sind, cosd, tand, secd, cscd, cotd, asin, acos,
                atan, asec, acsc, acot, asind, acosd, atand, asecd, acscd, acotd, sinh, cosh,
                tanh, sech, csch, coth, asinh, acosh, asech, acsch, acoth, deg2rad, rad2deg,
                zero, one, isless, iszero, sign, signbit, isinf, isnan, isfinite, isreal, Real,
                muladd, eps, float, angle

import DualNumbers: Dual, realpart, epsilon, dual

export PowerNumber, LogNumber, ϵ

include("LogNumber.jl")


"""
    PowerNumber(A, B, α, β)

represents a power series of the form `A*ϵ^α + B*ϵ^β + o(ϵ^β)` where `α ≤ β`.
When α == β we impose B = 0.

`PowerNumber` is a `Real` number type, so its coefficients `A` and `B` are real.
Complex expansions are represented as `Complex{<:PowerNumber}`, that is, as a pair of
real expansions, and the four-argument constructor returns one when given complex
coefficients.
"""
struct PowerNumber{T<:Real,V<:Real} <: Real
    A::T
    B::T
    α::V
    β::V
    function PowerNumber{T,V}(A,B,α,β) where {T<:Real,V<:Real}
        α > β && error("Must have α ≤ β")
        α == -Inf && error("Must have α≤ β")
        if α == β
            new{T,V}(A+B,0,β,β)
        elseif A == 0
            new{T,V}(B,0,β,β)
        else
            new{T,V}(A,B,α,β)
        end
    end
end

# `_pn` builds either a `PowerNumber` or, for complex coefficients, the equivalent
# `Complex{<:PowerNumber}`.  All constructors and operations funnel through it so that
# complex data never has to be stored inside a `PowerNumber`.
function _pn(A::Real, B::Real, α::Real, β::Real)
    a, b = promote(A, B)
    c, d = promote(α, β)
    PowerNumber{typeof(a),typeof(c)}(a, b, c, d)
end
_pn(A::Complex, B::Complex, α::Real, β::Real) =
    complex(_pn(real(A), real(B), α, β), _pn(imag(A), imag(B), α, β))
_pn(A::Number, B::Number, α::Real, β::Real) = _pn(promote(A,B)..., α, β)

PowerNumber(A::Number, B::Number, α::Real, β::Real) = _pn(A, B, α, β)
PowerNumber(A::Number, α::Real) = PowerNumber(A, zero(A), α, α)
PowerNumber(A::Number) = PowerNumber(A, zero(A), 0, Inf)
PowerNumber(A::Complex) = PowerNumber(A, zero(A), 0, Inf) # disambiguate from Base's Real(::Complex)

PowerNumber{T,V}(z::PowerNumber) where {T<:Real,V<:Real} =
    PowerNumber{T,V}(convert(T, z.A), convert(T, z.B), convert(V, z.α), convert(V, z.β))
PowerNumber{T,V}(a::Real) where {T<:Real,V<:Real} =
    PowerNumber{T,V}(convert(T,a), zero(T), zero(V), convert(V,Inf))

promote_rule(::Type{PowerNumber{V,W}}, ::Type{T}) where {V,W,T<:Real} =
    PowerNumber{promote_type(T,V),W}
promote_rule(::Type{PowerNumber{T,S}}, ::Type{PowerNumber{V,W}}) where {T,S,V,W} =
    PowerNumber{promote_type(T,V),promote_type(W,S)}

const ϵ = PowerNumber(1.0,1.0)

"""
    AnyPowerNumber

either a real [`PowerNumber`](@ref) or a `Complex` one.
"""
const AnyPowerNumber = Union{PowerNumber,Complex{<:PowerNumber}}

"""
    terms(z)

return `(A, B, α, β)` with `z == A*ϵ^α + B*ϵ^β + o(ϵ^β)` and `α ≤ β`.  For a
`Complex{<:PowerNumber}` the coefficients are complex and the two leading exponents of
the real and imaginary parts are merged.
"""
terms(z::PowerNumber) = (z.A, z.B, z.α, z.β)

"""
    _twoleading(ts, δ)

merge the `(coefficient, exponent)` pairs `ts` into the canonical `(A, B, α, β)` of the
two leading terms, given that everything beyond `ϵ^δ` is unknown.
"""
function _twoleading(ts, δ)
    kept = sort!([t for t in ts if !iszero(t[1]) && t[2] ≤ δ]; by=last)
    merged = similar(kept, 0)
    for t in kept
        if !isempty(merged) && merged[end][2] == t[2]
            merged[end] = (merged[end][1] + t[1], t[2])
        else
            push!(merged, t)
        end
    end
    if isempty(merged)
        z = zero(ts[1][1])
        (z, z, δ, δ)
    elseif length(merged) == 1
        A, α = merged[1]
        (A, zero(A), α, δ)
    else
        (merged[1][1], merged[2][1], merged[1][2], merged[2][2])
    end
end

apart(z::AnyPowerNumber) = terms(z)[1]
bpart(z::AnyPowerNumber) = terms(z)[2]
alpha(z::AnyPowerNumber) = terms(z)[3]
beta(z::AnyPowerNumber) = terms(z)[4]

PowerNumber(x::Dual) = PowerNumber(realpart(x), epsilon(x), 0, 1)
Dual(x::PowerNumber) = (x.α == 0 && x.β == 1) ? (return Dual(x.A, x.B)) : (throw("α, β must equal 0, 1 to convert to dual."))
dual(x::PowerNumber) = Dual(x)

zero(x::PowerNumber) = PowerNumber(zero(x.A), zero(x.B), x.α, x.β)
zero(::Type{PowerNumber{T,V}}) where {T,V} = PowerNumber(zero(T))
one(::Type{PowerNumber{T,V}}) where {T,V} = PowerNumber(one(T))

eps(::Type{<:PowerNumber{T}}) where T = eps(T)
float(P::PowerNumber) = PowerNumber(float(P.A), float(P.B), P.α, P.β)
float(::Type{PowerNumber{T,V}}) where {T,V} = PowerNumber{float(T),float(V)}

function (x::PowerNumber)(ε)
    (; A,B,α,β) = x
    A*ε^α + B*ε^β
end

function +(x::PowerNumber, y::PowerNumber)
    a,b,α,β = x.A,x.B,x.α,x.β
    c,d,γ,δ = y.A,y.B,y.α,y.β
    if δ ≈ β && d != 0
        return +(PowerNumber(a,b+d,α,β), PowerNumber(c,0,γ,δ))
    elseif δ < β #we assume β < δ
        return +(y, x)
    elseif γ > β || c == 0
        return x
    elseif γ ≈ β
        return PowerNumber(a,b+c,α,β)
    elseif γ < β && γ > α
        return PowerNumber(a,c,α,γ)
    elseif γ ≈ α
        return PowerNumber(a+c,b,α,β)
    else
        return PowerNumber(c,a,γ,α)
    end
end


function +(x::PowerNumber, y::Real)
    a,b,α,β = x.A,x.B,x.α,x.β
    if iszero(α)
        PowerNumber(a+y, b, α, β)
    elseif iszero(β)
        PowerNumber(a, b+y, α, β)
    elseif α > 0
        PowerNumber(y, a, 0, α)
    elseif β > 0
        PowerNumber(a, y, α, 0)
    else # both parts blow up so constants disapper
        x
    end
end

+(y::Real, x::PowerNumber) = +(x, y)

function *(x::PowerNumber, y::PowerNumber)
    a,b,α,β = x.A,x.B,x.α,x.β
    c,d,γ,δ = y.A,y.B,y.α,y.β
    # the four products are exact monomials; multiplying the o(ϵ^β) of one by the leading
    # term of the other is what sets the error order of the product
    _pn(_twoleading(((a*c,α+γ), (b*c,β+γ), (a*d,α+δ), (b*d,β+δ)), min(β+γ, α+δ))...)
end

*(x::PowerNumber, y::Real) = PowerNumber(y*x.A,y*x.B,x.α,x.β)
*(y::Real, x::PowerNumber) = x*y

muladd(x::Real, y::PowerNumber, z::Real) = x*y + z

function *(a::PowerNumber, l::LogNumber)
    @assert a.α == 0 && a.β == 1
    a.A * l
end

function *(l::LogNumber, a::PowerNumber)
    @assert a.α == 0 && a.β == 1
    l * a.A
end

LogNumber(a::PowerNumber{T}) where T = LogNumber{T}(a)

function LogNumber{T}(a::PowerNumber) where T
    if a.α == 0 && a.β > 0
        LogNumber{T}(zero(T), a.A)
    elseif a.α == 0 && a.β == 0
        LogNumber{T}(zero(T), a.A + a.B)
    else
        error("not implemented")
    end
end


-(x::PowerNumber) = PowerNumber(-x.A,-x.B,x.α,x.β)
-(x::PowerNumber, y::PowerNumber) = x + (-y)
-(x::PowerNumber, y::Real) = x + (-y)
-(y::Real, x::PowerNumber) = (-x) + y

# the shared kernels below take the canonical data `(A,B,α,β)` of a real or complex
# power number and return the canonical data of the result.

function _inv(A, B, α, β)
    α == Inf && error("Not defined for α = Inf")
    α != β ? (1/A, -B*A^(-2), -α, β-2*α) : (1/A, zero(A), -α, -α)
end

inv(x::PowerNumber) = _pn(_inv(terms(x)...)...)

/(z::PowerNumber, x::PowerNumber) = z*inv(x)
/(z::PowerNumber, x::Real) = z*inv(x)
/(x::Real, z::PowerNumber) = x*inv(z)

_pow(A, B, α, β, p) =
    α == β ? (A^p, zero(A), α*p, α*p) : (A^p, (A^(p-1))*B*p, p*α, β+(p-1)*α)

^(z::PowerNumber, p::Integer) = invoke(^, Tuple{Number,Integer}, z, p)
^(z::PowerNumber, p::Number) = _pn(_pow(terms(z)..., p)...)
^(z::PowerNumber, p::Rational) = z^float(p) # disambiguate from Base

sqrt(z::PowerNumber) = z^0.5
cbrt(z::PowerNumber) = z^(1/3)

iszero(z::PowerNumber) = z.α > 0 || (iszero(z.A) && (z.β > 0 || iszero(z.B)))
isnan(z::PowerNumber) = isnan(z.A)
isinf(z::PowerNumber) = z.α < 0 || isinf(z.A)
isfinite(z::PowerNumber) = !isinf(z) && !isnan(z)

==(a::PowerNumber, b::PowerNumber) = iszero(a - b)
==(a::Real, b::PowerNumber) = PowerNumber(a) == b
==(a::PowerNumber, b::Real) = a == PowerNumber(b)
==(a::AbstractIrrational, b::PowerNumber) = PowerNumber(a) == b # disambiguate from Base
==(a::PowerNumber, b::AbstractIrrational) = a == PowerNumber(b)

isapprox(a::PowerNumber, b::PowerNumber; opts...) = ≈(a.A, b.A; opts...) && ≈(a.B, b.B; opts...) &&
                                                ≈(a.α, b.α; opts...) && ≈(a.β, b.β; opts...)

function isapprox(a::PowerNumber, b::Real; opts...)
    a.α > 0 && return false
    iszero(a.α) && return isapprox(a.A, b; opts...)
    isapprox(zero(a.A), b; opts...)
end
isapprox(b::Real, a::PowerNumber; opts...) = isapprox(a, b; opts...)

function _log(A, B, α, β)
    iszero(A) && error("Cannot evaluate log at 0")
    LogNumber(promote(α, log(A))...)
end

log(z::PowerNumber) = _log(terms(z)...)
log1p(z::PowerNumber) = log(z+1)

atanh(z::PowerNumber) = (log(1+z)-log(1-z))/2

functionlist = (:log10, :log2, :exp, :exp2, :expm1,
                :sin, :cos, :tan, :sec, :csc, :cot, :sind, :cosd, :tand, :secd, :cscd, :cotd, :asin, :acos,
                :atan, :asec, :acsc, :acot, :asind, :acosd, :atand, :asecd, :acscd, :acotd, :sinh, :cosh,
                :tanh, :sech, :csch, :coth, :asinh, :acosh, :asech, :acsch, :acoth, :deg2rad, :rad2deg)

# applies `f` to the two leading terms via dual-number arithmetic
function _applyfun(f, A, B, α, β)
    if iszero(α)
        fz = f(dual(A,B))
        p = β
    elseif α > 0
        fz = f(dual(zero(A),A))
        p = α
    else
        error("Alpha must be non-negative")
    end
    (realpart(fz), epsilon(fz), zero(p), p)
end

for op in functionlist
    @eval $op(z::PowerNumber) = _pn(_applyfun($op, terms(z)...)...)
end

sign(z::PowerNumber) = sign(z.A)
signbit(z::PowerNumber) = signbit(z.A)
abs(z::PowerNumber) = sign(z) * z

# Comparison with any other `Real` promotes to a `PowerNumber` first, so ordering only
# needs to be defined between two of them.  `Real` requires `<`, whose Base fallback for
# two values of the same type just errors out.
<(x::PowerNumber, y::PowerNumber) = isless(x, y)
<=(x::PowerNumber, y::PowerNumber) = !isless(y, x)

function isless(z::PowerNumber, w::PowerNumber)
    d = z - w
    iszero(d) ? false : signbit(d.A)
end

function Base.show(io::IO, x::PowerNumber)
    if x.α == x.β || iszero(x.B)
        @assert iszero(x.B)
        if iszero(x.α)
            print(io, "$(x.A) + o(ϵ^$(x.β))")
        else
            print(io, "($(x.A))ϵ^$(x.α) + o(ϵ^$(x.β))")
        end
    else
        if iszero(x.α)
            print(io, "$(x.A) + ($(x.B))ϵ^$(x.β) + o(ϵ^$(x.β))")
        else
            print(io, "($(x.A))ϵ^$(x.α) + ($(x.B))ϵ^$(x.β) + o(ϵ^$(x.β))")
        end
    end
end

include("complex.jl")

end # module
