module PowerNumbers

using Base, DualNumbers, LinearAlgebra, Infinities

import Base: convert, *, +, -, ==, <, <=, >, !, !=, >=, /, ^, \, in, isapprox, promote_rule
import Base: exp, atanh, log1p, abs, log, inv, real, imag, conj, sqrt,
                sin, cos, cbrt, abs2, log10, log2, exp, exp2, expm1,
                sin, cos, tan, sec, csc, cot, sind, cosd, tand, secd, cscd, cotd, asin, acos,
                atan, asec, acsc, acot, asind, acosd, atand, asecd, acscd, acotd, sinh, cosh,
                tanh, sech, csch, coth, asinh, acosh, asech, acsch, acoth, deg2rad, rad2deg,
                zero, one, isless, iszero, sign, signbit, isinf, isnan, isfinite, isreal, Real,
                muladd, eps, float, angle

import DualNumbers: Dual, realpart, epsilon, dual
import Infinities: InfiniteCardinal, ℵ₀

# `ℵ₀` is part of the type's vocabulary now, so it is re-exported
export PowerNumber, LogNumber, ϵ, ℵ₀, InfiniteCardinal

include("LogNumber.jl")


"""
    PowerNumber(A, B, α, β)

represents a power series of the form `A*ϵ^α + B*ϵ^β + o(ϵ^β)` where `α ≤ β`.
When α == β we impose B = 0.

`PowerNumber` is a `Real` number type, so its coefficients `A` and `B` are real.
Complex expansions are represented as `Complex{<:PowerNumber}`, that is, as a pair of
real expansions, and the four-argument constructor returns one when given complex
coefficients.

The two exponents are typed independently, so that a value known exactly can carry
`β == ℵ₀` (`Infinities.InfiniteCardinal{0}`) without forcing the exponents to be floats:
converting an integer gives an expansion whose coefficients and orders are all integers.
"""
struct PowerNumber{T<:Real,V<:Real,W<:Real} <: Real
    A::T
    B::T
    α::V
    β::W
    function PowerNumber{T,V,W}(A, B, α, β) where {T<:Real,V<:Real,W<:Real}
        α > β && error("Must have α ≤ β")
        α == -Inf && error("Must have α≤ β")
        new{T,V,W}(A, B, α, β)
    end
end

## exponent (order) arithmetic
#
# `ℵ₀` and `Inf` both mark an order beyond which nothing is known, but they are different
# types, so the few places that have to move between them go through these helpers.

"""
    exactorder(W)

the order marking a value that is known exactly, as a `W`: `Inf` for a float exponent and
`ℵ₀` for an integer one, which `Infinities` compares and adds like any other integer.
"""
exactorder(::Type{W}) where W<:AbstractFloat = W(Inf)
exactorder(::Type{<:Real}) = ℵ₀

# the order type able to hold both a `W` and an exact marker, and both a `W` and a finite
# order (a type parameter of `ℵ₀` alone cannot represent `0`)
exactordertype(::Type{W}) where W<:AbstractFloat = W
exactordertype(::Type{W}) where W<:Real = promote_type(W, InfiniteCardinal{0})
finiteordertype(::Type{<:InfiniteCardinal}) = Integer
finiteordertype(::Type{W}) where W<:Real = W

ordertype(::Type{V}, ::Type{W}) where {V<:Real,W<:Real} = promote_type(V, W)
ordertype(::Type{V}, ::Type{<:InfiniteCardinal}) where V<:AbstractFloat = V
ordertype(::Type{<:InfiniteCardinal}, ::Type{W}) where W<:AbstractFloat = W

convertorder(::Type{W}, α) where W<:Real = convert(W, α)
convertorder(::Type{W}, ::InfiniteCardinal) where W<:AbstractFloat = W(Inf)

# `≈` on orders: `ℵ₀ ≈ Inf` has no common type to promote to, so infinite orders (which
# are exact markers, not measurements) are compared with `==` instead
const InfiniteOrder = Union{InfiniteCardinal,Infinities.RealInfinity,Infinities.Infinity}
approxorder(α, β; opts...) = ≈(α, β; opts...)
approxorder(α::InfiniteOrder, β; opts...) = α == β
approxorder(α, β::InfiniteOrder; opts...) = α == β
approxorder(α::InfiniteOrder, β::InfiniteOrder; opts...) = α == β

# `_pn` builds either a `PowerNumber` or, for complex coefficients, the equivalent
# `Complex{<:PowerNumber}`.  All constructors and operations funnel through it so that
# complex data never has to be stored inside a `PowerNumber`, and so that the canonical
# form (no `B` when `α == β`, no vanishing leading coefficient) is imposed in one place.
# The exponents keep their own types and are never promoted against each other.
function _pn(A::Real, B::Real, α::Real, β::Real)
    a, b = promote(A, B)
    if α == β
        PowerNumber{typeof(a),typeof(β),typeof(β)}(a+b, zero(a), β, β)
    elseif iszero(a) && iszero(b) && isfinite(β)
        # a vanishing expansion is o(ϵ^β), so record it at order β rather than as a zero leading term
        PowerNumber{typeof(a),typeof(β),typeof(β)}(a, b, β, β)
    elseif iszero(a) && !iszero(b)
        PowerNumber{typeof(b),typeof(β),typeof(β)}(b, zero(b), β, β)
    else
        PowerNumber{typeof(a),typeof(α),typeof(β)}(a, b, α, β)
    end
end
_pn(A::Complex, B::Complex, α::Real, β::Real) =
    complex(_pn(real(A), real(B), α, β), _pn(imag(A), imag(B), α, β))
_pn(A::Number, B::Number, α::Real, β::Real) = _pn(promote(A,B)..., α, β)

PowerNumber(A::Number, B::Number, α::Real, β::Real) = _pn(A, B, α, β)
PowerNumber(A::Number, α::Real) = PowerNumber(A, zero(A), α, α)
PowerNumber(A::Number) = PowerNumber(A, zero(A), 0, ℵ₀)
PowerNumber(A::Complex) = PowerNumber(A, zero(A), 0, ℵ₀) # disambiguate from Base's Real(::Complex)

PowerNumber{T,V,W}(z::PowerNumber) where {T<:Real,V<:Real,W<:Real} =
    PowerNumber{T,V,W}(convert(T, z.A), convert(T, z.B), convertorder(V, z.α), convertorder(W, z.β))
PowerNumber{T,V,W}(a::Real) where {T<:Real,V<:Real,W<:Real} =
    PowerNumber{T,V,W}(convert(T,a), zero(T), convert(V, 0), exactorder(W))

# `PowerNumber{T,V}` is taken to mean both exponents of type `V`
PowerNumber{T,V}(A, B, α, β) where {T<:Real,V<:Real} = PowerNumber{T,V,V}(A, B, α, β)
PowerNumber{T,V}(z::PowerNumber) where {T<:Real,V<:Real} = PowerNumber{T,V,V}(z)
PowerNumber{T,V}(a::Real) where {T<:Real,V<:Real} = PowerNumber{T,V,V}(a)

# `Infinities` claims `(::Type{T})(::Infinity) where T<:Real`, so spell out that an
# infinite argument is just a coefficient like any other
for Inf∞ in (:(Infinities.Infinity), :(Infinities.RealInfinity))
    @eval begin
        PowerNumber{T,V,W}(a::$Inf∞) where {T<:Real,V<:Real,W<:Real} = PowerNumber{T,V,W}(convert(T, a))
        PowerNumber{T,V}(a::$Inf∞) where {T<:Real,V<:Real} = PowerNumber{T,V,V}(convert(T, a))
    end
end

promote_rule(::Type{PowerNumber{T,V,W}}, ::Type{S}) where {T,V,W,S<:Real} =
    PowerNumber{promote_type(S,T),finiteordertype(V),exactordertype(W)}
promote_rule(::Type{PowerNumber{T1,V1,W1}}, ::Type{PowerNumber{T2,V2,W2}}) where {T1,V1,W1,T2,V2,W2} =
    PowerNumber{promote_type(T1,T2),ordertype(V1,V2),ordertype(W1,W2)}

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
# `zero`/`one` have to preserve the type, including the order types: handing back a
# `ℵ₀`-ordered value from a float-ordered type drags the whole computation over with it
zero(::Type{PowerNumber{T,V,W}}) where {T,V,W} = PowerNumber{T,V,W}(zero(T))
one(::Type{PowerNumber{T,V,W}}) where {T,V,W} = PowerNumber{T,V,W}(one(T))
zero(::Type{PowerNumber{T,V}}) where {T,V} = PowerNumber{T,V,V}(zero(T))
one(::Type{PowerNumber{T,V}}) where {T,V} = PowerNumber{T,V,V}(one(T))

eps(::Type{<:PowerNumber{T}}) where T = eps(T)
float(P::PowerNumber) = PowerNumber(float(P.A), float(P.B), P.α, P.β)
float(::Type{PowerNumber{T,V,W}}) where {T,V,W} = PowerNumber{float(T),V,W}

function (x::PowerNumber)(ε)
    (; A,B,α,β) = x
    A*ε^α + B*ε^β
end

function +(x::PowerNumber, y::PowerNumber)
    a,b,α,β = x.A,x.B,x.α,x.β
    c,d,γ,δ = y.A,y.B,y.α,y.β
    if approxorder(δ, β) && d != 0
        return +(PowerNumber(a,b+d,α,β), PowerNumber(c,0,γ,δ))
    elseif δ < β #we assume β < δ
        return +(y, x)
    elseif γ > β || c == 0
        return x
    elseif approxorder(γ, β)
        return PowerNumber(a,b+c,α,β)
    elseif γ < β && γ > α
        return PowerNumber(a,c,α,γ)
    elseif approxorder(γ, α)
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

# A `LogNumber` has no room for a power of `ϵ`, so only the `ϵ^0` coefficient of a power
# number survives a combination with one: the `ϵ^β` term times `log ϵ` is `o(1)`, and
# `ϵ^β` added to `s*log ϵ + c` is likewise beyond what the result can hold.
function _logconst(a::AnyPowerNumber)
    A, B, α, β = terms(a)
    α > 0 && return zero(A) # smaller than anything a LogNumber can resolve (ϵ*log ϵ → 0)
    α < 0 && throw(DomainError(a, "diverges faster than log ϵ, so cannot combine with a LogNumber"))
    A
end

# These have to take the whole complex number rather than being left to Base's
# componentwise `Complex` arithmetic: the real and imaginary parts of a
# `Complex{<:PowerNumber}` may carry different orders, and only their merged leading term
# is at `ϵ^0`. They also settle which of the generic `PowerNumber`-and-`Real` and
# `LogNumber`-and-`Real` methods applies, now that both types are `Real`.
for P in (:PowerNumber, :(Complex{<:PowerNumber})), L in (:LogNumber, :(Complex{<:LogNumber}))
    @eval begin
        *(a::$P, l::$L) = _logconst(a) * l
        *(l::$L, a::$P) = l * _logconst(a)
        +(a::$P, l::$L) = _logconst(a) + l
        +(l::$L, a::$P) = l + _logconst(a)
        -(a::$P, l::$L) = _logconst(a) - l
        -(l::$L, a::$P) = l - _logconst(a)
        /(a::$P, l::$L) = _logconst(a) / l
        /(l::$L, a::$P) = l / _logconst(a)
    end
end

# Combining the two always collapses to a `LogNumber`, so that is what they promote to.
# Without this the greedy `PowerNumber`-with-any-`Real` rule would nest a log number
# inside a power number wherever Base promotes before operating — `Base._mulsub`, which
# complex `muladd` goes through, is one such place.
promote_rule(::Type{PowerNumber{T,V,W}}, ::Type{LogNumber{S}}) where {T,V,W,S} = LogNumber{promote_type(T,S)}
promote_rule(::Type{LogNumber{S}}, ::Type{PowerNumber{T,V,W}}) where {T,V,W,S} = LogNumber{promote_type(T,S)}

isapprox(a::LogNumber, b::PowerNumber; opts...) = ≈(a, PowerNumber(b)(1); opts...)
isapprox(a::PowerNumber, b::LogNumber; opts...) = ≈(a(1), b; opts...)

LogNumber(a::PowerNumber{T}) where T = LogNumber{T}(a)
LogNumber{T}(a::PowerNumber) where T<:Real = LogNumber{T}(zero(T), convert(T, _logconst(a)))


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

# for α == 0 the error order is untouched, which is worth spelling out so that an exact
# value keeps its `ℵ₀` rather than picking up a `-0.0` and turning into a plain `∞`
_pow(A, B, α, β, p) =
    α == β ? (A^p, zero(A), α*p, α*p) :
             (A^p, (A^(p-1))*B*p, p*α, iszero(α) ? β : β+(p-1)*α)

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
                        approxorder(a.α, b.α; opts...) && approxorder(a.β, b.β; opts...)

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
