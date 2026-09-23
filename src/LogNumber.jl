"""
    LogNumber(s, c)

represents `s*log(ε) + c` as `ε → 0`.

`LogNumber` is a `Real` number type, so `s` and `c` are real.  Complex log expansions are
represented as `Complex{<:LogNumber}`, that is, as a pair of real ones, and the
two-argument constructor returns one when given complex parts.
"""
struct LogNumber{T<:Real} <: Real
    s::T
    c::T
end

# mirrors `_pn`: complex data lives in a `Complex{<:LogNumber}`, never inside a `LogNumber`
_ln(s::T, c::T) where T<:Real = LogNumber{T}(s, c)
_ln(s::Real, c::Real) = _ln(promote(s, c)...)
_ln(s::Complex, c::Complex) = complex(_ln(real(s), real(c)), _ln(imag(s), imag(c)))
_ln(s::Number, c::Number) = _ln(promote(s, c)...)

LogNumber(s::Number, c::Number) = _ln(s, c)
LogNumber(x::Real) = _ln(zero(x), x)
LogNumber(x::Complex) = _ln(zero(x), x) # disambiguate from Base's (::Type{<:Real})(::Complex)
LogNumber(x::LogNumber) = x
LogNumber{T}(x::Real) where T<:Real = LogNumber{T}(zero(T), convert(T, x))
LogNumber{T}(x::LogNumber) where T<:Real = LogNumber{T}(convert(T, x.s), convert(T, x.c))

# `Infinities` claims `(::Type{T})(::Infinity) where T<:Real`, so spell out that an
# infinite argument is the finite part
for Inf∞ in (:(Infinities.Infinity), :(Infinities.RealInfinity))
    @eval begin
        LogNumber(x::$Inf∞) = LogNumber(0.0, convert(Float64, x))
        LogNumber{T}(x::$Inf∞) where T<:Real = LogNumber{T}(convert(T, x))
    end
end

"""
    AnyLogNumber

either a real [`LogNumber`](@ref) or a `Complex` one.
"""
const AnyLogNumber = Union{LogNumber,Complex{<:LogNumber}}

@inline logpart(l::LogNumber) = l.s
@inline realpart(l::LogNumber) = l.c
logpart(l::Complex{<:LogNumber}) = complex(logpart(real(l)), logpart(imag(l)))
realpart(l::Complex{<:LogNumber}) = complex(realpart(real(l)), realpart(imag(l)))

# the ordinary number a finite `LogNumber` stands for
_value(l::LogNumber) = l.c

Base.promote_rule(::Type{LogNumber{T}}, ::Type{S}) where {T,S<:Real} = LogNumber{promote_type(T,S)}
Base.promote_rule(::Type{LogNumber{T}}, ::Type{LogNumber{S}}) where {T,S} = LogNumber{promote_type(T,S)}

==(a::LogNumber, b::LogNumber) = logpart(a) == logpart(b) && realpart(a) == realpart(b)
Base.isapprox(a::LogNumber, b::LogNumber; opts...) = ≈(logpart(a), logpart(b); opts...) && ≈(realpart(a), realpart(b); opts...)
Base.isapprox(a::LogNumber, b::Real; opts...) = ≈(a(1), b; opts...)
Base.isapprox(a::Real, b::LogNumber; opts...) = ≈(a, b(1); opts...)

(l::LogNumber)(ε) = logpart(l)*log(ε) + realpart(l)

# `zero`/`one` must preserve the type, or a `LogNumber`-valued computation silently drops
# back to plain numbers the first time one of them is called on the element type
zero(::Type{LogNumber{T}}) where T = LogNumber{T}(zero(T), zero(T))
one(::Type{LogNumber{T}}) where T = LogNumber{T}(zero(T), one(T))
zero(::LogNumber{T}) where T = zero(LogNumber{T})
one(::LogNumber{T}) where T = one(LogNumber{T})
iszero(l::LogNumber) = iszero(l.c) && iszero(l.s)

float(l::LogNumber) = LogNumber(float(l.s), float(l.c))
float(::Type{LogNumber{T}}) where T = LogNumber{float(T)}
eps(::Type{LogNumber{T}}) where T = eps(T)

for f in (:+, :-)
    @eval begin
        $f(a::LogNumber, b::LogNumber) = LogNumber($f(a.s, b.s), $f(a.c, b.c))
        $f(l::LogNumber, b::Real) = LogNumber(l.s, $f(l.c, b))
        $f(a::Real, l::LogNumber) = LogNumber($f(l.s), $f(a, l.c))
    end
end

-(l::LogNumber) = LogNumber(-l.s, -l.c)

*(l::LogNumber, b::Real) = LogNumber(l.s*b, l.c*b)
*(a::Real, l::LogNumber) = LogNumber(a*l.s, a*l.c)

# a product of two genuine log numbers has a `log(ε)^2` term, which this type cannot hold;
# it is only defined when one of them is finite
function *(a::LogNumber, b::LogNumber)
    iszero(a.s) && return _value(a) * b
    iszero(b.s) && return a * _value(b)
    throw(ArgumentError("$a * $b has a log(ε)^2 term, which a LogNumber cannot represent"))
end

muladd(x::Real, y::LogNumber, z::LogNumber) = x*y + z

function /(a::LogNumber, b::LogNumber)
    if !isinf(b)
        a / _value(b)
    elseif !isinf(a)
        _value(a) / b
    else
        a.s / b.s
    end
end
/(l::LogNumber, b::Real) = LogNumber(l.s/b, l.c/b)
# Base would promote to `Complex{<:LogNumber}` and run its complex division, which forms
# products of two log numbers; dividing the two parts is both correct and representable
/(l::LogNumber, b::Complex) = LogNumber(l.s/b, l.c/b)

function /(a::Real, l::LogNumber)
    isinf(l) && throw(DomainError(l))
    a / _value(l)
end

exp(l::LogNumber) = PowerNumber(exp(l.c), 0, l.s, l.s)
expm1(l::LogNumber) = exp(l) - 1

function abs(l::LogNumber)
    isinf(l) && throw(DomainError(l))
    abs(_value(l))
end

# as ε → 0⁺ the log part dominates, and `log ε → -∞`, so a larger `s` is a smaller number.
# Comparison with any other `Real` promotes to a `LogNumber` first.
isless(a::LogNumber, b::LogNumber) = a.s == b.s ? isless(a.c, b.c) : isless(b.s, a.s)
<(a::LogNumber, b::LogNumber) = isless(a, b)
<=(a::LogNumber, b::LogNumber) = !isless(b, a)

isinf(a::LogNumber) = !iszero(a.s)
isnan(a::LogNumber) = isnan(a.s) || isnan(a.c)
isfinite(a::LogNumber) = !isinf(a) && !isnan(a)

sign(a::LogNumber) = isinf(a) ? -sign(a.s) : sign(a.c)
signbit(a::LogNumber) = isinf(a) ? !signbit(a.s) : signbit(a.c)

Base.show(io::IO, x::LogNumber) = print(io, "($(logpart(x)))log ε + $(realpart(x))")
