"""
    LogNumber(s, c)

represents s*log(ε) + c as ε -> 0
"""
struct LogNumber{T} <: Number
    s::T
    c::T
end

LogNumber(s, t) = LogNumber(promote(s, t)...)
LogNumber(x::Number) = LogNumber{typeof(x)}(x)
LogNumber{T}(x::Number) where T = LogNumber{T}(zero(T), x)
LogNumber(x::LogNumber) = x
LogNumber{T}(x::LogNumber) where T = LogNumber{T}(x.s, x.c)

const RealOrComplex = Union{Real, Complex}

function Real(x::LogNumber{T}) where T
    isreal(x) || throw(InexactError(:Real, x))
    isinf(x) && convert(T, sign(x)*Inf)
    x.c
end

function RealOrComplex(x::LogNumber{T}) where T
    isinf(x) && convert(T, sign(x)*Inf)
    x.c
end


#@inline logpart(z::Number) = zero(z)
#@inline realpart(z::Number) = z

@inline logpart(l::LogNumber) = l.s
@inline realpart(l::LogNumber) = l.c


Base.promote_rule(::Type{LogNumber{T}}, ::Type{R}) where {T, R<:RealOrComplex} = LogNumber{promote_type(T,R)}
Base.convert(::Type{LogNumber}, z::LogNumber) = z
Base.convert(::Type{LogNumber}, z::RealOrComplex) = LogNumber(0, z)

==(a::LogNumber, b::LogNumber) = logpart(a) == logpart(b) && realpart(a) == realpart(b)
Base.isapprox(a::LogNumber, b::LogNumber; opts...) = ≈(logpart(a), logpart(b); opts...) && ≈(realpart(a), realpart(b); opts...)
Base.isapprox(a::LogNumber, b::RealOrComplex; opts...) = ≈(a(1),b; opts...)
Base.isapprox(a::RealOrComplex, b::LogNumber; opts...) = ≈(a,b(1); opts...)

(l::LogNumber)(ε) = logpart(l)*log(ε) + realpart(l)

zero(l::LogNumber{T}) where T = LogNumber(zero(T), zero(T))
one(l::LogNumber{T}) where T = one(T)
iszero(l::LogNumber{T}) where T = iszero(l.c) && iszero(l.s)


for op in (:real, :imag, :conj)
    @eval $op(l::LogNumber) = LogNumber($op(logpart(l)), $op(realpart(l)))
end

for f in (:+, :-)
    @eval begin
        $f(a::LogNumber, b::LogNumber) = LogNumber($f(a.s, b.s), $f(a.c, b.c))
        $f(l::LogNumber, b::RealOrComplex) = LogNumber(l.s, $f(l.c, b))
        $f(a::RealOrComplex, l::LogNumber) = LogNumber($f(l.s), $f(a, l.c))
    end
end

-(l::LogNumber) = LogNumber(-l.s, -l.c)

for Typ in (:Bool, :Real, :Complex)
    @eval begin
        *(l::LogNumber, b::$Typ) = LogNumber(l.s*b, l.c*b)
        *(a::$Typ, l::LogNumber) = LogNumber(a*l.s, a*l.c)
    end
end

muladd(x::RealOrComplex, y::LogNumber, z::LogNumber) = x*y + z

function /(a::LogNumber, b::LogNumber)
    if !isinf(b)
        a / RealOrComplex(b)
    elseif !isinf(a)
        RealOrComplex(a) / b
    else
        a.s / b.s
    end
end
/(l::LogNumber, b::RealOrComplex) = LogNumber(l.s/b, l.c/b)

function /(a::RealOrComplex, l::LogNumber)
    isinf(l) && throw(DomainError(l))
    a / RealOrComplex(l)
end

exp(l::LogNumber) = PowerNumber(exp(l.c), 0, l.s, l.s)
expm1(l::LogNumber) = exp(l) - 1

function abs(l::LogNumber)
    isinf(l) && throw(DomainError(l))
    abs(RealOrComplex(l))
end
isless(a::LogNumber, b::LogNumber) = a.s == b.s ? isless(a.c, b.c) : isless(b.s, a.s)
function isless(a::LogNumber, b::Number)
    if isinf(a)
        signbit(a) && (!isinf(b) || !signbit(b))
    else
        isless(RealOrComplex(a), b)
    end
end

isinf(a::LogNumber) = !iszero(a.s)
isreal(a::LogNumber) = isreal(a.s) && isreal(a.c)

sign(a::LogNumber) = isinf(a) ? -sign(a.s) : sign(a.c)
signbit(a::LogNumber) = isinf(a) ? !signbit(a.s) : sign(a.c)


Base.show(io::IO, x::LogNumber) = print(io, "($(logpart(x)))log ε + $(realpart(x))")
