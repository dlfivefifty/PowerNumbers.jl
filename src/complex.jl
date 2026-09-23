## Complex power numbers
#
# A complex expansion is stored as a `Complex{<:PowerNumber}`, i.e. as a pair of real
# expansions.  Since the real and imaginary parts may carry different exponents, most
# operations first put the pair back into the canonical `A*ϵ^α + B*ϵ^β + o(ϵ^β)` form
# with complex `A` and `B`; this is what `terms` returns.

function terms(z::Complex{<:PowerNumber})
    x, y = reim(z)
    # the error term of the pair is the worse of the two
    δ = min(x.β, y.β)
    _twoleading(((complex(x.A, zero(x.A)), x.α), (complex(x.B, zero(x.B)), x.β),
                 (complex(zero(y.A), y.A), y.α), (complex(zero(y.B), y.B), y.β)), δ)
end

(z::Complex{<:PowerNumber})(ε) = complex(real(z)(ε), imag(z)(ε))

zero(x::Complex{<:PowerNumber}) = complex(zero(real(x)), zero(imag(x)))

# Base's complex kernels are written for floating point numbers and would both lose the
# asymptotics and recurse through `float`, so the ones we need are given here directly.

inv(z::Complex{<:PowerNumber}) = _pn(_inv(terms(z)...)...)
/(z::Complex{<:PowerNumber}, w::Complex{<:PowerNumber}) = z*inv(w)
/(z::Complex{<:PowerNumber}, w::PowerNumber) = complex(real(z)/w, imag(z)/w)
/(z::PowerNumber, w::Complex{<:PowerNumber}) = z*inv(w)
/(z::Complex{<:PowerNumber}, w::Real) = complex(real(z)/w, imag(z)/w)
/(z::Real, w::Complex{<:PowerNumber}) = z*inv(w)
/(z::Complex{<:PowerNumber}, w::Complex) = z*inv(w)
/(z::Complex, w::Complex{<:PowerNumber}) = z*inv(w)

^(z::Complex{<:PowerNumber}, p::Integer) = invoke(^, Tuple{Number,Integer}, z, p)
^(z::Complex{<:PowerNumber}, p::Number) = _pn(_pow(terms(z)..., p)...)
^(z::Complex{<:PowerNumber}, p::Real) = _pn(_pow(terms(z)..., p)...) # disambiguate from Base
^(z::Complex{<:PowerNumber}, p::Rational) = z^float(p)
^(z::Complex{<:PowerNumber}, p::Bool) = p ? z : one(z)

sqrt(z::Complex{<:PowerNumber}) = z^0.5 # Base leaves `cbrt` undefined for `Complex`

abs2(z::Complex{<:PowerNumber}) = real(z)^2 + imag(z)^2
abs(z::Complex{<:PowerNumber}) = sqrt(abs2(z))
angle(z::Complex{<:PowerNumber}) = (A = terms(z)[1]; angle(A))

log(z::Complex{<:PowerNumber}) = _log(terms(z)...)
log1p(z::Complex{<:PowerNumber}) = log(z+1)
atanh(z::Complex{<:PowerNumber}) = (log(1+z)-log(1-z))/2

for op in functionlist
    @eval $op(z::Complex{<:PowerNumber}) = _pn(_applyfun($op, terms(z)...)...)
end

isapprox(a::Complex{<:PowerNumber}, b::Complex{<:PowerNumber}; opts...) =
    ≈(real(a), real(b); opts...) && ≈(imag(a), imag(b); opts...)
isapprox(a::Complex{<:PowerNumber}, b::Number; opts...) =
    ≈(real(a), real(b); opts...) && ≈(imag(a), imag(b); opts...)
isapprox(a::Number, b::Complex{<:PowerNumber}; opts...) = isapprox(b, a; opts...)
isapprox(a::LogNumber, b::Complex{<:PowerNumber}; opts...) = ≈(a, b(1); opts...)
isapprox(a::Complex{<:PowerNumber}, b::LogNumber; opts...) = ≈(a(1), b; opts...)

function Base.show(io::IO, z::Complex{<:PowerNumber})
    A, B, α, β = terms(z)
    if α == β || iszero(B)
        print(io, "($A)ϵ^$α + o(ϵ^$β)")
    else
        print(io, "($A)ϵ^$α + ($B)ϵ^$β + o(ϵ^$β)")
    end
end


## Complex log numbers

(l::Complex{<:LogNumber})(ε) = complex(real(l)(ε), imag(l)(ε))

# `exp(s*log ε + c) = exp(c)*ϵ^s`, which is only a power number when the log part is real
function exp(l::Complex{<:LogNumber})
    s, c = logpart(l), realpart(l)
    iszero(imag(s)) || throw(DomainError(l, "exponent of ϵ must be real"))
    _pn(exp(c), zero(c), real(s), real(s))
end
expm1(l::Complex{<:LogNumber}) = exp(l) - 1

# as above: divide the log and finite parts rather than letting Base run a complex
# division over `LogNumber` components
# a `Real` divisor is fine componentwise, so only the complex case needs taking over
/(l::Complex{<:LogNumber}, b::Complex) = LogNumber(logpart(l)/b, realpart(l)/b)

# compare the log and finite parts as complex numbers rather than componentwise, so that
# a rounding-level real part is weighed against the magnitude of the whole part
isapprox(a::Complex{<:LogNumber}, b::Complex{<:LogNumber}; opts...) =
    ≈(logpart(a), logpart(b); opts...) && ≈(realpart(a), realpart(b); opts...)
isapprox(a::Complex{<:LogNumber}, b::Number; opts...) = ≈(a(1), b; opts...)
isapprox(a::Number, b::Complex{<:LogNumber}; opts...) = ≈(a, b(1); opts...)
isapprox(a::Complex{<:LogNumber}, b::Complex{<:PowerNumber}; opts...) = ≈(a(1), b(1); opts...)
isapprox(a::Complex{<:PowerNumber}, b::Complex{<:LogNumber}; opts...) = ≈(a(1), b(1); opts...)

Base.show(io::IO, l::Complex{<:LogNumber}) = print(io, "($(logpart(l)))log ε + $(realpart(l))")
