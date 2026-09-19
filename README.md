# PowerNumbers.jl

[![CI](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl/graph/badge.svg)](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl)

`PowerNumbers.jl` provides number-like types for asymptotic expressions involving powers of an infinitesimal parameter (and logarithmic terms) as it approaches `0`.


## Quick start

```julia
julia> using PowerNumbers

julia> ϵ # analoguous to a dual number
(1.0)ϵ^1.0 + o(ϵ^1.0)

julia> x = 2 + 3sqrt(ϵ) # but we support fractional powers
2.0 + (3.0)ϵ^0.5 + o(ϵ^0.5)

julia> y = 1 - ϵ # simple algebraic relationships work
1.0 + (-1.0)ϵ^1.0 + o(ϵ^1.0)

julia> x * y
2.0 + (3.0)ϵ^0.5 + o(ϵ^0.5)

julia> x + y
3.0 + (3.0)ϵ^0.5 + o(ϵ^0.5)

julia> log(ϵ) + 5 # we also support logarithms
(1.0)log ε + 5.0
```

## Complex expansions

`PowerNumber <: Real`, so a complex expansion is an ordinary `Complex` whose real and
imaginary parts are each a `PowerNumber`. Complex arithmetic therefore comes from `Base`,
and `Complex{<:PowerNumber}` works wherever a complex number is expected:

```julia
julia> z = 2 + 3im*ϵ
(2.0 + 0.0im)ϵ^0.0 + (0.0 + 3.0im)ϵ^1.0 + o(ϵ^1.0)

julia> typeof(z)
Complex{PowerNumber{Float64, Float64}}

julia> 1/z
(0.5 + 0.0im)ϵ^-0.0 + (0.0 - 0.75im)ϵ^1.0 + o(ϵ^1.0)

julia> sqrt(-1 + ϵ + 0im)
(0.0 + 1.0im)ϵ^0.0 + (0.0 - 0.5im)ϵ^1.0 + o(ϵ^1.0)
```

The real and imaginary parts may carry different exponents; operations that need the
canonical `A*ϵ^α + B*ϵ^β + o(ϵ^β)` form merge them back together and truncate to the two
leading terms. Passing complex coefficients to the `PowerNumber` constructor builds the
corresponding `Complex{<:PowerNumber}`.
