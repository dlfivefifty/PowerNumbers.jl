# PowerNumbers.jl

[![CI](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl/graph/badge.svg)](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl)

`PowerNumbers.jl` provides number-like types for asymptotic expressions involving powers of an infinitesimal parameter (and logarithmic terms) as it approaches `0`.


## Quick start

```julia
julia> using PowerNumbers

julia> ϵ # analoguous to a dual number
(1)ϵ^1 + o(ϵ^1)

julia> x = 2 + 3sqrt(ϵ) # but we support fractional powers
2.0 + (3.0)ϵ^0.5 + o(ϵ^0.5)

julia> y = 1 - ϵ # simple algebraic relationships work
1 + (-1)ϵ^1 + o(ϵ^1)

julia> x * y
2.0 + o(ϵ^0.0)

julia> x + y
3.0 + (3.0)ϵ^0.5 + o(ϵ^0.5)

julia> log(ϵ) + 5 # we also support logarithms
(1.0)log ε + 5.0
```
