# PowerNumbers.jl

[![CI](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl)

`PowerNumbers.jl` provides number-like types for asymptotic expressions involving powers of `ϵ` (and logarithmic terms) as `ϵ → 0`.

## Installation

```julia
using Pkg
Pkg.add("PowerNumbers")
```

## Quick start

```julia
using PowerNumbers

ϵ = PowerNumber(1.0, 1.0)
x = 2 + 3*ϵ^0.5
y = 1 - ϵ

x * y
log1p(-y)
```

## Running tests

```julia
using Pkg
Pkg.test("PowerNumbers"; coverage=true)
```
