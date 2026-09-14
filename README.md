# PowerNumbers.jl

[![CI](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/dlfivefifty/PowerNumbers.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl/graph/badge.svg)](https://codecov.io/gh/dlfivefifty/PowerNumbers.jl)

`PowerNumbers.jl` provides number-like types for asymptotic expressions involving powers of an infinitesimal parameter (and logarithmic terms) as it approaches `0`.

## Installation

```julia
using Pkg
Pkg.add("PowerNumbers")
```

## Quick start

```julia
using PowerNumbers

eps = PowerNumber(1.0, 1.0)
x = 2 + 3*eps^0.5
y = 1 - eps

x * y
x + y
```

## Running tests

From the package repository root:

```julia
using Pkg
Pkg.activate(; temp=true)
Pkg.develop(path=".")
Pkg.test("PowerNumbers"; coverage=true)
```
