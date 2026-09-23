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

## Exact values

The two orders `α` and `β` are separate type parameters, so a value that is known exactly
can mark its error order with `ℵ₀` (`Infinities.InfiniteCardinal{0}`) instead of `Inf`.
An integer then converts to an expansion made entirely of integers, and stays that way
through arithmetic:

```julia
julia> PowerNumber(5)
5 + o(ϵ^ℵ₀)

julia> typeof(PowerNumber(5))
PowerNumber{Int64, Int64, InfiniteCardinal{0}}

julia> PowerNumber(1,2,0,1) * PowerNumber(3) # (1 + 2ϵ) * 3, all integers
3 + (6)ϵ^1 + o(ϵ^1)

julia> 2 + ϵ # a float ϵ still uses float orders, with Inf for an exact value
2.0 + (1.0)ϵ^1.0 + o(ϵ^1.0)
```

`ℵ₀` compares and adds like any other integer (`ℵ₀ + 1 == ℵ₀`, `1 < ℵ₀`, `ℵ₀ == Inf`), so
it slots into the order arithmetic without special cases. It is what lets an integer be
promoted into an expansion with integer orders at all: an `Int` field could not have held
`Inf`.

## Complex expansions

`PowerNumber <: Real`, so a complex expansion is an ordinary `Complex` whose real and
imaginary parts are each a `PowerNumber`. Complex arithmetic therefore comes from `Base`,
and `Complex{<:PowerNumber}` works wherever a complex number is expected:

```julia
julia> z = 2 + 3im*ϵ
(2.0 + 0.0im)ϵ^0.0 + (0.0 + 3.0im)ϵ^1.0 + o(ϵ^1.0)

julia> typeof(z)
Complex{PowerNumber{Float64, Float64, Float64}}

julia> 1/z
(0.5 + 0.0im)ϵ^-0.0 + (0.0 - 0.75im)ϵ^1.0 + o(ϵ^1.0)

julia> sqrt(-1 + ϵ + 0im)
(0.0 + 1.0im)ϵ^0.0 + (0.0 - 0.5im)ϵ^1.0 + o(ϵ^1.0)
```

The real and imaginary parts may carry different exponents; operations that need the
canonical `A*ϵ^α + B*ϵ^β + o(ϵ^β)` form merge them back together and truncate to the two
leading terms. Passing complex coefficients to the `PowerNumber` constructor builds the
corresponding `Complex{<:PowerNumber}`.

`LogNumber` works the same way — it is `Real`, and a complex log expansion is a
`Complex{<:LogNumber}`:

```julia
julia> LogNumber(2im, im+1)
(0 + 2im)log ε + 1 + 1im

julia> typeof(LogNumber(2im, im+1))
Complex{LogNumber{Int64}}

julia> log((1+im)*ϵ)
(1.0 + 0.0im)log ε + 0.34657359027997264 + 0.7853981633974483im
```

Where the two types meet, only the `ϵ^0` coefficient of the power number survives, since
a `LogNumber` has no room for a power of `ϵ` (and `ϵ*log ϵ → 0`). That is also what they
promote to, so code that promotes before operating agrees:

```julia
julia> (2+ϵ)*LogNumber(1.0,2.0)
(2.0)log ε + 4.0

julia> promote_type(typeof(ϵ), LogNumber{Float64})
LogNumber{Float64}
```

## Comparison with DualNumbers.jl, ForwardDiff.jl and TaylorSeries.jl

All four types push a truncated expansion in a small parameter through ordinary Julia
code. They differ in which expansions they can hold.

`DualNumbers.Dual`, `ForwardDiff.Dual` and `TaylorSeries.Taylor1` expand about a
**regular** point. Their exponents are non-negative integers fixed by the type rather
than by the value — `a + bε` with `ε² ≡ 0` for the duals, and coefficients of
`1, t, …, t^N` for a `Taylor1` of declared order `N`. No error order has to be tracked,
because everything dropped is of higher order than everything kept.

A `PowerNumber` stores its two exponents as data, so they can be any reals `α ≤ β` —
negative, fractional, and worked out by the arithmetic rather than declared up front. It
also carries the order `β` of what it does *not* know, so a term smaller than the error
is dropped instead of being silently kept. And `log ϵ`, which is not a power of `ϵ` at
all, becomes a `LogNumber`.

That is what makes it usable at a **singular** point. Starting from the infinitesimal
itself (`ϵ`, `Dual(0.0,1.0)`, `ForwardDiff.Dual(0.0,1.0)`, `Taylor1([0.0,1.0],4)`):

|                | `sqrt`                     | `inv`                      | `log`                      |
| -------------- | -------------------------- | -------------------------- | -------------------------- |
| `PowerNumbers` | `(1.0)ϵ^0.5 + o(ϵ^0.5)`    | `(1.0)ϵ^-1.0 + o(ϵ^-1.0)`  | `(1.0)log ε + 0.0`         |
| `DualNumbers`  | `Dual{Float64}(0.0,Inf)`   | `Dual{Float64}(Inf,-Inf)`  | `Dual{Float64}(-Inf,Inf)`  |
| `ForwardDiff`  | `Dual{Nothing}(0.0,Inf)`   | `Dual{Nothing}(Inf,-Inf)`  | `Dual{Nothing}(-Inf,Inf)`  |
| `TaylorSeries` | `DomainError`              | `ArgumentError`            | `DomainError`              |

The other three are not wrong here; a half-integer power, a pole and a logarithm are
simply not expressible in their expansions. Away from the singularity they all agree —
each of these computes the derivative of `exp` at `1`:

```julia
julia> exp(1+ϵ)
2.718281828459045 + (2.718281828459045)ϵ^1.0 + o(ϵ^1.0)

julia> exp(Dual(1.0,1.0))
2.718281828459045 + 2.718281828459045ɛ

julia> exp(ForwardDiff.Dual(1.0,1.0))
Dual{Nothing}(2.718281828459045,2.718281828459045)

julia> exp(Taylor1([1.0,1.0],3))
 2.718281828459045 + 2.718281828459045 t + 1.3591409142295225 t² + 0.45304697140984085 t³ + 𝒪(t⁴)
```

### Which to reach for

- **ForwardDiff.jl** for derivatives. It carries `N` partials at once, gives you
  gradients, Jacobians and (by nesting) Hessians, and is the standard choice for
  automatic differentiation in Julia. `PowerNumbers` does one variable and no more.
- **TaylorSeries.jl** for many orders, or for several variables via `TaylorN`. A
  `PowerNumber` keeps two terms and no more, so it cannot give you a fourth derivative.
  `Taylor1` handles fractional powers perfectly well as long as the expansion point is
  regular: `Taylor1([1.0,1.0],3)^0.5` is `1.0 + 0.5t - 0.125t² + 0.0625t³ + 𝒪(t⁴)`.
- **DualNumbers.jl** for a plain first-order dual with no extra machinery.
- **PowerNumbers.jl** for a limit at a point where the function blows up or has a branch
  point, and the leading exponent is not known before the computation runs. It grew out
  of evaluating Cauchy and Stieltjes transforms at the endpoints of an interval, where
  the answer goes like `ϵ^(-1/2)` or `log ϵ` and the singular parts have to cancel
  between neighbouring pieces before a finite answer appears.
