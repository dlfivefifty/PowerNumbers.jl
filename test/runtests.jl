using PowerNumbers, Test, IntervalArithmetic, HypergeometricFunctions
import PowerNumbers: PowerNumber, LogNumber, logpart, realpart, apart, bpart, alpha, beta
using Infinities, LinearAlgebra

@testset "RiemannDual -> PowerNumber" begin
    for h in (0.1,0.01), a in (2exp(0.1im),1.1)
        @test log(PowerNumber(0,a,0,1))(h) ≈ log(h*a)
        @test log(PowerNumber(a,Inf,-1,0))(h) ≈ log(a/h)
    end

    for h in (0.1,0.01), a in (2exp(0.1im),1.1)
        @test log1p(PowerNumber(-1,a,0,1))(h) ≈ log(h*a)
        @test log1p(PowerNumber(a,Inf,-1,0))(h) ≈ log(a/h)
    end

    h=0.0001
    for z in (PowerNumber(1,3exp(0.2im),0,1), PowerNumber(1,0.5exp(-1.3im),0,1),
              PowerNumber(-1,3exp(0.2im),0,1), PowerNumber(-1,0.5exp(-1.3im),0,1),
              PowerNumber(-1,1,0,1), PowerNumber(1,-1,0,1))
        @test atanh(z)(h) ≈  atanh(apart(z)*h^alpha(z)+bpart(z)*h^beta(z)) atol = 1E-4
    end

    #z = PowerNumber(1,-0.25,1)
    h = 0.0000001
    #@test speciallog(z)(h) ≈ speciallog(realpart(z)+epsilon(z)*h^alpha(z)) atol=1E-4

    z = PowerNumber(1,2,0,1)
    @test 1/(1-z(h)) ≈ (1/(1-z))(h) rtol=0.0001

    @test real(LogNumber(2im,im+1)) == LogNumber(0,1)
    @test imag(LogNumber(2im,im+1)) == LogNumber(2,1)
    @test conj(LogNumber(2im,im+1)) == LogNumber(-2im,1-im)

end

@testset "PowerNumbers arithmetic" begin
    h = 0.000000001
    @test inv(PowerNumber(1.0,2.0,-1/2,0))(h) ≈ inv(h^(-1/2) + 2) rtol = 0.001
    @test inv(PowerNumber(0.0,1.0,0,1/2))(h) ≈ inv(h^(1/2)) rtol = 0.001
    @test inv(PowerNumber(2.0,1.0,0,1/2))(h) ≈ inv(h^(1/2)+2) rtol = 0.001

    @test real(PowerNumber(2im,im+1,0,0.5)) == PowerNumber(0,1,0,0.5)
    @test imag(PowerNumber(2im,im+1,0,0.5)) == PowerNumber(2,1,0,0.5)
    @test conj(PowerNumber(2im,im+1,0,0.5)) == PowerNumber(-2im,1-im,0,0.5)

    @test (3*PowerNumber(0.01+3im,0.2im+1,0,0.4) - 4*PowerNumber(0.6-1.5im,1.5-0.3im,0,0.4))/3 == PowerNumber(-0.79+5.0im,-1.0 + 0.6im,0,0.4)
end

@testset "sin" begin
    @test sin(sqrt(ϵ))^2 === PowerNumber(1.0,1.0)
    @test sin(sqrt(ϵ))^2.0 === PowerNumber(1.0,1.0)
    @test sin(ϵ)/ϵ === PowerNumber(1.0,0.0,0.0,0.0)
    @test sin(sqrt(ϵ))/sqrt(ϵ) === PowerNumber(1.0,0.0,0.0,0.0)
end

@testset "LogNumber" begin
    z = 1-ϵ
    @test log1p(-z) isa LogNumber
    @test exp(LogNumber(2,3)) == exp(3)*ϵ^2
    HypergeometricFunctions.expm1(LogNumber(2,3))
end

@testset "Rational" begin
    @test PowerNumber(1.,0,0.,1.) + PowerNumber(-1.,0,2.,2.) == PowerNumber(1.,0.,0.,1)
    @test (1 + 1/ϵ + 1/ϵ^2) / (1 + 1/ϵ + 1/ϵ^2) == 1
end

@testset "HypergeometricFunctions" begin
    a,b,c = 1.154,1.2543,1.3543345
    z = 1-ϵ
    @test_throws MethodError _₂F₁(a,b,c,z)

    a,b,c = 1.1,1.2,1.3
    @test_throws MethodError _₂F₁(a,b,c,z)
    log1p(-z)
end

@testset "IntervalArithmetic" begin
    a = @interval(1.0)
    p = PowerNumber(a,a,0,1)
    @test (a^4+a^2-a) == (p^4+p^2-p)(0)
end

@testset "all 0" begin
    @test 0 == PowerNumber(0,1, -1, 2) == PowerNumber(1, 2, 3, 4)
    @test 1 == PowerNumber(1, 2, 0, 1)
    @test 1 ≠ PowerNumber(1, 2, -1, 0)
    @test PowerNumber(1, 2, -2, -1) == PowerNumber(1.0, 2, -2, -1)
    @test PowerNumber(1, 2, -2, -1) ≠ PowerNumber(1.0, 2, -2, 0)
    @test PowerNumber(1, 2, -2, 1) == PowerNumber(1.0, 3, -2, 2)
end

@testset "LogNumber arithmetic" begin
    @test LogNumber{Float64}(3.0) == LogNumber(0.0, 3.0)
    @test LogNumber(0.0, 3.0) == 3.0
    @test 3 + LogNumber(1.0, 2.0) == LogNumber(1.0, 5.0)
    @test 3 - LogNumber(1.0, 2.0) == LogNumber(-1.0, 1.0)
    @test -LogNumber(1.0, 2.0) == LogNumber(-1.0, -2.0)
    @test LogNumber(1.0, 2.0)*2.0 == LogNumber(2.0, 4.0)
    @test 2.0*LogNumber(1.0, 2.0) == LogNumber(2.0, 4.0)
end

@testset "exact values carry β == ℵ₀" begin
    # an integer converts to an expansion whose coefficients and orders are all integers:
    # `ℵ₀` marks "known exactly" where a float exponent would use `Inf`
    @test PowerNumber(5) isa PowerNumber{Int,Int,InfiniteCardinal{0}}
    @test alpha(PowerNumber(5)) === 0
    @test beta(PowerNumber(5)) === ℵ₀
    @test beta(PowerNumber(5)) == Inf
    @test PowerNumber(5)(0.5) == 5
    @test PowerNumber(5) == 5

    # exactness is preserved by arithmetic on exact values
    @test PowerNumber(5) * PowerNumber(3) == 15
    @test beta(PowerNumber(5) * PowerNumber(3)) === ℵ₀
    @test beta(PowerNumber(5) + PowerNumber(3)) === ℵ₀
    @test beta(inv(PowerNumber(4))) === ℵ₀
    @test inv(PowerNumber(4)) == 0.25
    @test beta(sqrt(PowerNumber(4))) === ℵ₀
    @test sqrt(PowerNumber(4)) == 2
    @test beta(exp(PowerNumber(1))) === ℵ₀
    @test PowerNumber(2)^3 === PowerNumber(8)

    # integer data stays integer all the way through
    @test PowerNumber(1,2,0,1) * PowerNumber(3) === PowerNumber(3,6,0,1)
    @test PowerNumber(1,2,0,1) * PowerNumber(3) isa PowerNumber{Int,Int,Int}

    # each exponent keeps its own type; they are never promoted against each other
    @test PowerNumber(1,1) isa PowerNumber{Int,Int,Int}
    @test PowerNumber(1,2,0,1) isa PowerNumber{Int,Int,Int}
    @test PowerNumber(1,2,0,1.5) isa PowerNumber{Int,Int,Float64}
    @test ϵ isa PowerNumber{Float64,Float64,Float64}

    # promoting a plain number into an expansion with integer exponents needs an order
    # that means "exact"; `Inf` could not be stored in an `Int` field, `ℵ₀` can
    e = PowerNumber(1,1)
    @test convert(PowerNumber{Int,Int,InfiniteCardinal{0}}, 3) == PowerNumber(3)
    @test 2 + im + e isa Complex{<:PowerNumber}
    @test (2 + im + e)(0.5) == 2.5 + im
    @test (1 + e)(0.5) == 1.5

    # `Infinities` also claims `(::Type{<:Real})(::Infinity)`; an infinite argument is
    # a coefficient, not an order
    @test apart(PowerNumber{Float64,Float64,Float64}(∞)) === Inf
    @test PowerNumber{Float64,Float64,Float64}(∞) isa PowerNumber{Float64,Float64,Float64}
end

@testset "PowerNumber constructors and conversions" begin
    @test PowerNumber{Float64,Float64}(PowerNumber(1,2,0,1)) == PowerNumber(1.0,2.0,0.0,1.0)
    @test zero(PowerNumber(3.0,4.0,0.0,1.0)) == PowerNumber(0.0,0.0,0.0,1.0)
    @test zero(PowerNumber{Float64,Float64}) == PowerNumber(0.0)
    @test one(PowerNumber{Float64,Float64}) == PowerNumber(1.0)
    @test eps(PowerNumber{Float64,Float64}) == eps(Float64)
end

@testset "PowerNumber addition merging" begin
    # y's leading order strictly between x's leading and subleading orders
    x = PowerNumber(1.0,2.0,0.0,2.0)
    y = PowerNumber(3.0,4.0,1.0,3.0)
    @test x + y == PowerNumber(1.0,3.0,0.0,1.0)

    # y's leading order below x's leading order, with y's subleading order truncated away
    x2 = PowerNumber(5.0,6.0,2.0,4.0)
    y2 = PowerNumber(3.0,4.0,0.0,5.0)
    @test x2 + y2 == PowerNumber(3.0,5.0,0.0,2.0)

    # adding a number where α < 0 < β
    @test PowerNumber(2.0,3.0,-1.0,1.0) + 5 == PowerNumber(2.0,5.0,-1.0,0.0)
end

@testset "PowerNumber multiplication tracks the error order" begin
    # the error order of a product is min(β+γ, α+δ), not the order of the leading product
    @test (1+ϵ)*(1+ϵ) == PowerNumber(1.0,2.0,0.0,1.0)
    @test (1+ϵ)^2 == PowerNumber(1.0,2.0,0.0,1.0)
    @test (1+ϵ)*(2-ϵ) == PowerNumber(2.0,1.0,0.0,1.0)
    @test ϵ*ϵ == PowerNumber(1.0,2.0)
    @test ϵ*inv(ϵ) === PowerNumber(1.0,0.0,0.0,0.0)
    @test sign(1+ϵ)*sqrt((1+ϵ)^2-1) ≈ sqrt(ϵ)*sqrt(2+ϵ)

    h = 1E-8
    for (x,y) in ((1+ϵ, 1+ϵ), (2+3ϵ, PowerNumber(0.5,-1.0,-1.0,0.0)), (sqrt(ϵ), 1-ϵ))
        @test (x*y)(h) ≈ x(h)*y(h) rtol=1E-6
    end
end

@testset "PowerNumber * LogNumber" begin
    # only the leading coefficient survives: (ϵ^β)*log(ϵ) is o(1)
    @test PowerNumber(2.0,3.0,0.0,1.0) * LogNumber(1.0,2.0) == LogNumber(1.0,2.0) * 2.0
    @test LogNumber(1.0,2.0) * PowerNumber(2.0,3.0,0.0,1.0) == LogNumber(1.0,2.0) * 2.0

    # the complex case must not fall through to LogNumber's generic `::Complex` method,
    # which would nest a power number inside the LogNumber
    z = (2.0+3.0im) + (4.0-1.0im)*ϵ
    @test z isa Complex{<:PowerNumber}
    @test z * LogNumber(1.0,2.0) === LogNumber(1.0,2.0) * (2.0+3.0im)
    @test LogNumber(1.0,2.0) * z === LogNumber(1.0,2.0) * (2.0+3.0im)
    @test z * LogNumber(1.0,2.0) isa Complex{<:LogNumber}
end

@testset "PowerNumber misc functions" begin
    @test_throws ErrorException inv(PowerNumber(1.0,2.0,Inf,Inf))
    @test sin(PowerNumber(0.3,0.5,0.0,1.0)) ≈ PowerNumber(sin(0.3), 0.5*cos(0.3), 0.0, 1.0)
    @test_throws ErrorException sin(PowerNumber(1.0,2.0,-1.0,0.0))

    @test sign(PowerNumber(-3.0,1.0,0.0,1.0)) == -1.0
    @test abs(PowerNumber(-3.0,2.0,0.0,1.0)) == PowerNumber(3.0,-2.0,0.0,1.0)

    @test isless(PowerNumber(0.0,1.0,1.0,2.0), 5.0)
    @test !isless(PowerNumber(2.0,1.0,-1.0,0.0), 5.0)
    @test isless(PowerNumber(3.0,1.0,0.0,1.0), 5.0)

    @test sprint(show, PowerNumber(2.0,0.0,1.0,1.0)) == "(2.0)ϵ^1.0 + o(ϵ^1.0)"
    @test sprint(show, PowerNumber(2.0,3.0,0.0,1.0)) == "2.0 + (3.0)ϵ^1.0 + o(ϵ^1.0)"
    @test sprint(show, PowerNumber(2.0,3.0,1.0,2.0)) == "(2.0)ϵ^1.0 + (3.0)ϵ^2.0 + o(ϵ^2.0)"
end


@testset "PowerNumber is Real" begin
    @test PowerNumber <: Real
    @test PowerNumber(1.0,2.0,0.0,1.0) isa Real
    @test ϵ isa Real
    @test real(ϵ) === ϵ
    @test imag(ϵ) == 0
    @test conj(ϵ) === ϵ
    @test isreal(ϵ)

    P = PowerNumber{Float64,Float64,Float64}
    @test promote_type(Float64, P) === P
    @test promote_type(Int, P) === P
    @test promote_type(ComplexF64, P) === Complex{P}

    @test ϵ < 1
    @test 1 > ϵ
    # ϵ and -ϵ agree to leading order, so neither is less than the other
    @test !(-ϵ < ϵ)
    @test !(ϵ < -ϵ)
    @test !(ϵ < ϵ)
    @test sort([1+ϵ, -1+ϵ, ϵ]) == [-1+ϵ, ϵ, 1+ϵ]

    @test signbit(-1+ϵ)
    @test !signbit(1+ϵ)
    @test isfinite(1+ϵ)
    @test isinf(1/ϵ)
    @test !isnan(ϵ)
end

@testset "Complex{PowerNumber}" begin
    h = 1E-8

    # complex coefficients produce a pair of real expansions
    z = (1+im)*ϵ
    @test z isa Complex{PowerNumber{Float64,Float64,Float64}}
    @test real(z) === ϵ
    @test imag(z) === ϵ
    @test PowerNumber(2im,im+1,0,0.5) isa Complex{<:PowerNumber}

    # the two leading terms of the pair are merged back together
    w = -1 + (1+im)*ϵ
    @test apart(w) == -1
    @test bpart(w) == 1+im
    @test alpha(w) == 0
    @test beta(w) == 1
    @test (apart(z), bpart(z), alpha(z), beta(z)) == (1+im, 0, 1, 1)

    # a constant keeps its `o(ϵ^Inf)` error term
    c = complex(PowerNumber(1.0), PowerNumber(2.0))
    @test (apart(c), bpart(c), alpha(c), beta(c)) == (1+2im, 0, 0, Inf)

    # real and imaginary parts may carry different orders
    m = complex(PowerNumber(1.0,2.0,0.0,1.0), PowerNumber(3.0,4.0,2.0,3.0))
    @test (apart(m), bpart(m), alpha(m), beta(m)) == (1, 2, 0, 1)

    @test w(h) ≈ -1 + (1+im)*h
    for f in (inv, sqrt, exp, sin, cos, tanh)
        @test f(w)(h) ≈ f(w(h)) rtol=1E-6
    end
    @test (w^3)(h) ≈ w(h)^3 rtol=1E-6
    @test (w^0.3)(h) ≈ w(h)^0.3 rtol=1E-6

    v = 2+3im+ϵ
    @test (w*v)(h) ≈ w(h)*v(h) rtol=1E-6
    @test (w/v)(h) ≈ w(h)/v(h) rtol=1E-6
    @test (v/w)(h) ≈ v(h)/w(h) rtol=1E-6
    @test (2/w)(h) ≈ 2/w(h) rtol=1E-6
    @test (w/2)(h) ≈ w(h)/2 rtol=1E-6
    @test (w+v)(h) ≈ w(h)+v(h) rtol=1E-6

    @test log(z) == LogNumber(1, log(1+im))
    @test log1p(-1+z) == LogNumber(1, log(1+im))

    @test sprint(show, z) == "(1.0 + 1.0im)ϵ^1.0 + o(ϵ^1.0)"
    @test sprint(show, w) == "(-1.0 + 0.0im)ϵ^0.0 + (1.0 + 1.0im)ϵ^1.0 + o(ϵ^1.0)"
end




@testset "LogNumber is Real" begin
    @test LogNumber <: Real
    @test LogNumber(1.0,2.0) isa Real
    @test real(LogNumber(1.0,2.0)) === LogNumber(1.0,2.0)
    @test imag(LogNumber(1.0,2.0)) == 0
    @test conj(LogNumber(1.0,2.0)) === LogNumber(1.0,2.0)

    @test zero(LogNumber{Float64}) === LogNumber(0.0,0.0)
    @test one(LogNumber{Float64}) === LogNumber(0.0,1.0)
    @test float(LogNumber(1,2)) === LogNumber(1.0,2.0)
    @test promote_type(Float64, LogNumber{Float64}) === LogNumber{Float64}
    @test promote_type(ComplexF64, LogNumber{Float64}) === Complex{LogNumber{Float64}}

    # `s*log ε → -∞`, so a larger log part is a smaller number
    @test LogNumber(1.0,0.0) < 5.0
    @test LogNumber(1.0,0.0) < LogNumber(0.0,-1E6)
    @test signbit(LogNumber(1.0,0.0))
    @test !signbit(LogNumber(-1.0,0.0))
    @test isinf(LogNumber(1.0,0.0))
    @test isfinite(LogNumber(0.0,1.0))

    # a product of two genuine log numbers is not representable
    @test LogNumber(0.0,3.0) * LogNumber(1.0,2.0) === LogNumber(3.0,6.0)
    @test LogNumber(1.0,2.0) * LogNumber(0.0,3.0) === LogNumber(3.0,6.0)
    @test_throws ArgumentError LogNumber(1.0,2.0) * LogNumber(1.0,2.0)
end

@testset "Complex{LogNumber}" begin
    l = LogNumber(2im, im+1)
    @test l isa Complex{<:LogNumber}
    @test real(l) === LogNumber(0,1)
    @test imag(l) === LogNumber(2,1)
    @test conj(l) === LogNumber(-2im, 1-im)
    @test logpart(l) == 2im
    @test realpart(l) == im+1
    @test l(ℯ) ≈ 2im + im + 1

    # log of a complex power number is a complex log number
    @test log((1+im)*ϵ) isa Complex{<:LogNumber}
    @test log((1+im)*ϵ) == LogNumber(1, log(1+im))
    @test exp(log((1+im)*ϵ)) ≈ (1+im)*ϵ

    # dividing by a complex scalar must not run Base's complex division over the parts
    @test LogNumber(1.0,2.0)/(2im) === LogNumber(1.0/(2im), 2.0/(2im))
    @test (l/(2im))(0.5) ≈ l(0.5)/(2im)

    # a rounding-level real part is weighed against the whole part, not on its own
    @test LogNumber(1, π*im + eps()^2) ≈ LogNumber(1, π*im)
end

@testset "PowerNumber and LogNumber mix" begin
    # only the ϵ^0 coefficient of a power number reaches a log number
    @test (2+ϵ) * LogNumber(1.0,2.0) === 2.0 * LogNumber(1.0,2.0)
    @test (2+ϵ) + LogNumber(1.0,2.0) === LogNumber(1.0,4.0)
    @test ϵ * LogNumber(1.0,2.0) === zero(LogNumber{Float64})   # ϵ*log ϵ → 0
    @test ϵ + LogNumber(1.0,2.0) === LogNumber(1.0,2.0)
    @test_throws DomainError inv(ϵ) * LogNumber(1.0,2.0)

    # and that is what they promote to, so Base code that promotes first agrees
    @test promote_type(PowerNumber{Float64,Float64,Float64}, LogNumber{Float64}) === LogNumber{Float64}
    @test convert(LogNumber{Float64}, 2+ϵ) === LogNumber(0.0,2.0)
    @test muladd(2+ϵ, LogNumber(1.0,2.0), LogNumber(0.0,1.0)) === LogNumber(2.0,5.0)
end


#0.19999999999999996, 0.10000000000000009, 1.3, 1.0 + (-1.0)ϵ^1.0 + o(ϵ^1.0)
#(HypergeometricFunctions._₂F₁general)(0.10000000000000009, 0.19999999999999996, 1.3, 1.0 + (-1.0)ϵ^1.0 + o(ϵ^1.0))
#(HypergeometricFunctions._₂F₁Inf)(0.10000000000000009, 1.1, 1.3, (-1.0)ϵ^-1.0 + (1.0)ϵ^0.0 + o(ϵ^0.0))
#(HypergeometricFunctions.AInf)(0.10000000000000009, 1.1, 1.3, (-1.0)ϵ^1.0 + (-1.0)ϵ^2.0 + o(ϵ^2.0), 1, 0.0)
#z = (-1.0)ϵ^1.0 + (-1.0)ϵ^2.0 + o(ϵ^2.0)
# z = PowerNumber(-1.,-1.,1.,2.)
# @enter (HypergeometricFunctions.BInf)(0.10000000000000009, 1.1, 1.3, z, 1, 0.0)
# @enter (HypergeometricFunctions.recInfβ₀)(0.10000000000000009, 1.1, 1.3, z, 1, 0.0)
