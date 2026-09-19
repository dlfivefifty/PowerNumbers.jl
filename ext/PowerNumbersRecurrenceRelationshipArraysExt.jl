module PowerNumbersRecurrenceRelationshipArraysExt
using PowerNumbers, RecurrenceRelationshipArrays
using PowerNumbers: AnyPowerNumber, apart, alpha, beta


# used for singularintegrals w/ powernumbers
function RecurrenceRelationshipArrays.RecurrenceArray(z::AnyPowerNumber, (A,B,C), data::AbstractVector{<:LogNumber})
    @assert alpha(z) == 0 && beta(z) == 1
    RecurrenceArray(apart(z), (A,B,C), data)
end
end
