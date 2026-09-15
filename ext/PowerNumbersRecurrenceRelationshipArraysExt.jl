module PowerNumbersRecurrenceRelationshipArraysExt
using PowerNumbers, RecurrenceRelationshipArrays


# used for singularintegrals w/ powernumbers
function RecurrenceRelationshipArrays.RecurrenceArray(z::PowerNumber, (A,B,C), data::AbstractVector{<:LogNumber})
    @assert z.α == 0 && z.β == 1
    RecurrenceArray(z.A, (A,B,C), data)
end
end