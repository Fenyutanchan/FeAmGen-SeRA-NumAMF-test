### begin env
import  Pkg
Pkg.activate(
    (dirname ∘ dirname)(@__FILE__)
)
Pkg.instantiate()

using   SymEngine
### end env

### begin my scripts
include("generate_denominator_momenta.jl")
### end my scripts

function main() # resolve the ambiguous variables in soft scope

### begin random generation 
n_loop  =   rand(1:3)
n_den   =   10
println("loop: $n_loop")

@vars   q1, q2, q3
q_list  =   [q1, q2, q3][1:n_loop]
vacumm_mom_list =   (first ∘ generate_denominator_momenta)(
    0, n_loop, n_den
)
unique!(vacumm_mom_list)
map!(
    mom_ -> begin
        tmp_index   =   findfirst(!iszero, coeff.(mom_, q_list))
        @assert !isnothing(tmp_index)
        expand(mom_ * coeff(mom_, q_list[tmp_index]))
    end,
    vacumm_mom_list,
    vacumm_mom_list
) # just the thing `SeRA.normalize_loop_mom()` will do.
unique!(vacumm_mom_list)
sort!(vacumm_mom_list, by=string)
println("origin momenta list: $vacumm_mom_list")
### end random generation

### begin first part
has_single_q1   =   q1 ∈ vacumm_mom_list
if !has_single_q1
    first_has_q1_mom    =   vacumm_mom_list[
        findfirst(
            mom_ -> (!iszero ∘ coeff)(mom_, q1),
            vacumm_mom_list
        )
    ]
    q1_coeff    =   coeff(first_has_q1_mom, q1)
    q1_replace  =   (1 + q1_coeff) * q1 - q1_coeff * first_has_q1_mom
    
    q_list          =   (expand ∘ subs).(q_list, q1 => q1_replace)
    vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q1 => q1_replace)
end
println("single q₁: $has_single_q1")
println("q list: $q_list")
println("momenta list: $vacumm_mom_list")

if n_loop ≥ 2
    has_single_q2   =   q2 ∈ vacumm_mom_list
    if !has_single_q2
        first_has_q2_mom    =   vacumm_mom_list[
            findfirst(mom_ -> (!iszero ∘ coeff)(mom_, q2), vacumm_mom_list)
        ]
        q2_coeff    =   coeff(first_has_q2_mom, q2)
        q2_replace  =   (1 + q2_coeff) * q2 - q2_coeff * first_has_q2_mom

        q_list          =   (expand ∘ subs).(q_list, q2 => q2_replace)
        vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q2 => q2_replace)
    end

    println("single q₂: $has_single_q2")
    println("q list: $q_list")
    println("momenta list: $vacumm_mom_list")
end

if n_loop ≥ 3
    has_single_q3   =   q3 ∈ vacumm_mom_list
    if !has_single_q3
        first_has_q3_mom    =   vacumm_mom_list[
            findfirst(mom_ -> (!iszero ∘ coeff)(mom_, q3), vacumm_mom_list)
        ]
        q3_coeff    =   coeff(first_has_q3_mom, q3)
        q3_replace  =   (1 + q3_coeff) * q3 - q3_coeff * first_has_q3_mom

        q_list          =   (expand ∘ subs).(q_list, q3 => q3_replace)
        vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q3 => q3_replace)
    end

    println("single q₃: $has_single_q3")
    println("q list: $q_list")
    println("momenta list: $vacumm_mom_list")
end

# Find a linear transformation to make single qᵢ appearing.
# q_list: how to make original q's from transformed q's.
# vacumm_mom_list: of transformed q's.
### end first part

### begin second part
has_both_q1_q2_pos  =   findfirst(
    mom_ -> !iszero(coeff(mom_, q1) * coeff(mom_, q2)),
    vacumm_mom_list
)
if !isnothing(has_both_q1_q2_pos)
    println("has both q1 and q2 @ $has_both_q1_q2_pos")

    has_both_q1_q2_mom      =   vacumm_mom_list[has_both_q1_q2_pos]
    q1_q2_diff_sign_flag    =   coeff(has_both_q1_q2_mom, q1) * coeff(has_both_q1_q2_mom, q2) == -1
    if q1_q2_diff_sign_flag
        q_list          =   (expand ∘ subs).(q_list, q2 => -q2)
        vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q2 => -q2)
        println("q1 and q2 have different sings.")
        println("q list: $q_list")
        println("momenta list: $vacumm_mom_list")
    end
end
# it seems that only flip q2 when it has different sign with q1.

if n_loop ≥ 3
    has_both_q1_q3_pos  =   findfirst(
        mom_ -> !iszero(
            coeff(mom_, q1) * coeff(mom_, q3)
        ),
        vacumm_mom_list
    )
    has_both_q2_q3_pos  =   findfirst(
        mom_ -> !iszero(
            coeff(mom_, q2) * coeff(mom_, q3)
        ),
        vacumm_mom_list
    )
    has_both_q1_q3_mom  =   isnothing(has_both_q1_q3_pos) ? nothing : vacumm_mom_list[has_both_q1_q3_pos]
    has_both_q2_q3_mom  =   isnothing(has_both_q2_q3_pos) ? nothing : vacumm_mom_list[has_both_q2_q3_pos]
    if !isnothing(has_both_q1_q3_mom) && isnothing(has_both_q2_q3_mom)
        q1_q3_diff_sign_flag    =   coeff(has_both_q1_q3_mom, q1) * coeff(has_both_q1_q3_mom, q3) == -1
        if q1_q3_diff_sign_flag
            q_list          =   (expand ∘ subs).(q_list, q3 => - q3)
            vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q3 => -q3)
        end
    elseif isnothing(has_both_q1_q3_mom) && !isnothing(has_both_q2_q3_mom)
        q2_q3_diff_sign_flag    =   coeff(has_both_q2_q3_mom, q2) * coeff(has_both_q2_q3_mom, q3) == -1
        if q2_q3_diff_sign_flag
            q_list          =   (expand ∘ subs).(q_list, q3 => -q3)
            vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q3 => -q3)
        end
    elseif !isnothing(has_both_q1_q3_mom) && !isnothing(has_both_q2_q3_mom)
        q1_q3_diff_sign_flag    =   coeff(has_both_q1_q3_mom, q1) * coeff(has_both_q1_q3_mom, q3) == -1
        q2_q3_diff_sign_flag    =   coeff(has_both_q2_q3_mom, q2) * coeff(has_both_q2_q3_mom, q3) == -1
        if q1_q3_diff_sign_flag && !q2_q3_diff_sign_flag
            q_list          =   (expand ∘ subs).(q_list, q1 => -q1)
            vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q1 => -q1)
        elseif !q1_q3_diff_sign_flag && q2_q3_diff_sign_flag
            q_list          =   (expand ∘ subs).(q_list, q2 => -q2)
            vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q2 => -q2)
        elseif q1_q3_diff_sign_flag && q2_q3_diff_sign_flag
            q_list          =   (expand ∘ subs).(q_list, q3 => -q3)
            vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q3 => -q3)
        else
            nothing
        end
    end
end
# add the case of three loops.

for ii ∈ eachindex(vacumm_mom_list)
    one_mom             =   vacumm_mom_list[ii]
    coeff_list          =   coeff.(one_mom, [q1, q2, q3])
    first_nonzero_pos   =   findfirst(!iszero, coeff_list)
    first_nonzero_coeff =   coeff_list[first_nonzero_pos]
    @assert abs(first_nonzero_coeff) == 1
    vacumm_mom_list[ii] =   expand(first_nonzero_coeff * one_mom)
end # just SeRA.normalize_loop_mom() does.

println("after the second part:")
println("q list: $q_list")
println("momenta list: $vacumm_mom_list")

# bug generation
# loop: 3
# origin momenta list: Basic[q1, q1 + q3, q1 - q2, q1 - q2 + q3, q2 - q3]
# single q₁: true
# q list: Basic[q1, q2, q3]
# momenta list: Basic[q1, q1 + q3, q1 - q2, q1 - q2 + q3, q2 - q3]
# single q₂: false
# q list: Basic[q1, q1 - q2, q3]
# momenta list: Basic[q1, q1 + q3, q2, q2 + q3, q1 - q2 - q3]
# single q₃: false
# q list: Basic[q1, q1 - q2, -q1 + q3]
# momenta list: Basic[q1, q3, q2, -q1 + q2 + q3, 2*q1 - q2 - q3]
# has both q1 and q2 @ 4
# q1 and q2 have different sings.
# q list: Basic[q1, q1 + q2, -q1 + q3]
# momenta list: Basic[q1, q3, -q2, -q1 - q2 + q3, 2*q1 + q2 - q3]
### end second part

end # main

main()
