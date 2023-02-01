### begin env
import  Pkg
Pkg.activate(
    (dirname ∘ dirname)(@__FILE__)
)
Pkg.instantiate()

using   SymEngine
### end env

### begin random generation 
n_loop  =   rand(1:3)
n_den   =   8

@vars   q1, q2, q3
@funs   DenseArray
q_list  =   [q1, q2, q3]
vanish_qi_dict  =   Dict{Basic, Basic}(
    q_list .=> 0
)

vacumm_mom_list =   Vector{Basic}(undef, n_den)
while true
    for ii ∈ 1:n_den
        q_sign_list =   rand([-1, 0, 1], 3)
        while (iszero ∘ sum)(q_sign_list)
            q_sign_list =   rand([-1, 0, 1], 3)
        end
        # println(q_sign_list)
        vacumm_mom_list[ii] =   (expand ∘ sum)(q_sign_list .* q_list)
    end
    if all(
        q_ -> any(
            mom_ -> (!iszero ∘ coeff)(mom_, q_),
            vacumm_mom_list
        ),
        q_list
    )
        break
    end
end
unique!(vacumm_mom_list)
map!(
    mom_ -> begin
        tmp_index   =   findfirst(!iszero, coeff.(mom_, q_list))
        if isnothing(tmp_index)
            mom_
        else
            expand(mom_ * coeff(mom_, q_list[tmp_index]))
        end
    end,
    vacumm_mom_list,
    vacumm_mom_list
)
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
end
println("single q₂: $has_single_q2")
println("q list: $q_list")
println("momenta list: $vacumm_mom_list")

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
end
println("single q₃: $has_single_q3")
println("q list: $q_list")
println("momenta list: $vacumm_mom_list")

# Find a linear transformation to make single qᵢ appearing.
# q_list: how to make original q's from transformed q's.
# vacumm_mom_list: of transformed q's.
### end first part
