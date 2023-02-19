### begin env
import  Pkg
Pkg.activate(
    (dirname ∘ dirname)(@__FILE__)
)
Pkg.instantiate()

using   SymEngine
### end env

### begin my scripts
include("disjoint_loop_momenta_partition.jl")
include("generate_denominator_momenta.jl")
### end my scripts

function main() # resolve the ambiguous variables in soft scope

### begin random generation 
n_loop  =   4
n_den   =   10

# @vars   q1, q2, q3
qi_list =   [Basic("q$ii") for ii ∈ 1:n_loop]
q_list  =   deepcopy(qi_list)
# vanish_qi_dict  =   Dict{Basic, Basic}(
#     q_list .=> 0
# )

# vacumm_mom_list =   Vector{Basic}(undef, n_den)
# while true
#     for ii ∈ 1:n_den
#         qi_sign_list    =   rand([-1, 0, 1], n_loop)
#         while (iszero ∘ sum)(qi_sign_list)
#             qi_sign_list    =   rand([-1, 0, 1], n_loop)
#         end
#         # println(q_sign_list)
#         vacumm_mom_list[ii] =   (expand ∘ sum)(qi_sign_list .* qi_list)
#     end
#     if all(
#         qi_ -> any(
#             mom_ -> (!iszero ∘ coeff)(mom_, qi_),
#             vacumm_mom_list
#         ),
#         qi_list
#     )
#         break
#     end
# end

# vacumm_mom_list =   (first ∘ generate_denominator_momenta)(
#     0, n_loop, n_den
# )
# unique!(vacumm_mom_list)
# map!(
#     mom_ -> begin
#         tmp_index   =   findfirst(!iszero, coeff.(mom_, qi_list))
#         @assert !isnothing(tmp_index)
#         expand(mom_ * coeff(mom_, qi_list[tmp_index]))
#     end,
#     vacumm_mom_list,
#     vacumm_mom_list
# ) # just the thing `SeRA.normalize_loop_mom()` will do.
# unique!(vacumm_mom_list)
# sort!(vacumm_mom_list, by=string)
# println("origin momenta list: $vacumm_mom_list")

vacumm_mom_list =   Basic.(
    [
        "q1",
        "q1 - q2",
        "q1 - q2 - q3",
        "q1 - q2 - q4",
        "q1 - q3",
        "q1 - q4",
        "q2 + q3 + q4",
        "q2 + q4",
        "q4"
    ]
)
println("origin momenta list: $vacumm_mom_list")
### end random generation

### begin first part
# has_single_q1   =   q1 ∈ vacumm_mom_list
# if !has_single_q1
#     first_has_q1_mom    =   vacumm_mom_list[
#         findfirst(
#             mom_ -> (!iszero ∘ coeff)(mom_, q1),
#             vacumm_mom_list
#         )
#     ]
#     q1_coeff    =   coeff(first_has_q1_mom, q1)
#     q1_replace  =   (1 + q1_coeff) * q1 - q1_coeff * first_has_q1_mom
    
#     q_list          =   (expand ∘ subs).(q_list, q1 => q1_replace)
#     vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q1 => q1_replace)
# end
# println("single q₁: $has_single_q1")
# println("q list: $q_list")
# println("momenta list: $vacumm_mom_list")

# if n_loop ≥ 2
#     has_single_q2   =   q2 ∈ vacumm_mom_list
#     if !has_single_q2
#         first_has_q2_mom    =   vacumm_mom_list[
#             findfirst(mom_ -> (!iszero ∘ coeff)(mom_, q2), vacumm_mom_list)
#         ]
#         q2_coeff    =   coeff(first_has_q2_mom, q2)
#         q2_replace  =   (1 + q2_coeff) * q2 - q2_coeff * first_has_q2_mom

#         q_list          =   (expand ∘ subs).(q_list, q2 => q2_replace)
#         vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q2 => q2_replace)
#     end

#     println("single q₂: $has_single_q2")
#     println("q list: $q_list")
#     println("momenta list: $vacumm_mom_list")
# end

# if n_loop ≥ 3
#     has_single_q3   =   q3 ∈ vacumm_mom_list
#     if !has_single_q3
#         first_has_q3_mom    =   vacumm_mom_list[
#             findfirst(mom_ -> (!iszero ∘ coeff)(mom_, q3), vacumm_mom_list)
#         ]
#         q3_coeff    =   coeff(first_has_q3_mom, q3)
#         q3_replace  =   (1 + q3_coeff) * q3 - q3_coeff * first_has_q3_mom

#         q_list          =   (expand ∘ subs).(q_list, q3 => q3_replace)
#         vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, q3 => q3_replace)
#     end

#     println("single q₃: $has_single_q3")
#     println("q list: $q_list")
#     println("momenta list: $vacumm_mom_list")
# end

for ii ∈ eachindex(qi_list)
    qi  =   qi_list[ii]
    has_single_qi   =   qi ∈ vacumm_mom_list
    if !has_single_qi
        first_has_qi_mom    =   vacumm_mom_list[
            findfirst(
                mom_ -> (!iszero ∘ coeff)(mom_, qi),
                vacumm_mom_list
            )
        ]
        qi_coeff    =   coeff(first_has_qi_mom, qi)
        # qi_replace  =   (1 + qi_coeff) * qi - qi_coeff * first_has_qi_mom
        qi_replace  =   (1 + (abs ∘ inv)(qi_coeff)) * qi - inv(qi_coeff) * first_has_qi_mom
        
        q_list          =   (expand ∘ subs).(q_list, qi => qi_replace)
        vacumm_mom_list =   (expand ∘ subs).(vacumm_mom_list, qi => qi_replace)
    end
    println("single q$ii: $has_single_qi")
    println("q list: $q_list")
    println("momenta list: $vacumm_mom_list")
end

# Find a linear transformation to make single qᵢ appearing.
# q_list: how to make original q's from transformed q's.
# vacumm_mom_list: of transformed q's.
### end first part

### begin second part
# the second part is just analyze the relative signs of qᵢ and qⱼ.
# I consider import the v2 algorithm into this.
qi_indices_partition, mom_indices_partition =   disjoint_loop_momenta_partition(qi_list, vacumm_mom_list)
println("qi_indices_partition: $qi_indices_partition")
println("mom_indices_partition: $mom_indices_partition")

qi_to_be_flipped_list   =   Vector{Basic}()

for (qi_partition, mom_partition) ∈ zip(qi_indices_partition, mom_indices_partition)
    this_qi_list    =   qi_list[qi_partition]
    this_mom_list   =   vacumm_mom_list[mom_partition]

    same_sign_pair_list     =   Vector{Tuple{Basic, Basic}}()
    opposite_sign_pair_list =   Vector{Tuple{Basic, Basic}}()
    for qi_index ∈ eachindex(this_qi_list)
        qi  =   this_qi_list[qi_index]
        for qj ∈ this_qi_list[qi_index+1:end]

            # find the momentum where both qi and qj have non-zero coefficients
            same_sign_qiqj_mom_list     =   Vector{Basic}()
            opposite_sign_qiqj_mom_list =   Vector{Basic}()
            for mom ∈ this_mom_list
                qi_coeff    =   coeff(mom, qi)
                iszero(qi_coeff) && continue
                qj_coeff    =   coeff(mom, qj)
                iszero(qj_coeff) && continue

                if qi_coeff * qj_coeff > 0
                    push!(same_sign_qiqj_mom_list, mom)
                else
                    push!(opposite_sign_qiqj_mom_list, mom)
                end # if
            end # for mom
            @assert isempty(same_sign_qiqj_mom_list) || isempty(opposite_sign_qiqj_mom_list)

            !isempty(same_sign_qiqj_mom_list) && push!(same_sign_pair_list, (qi, qj))
            !isempty(opposite_sign_qiqj_mom_list) && push!(opposite_sign_pair_list, (qi, qj))
        end # for qj
    end # for qi_index

    first_same_sign_qi_list     =   Basic[first(this_qi_list)]
    first_opposite_sign_qi_list =   Basic[]
    while !isempty(same_sign_pair_list) || !isempty(opposite_sign_pair_list)
        first_same_sign_pair_list   =   filter(
            one_pair -> !(isempty ∘ intersect)(one_pair, first_same_sign_qi_list),
            same_sign_pair_list
        )
        union!(first_same_sign_qi_list, first_same_sign_pair_list...)
        setdiff!(same_sign_pair_list, first_same_sign_pair_list)

        first_opposite_sign_pair_list   =   filter(
            one_pair -> !(isempty ∘ intersect)(one_pair, first_same_sign_qi_list),
            opposite_sign_pair_list
        )
        union!(
            first_opposite_sign_qi_list,
            [
                (first ∘ setdiff)(one_pair, first_same_sign_qi_list)
                for one_pair ∈ first_opposite_sign_pair_list
            ]
        )
        setdiff!(opposite_sign_pair_list, first_opposite_sign_pair_list)

        first_opposite_sign_pair_list   =   filter(
            one_pair -> !(isempty ∘ intersect)(one_pair, first_opposite_sign_qi_list),
            same_sign_pair_list
        )
        union!(first_opposite_sign_qi_list, first_opposite_sign_pair_list...)
        setdiff!(same_sign_pair_list, first_opposite_sign_pair_list)

        first_same_sign_pair_list   =   filter(
            one_pair -> !(isempty ∘ intersect)(one_pair, first_opposite_sign_qi_list),
            opposite_sign_pair_list
        )
        union!(
            first_same_sign_qi_list,
            [
                (first ∘ setdiff)(one_pair, first_opposite_sign_qi_list)
                for one_pair ∈ first_same_sign_pair_list
            ]
        )
        setdiff!(opposite_sign_pair_list, first_same_sign_pair_list)
    end # while

    union!(qi_to_be_flipped_list, first_opposite_sign_qi_list)
end # for
### end second part

end # main

main()

### Found error: there are confilicts between q1 and q4 signs.
# julia> include("test-generate_loop_mom_canonicalization_map_v1.jl")
#   Activating project at `~/workspace/FeAmGen-SeRA-NumAMF-test`
# origin momenta list: Basic[q1, q1 - q2, q1 - q2 - q3, q1 - q2 - q4, q1 - q3, q1 - q4, q2 + q3 + q4, q2 + q4, q4]
# single q1: true
# q list: Basic[q1, q2, q3, q4]
# momenta list: Basic[q1, q1 - q2, q1 - q2 - q3, q1 - q2 - q4, q1 - q3, q1 - q4, q2 + q3 + q4, q2 + q4, q4]
# single q2: false
# q list: Basic[q1, q1 - q2, q3, q4]
# momenta list: Basic[q1, q2, q2 - q3, q2 - q4, q1 - q3, q1 - q4, q1 - q2 + q3 + q4, q1 - q2 + q4, q4]
# single q3: false
# q list: Basic[q1, q1 - q2, q2 - q3, q4]
# momenta list: Basic[q1, q2, q3, q2 - q4, q1 - q2 + q3, q1 - q4, q1 - q3 + q4, q1 - q2 + q4, q4]
# single q4: true
# q list: Basic[q1, q1 - q2, q2 - q3, q4]
# momenta list: Basic[q1, q2, q3, q2 - q4, q1 - q2 + q3, q1 - q4, q1 - q3 + q4, q1 - q2 + q4, q4]
# qi_indices_partition: [[1, 2, 3, 4]]
# mom_indices_partition: [[1, 2, 3, 4, 5, 6, 7, 8, 9]]

### There is no way for keeping the relative signs between q's.
# julia> include("test-generate_loop_mom_canonicalization_map_v1.jl")
#   Activating project at `~/workspace/FeAmGen-SeRA-NumAMF-test`
# origin momenta list: Basic[q1, q1 - q2, q1 - q2 - q3, q1 - q2 - q4, q1 - q3, q1 - q4, q2 + q3 + q4, q2 + q4, q4]
# single q1: true
# q list: Basic[q1, q2, q3, q4]
# momenta list: Basic[q1, q1 - q2, q1 - q2 - q3, q1 - q2 - q4, q1 - q3, q1 - q4, q2 + q3 + q4, q2 + q4, q4]
# single q2: false
# q list: Basic[q1, q1 + q2, q3, q4]
# momenta list: Basic[q1, -q2, -q2 - q3, -q2 - q4, q1 - q3, q1 - q4, q1 + q2 + q3 + q4, q1 + q2 + q4, q4]
# single q3: false
# q list: Basic[q1, q1 + q2, -q2 + q3, q4]
# momenta list: Basic[q1, -q2, -q3, -q2 - q4, q1 + q2 - q3, q1 - q4, q1 + q3 + q4, q1 + q2 + q4, q4]
# single q4: true
# q list: Basic[q1, q1 + q2, -q2 + q3, q4]
# momenta list: Basic[q1, -q2, -q3, -q2 - q4, q1 + q2 - q3, q1 - q4, q1 + q3 + q4, q1 + q2 + q4, q4]
# qi_indices_partition: [[1, 2, 3, 4]]
# mom_indices_partition: [[1, 2, 3, 4, 5, 6, 7, 8, 9]]
