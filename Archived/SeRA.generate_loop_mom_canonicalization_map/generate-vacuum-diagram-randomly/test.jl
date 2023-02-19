### begin env
import  Pkg
Pkg.activate(
    (dirname ∘ dirname)(@__FILE__)
)
Pkg.instantiate()

using   Graphs
using   LinearAlgebra
using   StatsBase
using   SymEngine
### end env

### begin my scripts
include("disjoint_loop_momenta_partition.jl")
include("generate_denominator_momenta.jl")
### end my scripts

function main() # resolve the ambiguous variables in soft scope

### begin random generation 
# n_loop  =   rand(2:10)
n_loop  =   4

qi_list =   [Basic("q$ii") for ii ∈ 1:n_loop]
q_list  =   deepcopy(qi_list)

_, incidence_mat, _ =   generate_random_connected_graph(n_loop)
vacumm_mom_list     =   generate_denominator_momentum_list(incidence_mat)
@show vacumm_mom_list

filter!(!iszero, vacumm_mom_list)
unique!(vacumm_mom_list)
map!(
    mom_ -> begin
        tmp_index   =   findfirst(!iszero, coeff.(mom_, qi_list))
        if isnothing(tmp_index)
            mom_
        else
            expand(mom_ * coeff(mom_, qi_list[tmp_index]))
        end
    end,
    vacumm_mom_list,
    vacumm_mom_list
) # just the thing `SeRA.normalize_loop_mom()` will do.
unique!(vacumm_mom_list)
sort!(vacumm_mom_list, by=string)
# vacumm_mom_list =   Basic.(
#     [
#         "q1",
#         "q2",
#         "q2 + q4 + q5",
#         "q2 - q3",
#         "q2 - q6",
#         "q3",
#         "q3 + q4",
#         "q4",
#         "q4 + q5",
#         "q5",
#         "q6"
#     ]
# )

# n_den  =   length(vacumm_mom_list)

println("origin momenta list: $vacumm_mom_list")
### end random generation

### begin first part
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
println("To be flipped qᵢ's: $qi_to_be_flipped_list")
### end second part

end # main

main()
