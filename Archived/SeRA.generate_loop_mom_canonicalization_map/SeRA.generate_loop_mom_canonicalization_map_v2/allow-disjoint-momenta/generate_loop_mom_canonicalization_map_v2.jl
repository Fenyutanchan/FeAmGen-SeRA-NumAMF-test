function my_generate_loop_mom_canonicalization_map_v2(
    n_loop::Int64, 
    loop_den_list::Vector{Basic}
)::Dict{Basic, Basic}
    # In the beginning, we suppose normalize_loop_mom has been used.

    @funs Den
    qi_list         =   [Basic("q$ii") for ii ∈ 1:n_loop]
    vanish_qi_dict  =   Dict{Basic, Basic}(qi_list .=> zero(Basic))

    @assert all(den_ -> get_name(den_) == "Den", loop_den_list)

    den_mom_list        =   map(den_ -> (first ∘ get_args)(den_), loop_den_list)
    den_ext_mom_list    =   subs.(den_mom_list, Ref(vanish_qi_dict))
    vac_den_mom_list    =   expand.(den_mom_list - den_ext_mom_list)
    unique!(vac_den_mom_list)

    qi_indices_partition, mom_indices_partition =   disjoint_loop_momenta_partition(
        qi_list, vac_den_mom_list
    )

    qi_to_be_flipped_list   =   Vector{Basic}()

    for (qi_partition, mom_partition) ∈ zip(qi_indices_partition, mom_indices_partition)
        this_qi_list    =   qi_list[qi_partition]
        this_mom_list   =   vac_den_mom_list[mom_partition]

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

    norm_dict   =   Dict{Basic, Basic}(
        map(qi -> (qi => -qi), qi_to_be_flipped_list) 
    )

    return norm_dict
end # function 