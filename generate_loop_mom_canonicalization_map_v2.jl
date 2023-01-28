function my_generate_loop_mom_canonicalization_map_v2(
    n_loop::Int64, 
    loop_den_list::Vector{Basic}
)::Dict{Basic, Basic}
    # In the beginning, we suppose normalize_loop_mom has been used.

    @vars q1
    @funs Den
    qi_list         =   [Basic("q$ii") for ii ∈ 1:n_loop]
    vanish_qi_dict  =   Dict{Basic, Basic}(qi_list .=> zero(Basic))

    @assert all(den_ -> get_name(den_) == "Den", loop_den_list)

    den_mom_list        =   map(den_ -> (first ∘ get_args)(den_), loop_den_list)
    den_ext_mom_list    =   subs.(den_mom_list, Ref(vanish_qi_dict))
    vac_den_mom_list    =   expand.(den_mom_list - den_ext_mom_list)
    unique!(vac_den_mom_list)

    opposite_sign_pair_list =   Vector{Tuple{Basic, Basic}}()
    same_sign_pair_list     =   Vector{Tuple{Basic, Basic}}()
    for qi_index ∈ 1:(n_loop-1)
        qi  =   qi_list[qi_index]
        for qj_index ∈ (qi_index+1):n_loop
        qj = qi_list[qj_index]

        # find the momentum where both qi and qj have non-zero coefficients
        same_sign_qiqj_mom_list     =   Vector{Basic}()
        opposite_sign_qiqj_mom_list =   Vector{Basic}()
        for mom ∈ vac_den_mom_list 
            qi_coeff    =   SymEngine.coeff(mom, qi)
            if iszero(qi_coeff)
                continue
            end # if
            qj_coeff    =   SymEngine.coeff(mom, qj)
            if iszero(qj_coeff)
                continue
            end # if

            if qi_coeff * qj_coeff > 0
                push!(same_sign_qiqj_mom_list, mom)
            else
                push!(opposite_sign_qiqj_mom_list, mom)
            end # if
        end # for mom 

        @assert isempty(same_sign_qiqj_mom_list) || isempty(opposite_sign_qiqj_mom_list)

        if isempty(same_sign_qiqj_mom_list) && isempty(opposite_sign_qiqj_mom_list)
            continue
        elseif !isempty(same_sign_qiqj_mom_list)
            push!(same_sign_pair_list, (qi, qj))
        elseif !isempty(opposite_sign_qiqj_mom_list)
            push!(opposite_sign_pair_list, (qi, qj))
        end # if

        end # for qj_index
    end # for qi_index

    q1_same_sign_qi_list = Basic[q1]
    while true
        pair_list   =   filter(
            one_pair -> !(isempty ∘ intersect)(q1_same_sign_qi_list, one_pair),
            same_sign_pair_list
        )
        if isempty(pair_list)
            break # none in same_sign_pair_list is related to q1_same_sign_qi_list
        end # if
        union!(q1_same_sign_qi_list, pair_list...)
        setdiff!(same_sign_pair_list, pair_list)
    end # while

    pos_list    =   findall(
        one_pair -> !(isempty ∘ intersect)(one_pair, q1_same_sign_qi_list),
        opposite_sign_pair_list
    )
    q1_opposite_sign_qi_list    =   (
        isempty(pos_list) ?
        Vector{Basic}() : 
        (first ∘ setdiff).(opposite_sign_pair_list[pos_list], Ref(q1_same_sign_qi_list))
    )

    while !isempty(q1_opposite_sign_qi_list)
        pair_list   =   filter(
            one_pair -> !(isempty ∘ intersect)(q1_opposite_sign_qi_list, one_pair),
            same_sign_pair_list
        )
        if isempty(pair_list)
            break # none in same_sign_pair_list is related to q1_opposite_sign_qi_list
        end # if
        union!(q1_opposite_sign_qi_list, pair_list...)
        setdiff!(same_sign_pair_list, pair_list)
    end # while
    
    norm_dict   =   Dict{Basic, Basic}(
        map(qi -> (qi => -qi), q1_opposite_sign_qi_list) 
    )

    return norm_dict
end # function 