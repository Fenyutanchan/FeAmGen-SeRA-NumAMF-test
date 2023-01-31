function my_normalize_loop_mom(
    loop_den_list::Vector{Basic}
)::Vector{Basic}
    # qi_list =   begin
    #     variable_str_list   =   string.(free_symbols(loop_den_list))
    #     filter!(v_str -> first(v_str) == 'q', variable_str_list)
    #     n_loop  =   (first ∘ findmax)(
    #         map(
    #             v_str -> Meta.parse(v_str[2:end]),
    #             variable_str_list
    #         )
    #     )
    #     [Basic("q$ii") for ii ∈ 1:n_loop]
    # end
    qi_list =   begin
        variable_str_list   =   string.(free_symbols(loop_den_list))
        filter!(v_str -> first(v_str) == 'q', variable_str_list)
        index_list  =   map(
            v_str -> Meta.parse(v_str[2:end]),
            variable_str_list
        )
        sort!(index_list)
        qi_list =   [Basic("q$ii") for ii ∈ index_list]
    end
    @funs   Den
    
    n_den   =   length(loop_den_list)
    
    new_loop_den_list   =   Vector{Basic}(undef, n_den)
    
    # Correct the sign locally for each loop momentum. 
    for index ∈ 1:n_den
        one_den =   loop_den_list[index]
        den_mom, den_mass, ieta_term    =   get_args(one_den)
        
        # den_width is supposed to be zero
        @assert ieta_term ∈ [zero(Basic), Basic("ieta")]
        coeff_list          =   map(qi_ -> SymEngine.coeff(den_mom, qi_), qi_list)
        first_nonzero_pos   =   findfirst(x_ -> x_ != 0, coeff_list)
        first_nonzero_coeff =   coeff_list[first_nonzero_pos]
        @assert abs(first_nonzero_coeff) == 1
        new_loop_den_list[index]    =   Den(expand(first_nonzero_coeff * den_mom), den_mass, ieta_term)
    end # for index

    return new_loop_den_list
    
end # function normalize_loop_mom
