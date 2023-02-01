##########################################################
function my_generate_loop_mom_canonicalization_map_v1( 
    n_loop::Int64, 
    loop_den_list::Vector{Basic} 
)::Dict{Basic,Basic}
##########################################################

  # In the beginning, we suppose normalize_loop_mom has been used.

  @vars q1,q2,q3
  @funs Den
  qi_list = Basic[q1,q2,q3]
  vanish_qi_dict = Dict{Basic,Basic}( q1=>0, q2=>0, q3=>0 )

  @assert all( den_ -> get_name(den_) == "Den", loop_den_list )

  den_mom_list = map( den_ -> get_args(den_)[1], loop_den_list )
  vac_den_mom_list = map( mom_ -> expand( mom_ - subs( mom_, vanish_qi_dict... ) ), den_mom_list )
  unique!( vac_den_mom_list )

  has_single_q1 = q1 ∈ vac_den_mom_list

  if has_single_q1 == false
    first_has_q1_mom = vac_den_mom_list[ findfirst( m_ -> SymEngine.coeff(m_,q1) != 0, vac_den_mom_list ) ]
    q1_coeff = SymEngine.coeff(first_has_q1_mom,q1)
    q1_replace = (1+q1_coeff)*q1-q1_coeff*first_has_q1_mom 
    #------
    qi_list = map( qi_ -> subs( qi_, q1 => q1_replace ), qi_list )
    vac_den_mom_list = map( qi_ -> subs( qi_, q1 => q1_replace ), vac_den_mom_list )
  end # if

  if n_loop ≥ 2
    has_single_q2 = q2 ∈ vac_den_mom_list
    if has_single_q2 == false
      first_has_q2_mom = vac_den_mom_list[ findfirst( m_ -> SymEngine.coeff(m_,q2) != 0, vac_den_mom_list ) ]
      q2_coeff = SymEngine.coeff(first_has_q2_mom,q2)
      q2_replace = (1+q2_coeff)*q2-q2_coeff*first_has_q2_mom 
      #------
      qi_list = map( qi_ -> subs( qi_, q2 => q2_replace ), qi_list )
      vac_den_mom_list = map( qi_ -> subs( qi_, q2 => q2_replace ), vac_den_mom_list )
    end # if
  end # if

  if n_loop == 3 
    has_single_q3 = q3 in vac_den_mom_list
    if has_single_q3 == false
      first_has_q3_mom = vac_den_mom_list[ findfirst( m_ -> SymEngine.coeff(m_,q3) != 0, vac_den_mom_list ) ]
      q3_coeff = SymEngine.coeff(first_has_q3_mom,q3)
      q3_replace = (1+q3_coeff)*q3-q3_coeff*first_has_q3_mom 
      #------
      qi_list = map( qi_ -> subs( qi_, q3 => q3_replace ), qi_list )
      vac_den_mom_list = map( qi_ -> subs( qi_, q3 => q3_replace ), vac_den_mom_list )
    end # if
  end # if
 
  has_both_q1_q2_pos = findfirst( m_ -> SymEngine.coeff(m_,q1) * SymEngine.coeff(m_,q2) != 0, vac_den_mom_list )
  if has_both_q1_q2_pos ≠ nothing
    has_both_q1_q2_mom = vac_den_mom_list[ has_both_q1_q2_pos ]
    q1_q2_diff_sign = SymEngine.coeff(has_both_q1_q2_mom,q1)*SymEngine.coeff(has_both_q1_q2_mom,q2) == (-1)
    #------
    if q1_q2_diff_sign == true 
      qi_list = map( qi_ -> subs( qi_, q2 => -q2 ), qi_list )
      vac_den_mom_list = map( qi_ -> subs( qi_, q2 => -q2 ), vac_den_mom_list )
    end # if
  end # if

  if n_loop == 3
    has_both_q1_q3_pos = findfirst( m_ -> SymEngine.coeff(m_,q1)*SymEngine.coeff(m_,q3) != 0, vac_den_mom_list )
    has_both_q2_q3_pos = findfirst( m_ -> SymEngine.coeff(m_,q2)*SymEngine.coeff(m_,q3) != 0, vac_den_mom_list )
    has_both_q1_q3_mom = has_both_q1_q3_pos == nothing ? nothing : vac_den_mom_list[ has_both_q1_q3_pos ]
    has_both_q2_q3_mom = has_both_q2_q3_pos == nothing ? nothing : vac_den_mom_list[ has_both_q2_q3_pos ]
    if has_both_q1_q3_mom ≠ nothing && has_both_q2_q3_mom == nothing
      q1_q3_diff_sign = SymEngine.coeff(has_both_q1_q3_mom,q1)*SymEngine.coeff(has_both_q1_q3_mom,q3) == (-1)
      if q1_q3_diff_sign == true
        qi_list = map( qi_ -> subs( qi_, q3 => -q3 ), qi_list )
        vac_den_mom_list = map( qi_ -> subs( qi_, q3 => -q3 ), vac_den_mom_list )
      end # if
    elseif has_both_q1_q3_mom == nothing && has_both_q2_q3_mom ≠ nothing
      q2_q3_diff_sign = SymEngine.coeff(has_both_q2_q3_mom,q2)*SymEngine.coeff(has_both_q2_q3_mom,q3) == (-1)
      if q2_q3_diff_sign == true
        qi_list = map( qi_ -> subs( qi_, q3 => -q3 ), qi_list )
        vac_den_mom_list = map( qi_ -> subs( qi_, q3 => -q3 ), vac_den_mom_list )
      end # if
    elseif has_both_q1_q3_mom ≠ nothing && has_both_q2_q3_mom ≠ nothing
      q1_q3_diff_sign = SymEngine.coeff(has_both_q1_q3_mom,q1)*SymEngine.coeff(has_both_q1_q3_mom,q3) == (-1)
      q2_q3_diff_sign = SymEngine.coeff(has_both_q2_q3_mom,q2)*SymEngine.coeff(has_both_q2_q3_mom,q3) == (-1)
      if q1_q3_diff_sign == true && q2_q3_diff_sign == false
        qi_list = map( qi_ -> subs( qi_, q1 => -q1 ), qi_list )
        vac_den_mom_list = map( qi_ -> subs( qi_, q1 => -q1 ), vac_den_mom_list )
      elseif q1_q3_diff_sign == false && q2_q3_diff_sign == true
        qi_list = map( qi_ -> subs( qi_, q2 => -q2 ), qi_list )
        vac_den_mom_list = map( qi_ -> subs( qi_, q2 => -q2 ), vac_den_mom_list )
      elseif q1_q3_diff_sign == true && q2_q3_diff_sign == true
        qi_list = map( qi_ -> subs( qi_, q3 => -q3 ), qi_list )
        vac_den_mom_list = map( qi_ -> subs( qi_, q3 => -q3 ), vac_den_mom_list )
      else # q1_q3_diff_sign == false && q2_q3_diff_sign == false
        # do nothing for e.g. --θ--O--
      end # if
    else # has_both_q1_q3_pos == nothing && has_both_q2_q3_pos == nothing
      # do nothing
    end # if
        
  end # if


  # Correct the sign locally for each loop momentum. 
  for index ∈ 1:length(vac_den_mom_list)
    one_mom = vac_den_mom_list[index]
    coeff_list = map( qi_ -> SymEngine.coeff(one_mom,qi_), Basic[q1,q2,q3] )
    first_nonzero_pos = findfirst( x_->x_!=0, coeff_list )
    first_nonzero_coeff = coeff_list[first_nonzero_pos]
    @assert abs( first_nonzero_coeff ) == 1
    vac_den_mom_list[index] = first_nonzero_coeff*one_mom
  end # for index

  @vars q1x,q2x,q3x
  has_q12 = q1+q2 ∈ vac_den_mom_list
  has_q13 = q1+q3 ∈ vac_den_mom_list
  has_q23 = q2+q3 ∈ vac_den_mom_list
  if has_q23 == true && has_q13 == false && has_q12 == true 
    # q1 <=> q2
    qi_list = map( qi_ -> subs( qi_, q1=>q2x, q2=>q1x ), qi_list )
    vac_den_mom_list = map( qi_ -> subs( qi_, q1=>q2x, q2=>q1x ), vac_den_mom_list )
  elseif has_q23 == true && has_q13 == true && has_q12 == false
    # q1, q2, q3 <=> q2, q3, q1
    qi_list = map( qi_ -> subs( qi_, q1=>q2x, q2=>q3x, q3=>q1x ), qi_list )
    vac_den_mom_list = map( qi_ -> subs( qi_, q1=>q2x, q2=>q3x, q3=>q1x ), vac_den_mom_list )
  elseif has_q13 == true && has_q12 == false && has_q23 == false
    # q2 <=> q3
    qi_list = map( qi_ -> subs( qi_, q2=>q3x, q3=>q2x ), qi_list )
    vac_den_mom_list = map( qi_ -> subs( qi_, q2=>q3x, q3=>q2x ), vac_den_mom_list )
  elseif has_q23 == true && has_q12 == false && has_q13 == false 
    # q1 <=> q3
    qi_list = map( qi_ -> subs( qi_, q1=>q3x, q3=>q1x ), qi_list )
    vac_den_mom_list = map( qi_ -> subs( qi_, q1=>q3x, q3=>q1x ), vac_den_mom_list )
  end # if


  # Now vac_den_mom_list is useless. Finalizing
  qi_list = map( qi_ -> subs( qi_, q1=>q1x, q2=>q2x, q3=>q3x ), qi_list )
  return Dict{Basic,Basic}( q1 => expand(qi_list[1]), q2 => expand(qi_list[2]), q3 => expand(qi_list[3]) )

end # function generate_loop_mom_canonicalization_map
