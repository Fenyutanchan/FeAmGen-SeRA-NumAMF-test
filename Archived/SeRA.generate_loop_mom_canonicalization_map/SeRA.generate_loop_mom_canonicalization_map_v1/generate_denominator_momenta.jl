function generate_denominator_momenta(
    n_leg::Int,
    n_loop::Int,
    n_den::Int
)::Tuple{
    Vector{Basic}, # denominator momenta list
    Vector{Int} # q sign list
}

    q_list  =   [Basic("q$ii") for ii ∈ 1:n_loop]
    k_list  =   [Basic("k$ii") for ii ∈ 1:n_leg]

    q_sign_list =   rand([-1, 1], n_loop)
    
    mom_list    =   Vector{Basic}(undef, n_den)
    while true
        for ii ∈ 1:n_den
            q_sum   =   rand([-1, 1]) * sum(
                rand(
                    [zero(Basic), one(Basic)],
                    n_loop
                ) .* q_sign_list .* q_list
            )
            k_sum   =   sum(rand(-1:1, n_leg) .* k_list)

            if iszero(q_sum)
                q_sum   +=  rand(q_list[1:n_loop])
            end

            mom_list[ii]    =   expand(q_sum + k_sum)
        end
        if all(
            qi_ -> any(
                mom_ -> !(iszero ∘ coeff)(mom_, qi_),
                mom_list
            ),
            q_list[1:n_loop]
        )
            break
        end
    end
    sort!(mom_list, by=string)

    return mom_list, q_sign_list
end