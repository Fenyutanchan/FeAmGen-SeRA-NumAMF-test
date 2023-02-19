function disjoint_loop_momenta_partition(
    q_list::Vector{Basic},      # loop momenta list
    mom_list::Vector{Basic},    # linear combination of any momenta list
)::Tuple{
    Vector{Vector{Int}},    # partition of indices for input q_list
    Vector{Vector{Int}},    # partition of indices for input mom_list
}
    q_indices_partition     =   Vector{Int}[]
    mom_indices_partition   =   Vector{Int}[]

    q_indices   =   [ii for ii ∈ eachindex(q_list)]
    while !isempty(q_indices)
        this_q_partition    =   [first(q_indices)]
        this_mom_partition  =   Int[]

        setdiff!(q_indices, this_q_partition)

        last_length =   0
        this_length =   1
        while true
            for index ∈ this_q_partition[last_length+1:end]
                match_mom_indices   =   findall(
                    !iszero, coeff.(mom_list, q_list[index])
                )
                match_q_indices     =   filter(
                    ii -> any(
                        (!iszero ∘ coeff).(
                            mom_list[match_mom_indices],
                            q_list[ii]
                        )
                    ),
                    q_indices
                )
                union!(this_q_partition, match_q_indices)
                union!(this_mom_partition, match_mom_indices)
                setdiff!(q_indices, this_q_partition)
            end

            last_length =   this_length
            this_length =   length(this_q_partition)
            if last_length == this_length
                sort!(this_q_partition)
                sort!(this_mom_partition)
                break
            end
        end

        push!(q_indices_partition, this_q_partition)
        push!(mom_indices_partition, this_mom_partition)
    end

    return  q_indices_partition, mom_indices_partition
end # function
