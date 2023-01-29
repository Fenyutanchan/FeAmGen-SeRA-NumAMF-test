import  Pkg
Pkg.activate(".")
Pkg.instantiate()

using   SymEngine

@vars   q1, q2, q3, q4, q5
@vars   k1
@vars   m2, m4
@vars   ieta

n_loop  =   5

q_list  =   [Basic("q$ii") for ii ∈ 1:5]

mom_list    =   [
    -k1 + q2 - q4,
    k1 - q1 + q3,
    -q4 + q5,
    -k1 - q2 - q3,
    k1 - q1 + q2
]

disjoint_partition_list =   begin
    indices =   [ii for ii ∈ 1:n_loop]
    output  =   Vector{Int}[]
    while !isempty(indices)
        this_partition  =   [first(indices)]
        setdiff!(indices, this_partition)

        last_length =   0
        this_length =   1
        while true
            for index ∈ this_partition[last_length+1:end]
                mom_indices     =   findall(
                    !iszero, coeff.(mom_list, q_list[index])
                )
                indices_indices =   findall(
                    ii -> any(
                        (!iszero ∘ coeff).(
                            mom_list[mom_indices],
                            q_list[ii]
                        )
                    ),
                    indices
                )
                union!(this_partition, indices[indices_indices])
                setdiff!(indices, this_partition)
            end

            last_length =   this_length
            this_length =   length(this_partition)
            if last_length == this_length
                break
            end
        end

        push!(output, this_partition)
    end
    output
end
