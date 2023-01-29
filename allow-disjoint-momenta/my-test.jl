import  Pkg
Pkg.activate(".")
Pkg.instantiate()

using   SymEngine

include("disjoint_loop_momenta_partition.jl")
include("generate_loop_mom_canonicalization_map.jl")
include("generate_loop_mom_canonicalization_map_v2.jl")
include("normalize_loop_mom.jl")

n_den_tot   =   10
n_leg_tot   =   10
n_loop_tot  =   10

q_str_list  =   ["q$ii" for ii ∈ 1:n_loop_tot]
k_str_list  =   ["k$ii" for ii ∈ 1:n_leg_tot]
m_str_list  =   ["m$ii" for ii ∈ 1:n_den_tot]
(eval ∘ Meta.parse)(
    "@vars " * join(
        vcat(
            q_str_list,
            k_str_list,
            m_str_list
        ),
        ", "
    )
)
@funs   Den

q_list  =   [Basic(q_str) for q_str ∈ q_str_list]
k_list  =   [Basic(k_str) for k_str ∈ k_str_list]

for n_den ∈ 1:n_den_tot, n_leg ∈ 1:n_leg_tot, n_loop ∈ 1:n_loop_tot
    m_list      =   [
        rand([Basic(m_str), zero(Basic)])
        for m_str ∈ m_str_list[1:n_den]
    ]
    ieta_list   =   rand([Basic("ieta"), zero(Basic)], n_den)
    q_sign_list =   rand([-1, 1], n_loop)

    mom_list    =   Vector{Basic}(undef, n_den)
    while true
        for ii ∈ 1:n_den
            q_sum   =   rand([-1, 1]) * sum(
                rand(
                    [zero(Basic), one(Basic)],
                    n_loop
                ) .* q_sign_list .* q_list[1:n_loop]
            )
            k_sum   =   sum(rand(-1:1, n_leg) .* k_list[1:n_leg])
    
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

    mom_list_bench  =   Vector{Basic}(undef, n_den)

    q1_same_sign_qi_list            =   q_list[
        findall(
            q_sign -> q_sign * first(q_sign_list) > 0,
            q_sign_list
        )
    ]
    q1_opposite_sign_qi_list        =   setdiff(
        q_list[1:n_loop],
        q1_same_sign_qi_list
    )
    
    # non_trivial_no_intersection_flag    =   if isempty(q1_opposite_sign_qi_list)
    #     false
    # else
    #     have_q1_same_sign_qi_list       =   [
    #         (isempty ∘ findall)(!iszero, coeff.(mom, q1_same_sign_qi_list))
    #         for mom ∈ mom_list
    #     ]
    #     have_q1_opposite_sign_qi_list   =   [
    #         (isempty ∘ findall)(!iszero, coeff.(mom, q1_opposite_sign_qi_list))
    #         for mom ∈ mom_list
    #     ]

    #     all(
    #         xor.(
    #             have_q1_same_sign_qi_list,
    #             have_q1_opposite_sign_qi_list
    #         )
    #     )
    # end

    q_indices_partition, mom_indices_partition  =   disjoint_loop_momenta_partition(
        q_list[1:n_loop], mom_list
    )

    for (this_q_partition, this_mom_partition) ∈ zip(q_indices_partition, mom_indices_partition)
        this_q_list     =   q_list[this_q_partition]
        this_mom_list   =   mom_list[this_mom_partition]

        this_base_sign              =   q_sign_list[first(this_q_partition)]
        this_same_sign_qi_list      =   this_q_list[
            findall(
                q_sign -> q_sign * this_base_sign > 0,
                q_sign_list[this_q_partition]
            )
        ]
        this_opposite_sign_qi_list  =   setdiff(
            this_q_list,
            this_same_sign_qi_list
        )

        for this_mom_index ∈ eachindex(this_mom_list)
            mom     =   this_mom_list[this_mom_index]
            k_sum   =   subs(
                mom,
                Dict{Basic, Basic}(
                    [qi_ => zero(Basic) for qi_ ∈ this_q_list]
                )
            )

            bench_sign, q_sum_bench =   begin
                q_coeff_list        =   coeff.(mom, this_q_list)
                this_first_qi_index =   findfirst(!iszero, q_coeff_list)
                @assert !isnothing(this_first_qi_index)
                first_qi            =   this_q_list[this_first_qi_index]
                (first_qi ∈ this_same_sign_qi_list ? 1 : -1) * q_coeff_list[this_first_qi_index], sum(
                    this_q_list[findall(!iszero, q_coeff_list)]
                )
            end

            mom_list_bench[this_mom_partition[this_mom_index]]  =   expand(
                q_sum_bench + bench_sign * k_sum
            )
        end
    end

    loop_den_list       =   Den.(mom_list, m_list, ieta_list)
    norm_dict_v1        =   my_generate_loop_mom_canonicalization_map(n_loop, loop_den_list)
    norm_dict_v2        =   my_generate_loop_mom_canonicalization_map_v2(n_loop, loop_den_list)
    loop_den_list_v1    =   my_normalize_loop_mom(subs.(loop_den_list, Ref(norm_dict_v1)))
    loop_den_list_v2    =   my_normalize_loop_mom(subs.(loop_den_list, Ref(norm_dict_v2)))
    loop_den_list_bench =   Den.(mom_list_bench, m_list, ieta_list)

    try
        @assert loop_den_list_v2 == loop_den_list_bench
    catch
        println()
        println(n_loop)
        # println(disjoint_momenta_indices_partition_list)
        println(q_indices_partition)
        println(mom_indices_partition)
        println(norm_dict_v1)
        println(norm_dict_v2)
        println(q1_same_sign_qi_list)
        println(Den.(mom_list, m_list, ieta_list))
        println()
        println(loop_den_list_v1)
        println()
        println(loop_den_list_v2)
        println()
        # println(non_trivial_no_intersection_flag)
        println(loop_den_list_bench)
    end
end
