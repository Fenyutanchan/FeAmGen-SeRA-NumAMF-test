import  Pkg
Pkg.activate(".")
Pkg.instantiate()

using   SeRA
using   SymEngine

include("normalize_loop_mom.jl")
include("generate_loop_mom_canonicalization_map.jl")

n_den_tot   =   10
n_leg_tot   =   4
n_loop_tot  =   6

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

    q_sign_list                 =   rand([-1, 1], n_loop)

    mom_list        =   Vector{Basic}(undef, n_den)
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
    
    non_trivial_no_intersection_flag    =   if isempty(q1_opposite_sign_qi_list)
        false
    else
        have_q1_same_sign_qi_list       =   [
            (isempty ∘ findall)(!iszero, coeff.(mom, q1_same_sign_qi_list))
            for mom ∈ mom_list
        ]
        have_q1_opposite_sign_qi_list   =   [
            (isempty ∘ findall)(!iszero, coeff.(mom, q1_opposite_sign_qi_list))
            for mom ∈ mom_list
        ]

        all(
            xor.(
                have_q1_same_sign_qi_list,
                have_q1_opposite_sign_qi_list
            )
        )
    end

    # if non_trivial_no_intersection_flag
    #     for ii ∈ 1:n_den
    #         mom     =   mom_list[ii]
    #         k_sum   =   subs(
    #             mom,
    #             Dict{Basic, Basic}(
    #                 [Basic("q$jj") => zero(Basic) for jj ∈ 1:n_loop]
    #             )
    #         )

    #         bench_sign, q_sum_bench =   begin
    #             q_coeff_list    =   coeff.(mom, q_list[1:n_loop])
    #             first_qi_index  =   findfirst(!iszero, q_coeff_list)
    #             @assert !isnothing(first_qi_index)
    #             # first_qi        =   q_list[first_qi_index]
    #             q_coeff_list[first_qi_index], sum(
    #                 q_list[findall(!iszero, q_coeff_list)]
    #             )
    #         end
    #     end
    # else
    #     for ii ∈ 1:n_den
    #         mom     =   mom_list[ii]
    #         k_sum   =   subs(
    #             mom,
    #             Dict{Basic, Basic}(
    #                 [Basic("q$jj") => zero(Basic) for jj ∈ 1:n_loop]
    #             )
    #         )

    #         bench_sign, q_sum_bench =   begin
    #             q_coeff_list    =   coeff.(mom, q_list[1:n_loop])
    #             first_qi_index  =   findfirst(!iszero, q_coeff_list)
    #             @assert !isnothing(first_qi_index)
    #             first_qi        =   q_list[first_qi_index]
    #             (first_qi ∈ q1_same_sign_qi_list ? 1 : -1) * q_coeff_list[first_qi_index], sum(
    #                 q_list[findall(!iszero, q_coeff_list)]
    #             )
    #         end

    #         mom_list_bench[ii]  =   expand(q_sum_bench + bench_sign * k_sum)
    #     end
    # end

    for ii ∈ 1:n_den
        mom     =   mom_list[ii]
        k_sum   =   subs(
            mom,
            Dict{Basic, Basic}(
                [Basic("q$jj") => zero(Basic) for jj ∈ 1:n_loop]
            )
        )

        bench_sign, q_sum_bench =   begin
            q_coeff_list    =   coeff.(mom, q_list[1:n_loop])
            first_qi_index  =   findfirst(!iszero, q_coeff_list)
            @assert !isnothing(first_qi_index)
            first_qi        =   q_list[first_qi_index]
            (
                non_trivial_no_intersection_flag ? 
                q_coeff_list[first_qi_index] :
                (first_qi ∈ q1_same_sign_qi_list ? 1 : -1) * q_coeff_list[first_qi_index]
            ), sum(
                q_list[findall(!iszero, q_coeff_list)]
            )
        end

        mom_list_bench[ii]  =   expand(q_sum_bench + bench_sign * k_sum)
    end

    loop_den_list   =   Den.(mom_list, m_list, ieta_list)
    norm_dict       =   my_generate_loop_mom_canonicalization_map(n_loop, loop_den_list)
    loop_den_list   =   my_normalize_loop_mom(subs.(loop_den_list, Ref(norm_dict)))

    loop_den_list_bench =   Den.(mom_list_bench, m_list, ieta_list)

    try
        @assert loop_den_list == loop_den_list_bench
    catch
        println()
        println(n_loop)
        println(norm_dict)
        println(q1_same_sign_qi_list)
        println(Den.(mom_list, m_list, ieta_list))
        println(loop_den_list)
        println(loop_den_list_bench)
    end
end
