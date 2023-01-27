import  Pkg
Pkg.activate(".")
Pkg.instantiate()

using   SeRA
using   SymEngine

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

    mom_list    =   Vector{Basic}(undef, n_den)
    q_sign_list =   rand([-1, 1], n_loop)
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

    loop_den_list   =   Den.(mom_list, m_list, ieta_list)
    (println ∘ my_normalize_loop_mom)(loop_den_list)
end

# norm_dict   =   SeRA.generate_loop_mom_canonicalization_map_v2(3, loop_den_list)
