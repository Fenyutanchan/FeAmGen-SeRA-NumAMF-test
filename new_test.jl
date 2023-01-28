import  Pkg
Pkg.activate(".")
Pkg.instantiate()

# using   SeRA
using   SymEngine

include("normalize_loop_mom.jl")

n_den_tot   =   10
n_leg_tot   =   4
n_loop      =   3

q_str_list  =   ["q$ii" for ii ∈ 1:n_loop]
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
@vars   ieta
@funs   Den

qi_list =   [Basic(q_str) for q_str ∈ q_str_list]
k_list  =   [Basic(k_str) for k_str ∈ k_str_list]

loop_den_list   =   [
    Den(k1 - k2 + k3 - q1 + q2, m1, ieta),
    Den(-k1 - k2 - k3 + q1, 0, 0),
    Den(k3 + q3, 0, 0),
    Den(-k1 - k2 - k3 + q3, m4, ieta),
    Den(-k2 + k3 - q1 + q2, 0, ieta),
    Den(-k1 + k2 - k3 + q3, 0, ieta),
    Den(k1 - k2 + k3 + q2, m7, 0)
]

vanish_qi_dict  =   Dict{Basic, Basic}(qi_list .=> zero(Basic))

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
