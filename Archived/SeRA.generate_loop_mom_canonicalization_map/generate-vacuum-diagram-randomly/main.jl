### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 9c6deb19-72f7-4c65-a560-685edb411f8c
using	Combinatorics, Graphs, GraphPlot, LinearAlgebra, StatsBase, SymEngine

# ╔═╡ 052e7e7c-ec0e-456a-8a2e-7479742eade6
calc_vertex_weights(incidence_mat::Matrix{<:Integer})::Vector{<:Integer}	=	Int[
	4 - sum(vertex_row)
	for vertex_row ∈ eachrow(incidence_mat)
]

# ╔═╡ bb8e2540-6a49-4869-988d-924447f10014
function generate_weighted_vertex_list(incidence_mat::Matrix{<:Integer})::Vector{<:Integer}
	weight_list				=	calc_vertex_weights(incidence_mat)
	weighted_vertex_list	=	Int[]
	for vertex_index ∈ eachindex(weight_list)
		append!(
			weighted_vertex_list,
			[vertex_index for ii ∈ 1:weight_list[vertex_index]]
		)
	end
	return	weighted_vertex_list
end

# ╔═╡ ab0d80a3-ef8a-4271-98f3-bb5d6ac09f36
function incidence_to_adjacency(incidence_mat::Matrix{<:Integer})::Matrix{<:Integer}
	num_vertex		=	size(incidence_mat, 1)
	adjacency_mat	=	zeros(Int, num_vertex, num_vertex)
	for edge ∈ eachcol(incidence_mat)
		vertex_indices	=	findall(index -> !iszero(edge[index]), 1:num_vertex)
		@assert sum(edge) == 2 && (length(vertex_indices) == 2 || length(vertex_indices) == 1)

		first_index	=	first(vertex_indices)
		last_index	=	last(vertex_indices)
		
		adjacency_mat[first_index, last_index]	+=	1
		adjacency_mat[last_index, first_index]	+=	1
	end
	@assert sum(adjacency_mat) == sum(incidence_mat)
	@show adjacency_mat

	return	adjacency_mat
end

# ╔═╡ 0ffac9a9-baf0-4166-9a80-3de5a4dae081
function generate_random_incidence_matrix(
	n_loop::Integer
)::Matrix{<:Integer}
	n_edge		=	rand((2 * (n_loop - 1)):(3 * (n_loop - 1)))
	n_vertex	=	n_edge - n_loop + 1
	avg_degree	=	2 * n_edge / n_vertex

	println("loop: $n_loop")
	println("edge: $n_edge")
	println("vertex: $n_vertex")
	println("average degree of graph: $avg_degree")

	incidence_mat	=	zeros(Int, n_vertex, n_edge)
	while true
		for edge_index ∈ 1:n_edge
			weighted_vetex_list		=	generate_weighted_vertex_list(incidence_mat)
			selected_vetex_indices	=	sample(weighted_vetex_list, 2)
			for vertex_index ∈ selected_vetex_indices
				incidence_mat[vertex_index, edge_index]	+=	1
			end
		end

		degree_check_flag			=	all(degree -> degree ∈ [3, 4], sum(incidence_mat, dims=2))
		# connectedness_check_flag	=	check_connectedness(incidence_mat)
		if degree_check_flag
			break
		else
			incidence_mat	=	zero(incidence_mat)
		end
	end
	
	return incidence_mat
end

# ╔═╡ d149b057-12e4-4bb8-86a2-ef3390868a2e
# function generate_random_connected_graph(n_loop::Integer)::SimpleGraph{<:Integer}
function generate_random_connected_graph(n_loop::Integer)::Tuple{
	SimpleGraph{<:Integer},
	Matrix{<:Integer},
	Matrix{<:Integer}
}
	incidence_mat	=	generate_random_incidence_matrix(n_loop)
	adjacency_mat	=	incidence_to_adjacency(incidence_mat)
	the_graph		=	Graph(adjacency_mat)
	while !is_connected(the_graph)
		incidence_mat	=	generate_random_incidence_matrix(n_loop)
		adjacency_mat	=	incidence_to_adjacency(incidence_mat)
		the_graph		=	Graph(adjacency_mat)
	end

	return the_graph, incidence_mat, adjacency_mat
end

# ╔═╡ 7241a867-2b7d-4a36-8b30-4597b38e0ae8
function generate_denominator_momentum_list(incidence_mat::Matrix{<:Integer})::Vector{Basic}

	n_edge		=	size(incidence_mat, 2)
	n_vertex	=	size(incidence_mat, 1)
	n_loop		=	n_edge - n_vertex + 1

	signed_incidence_mat	=	deepcopy(incidence_mat)
	self_loop_edge_indices	=	Int[]
	for edge_index ∈ 1:n_edge
		vertex_index		=	findlast(
			!iszero,
			signed_incidence_mat[:, edge_index]	
		)
		signed_incidence_mat[vertex_index, edge_index]	-=	2
		if iszero(signed_incidence_mat[:, edge_index])
			push!(self_loop_edge_indices, edge_index)
		end
	end
	n_rank	=	rank(signed_incidence_mat)
	@assert	(iszero ∘ sum)(signed_incidence_mat, dims=1)
	@assert	n_rank == n_edge - n_loop

	qi_list 				=	[Basic("q$ii") for ii ∈ 1:n_loop]
	mom_list 				=	zeros(Basic, n_edge)
	# mom_determination_flags =	zeros(Bool, n_edge)

	redundant_vertex_index	=	findlast(
		vertex_index -> n_rank == rank(
			signed_incidence_mat[
				setdiff(1:n_vertex, vertex_index),
				begin:end
			]
		),
		1:n_vertex
	)
	selected_vertex_indices	=	setdiff(1:n_vertex, redundant_vertex_index)
	# @show selected_vertex_indices

	# selected_edge_indices		=	Vector{Int}(undef, n_rank)
	# for comb ∈ combinations(
	# 	setdiff(
	# 		n_edge:-1:1, self_loop_edge_indices
	# 	),
	# 	n_rank
	# )
	# 	# @show comb
	# 	if rank(
	# 		signed_incidence_mat[
	# 			setdiff(1:n_vertex, redundant_vertex_index),
	# 			reverse(comb)
	# 		]
	# 	) == n_rank
	# 		selected_edge_indices	=	reverse(comb)
	# 	end
	# end
	selected_edge_indices		=	Int[]
	tmp_edge_indices			=	setdiff(n_edge:-1:1, self_loop_edge_indices)
	while length(selected_edge_indices) < n_rank
		for edge_index ∈ tmp_edge_indices
			last_matrix	=	signed_incidence_mat[:, selected_edge_indices]
			this_matrix	=	hcat(
				last_matrix,
				signed_incidence_mat[:, edge_index]
			)
			if rank(last_matrix) + 1 == rank(this_matrix)
				push!(selected_edge_indices, edge_index)
			end
		end
	end
	sort!(selected_edge_indices)
	# @show selected_edge_indices
	independent_edge_indices	=	setdiff(1:n_edge, selected_edge_indices)
	
	mom_list[independent_edge_indices]	=	qi_list
	# mom_determination_flags[independet_edge_indices]	.=	true

	selected_incidence_mat	=	Basic.(
		signed_incidence_mat[
			selected_vertex_indices,
			selected_edge_indices
		]
	)
	RHS_vector	=	- signed_incidence_mat[
		selected_vertex_indices,
		independent_edge_indices
	] * qi_list
	
	mom_list[selected_edge_indices]	=	inv(selected_incidence_mat) * RHS_vector

	mom_list	=	expand.(mom_list)
	
	@assert	iszero((expand).(signed_incidence_mat * mom_list))
	
	return	mom_list
end

# ╔═╡ 853fb92e-d77a-48ce-a802-48d3b22ff994
function disjoint_loop_momenta_partition(
	mom_list::Vector{Basic}	# linear combination of any momenta list
)::Vector{Vector{Basic}}	# partition of indices for input mom_list
	q_list			=	free_symbols(mom_list)
	sort!(q_list,  by=(qi -> Meta.parse(string(qi)[2:end])))
	
	rest_mom_list 	=	deepcopy(mom_list)
	mom_list_partition	=	Vector{Basic}[]

	while !isempty(q_list)
		this_q_partition	=	[first(q_list)]
		this_mom_partition	=	Basic[]

		setdiff!(q_list, this_q_partition)

		match_mom_list	=	filter(
			mom_ -> (!iszero ∘ coeff)(mom_, first(this_q_partition)),
			rest_mom_list
		)
		match_q_list	=	filter(
			q_ -> any((!iszero ∘ coeff).(match_mom_list, q_)),
			q_list
		)
		union!(this_q_partition, match_q_list)
		union!(this_mom_partition, match_mom_list)
		setdiff!(q_list, match_q_list)
		setdiff!(rest_mom_list, match_mom_list)

		while !isempty(match_q_list)
			match_mom_list	=	filter(
				mom_ -> any((!iszero ∘ coeff).(mom_, match_q_list)),
				rest_mom_list
			)
			match_q_list 	=	filter(
				q_ -> any((!iszero ∘ coeff).(match_mom_list, q_)),
				q_list
			)
			union!(this_q_partition, match_q_list)
			union!(this_mom_partition, match_mom_list)
			setdiff!(q_list, match_q_list)
			setdiff!(rest_mom_list, match_mom_list)
		end

		# push!(q_indices_partition, this_q_partition)
		push!(mom_list_partition, this_mom_partition)
	end
	@assert all(iszero, rest_mom_list)

	# return	q_indices_partition, mom_indices_partition
	return	mom_list_partition
end	# function


# ╔═╡ 20a6f040-a808-4826-95f9-9ab70b72fb89
function make_same_opposite_sign_pair_list(
	mom_list::Vector{Basic}
)::Tuple{
	Vector{Tuple{Basic, Basic}},	# :same
	Vector{Tuple{Basic, Basic}}		# :opposite
}
	q_list	=	free_symbols(mom_list)
	sort!(q_list, by=(qi -> Meta.parse(string(qi)[2:end])))

	same_sign_pair_list		=	Vector{Tuple{Basic, Basic}}()
    opposite_sign_pair_list	=	Vector{Tuple{Basic, Basic}}()
	for qi_index ∈ eachindex(q_list)
		qi	=	q_list[qi_index]
		for qj ∈ q_list[qi_index+1:end]
			# find the momentum where both qi and qj have non-zero coefficients
			same_sign_qiqj_mom_list		=	Vector{Basic}()
			opposite_sign_qiqj_mom_list	=	Vector{Basic}()
			for mom ∈ mom_list
				qi_coeff	=	coeff(mom, qi)
				iszero(qi_coeff) && continue
				qj_coeff	=	coeff(mom, qj)
				iszero(qj_coeff) && continue

				if qi_coeff * qj_coeff > 0
					push!(same_sign_qiqj_mom_list, mom)
				else
					push!(opposite_sign_qiqj_mom_list, mom)
				end # if
			end # for mom
			# @assert isempty(same_sign_qiqj_mom_list) || isempty(opposite_sign_qiqj_mom_list)

			!isempty(same_sign_qiqj_mom_list) && push!(same_sign_pair_list, (qi, qj))
			!isempty(opposite_sign_qiqj_mom_list) && push!(opposite_sign_pair_list, (qi, qj))
		end # for qj
	end # for qi_index

	return	same_sign_pair_list, opposite_sign_pair_list
	
end

# ╔═╡ dcb7ba6c-18e4-498a-8d31-16242c238aee
function has_loop_momenta_sign_conflict(
	same_sign_pair_list::Vector{Tuple{Basic, Basic}},
	opposite_sign_pair_list::Vector{Tuple{Basic, Basic}}
)::Tuple{
	Bool,
	Union{
		Tuple{Basic, Basic},
		Missing
	}
}
	# for (qi, qj) ∈ combinations(q_list, 2)
	# 	has_qi_qj_mom_list	=	filter(
	# 		mom_ -> !iszero(coeff(mom_, qi) * coeff(mom_, qj)),
	# 		mom_list
	# 	)
	# 	if isempty(has_qi_qj_mom_list)
	# 		continue
	# 	end
		
	# 	relative_sign_list	=	coeff.(has_qi_qj_mom_list, qi) .* coeff.(has_qi_qj_mom_list, qj)
	# 	@assert	abs.(relative_sign_list) == one.(relative_sign_list)
	# 	first_sign	=	first(relative_sign_list)
	# 	for sign ∈ Base.rest(relative_sign_list)
	# 		if first_sign != sign
	# 			return	true, qi, qj
	# 		end
	# 	end
	# end
	# return	false, missing, missing

	if isempty(same_sign_pair_list) && isempty(opposite_sign_pair_list)
		return	false, missing
	end

	q_list	=	union(same_sign_pair_list..., opposite_sign_pair_list...)
	sort!(q_list, by=(qi -> Meta.parse(string(qi)[2:end])))
	
	first_same_sign_qi_list		=	Basic[first(q_list)]
    first_opposite_sign_qi_list	=	Basic[]
	while !isempty(same_sign_pair_list) || !isempty(opposite_sign_pair_list)
		first_same_sign_pair_list   =   filter(
			one_pair -> !(isempty ∘ intersect)(one_pair, first_same_sign_qi_list),
			same_sign_pair_list
		)
        union!(first_same_sign_qi_list, first_same_sign_pair_list...)
        setdiff!(same_sign_pair_list, first_same_sign_pair_list)

		first_opposite_sign_pair_list   =   filter(
			one_pair -> !(isempty ∘ intersect)(one_pair, first_same_sign_qi_list),
			opposite_sign_pair_list
		)
		for one_pair ∈ first_opposite_sign_pair_list
			tmp	=	setdiff(one_pair, first_same_sign_qi_list)
			if !isempty(tmp)
				union!(first_opposite_sign_qi_list, tmp)
			else
				return true, one_pair
			end
		end
		# union!(
		# 	first_opposite_sign_qi_list,
		# 	[
		# 		setdiff(one_pair, first_same_sign_qi_list)
		# 		for one_pair ∈ first_opposite_sign_pair_list
		# 	]...
		# )
		setdiff!(opposite_sign_pair_list, first_opposite_sign_pair_list)

		first_opposite_sign_pair_list   =   filter(
			one_pair -> !(isempty ∘ intersect)(one_pair, first_opposite_sign_qi_list),
			same_sign_pair_list
		)
		union!(first_opposite_sign_qi_list, first_opposite_sign_pair_list...)
		setdiff!(same_sign_pair_list, first_opposite_sign_pair_list)

		first_same_sign_pair_list   =   filter(
			one_pair -> !(isempty ∘ intersect)(one_pair, first_opposite_sign_qi_list),
			opposite_sign_pair_list
		)
		for one_pair ∈ first_same_sign_pair_list
			tmp	=	setdiff(one_pair, first_opposite_sign_qi_list)
			if !isempty(tmp)
				union!(first_same_sign_qi_list, tmp)
			else
				return	true, one_pair
			end
		end
		# union!(
		# 	first_same_sign_qi_list,
		# 	[
		# 		setdiff(one_pair, first_opposite_sign_qi_list)
		# 		for one_pair ∈ first_same_sign_pair_list
		# 	]...
		# )
		setdiff!(opposite_sign_pair_list, first_same_sign_pair_list)
	end # while

	for (qi, qj) ∈ Iterators.product(first_same_sign_qi_list, first_opposite_sign_qi_list)
		if qi == qj
			return	true, (qi, qj)
		end
	end
	
	return	false, missing
end

# ╔═╡ 580e5cb1-a991-4efc-935b-b59216a6e3a1
function main(
	vacuum_momenta_list::Vector{Basic}
)
	mom_list	=	deepcopy(vacuum_momenta_list)
	q_list		=	free_symbols(vacuum_momenta_list)
	sort!(q_list,  by=(qi -> Meta.parse(string(qi)[2:end])))
	for (mom_index, mom_) ∈ enumerate(mom_list)
		coeff_list					=	coeff.(mom_, q_list)
		first_non_zero_coeff_index	=	findfirst(!iszero, coeff_list)
		if !isnothing(first_non_zero_coeff_index)
			first_coeff			=	coeff_list[first_non_zero_coeff_index]
			mom_				*=	inv(first_coeff)
			mom_list[mom_index]	= 	expand(mom_)
		end
	end
	unique!(mom_list)
	filter!(!iszero, mom_list)
	sort!(
		mom_list,
		by=(
			mom_ -> (
				(length ∘ free_symbols)(mom_),
				(first ∘ findmin)((Meta.parse ∘ (str -> str[2:end]) ∘ string).(free_symbols(mom_)))
			)
		)
	)
	@show mom_list
	
	num_edge 	=	size(mom_list, 1)
	num_q		=	length(q_list)

	coeff_mat	=	reduce(hcat, [coeff.(mom_list, q_) for q_ ∈ q_list])
	@show coeff_mat

	target_coeff_mat	=	zero(coeff_mat)
	if !iszero(num_edge - num_q)
		for comb ∈ combinations(
			sort(
				[
					[parse(Int, num_str, base=10) for num_str ∈ str]
					for str ∈ string.(setdiff((2^num_q):-1:1, 2 .^ (0:num_q)), base=2, pad=num_q)
				],
				by=(v -> (sum(v), -v))
			),
			num_edge - num_q
		)
			tmp_coeff_mat	=	vcat(
				I(num_q),
				(transpose ∘ reduce)(hcat, comb)
			)
			@show tmp_coeff_mat
			@show (rank ∘ float)(coeff_mat)
			@show (rank ∘ float ∘ hcat)(coeff_mat, tmp_coeff_mat)
			if (rank ∘ float)(coeff_mat) == (rank ∘ float ∘ hcat)(coeff_mat, tmp_coeff_mat)
				target_coeff_mat	=	tmp_coeff_mat
				break
			end
		end
	else
		target_coeff_mat	=	coeff_mat
	end
	target_coeff_mat
end

# ╔═╡ fe8dc9f0-461e-4aaa-b926-98ed7c8fd9b7
the_graph, incidence_mat, adjacency_mat	=	generate_random_connected_graph(4)

# ╔═╡ 09547c5f-a7ef-4c42-be8a-edc83153e21a
gplot(the_graph, nodelabel=1:nv(the_graph), edgelabel=1:ne(the_graph))

# ╔═╡ 15399315-0046-4142-aa83-c3358b261123
# ╠═╡ disabled = true
#=╠═╡
incidence_mat = [0 1 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 0; 0 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1; 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0; 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 0 1 0 1 1 0 0 0 0 0 0 0 0 0; 1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1; 0 1 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 0 1 0 0 0 0]
  ╠═╡ =#

# ╔═╡ 4ead37cd-2732-47bb-ad0f-2493678857b0
original_vacuum_momenta_list	=	generate_denominator_momentum_list(incidence_mat)

# ╔═╡ 3a650f45-7832-42ac-b27a-f0b2e0b5c432
main(original_vacuum_momenta_list)

# ╔═╡ cfd75835-e569-487a-83e7-c6bbe907e3ad
begin
	test_mom_matrix	=	Basic.(
		[
			"q1 + q2 - q3",
			"q2 - q3 - q4",
			"q1 - q4",
			"q1 + q2",
			"q2 - q4",
			"q3 - q4"
		]
	)
	q_list		=	free_symbols(test_mom_matrix)
	sort!(q_list,  by=(qi -> Meta.parse(string(qi)[2:end])))
	for (mom_index, mom_) ∈ enumerate(test_mom_matrix)
		coeff_list					=	coeff.(mom_, q_list)
		first_non_zero_coeff_index	=	findfirst(!iszero, coeff_list)
		if !isnothing(first_non_zero_coeff_index)
			first_coeff					=	coeff_list[first_non_zero_coeff_index]
			mom_						*=	inv(first_coeff)
			test_mom_matrix[mom_index]	= 	expand(mom_)
		end
	end

	unique!(test_mom_matrix)
	filter!(!iszero, test_mom_matrix)
	sort!(
		test_mom_matrix,
		by=(
			mom_ -> (
				(length ∘ free_symbols)(mom_),
				(first ∘ findmin)((Meta.parse ∘ (str -> str[2:end]) ∘ string).(free_symbols(mom_)))
			)
		)
	)

	num_edge 	=	size(test_mom_matrix, 1)
	num_q		=	length(q_list)

	coeff_mat	=	hcat([coeff.(test_mom_matrix, q_) for q_ ∈ q_list]...)
	@show tmp_mat		=	vcat(coeff_mat, Basic.(I(num_q)))
	
	for qi_index ∈ 1:num_q
		
		q_indices_perm	=	sort(
			qi_index:num_q,
			by=q_index -> findfirst(!iszero, tmp_mat[:, q_index])
		)
		tmp_mat[:, qi_index:num_q]	=	tmp_mat[:, q_indices_perm]
		@assert	!iszero(tmp_mat[qi_index, qi_index])
		tmp_mat[:, qi_index]		/=	tmp_mat[qi_index, qi_index]
		
		for qj_index ∈ (qi_index + 1):num_q
			tmp_factor				=	tmp_mat[qi_index, qj_index] / tmp_mat[qi_index, qi_index]
			tmp_mat[:, qj_index]	-=	tmp_factor * tmp_mat[:, qi_index]
		end
	end

	for qi_index ∈ num_q:-1:1
		base_input	=	tmp_mat[:, qi_index]
	end
	tmp_mat
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[compat]
Combinatorics = "~1.0.2"
GraphPlot = "~0.5.2"
Graphs = "~1.7.4"
StatsBase = "~0.33.21"
SymEngine = "~0.8.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "99a8291ea9b215150fe22bfa2671e4c000f4fa1c"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "844b061c104c408b24537482469400af6075aae4"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.5"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "d853e57661ba3a57abcdaa201f4c9917a93487a2"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.4"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"
version = "6.2.1+2"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5cd479730a0cb01f880eff119e9803c13f214cab"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.5.2"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "ba2d094a88b6b287bd25cfa86f301e7693ffae2f"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.4"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "45b288af6956e67e621c5cbb2d75a261ab58300b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.20"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MPC_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPFR_jll", "Pkg"]
git-tree-sha1 = "9618bed470dcb869f944f4fe4a9e76c4c8bf9a11"
uuid = "2ce0c516-f11f-5db3-98ad-e0e1048fbd70"
version = "1.2.1+0"

[[deps.MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"
version = "4.1.1+1"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "151d91d63d8d6c1a5789ecb7de51547e00480f1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.SymEngine]]
deps = ["Compat", "Libdl", "LinearAlgebra", "RecipesBase", "SpecialFunctions", "SymEngine_jll"]
git-tree-sha1 = "6cf88a0b98c758a36e6e978a41e8a12f6f5cdacc"
uuid = "123dc426-2d89-5057-bbad-38513e3affd8"
version = "0.8.7"

[[deps.SymEngine_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPC_jll", "MPFR_jll", "Pkg"]
git-tree-sha1 = "3cd0f249ae20a0093f839738a2f2c1476d5581fe"
uuid = "3428059b-622b-5399-b16f-d347a77089a4"
version = "0.8.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═9c6deb19-72f7-4c65-a560-685edb411f8c
# ╟─052e7e7c-ec0e-456a-8a2e-7479742eade6
# ╟─bb8e2540-6a49-4869-988d-924447f10014
# ╟─ab0d80a3-ef8a-4271-98f3-bb5d6ac09f36
# ╟─0ffac9a9-baf0-4166-9a80-3de5a4dae081
# ╟─d149b057-12e4-4bb8-86a2-ef3390868a2e
# ╠═7241a867-2b7d-4a36-8b30-4597b38e0ae8
# ╟─853fb92e-d77a-48ce-a802-48d3b22ff994
# ╟─20a6f040-a808-4826-95f9-9ab70b72fb89
# ╟─dcb7ba6c-18e4-498a-8d31-16242c238aee
# ╠═580e5cb1-a991-4efc-935b-b59216a6e3a1
# ╠═fe8dc9f0-461e-4aaa-b926-98ed7c8fd9b7
# ╠═09547c5f-a7ef-4c42-be8a-edc83153e21a
# ╠═15399315-0046-4142-aa83-c3358b261123
# ╠═4ead37cd-2732-47bb-ad0f-2493678857b0
# ╠═3a650f45-7832-42ac-b27a-f0b2e0b5c432
# ╠═cfd75835-e569-487a-83e7-c6bbe907e3ad
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
