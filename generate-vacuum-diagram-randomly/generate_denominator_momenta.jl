calc_vertex_weights(incidence_mat::Matrix{<:Integer})::Vector{<:Integer}	=	Int8[
	4 - sum(vertex_row)
	for vertex_row ∈ eachrow(incidence_mat)
]

function generate_weighted_vertex_list(incidence_mat::Matrix{<:Integer})::Vector{Int8}
	weight_list				=	calc_vertex_weights(incidence_mat)
	weighted_vertex_list	=	Int8[]
	for vertex_index ∈ eachindex(weight_list)
		append!(
			weighted_vertex_list,
			[vertex_index for ii ∈ 1:weight_list[vertex_index]]
		)
	end
	return	weighted_vertex_list
end

function incidence_to_adjacency(incidence_mat::Matrix{<:Integer})::Matrix{Int8}
	num_vertex		=	size(incidence_mat, 1)
	adjacency_mat	=	zeros(Int8, num_vertex, num_vertex)
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

function generate_random_incidence_matrix(
	n_loop::Integer
)::Matrix{Int8}
	n_edge		=	rand((2 * (n_loop - 1)):(3 * (n_loop - 1)))
	n_vertex	=	n_edge - n_loop + 1
	avg_degree	=	2 * n_edge / n_vertex

	println("loop: $n_loop")
	println("edge: $n_edge")
	println("vertex: $n_vertex")
	println("average degree of graph: $avg_degree")

	incidence_mat	=	zeros(Int8, n_vertex, n_edge)
	while true
		for edge_index ∈ 1:n_edge
			weighted_vetex_list		=	generate_weighted_vertex_list(incidence_mat)
			selected_vetex_indices	=	sample(weighted_vetex_list, 2)
			for vertex_index ∈ selected_vetex_indices
				incidence_mat[vertex_index, edge_index]	+=	1
			end
		end

		degree_check_flag			=	all(degree -> degree ∈ [3, 4], sum(incidence_mat, dims=2))
		if degree_check_flag
			break
		else
			incidence_mat	=	zero(incidence_mat)
		end
	end
	
	return incidence_mat
end

function generate_random_connected_graph(n_loop::Integer)::Tuple{
	SimpleGraph{<:Integer},
	Matrix{Int8},
	Matrix{Int8}
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

function generate_denominator_momentum_list(incidence_mat::Matrix{<:Integer})::Vector{Basic}
	signed_incidence_mat	=	deepcopy(incidence_mat)
	
	for col ∈ eachcol(signed_incidence_mat)
		vertex_index		=	findlast(!iszero, col)
		col[vertex_index]	-=	2
		@assert sum(col) == 0
	end
	
	n_edge		=	size(signed_incidence_mat, 2)
	n_vertex	=	size(signed_incidence_mat, 1)
	n_loop		=	n_edge - n_vertex + 1
	n_rank		=	rank(signed_incidence_mat)
	@assert n_rank == n_edge - n_loop

	qi_list 	=	[Basic("q$ii") for ii ∈ 1:n_loop]
	mom_list	=	Vector{Basic}(undef, n_edge)

	self_loop_indices	=	Integer[]

	# Gaussian elimination
	row_index	=	1
	for edge_index ∈ 1:n_edge
		mat_col	=	signed_incidence_mat[:, edge_index]

		vertex_indices	=	findall(!iszero, mat_col)
		if isempty(vertex_indices)
			push!(self_loop_indices, edge_index)
		else
			setdiff!(vertex_indices, 1:(row_index-1))
			if !isempty(vertex_indices)
				row_main_element, the_vertex_index	=	findmax(mat_col[vertex_indices])
				the_vertex_index 					=	vertex_indices[the_vertex_index]
				setdiff!(vertex_indices, the_vertex_index)
				for vertex_index ∈ vertex_indices
					signed_incidence_mat[vertex_index, :]	-=	(mat_col[vertex_index] // row_main_element) * signed_incidence_mat[the_vertex_index, :]
				end
				signed_incidence_mat[row_index, :], signed_incidence_mat[the_vertex_index, :] = signed_incidence_mat[the_vertex_index, :], signed_incidence_mat[row_index, :]
			end
			row_index	+=	1
		end
	end
	# end Gaussian elimination

	mom_list[self_loop_indices]		=	qi_list[1:length(self_loop_indices)]
	num_not_self_loop				=	n_loop - length(self_loop_indices)
	not_self_loop_indices			=	setdiff(1:n_edge, self_loop_indices)
	non_simple_mom_indices			=	Integer[]
	for mom_index ∈ not_self_loop_indices
		if rank(signed_incidence_mat[:, non_simple_mom_indices]) == n_rank
			break
		end
		if rank(signed_incidence_mat[:, non_simple_mom_indices]) != rank(signed_incidence_mat[:, union(non_simple_mom_indices, mom_index)])
			push!(non_simple_mom_indices, mom_index)
		end
	end
	setdiff!(not_self_loop_indices, non_simple_mom_indices)
	mom_list[not_self_loop_indices]	=	qi_list[length(self_loop_indices)+1:end]

	@assert num_not_self_loop == length(not_self_loop_indices)
	@assert	n_rank == length(non_simple_mom_indices)

	signed_incidence_mat	=	Basic.(signed_incidence_mat)
	
	mom_list[non_simple_mom_indices]	=	- inv(
		signed_incidence_mat[1:n_rank, non_simple_mom_indices]
	) * (
		signed_incidence_mat[1:n_rank, not_self_loop_indices] * mom_list[not_self_loop_indices]
	)
	return	expand.(mom_list)
end
