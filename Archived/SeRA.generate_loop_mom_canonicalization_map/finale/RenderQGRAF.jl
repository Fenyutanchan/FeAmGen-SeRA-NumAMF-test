import  Graphs
import	SymEngine
import  YAML

function eliminate_degree_two_vertex(
	in_incidence_mat::Matrix{<:Integer}
)::Matrix{<:Integer}
	incidence_mat	=	deepcopy(in_incidence_mat)
	while true
		edge_list	=	[
			[findfirst(!iszero, col), findlast(!iszero, col)]
			for col ∈ eachcol(incidence_mat)
		]
	
		degree_two_vertices	=	findall(row -> sum(row) == 2, eachrow(incidence_mat))
		if isempty(degree_two_vertices)
			break
		end
		vertex_list			=	setdiff(1:size(incidence_mat, 1), degree_two_vertices)
	
		for vertex ∈ degree_two_vertices
			relavant_edge_indices		=	findall(!iszero, incidence_mat[vertex, :])
			@assert length(relavant_edge_indices) == 2
	
			deleteat!(
				edge_list,
				findall(
					edge -> edge ∈ [
						findall(!iszero, incidence_mat[:, first(relavant_edge_indices)]), 
						findall(!iszero, incidence_mat[:, last(relavant_edge_indices)])
					],
					edge_list
				)
			)
			to_be_connected_vertices	=	setdiff(
				union(
					findall(!iszero, incidence_mat[:, first(relavant_edge_indices)]),
					findall(!iszero, incidence_mat[:, last(relavant_edge_indices)])
				),
				vertex
			)
			push!(edge_list, sort(to_be_connected_vertices))
		end
	
		incidence_mat	=	zeros(Int, length(vertex_list), length(edge_list))
		for (edge_index, edge) ∈ enumerate(edge_list)
			first_vertex_index	=	findfirst(v -> v == first(edge), vertex_list)
			last_vertex_index	=	findfirst(v -> v == last(edge), vertex_list)
	
			incidence_mat[first_vertex_index, edge_index]	+= 1
			incidence_mat[last_vertex_index, edge_index]	+= 1
		end
		
		@assert sum(incidence_mat, dims=1) == 2 * ones(Int, 1, length(edge_list))
	end

	return	incidence_mat
end

function convert_simple_incidence_mat(
	incidence_mat::Matrix{<:Integer}
)::Matrix{<:Integer}
	num_vertex	=	size(incidence_mat, 1)

	edge_list	=	[
		(findfirst(!iszero, col), findlast(!iszero, col))
		for col ∈ eachcol(incidence_mat)
	]
	sort!(edge_list)

	unique_edge_list	=	unique(edge_list)

	for (unique_edge_index, unique_edge) ∈ enumerate(unique_edge_list)
		edge_indices	=	findall(edge -> edge == unique_edge, edge_list)

		if isempty(edge_indices)
			continue
		elseif first(unique_edge) == last(unique_edge)
			the_vertex		=	first(unique_edge)
			deleteat!(edge_list, edge_indices)
			for vertex ∈ (num_vertex + 1):2:(num_vertex + 2 * length(edge_indices))
				push!(
					edge_list,
					(the_vertex, vertex),
					(vertex, vertex + 1),
					(the_vertex, vertex + 1)
				)
			end
			num_vertex	+=	2 * length(edge_indices)
		elseif length(edge_indices) ≥ 2
			the_first_vertex, the_second_vertex	=	unique_edge
			deleteat!(edge_list, edge_indices[2:end])
			for vertex ∈ (num_vertex + 1):(num_vertex + length(edge_indices) - 1)
				push!(
					edge_list,
					(the_first_vertex, vertex),
					(the_second_vertex, vertex)
				)
			end
			num_vertex	+=	length(edge_indices) - 1
		else
			continue
		end
	end
	sort!(edge_list)
	
	num_edge	=	length(edge_list)

	out_incidence_mat	=	zeros(Int, num_vertex, num_edge)

	for (edge_index, edge) ∈ enumerate(edge_list)
		out_incidence_mat[first(edge), edge_index]	+=	1
		out_incidence_mat[last(edge), edge_index]	+=	1
	end
	
	return out_incidence_mat
end

function render_Feynman_diagram_yaml(
	input_Dict::Dict{Any, Any}
)::Tuple{
	Matrix{<:Integer},		# incidence matrix
	Vector{SymEngine.Basic}	# momenta list
}

	# key_list	=	keys(input_Dict)
	# @assert union(string.(key_list)) == union(["n_loop", "outgoing_propagators", "sign", "remnant_propagators", "diagram_index", "Diagram", "incoming_propagators", "vertices", "symmetry_factor"])

	incidence_mat	=	zeros(
		Int,
		length(input_Dict["vertices"]),
		length(input_Dict["remnant_propagators"])
	)
	mom_list		=	SymEngine.Basic[]

	for vertex ∈ input_Dict["vertices"]
		for edge ∈ filter(index -> index > 0, vertex["propagator_index_list"])
			incidence_mat[vertex["vertex"], edge]	+= 	1
		end
	end
	
	@assert sum(incidence_mat, dims=1) == 2 * ones(Int, 1, length(input_Dict["remnant_propagators"]))
	for edge ∈ input_Dict["remnant_propagators"]
		if edge["birth_index"] == edge["death_index"]
			@assert incidence_mat[edge["birth_index"], edge["propagator_index"]] == 2
		else
			@assert incidence_mat[edge["birth_index"], edge["propagator_index"]] == 1
			@assert incidence_mat[edge["death_index"], edge["propagator_index"]] == 1
		end
		push!(mom_list, SymEngine.Basic(edge["momentum"]))
	end
	@assert	all(col -> sum(col) == 2, eachcol(incidence_mat))
	return	incidence_mat, mom_list
end

function construct_graph_from_incidence_matrix(
	incidence_mat::Matrix{<:Integer}
)::Graphs.SimpleGraphs.SimpleGraph
	g = Graphs.Graph(size(incidence_mat, 1))

	for col ∈ eachcol(incidence_mat)
		vertex_list = findall(!iszero, col)
		@assert length(vertex_list) == 2
		Graphs.add_edge!(g, vertex_list...)
	end
	
	return g
end

function render_qgraf_yaml(
	file_name::String
)::Tuple{
	Vector{Graphs.SimpleGraphs.SimpleGraph},	# graph list
	Vector{Vector{SymEngine.Basic}}						# momenta list
}
	Feynman_diagram_yaml_list			=	YAML.load_file(file_name)["FeynmanDiagrams"]
	
	incidence_mat_list	=	Matrix{Int}[]
	mom_list_list		=	Vector{SymEngine.Basic}[]
	for Feynman_diagram_yaml ∈ Feynman_diagram_yaml_list
		incidence_mat, mom_list	=	render_Feynman_diagram_yaml(Feynman_diagram_yaml)
		push!(incidence_mat_list, incidence_mat)
		push!(mom_list_list, mom_list)
	end
	 
	incidence_mat_list	=	(convert_simple_incidence_mat ∘ eliminate_degree_two_vertex).(incidence_mat_list)
	unique!(incidence_mat_list)
	unique!.(mom_ -> SymEngine.expand(- mom_), mom_list_list)
	construct_graph_from_incidence_matrix.(incidence_mat_list), mom_list_list
end
