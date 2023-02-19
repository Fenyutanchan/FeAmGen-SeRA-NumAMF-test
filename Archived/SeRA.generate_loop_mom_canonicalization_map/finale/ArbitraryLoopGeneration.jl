import  Combinatorics
import  Graphs

export	construct_n_loop_graph
export	has_tadpole

function construct_n_loop_graph(n_loop::Integer)::Vector{Graphs.SimpleGraphs.SimpleGraph}
	graph_list				=	Graphs.SimpleGraphs.SimpleGraph[]
	num_graph				=   0
	vertex_edge_num_pair	=	[
		(v, n_loop - 1 + v)
		for v in (n_loop - 1):(2 * (n_loop - 1))
	]

	for (num_vertex, num_edge) ∈ vertex_edge_num_pair
		could_be_edge_list	=	Graphs.SimpleGraphs.SimpleEdge[]
		for vi ∈ 1:num_vertex
			for vj ∈ vi:num_vertex
				push!(could_be_edge_list, Graphs.Edge(vi, vj))
			end
		end

		for selected_edge_list ∈ Combinatorics.with_replacement_combinations(could_be_edge_list, num_edge)
			different_edge_list		=	unique(selected_edge_list)
			edge_multiplicity_list	=	[
				(length ∘ filter)(e -> e == edge, selected_edge_list)
				for edge in different_edge_list
			]
			vertex_for_adding	=	num_vertex + 1
			edges_for_adding	=	[]
			vertex_degree_list	=	zeros(Int, num_vertex)
			for (edge_index, edge) ∈ enumerate(different_edge_list)
				this_edge_multiplicity			=	edge_multiplicity_list[edge_index]
				vertex_degree_list[Graphs.src(edge)]	+=	this_edge_multiplicity
				vertex_degree_list[Graphs.dst(edge)]	+=	this_edge_multiplicity
	
				if any(d -> d > 4, vertex_degree_list[[Graphs.src(edge), Graphs.dst(edge)]])
					break
				end
	
				if Graphs.src(edge) == Graphs.dst(edge)
					while this_edge_multiplicity > 0
						push!(
							edges_for_adding,
							Graphs.Edge(Graphs.src(edge), vertex_for_adding)
						)
						push!(
							edges_for_adding,
							Graphs.Edge(Graphs.dst(edge), vertex_for_adding + 1)
						)
						push!(
							edges_for_adding,
							Graphs.Edge(vertex_for_adding, vertex_for_adding + 1)
						)
						
						this_edge_multiplicity	-=	1
						vertex_for_adding		+=	2
					end
				else
					while this_edge_multiplicity > 1
						push!(
							edges_for_adding,
							Graphs.Edge(Graphs.src(edge), vertex_for_adding)
						)
						push!(
							edges_for_adding,
							Graphs.Edge(Graphs.dst(edge), vertex_for_adding)
						)
	
						this_edge_multiplicity	-=	1
						vertex_for_adding		+=	1
					end
				end
			end

			degree_flag	=	all(d -> 3 ≤ d ≤ 4, vertex_degree_list)
			if !degree_flag
				continue
			else
				edge_list	=	deepcopy(different_edge_list)
				filter!(e -> Graphs.src(e) != Graphs.dst(e), edge_list)
				push!(edge_list, edges_for_adding...)
				
				could_be_graph  =   Graphs.Graph(edge_list)
	
				if Graphs.is_connected(could_be_graph)
					graph_index	=	findfirst(
						g -> Graphs.Experimental.has_isomorph(could_be_graph, g),
						graph_list
					)
					if isnothing(graph_index)
						push!(graph_list, could_be_graph)
						num_graph   +=  1
						print("Found $num_graph graphs.\r")
					end
				end
			end
		end
	end
	
	println("There are $num_graph graphs in total!")
	return	graph_list
end

function has_tadpole(g::Graphs.SimpleGraphs.SimpleGraph)::Bool
	degree_two_vertex_indices	=	findall(d -> d == 2, Graphs.degree(g))
	
	if isempty(degree_two_vertex_indices)
		return	false
	else
		degree_two_vertices	=	Graphs.vertices(g)[degree_two_vertex_indices]
		edge_list			=	Graphs.edges(g)
		could_be_edges		=	(vec ∘ collect)(
			Base.product(degree_two_vertices, degree_two_vertices)
		)
		unique!(x -> (x[2], x[1]), could_be_edges)
		for vertex_pair ∈ could_be_edges
			if Graphs.Edge(vertex_pair) ∈ edge_list
				return	true
			end
		end
		return	false
	end
end
