### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 3a22370c-4de7-40f2-a60b-8cea0cba9d0a
using Combinatorics, Graphs, GraphPlot, YAML

# ╔═╡ 7ef6708e-1533-4375-a97b-0a07c0fa92f4
include("./ArbitraryLoopGeneration.jl")

# ╔═╡ c0a4d659-c85d-4081-8dcd-25e52f722281
md"""
# Cross check with `qgraf`
"""

# ╔═╡ e30bbcb7-eb90-4210-9019-a1be783649d0
md"""
## Import packages
"""

# ╔═╡ be59e0fb-646d-4768-8b22-e2ce2696dcc7
md"""
`Combinatorics.jl`, `Graphs.jl`, `GraphPlot.jl` and `YAML.jl` are needed.
`Combinatorics.jl` for generation graphs.
`Graphs.jl` for basic manuplations on the graphs.
`GraphPlot.jl` for showing the graphs.
`YAML.jl` for reading output of `qgraf`.
"""

# ╔═╡ de9e0f94-20cf-41fc-a0fc-3091d32b4a81
md"""
## Include files
"""

# ╔═╡ bb406160-32b9-4554-b2f4-0e1f806e6bf9
md"""
Here, I import the function `construct_n_loop_graph`.
Please check the Pluto notebook `Generations_of_arbitrary_n-loop_vacuum_Feynman_diagrams.jl` for details.
"""

# ╔═╡ e9e81a8b-4ee5-46ad-bdc9-584563690263
md"""
## Useful functions
"""

# ╔═╡ 16d48688-09b1-4158-a0ae-cacbeae47c69
md"""
Function `render_Feynman_diagram_yaml(input_Dict::Dict{Any, Any})::Matrix{<:Integer}` resolves the information from the YAML output of `qgraf`.
Then, it returns the incidence matrix of the corresponding graph without the external legs, i.e., degree-1 vertices and the edges ended with them.
"""

# ╔═╡ 6811496b-952e-4c94-8722-a6359012fdba
function render_Feynman_diagram_yaml(
	input_Dict::Dict{Any, Any}
)::Matrix{<:Integer}
	key_list	=	keys(input_Dict)
	# @assert union(string.(key_list)) == union(["n_loop", "outgoing_propagators", "sign", "remnant_propagators", "diagram_index", "Diagram", "incoming_propagators", "vertices", "symmetry_factor"])

	incidence_mat	=	zeros(
		Int,
		length(input_Dict["vertices"]),
		length(input_Dict["remnant_propagators"])
	)

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
		@assert iszero(
			incidence_mat[
				setdiff(
					1:length(input_Dict["vertices"]),
					[edge["birth_index"], edge["death_index"]]
				),
				edge["propagator_index"]
			]
		)
	end

	

	return	incidence_mat
end

# ╔═╡ a5774171-487c-4821-8a8c-07f906d0f3b5
md"""
Function `eliminated_degree_two_vertex(in_incidence_mat::Matrix{<:Integer})::Matrix{<:Integer}` eliminate the degree-$2$ vertices on the incidence matrices.
"""

# ╔═╡ ed1f8beb-3979-4cc8-b564-8d23c7cb6fbc
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

# ╔═╡ e63f57f9-d191-48c2-b80e-d4409b721a58
md"""
Function `function convert_simple_incidence_mat(incidence_mat::Matrix{<:Integer})::Matrix{<:Integer}` convert graphs with multi-edges and self-loops to simple graphs via adding appropriate number of degree-$2$ vertices.
"""

# ╔═╡ 2fe8a8c9-01e9-4dea-8665-130b90e5ac9b
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

# ╔═╡ eb477a74-b599-4fed-85bb-80578c5fa203
md"""
Function `construct_graph_from_incidence_matrix(incidence_mat::Matrix{<:Integer})::Graphs.SimpleGraphs.SimpleGraph` constructs the graphs from incidence matrix directly
"""

# ╔═╡ a0f8fc6b-99a4-4544-97d0-64c4a2fb9064
function construct_graph_from_incidence_matrix(
	incidence_mat::Matrix{<:Integer}
)::Graphs.SimpleGraphs.SimpleGraph
	g = Graph(size(incidence_mat, 1))

	for col ∈ eachcol(incidence_mat)
		vertex_list = findall(!iszero, col)
		@assert length(vertex_list) == 2
		add_edge!(g, vertex_list...)
	end
	
	return g
end

# ╔═╡ a5f7c3a8-5d4d-44fc-84e2-77095f948cc8
md"""
## Main code
"""

# ╔═╡ 3e7b8101-73c1-4973-84eb-1a7a0d085097
md"""
Generate all irreducible loop-$4$ Feynman vacuum diagrams (loop-$4$ are also known as $4$ independent cycles in the graph theory).
"""

# ╔═╡ b1aba99e-fcf9-4298-91ea-ecc6d0fc200c
# ╠═╡ show_logs = false
irreducible_loop_4_graph_list	=	filter(
	g -> (isempty ∘ Graphs.bridges)(g),
	construct_n_loop_graph(4)
)

# ╔═╡ a21a5559-2669-44f2-9d80-7dfc8ac8d67c
md"""
Render all graphs from `qgraf`.
"""

# ╔═╡ bc112850-e851-4e05-8fd2-c2d8bfe4e18b
graph_list_from_qgraf	=	begin
	Feynman_diagram_yaml_list			=	YAML.load_file("qgraf_out.yaml")["FeynmanDiagrams"]
	incidence_matrix_list_from_qgraf	=	(convert_simple_incidence_mat ∘ eliminate_degree_two_vertex ∘ render_Feynman_diagram_yaml).(Feynman_diagram_yaml_list)
	unique!(incidence_matrix_list_from_qgraf)
	construct_graph_from_incidence_matrix.(incidence_matrix_list_from_qgraf)
end

# ╔═╡ 75ca6d40-192f-424b-824d-a2a3f997b37c
md"""
Find the corresponding graphs in all irreducible loop-$4$.
"""

# ╔═╡ a6165d80-5cfd-417f-90c0-f7cffa3b3434
md"""
Show them.
"""

# ╔═╡ fecad945-8a58-4e21-a764-b357ec03f895
irreducible_graph_indices	=	[
	findfirst(
		g_irreducible -> Graphs.Experimental.has_isomorph(g_qgraf, g_irreducible),
		irreducible_loop_4_graph_list
	)
	for g_qgraf ∈ graph_list_from_qgraf
]

# ╔═╡ 93d172f4-5ac1-44ea-85eb-b2d53ad889ab
[
	(
		gplot(irreducible_loop_4_graph_list[irreducible_graph_indices]),
		gplot(graph_from_qgraf)
	)
	for (irreducible_graph_indices, graph_from_qgraf) ∈ zip(irreducible_graph_indices, graph_list_from_qgraf)
]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
YAML = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"

[compat]
Combinatorics = "~1.0.2"
GraphPlot = "~0.5.2"
Graphs = "~1.7.4"
YAML = "~0.4.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "70b875d38dfe6cc6bff68bcd807b6c0cee86eda4"

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
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "d853e57661ba3a57abcdaa201f4c9917a93487a2"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.4"

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

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "18f84637e00b72ba6769034a4b50d79ee40c84a9"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.5"

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

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "129703d62117c374c4f2db6d13a027741c46eafd"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.13"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringEncodings]]
deps = ["Libiconv_jll"]
git-tree-sha1 = "33c0da881af3248dafefb939a21694b97cfece76"
uuid = "69024149-9ee7-55f6-a4c4-859efe599b68"
version = "0.3.6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.YAML]]
deps = ["Base64", "Dates", "Printf", "StringEncodings"]
git-tree-sha1 = "dbc7f1c0012a69486af79c8bcdb31be820670ba2"
uuid = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"
version = "0.4.8"

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
# ╟─c0a4d659-c85d-4081-8dcd-25e52f722281
# ╟─e30bbcb7-eb90-4210-9019-a1be783649d0
# ╟─be59e0fb-646d-4768-8b22-e2ce2696dcc7
# ╠═3a22370c-4de7-40f2-a60b-8cea0cba9d0a
# ╟─de9e0f94-20cf-41fc-a0fc-3091d32b4a81
# ╟─bb406160-32b9-4554-b2f4-0e1f806e6bf9
# ╟─7ef6708e-1533-4375-a97b-0a07c0fa92f4
# ╟─e9e81a8b-4ee5-46ad-bdc9-584563690263
# ╟─16d48688-09b1-4158-a0ae-cacbeae47c69
# ╟─6811496b-952e-4c94-8722-a6359012fdba
# ╟─a5774171-487c-4821-8a8c-07f906d0f3b5
# ╟─ed1f8beb-3979-4cc8-b564-8d23c7cb6fbc
# ╟─e63f57f9-d191-48c2-b80e-d4409b721a58
# ╟─2fe8a8c9-01e9-4dea-8665-130b90e5ac9b
# ╟─eb477a74-b599-4fed-85bb-80578c5fa203
# ╟─a0f8fc6b-99a4-4544-97d0-64c4a2fb9064
# ╟─a5f7c3a8-5d4d-44fc-84e2-77095f948cc8
# ╟─3e7b8101-73c1-4973-84eb-1a7a0d085097
# ╟─b1aba99e-fcf9-4298-91ea-ecc6d0fc200c
# ╟─a21a5559-2669-44f2-9d80-7dfc8ac8d67c
# ╟─bc112850-e851-4e05-8fd2-c2d8bfe4e18b
# ╟─75ca6d40-192f-424b-824d-a2a3f997b37c
# ╟─a6165d80-5cfd-417f-90c0-f7cffa3b3434
# ╟─fecad945-8a58-4e21-a764-b357ec03f895
# ╟─93d172f4-5ac1-44ea-85eb-b2d53ad889ab
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
