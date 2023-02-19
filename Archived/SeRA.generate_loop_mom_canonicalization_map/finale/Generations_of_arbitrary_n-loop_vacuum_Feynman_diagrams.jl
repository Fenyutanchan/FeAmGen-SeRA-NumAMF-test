### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 7eb4235f-646b-43a7-9427-ed413812e837
using	Combinatorics

# ╔═╡ a419b90c-622b-42ca-bbc7-1347ebe20cdf
using	Graphs, GraphPlot

# ╔═╡ 84e8b6c0-ab94-11ed-12f9-95a699ce219f
md"""
# Generations of arbitrary $n$-loop vacuum Feynman diagrams
"""

# ╔═╡ b2089e09-94e8-44c4-b70b-792431108356
md"""
## Main code
"""

# ╔═╡ 9c216a57-37f3-4e24-ba72-39c779deea13
md"""
### Import packages
"""

# ╔═╡ 32fff9ef-fe8b-43a5-9d8d-2406c8ee0a3e
md"""
Package `Combinatorics.jl` provides the method `Combinatorics.with_replacement_combinations(a, t)`: Generate all combinations with replacement of size t from an array a.
"""

# ╔═╡ c80b1f28-7f63-4a99-ba97-608ed0ec568b
md"""
Packages `Graphs.jl` and `GraphPlot` provide many useful function for manipulation on the graphs.
"""

# ╔═╡ dc0400d5-cf16-41ff-abb6-6abb1003c52d
md"""
### Function defined by myself
"""

# ╔═╡ 87ec2a49-3409-4524-97bb-7f809ce116d5
md"""
Let's construct the list of the number pair $(n_v, n_e)$ for an arbitrary $n$-loop vacuum Feynman diagram, where the diagram consists of $n_v$ vertices and $n_e$ edges.
Notice that the first Betti number for any graph is given as
```math
n_\ell = n_e - n_v + n_c,
```
where the graph has $n_c$ connected components and $n_\ell$ independent cycles.
We just constraints our consideration with $c = 1$.
In addition, there is a theorem says
```math
\sum_{v \in V} d(v) = 2 n_e,
```
where $V$ is the set of the vertices in the graph $G$, and $d(v)$ is the degree of vertex $v$, which means the number of edges end with this vertex.
For Feynman diagrams in the renormalization theory Feynman, we have that $\forall v \in V: d(v) = 3 \text{ or } 4$.
Therefore,
```math
n_\ell - 1 \le n_v \le 2 (n_\ell - 1).
```
Finally, the function are given as
"""

# ╔═╡ fdb4bb75-62a8-4293-b445-3d56dcb44d89
generate_vertex_edge_num_pair_list(n_loop::Integer)	=	[
	(v, n_loop - 1 + v)
	for v in (n_loop - 1):(2 * (n_loop - 1))
]

# ╔═╡ 01c021ff-dbb3-4bf2-ae2e-7765c90c03b7
md"""
Then, We could just enumerate all possible edges for generating all arbitrary $n$-loop vacuum Feynman diagrams.
It seems that the this algorithm is not the best one.
But for the case of $4$-loop, there are $10346677$ possible cases to be check, and it only takes $45~\mathrm{s}$ to run out all of them on one single kernel of my M1 MacBook Pro .
Finally, it shows that only 73 topologically non-isomorphic graphs.
For $5$-loop, $53485609654$ possible cases will be checked.
And for $6$-loop, there are $577217553349739$ possible cases to check.
The general formula to calculate the number of possible cases shows
```math
\sum_{n_v = n_l - 1}^{2 (n_l - 1)} \begin{pmatrix} n_e + n_v - 1 \\ n_v \end{pmatrix} = \sum_{n_v = n_l - 1}^{2 (n_l - 1)} \begin{pmatrix} n_l + 2 (n_v - 1) \\ n_v \end{pmatrix}.
```
"""

# ╔═╡ e77404b1-6d4c-42ad-b8ca-38a172aadb9a
function construct_n_loop_graph(n_loop::Integer)::Vector{Graphs.SimpleGraphs.SimpleGraph}
	graph_list	=	Graphs.SimpleGraphs.SimpleGraph[]
	num_graph   =   0

	for (num_vertex, num_edge) ∈ generate_vertex_edge_num_pair_list(n_loop)
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
				vertex_degree_list[src(edge)]	+=	this_edge_multiplicity
				vertex_degree_list[dst(edge)]	+=	this_edge_multiplicity
	
				if any(d -> d > 4, vertex_degree_list[[src(edge), dst(edge)]])
					break
				end
	
				if src(edge) == dst(edge)
					while this_edge_multiplicity > 0
						push!(
							edges_for_adding,
							Graphs.Edge(src(edge), vertex_for_adding)
						)
						push!(
							edges_for_adding,
							Graphs.Edge(dst(edge), vertex_for_adding + 1)
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
							Graphs.Edge(src(edge), vertex_for_adding)
						)
						push!(
							edges_for_adding,
							Graphs.Edge(dst(edge), vertex_for_adding)
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
				filter!(e -> src(e) != dst(e), edge_list)
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

# ╔═╡ 6a212916-ff3e-4530-9eed-84c7a6f31a61
md"""
Terminology "tadpole" in the Feynman diagrams means the terminology "self-loop" in the graphs.
Because the method to check graph isomorphism between non-simple graphs with multi-edges and self-loops is not implemented in `Graphs.jl` or other graph tools.
I'd like to add some degree-$2$ vertices to make them simple.
For "tadpole" or "self-loop", two degree-$2$ vertices will be added.
So I will find the edge ends with two degree-$2$ vertices for checking if it has "tadpole" or "self-loop".
"""

# ╔═╡ 76b19dfc-9ee1-4014-a2d9-8bdd506ebc75
function has_tadpole(g::Graphs.SimpleGraphs.SimpleGraph)::Bool
	degree_two_vertex_indices	=	findall(d -> d == 2, degree(g))
	
	if isempty(degree_two_vertex_indices)
		return	false
	else
		degree_two_vertices	=	vertices(g)[degree_two_vertex_indices]
		edge_list			=	edges(g)
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

# ╔═╡ 5671babd-29bd-4fca-bbc9-5e253539a2d8
md"""
## Main test
"""

# ╔═╡ c2238331-6e97-4f60-b91d-27b15259274c
md"""
First, $3$-loop cases show as below.
There are $12$ topologically non-isomorphic graphs generated.
"""

# ╔═╡ ac7adda9-846f-4c9f-8dd0-712f8b23d128
loop_3_graph_list	=	construct_n_loop_graph(3)

# ╔═╡ 5da0c20b-98f2-4b36-9c3e-9408fcfb3bef
md"""
Second, $4$-loop cases show as below.
There are $73$ topologically non-isomorphic graphs generated.
"""

# ╔═╡ 009b1d36-6abe-454d-b598-a79578ab0507
loop_4_graph_list	=	construct_n_loop_graph(4)

# ╔═╡ 40baede1-aef0-4d98-a178-55f238eac0fe
md"""
Now, I will find all irreducible and tadpole-free $3$-loop Feynman diagrams.
Terminology "irreducible" in the Feynman diagrams means that the graphs have no bridge edges.
Therefore, there are only $4$ graphs after this filtering.
"""

# ╔═╡ a91a3042-72b1-4a0e-9a55-3e9815034fd2
special_loop_3_graph_list	=	filter(
	g -> (isempty ∘ bridges)(g) && !has_tadpole(g),
	loop_3_graph_list
)

# ╔═╡ ff1587a6-6123-4d48-a117-4b7f0406f8dc
md"""
For irreducible and tadpole-free $4$-loop Feynman diagrams, there are only $15$ graphs.
"""

# ╔═╡ cfc35568-1535-41e2-b661-d89171b55965
special_loop_4_graph_list	=	filter(
	g -> (isempty ∘ bridges)(g) && !has_tadpole(g),
	loop_4_graph_list
)

# ╔═╡ 6d634666-432b-4a09-9074-002f672d55ed
md"""
Plot them:
"""

# ╔═╡ 08bc4501-97e6-40b4-aa0c-2ebaf6350eb9
gplot.(special_loop_3_graph_list)

# ╔═╡ f98b30ad-2410-4926-90d8-de6f64aa471b
gplot.(special_loop_4_graph_list)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"

[compat]
Combinatorics = "~1.0.2"
GraphPlot = "~0.5.2"
Graphs = "~1.7.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "1d7a71c2c96b918ea7c7336582220eb1ed192a8a"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

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

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╟─84e8b6c0-ab94-11ed-12f9-95a699ce219f
# ╟─b2089e09-94e8-44c4-b70b-792431108356
# ╟─9c216a57-37f3-4e24-ba72-39c779deea13
# ╟─32fff9ef-fe8b-43a5-9d8d-2406c8ee0a3e
# ╠═7eb4235f-646b-43a7-9427-ed413812e837
# ╟─c80b1f28-7f63-4a99-ba97-608ed0ec568b
# ╠═a419b90c-622b-42ca-bbc7-1347ebe20cdf
# ╟─dc0400d5-cf16-41ff-abb6-6abb1003c52d
# ╟─87ec2a49-3409-4524-97bb-7f809ce116d5
# ╟─fdb4bb75-62a8-4293-b445-3d56dcb44d89
# ╟─01c021ff-dbb3-4bf2-ae2e-7765c90c03b7
# ╟─e77404b1-6d4c-42ad-b8ca-38a172aadb9a
# ╟─6a212916-ff3e-4530-9eed-84c7a6f31a61
# ╟─76b19dfc-9ee1-4014-a2d9-8bdd506ebc75
# ╟─5671babd-29bd-4fca-bbc9-5e253539a2d8
# ╟─c2238331-6e97-4f60-b91d-27b15259274c
# ╠═ac7adda9-846f-4c9f-8dd0-712f8b23d128
# ╟─5da0c20b-98f2-4b36-9c3e-9408fcfb3bef
# ╠═009b1d36-6abe-454d-b598-a79578ab0507
# ╟─40baede1-aef0-4d98-a178-55f238eac0fe
# ╠═a91a3042-72b1-4a0e-9a55-3e9815034fd2
# ╟─ff1587a6-6123-4d48-a117-4b7f0406f8dc
# ╠═cfc35568-1535-41e2-b661-d89171b55965
# ╟─6d634666-432b-4a09-9074-002f672d55ed
# ╠═08bc4501-97e6-40b4-aa0c-2ebaf6350eb9
# ╠═f98b30ad-2410-4926-90d8-de6f64aa471b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
