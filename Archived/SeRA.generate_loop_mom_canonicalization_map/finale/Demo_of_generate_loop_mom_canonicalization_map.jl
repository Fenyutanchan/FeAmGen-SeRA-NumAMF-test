### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ e235112a-a392-11ed-2004-4d1f86b09e75
using Combinatorics, Graphs, GraphPlot, LinearAlgebra, SparseArrays, SymEngine, Test, YAML

# ╔═╡ cbc8faef-fd94-4229-86f3-50b887912f3c
md"""
# Demo of `generate_loop_mom_canonicalization_map`
"""

# ╔═╡ fff1ef78-35d2-455e-a07e-0781412185d8
md"""
## Main code
"""

# ╔═╡ 0fd330a4-798c-4367-a3af-dee368b95dde
md"""
### Import packages
"""

# ╔═╡ 1ade2b3e-a490-4190-a0d3-b70894ace130
md"""
### Include files
"""

# ╔═╡ 7273a0ba-e781-41b3-8dcc-b7f7e58c1a44
md"""
Function `ingredients` is for reading outer julia scripts.
"""

# ╔═╡ f07ce9d4-a7c1-4cd4-aabc-ab79265a49ce
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(
		m,
		Expr(
			:toplevel,
			:(eval(x) = $(Expr(:core, :eval))($name, x)),
			:(include(x) = $(Expr(:top, :include))($name, x)),
			:(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
			:(include($path))
		)
	)
	m
end

# ╔═╡ 4b8a966f-1d50-4ced-9618-880a668c52db
ArbitraryLoopGeneration	=	ingredients("ArbitraryLoopGeneration.jl")

# ╔═╡ bdc0c135-0940-4b72-949e-7083ba54bdce
RenderQGRAF	=	ingredients("RenderQGRAF.jl")

# ╔═╡ 2c12d1ff-3eb3-4653-9a7e-350d049daa28
md"""
### Useful functions
"""

# ╔═╡ 28adf81e-02e7-460d-bd7a-878ee3d0f371
md"""
For the Feynman diagrams with momenta, directed edges are needed.
We just flip the sign of the last end in any edge.
Then the function `signed_incidence_matrix` is presented.
"""

# ╔═╡ d4aca402-dec7-4ce7-862a-7f6ac6d1dd26
function signed_incidence_matrix(g::Graphs.SimpleGraphs.SimpleGraph)::SparseArrays.SparseMatrixCSC
	incidence_mat	=	incidence_matrix(g)
	for col ∈ eachcol(incidence_mat)
		col[findlast(!iszero, col)]	-=	2
	end
	return	incidence_mat
end

# ╔═╡ e0b343c4-3597-4a08-a2e1-93e579aa0569
md"""
The complicated function `generate_loop_mom_canonicalization` will generate the canonical loop momenta for any simple graph.
"""

# ╔═╡ b4eefbab-ee96-424f-95cd-93920973a273
function generate_loop_mom_canonicalization(
	g::Graphs.SimpleGraphs.SimpleGraph
)::Vector{SymEngine.Basic}
	@assert	Graphs.is_connected(g)
	vertex_list		=	Graphs.vertices(g)
	edge_list		=	Graphs.edges(g)
	n_loop			=	length(cycle_basis(g))
	@assert	n_loop	==	ne(g) - nv(g) + 1

	q_list	=	[Basic("q$ii") for ii ∈ 1:n_loop]

	signed_incidence_mat	=	(Matrix ∘ signed_incidence_matrix)(g)
	n_rank					=	LinearAlgebra.rank(signed_incidence_mat)
	@assert	n_rank == ne(g) - n_loop

	for selected_vertex_indices ∈ Combinatorics.combinations(1:nv(g), n_rank)
		if LinearAlgebra.rank(signed_incidence_mat[selected_vertex_indices, :]) < n_rank
			continue
		end
		for selected_edge_indices ∈ Combinatorics.combinations(1:ne(g), n_rank)
			if LinearAlgebra.rank(signed_incidence_mat[:, selected_edge_indices]) < n_rank
				continue
			end
			free_edge_indices	=	setdiff(1:ne(g), selected_edge_indices)
			mom_list			=	zeros(SymEngine.Basic, ne(g))
			for sign_list ∈ Base.product([(1, -1) for _ in q_list]...)
				mom_list[free_edge_indices]		=	sign_list .* q_list
				mom_list[selected_edge_indices]	=	inv(
					SymEngine.Basic.(
						signed_incidence_mat[
							selected_vertex_indices,
							selected_edge_indices
						]
					)
				) * (- signed_incidence_mat[selected_vertex_indices, free_edge_indices] * mom_list[free_edge_indices])

				check_flag	=	true

				mom_list	=	expand.(mom_list)
				for mom ∈ mom_list[selected_edge_indices]
					coeff_list	=	SymEngine.coeff.(mom, q_list)
					check_flag	=	all(the_coeff -> the_coeff ≥ 0, coeff_list) || all(the_coeff -> the_coeff ≤ 0, coeff_list)
					if !check_flag
						break
					end
				end
				
				if check_flag
					@assert iszero(SymEngine.expand.(signed_incidence_mat * mom_list))
					return mom_list
				else
					continue
				end
			end
		end
	end

	
	throw("Wrong!")
end

# ╔═╡ 258e7470-e4e3-4864-bb86-ed342b49ef2a
md"""
If we consider that input momenta list and make the canonicalization map for it, the function `generate_loop_mom_canonicalization_map` may works.
"""

# ╔═╡ 13d01e36-921d-4401-9eaa-fac0bd73b920
function generate_loop_mom_canonicalization_map(
	mom_list::Vector{SymEngine.Basic}
)::Dict{SymEngine.Basic, SymEngine.Basic}
	q_list		=	SymEngine.free_symbols(mom_list)
	filter!(q -> (first ∘ string)(q) == 'q', q_list)
	sort!(q_list, by=q->parse(Int, string(q)[2:end]))
	
	tmp_mom_list	=	mom_list - SymEngine.subs.(mom_list, Ref(Dict(q_ => 0 for q_ ∈ q_list)))
	tmp_mom_list	=	SymEngine.expand.(tmp_mom_list)
	filter!(!iszero, tmp_mom_list)
	unique!(abs, tmp_mom_list)
	sort!(
		tmp_mom_list,
		by=mom_->(findfirst(!iszero, SymEngine.coeff.(mom_, q_list)), length(SymEngine.free_symbols(mom_)))
	)

	if (isempty ∘ setdiff)(q_list, tmp_mom_list)
		check_flag	=	all(
			mom_ -> begin
				coeff_list	=	SymEngine.coeff.(mom_, q_list)
				all(the_coeff -> the_coeff ≥ 0, coeff_list) || all(the_coeff -> the_coeff ≤ 0, coeff_list)
			end,
			tmp_mom_list
		)
		if check_flag
			return Dict{SymEngine.Basic, SymEngine.Basic}()
		end
	end

	for selected_mom_indices ∈ Combinatorics.combinations(eachindex(tmp_mom_list), length(q_list))
		for sign_list ∈ Base.product([(1, -1) for _ in q_list]...)
			coeff_matrix	=	reduce(
				vcat,
				transpose(SymEngine.coeff.(mom_, q_list))
					for mom_ ∈ (sign_list .* tmp_mom_list[selected_mom_indices])
			)
			if (iszero ∘ LinearAlgebra.det)(coeff_matrix)
				break
			end
			replacement_rules	=	Dict(q_list .=> inv(coeff_matrix) * q_list)
			new_mom_list		=	(SymEngine.expand ∘ SymEngine.subs).(tmp_mom_list, Ref(replacement_rules))
			@assert	new_mom_list[selected_mom_indices] == sign_list .* q_list

			check_flag	=	all(
				mom_ -> begin
					coeff_list	=	SymEngine.coeff.(mom_, q_list)
					all(the_coeff -> the_coeff ≥ 0, coeff_list) || all(the_coeff -> the_coeff ≤ 0, coeff_list)
				end,
				new_mom_list
			)
			if check_flag
				return	replacement_rules
			end
		end
	end
end

# ╔═╡ 40427506-7327-4cbc-a45c-5ecedea4bfbd
md"""
## Test
"""

# ╔═╡ 7382e6e5-148f-46a1-86b3-3935c19c4dc5
md"""
### Loop $4$ case
"""

# ╔═╡ 10ca1f2c-16f1-4776-bdc7-2e117543193a
md"""
Construct all topologically non-isomorphic vacuum Feynman diagrams of $4$-loop (also known as $4$ independent cycles in the graph theory).
"""

# ╔═╡ 10214c5c-a891-490c-9400-769d5fef491f
# ╠═╡ show_logs = false
loop_4_graph_list	=	ArbitraryLoopGeneration.construct_n_loop_graph(4)

# ╔═╡ cf2913db-0b92-4cfb-9ddc-d1790fb72a2b
md"""
Filter all irreducible diagrams (without bridges in the graph theory) without tadpoles (without self-loops in the graph theory).
"""

# ╔═╡ 4773b99e-e409-4cf7-ad27-5fcbb3bf7896
# ╠═╡ show_logs = false
irreducible_loop_4_graph_without_tadpoles_list	=	filter(
	g -> isempty(Graphs.bridges(g)) && !ArbitraryLoopGeneration.has_tadpole(g),
	loop_4_graph_list
)

# ╔═╡ 224e1c28-263f-414f-acfc-aa0c930dec09
gplot.(irreducible_loop_4_graph_without_tadpoles_list)

# ╔═╡ 3bb32b80-93b9-42d7-91c9-aac98d066c58
md"""
Generate all canonical momenta for every propagators (edges in the graph theory).
"""

# ╔═╡ f5e9852d-1cd8-4f64-885f-c98625069173
generate_loop_mom_canonicalization.(loop_4_graph_list)

# ╔═╡ 247c4e14-2b49-4f20-8b0c-09e86be6770d
generate_loop_mom_canonicalization.(irreducible_loop_4_graph_without_tadpoles_list)

# ╔═╡ 9b0c0e38-a80e-49e7-969e-a0b238c1eb85
md"""
### Loop $3$ case
"""

# ╔═╡ 8cbceca5-eab6-4162-9b42-09a38083ac79
md"""
Same as the loop-$4$ case.
"""

# ╔═╡ 1406dc85-9493-4272-85cb-b5d29dd7d1b1
# ╠═╡ show_logs = false
loop_3_graph_list	=	ArbitraryLoopGeneration.construct_n_loop_graph(3)

# ╔═╡ 85f0d276-c13a-4002-97f7-760c7fe30161
generate_loop_mom_canonicalization.(loop_3_graph_list)

# ╔═╡ 81c46d9c-6de2-483f-a9bf-a51de72bf889
md"""
### For `qgraf`
"""

# ╔═╡ f4189f3e-fa8f-463a-9e5e-35af3ddddea1
md"""
Test for the `generate_loop_mom_canonicalization_map`.
Now, we need render the `.yaml` file from `qgraf`.
"""

# ╔═╡ 7784a8b0-6a01-43e2-b4f6-c9e1ae2ed8cd
_, mom_list_list	=	RenderQGRAF.render_qgraf_yaml("./qgraf_out.yaml")

# ╔═╡ 02280379-bb70-4e65-bcd0-70f7198b7c54
md"""
Then check the function of `generate_loop_mom_canonicalization_map`.
"""

# ╔═╡ aba9fcd0-78d2-4559-b5c8-feb452820f28
Test.@testset "Checking generate_loop_mom_canonicalization_map with qgraf" begin
	repl_rule_list	=	generate_loop_mom_canonicalization_map.(mom_list_list)
	for (mom_list, repl_rule) ∈ zip(mom_list_list, repl_rule_list)
		new_mom_list	=	(SymEngine.expand ∘ SymEngine.subs).(mom_list, Ref(repl_rule))
		
		q_list			=	SymEngine.free_symbols(new_mom_list)
		filter!(sym -> (first ∘ string)(sym) == 'q', q_list)
		sort!(q_list, by=q->parse(Int, string(q)[2:end]))
		unique!(abs, new_mom_list)

		for (index, mom) ∈ enumerate(new_mom_list)
			if sum(SymEngine.coeff.(mom, q_list)) < 0
				new_mom_list[index]	=	expand(-mom)
			end
		end
		
		@test (isempty ∘ setdiff)(q_list, new_mom_list)
		for mom ∈ new_mom_list
			coeff_list	=	SymEngine.coeff.(mom, q_list)
			@test all(the_coeff -> the_coeff ∈ (0, 1), coeff_list)
		end
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
GraphPlot = "a2cc645c-3eea-5389-862e-a155d0052231"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
YAML = "ddb6d928-2868-570f-bddf-ab3f9cf99eb6"

[compat]
Combinatorics = "~1.0.2"
GraphPlot = "~0.5.2"
Graphs = "~1.8.0"
SymEngine = "~0.8.7"
YAML = "~0.4.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "7a850c538bab26bede878cfd9d1403d6d62856b4"

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
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

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

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "680e733c3a0a9cea9e935c8c2184aea6a63fa0b5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.21"

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
git-tree-sha1 = "946b56b2135c6c10bbb93efad8a78b699b6383ab"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.6"

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
git-tree-sha1 = "cee507162ecbb677450f20058ca83bd559b6b752"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.14"

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
# ╟─cbc8faef-fd94-4229-86f3-50b887912f3c
# ╟─fff1ef78-35d2-455e-a07e-0781412185d8
# ╟─0fd330a4-798c-4367-a3af-dee368b95dde
# ╠═e235112a-a392-11ed-2004-4d1f86b09e75
# ╟─1ade2b3e-a490-4190-a0d3-b70894ace130
# ╟─7273a0ba-e781-41b3-8dcc-b7f7e58c1a44
# ╟─f07ce9d4-a7c1-4cd4-aabc-ab79265a49ce
# ╟─4b8a966f-1d50-4ced-9618-880a668c52db
# ╟─bdc0c135-0940-4b72-949e-7083ba54bdce
# ╟─2c12d1ff-3eb3-4653-9a7e-350d049daa28
# ╟─28adf81e-02e7-460d-bd7a-878ee3d0f371
# ╟─d4aca402-dec7-4ce7-862a-7f6ac6d1dd26
# ╟─e0b343c4-3597-4a08-a2e1-93e579aa0569
# ╠═b4eefbab-ee96-424f-95cd-93920973a273
# ╟─258e7470-e4e3-4864-bb86-ed342b49ef2a
# ╠═13d01e36-921d-4401-9eaa-fac0bd73b920
# ╟─40427506-7327-4cbc-a45c-5ecedea4bfbd
# ╟─7382e6e5-148f-46a1-86b3-3935c19c4dc5
# ╟─10ca1f2c-16f1-4776-bdc7-2e117543193a
# ╠═10214c5c-a891-490c-9400-769d5fef491f
# ╟─cf2913db-0b92-4cfb-9ddc-d1790fb72a2b
# ╟─4773b99e-e409-4cf7-ad27-5fcbb3bf7896
# ╟─224e1c28-263f-414f-acfc-aa0c930dec09
# ╟─3bb32b80-93b9-42d7-91c9-aac98d066c58
# ╠═f5e9852d-1cd8-4f64-885f-c98625069173
# ╠═247c4e14-2b49-4f20-8b0c-09e86be6770d
# ╟─9b0c0e38-a80e-49e7-969e-a0b238c1eb85
# ╟─8cbceca5-eab6-4162-9b42-09a38083ac79
# ╟─1406dc85-9493-4272-85cb-b5d29dd7d1b1
# ╠═85f0d276-c13a-4002-97f7-760c7fe30161
# ╟─81c46d9c-6de2-483f-a9bf-a51de72bf889
# ╟─f4189f3e-fa8f-463a-9e5e-35af3ddddea1
# ╠═7784a8b0-6a01-43e2-b4f6-c9e1ae2ed8cd
# ╟─02280379-bb70-4e65-bcd0-70f7198b7c54
# ╠═aba9fcd0-78d2-4559-b5c8-feb452820f28
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
