Test for [FeAmGen.jl](https://github.com/zhaoli-IHEP/FeAmGen.jl.git), [NumAMF.jl](https://github.com/zhaoli-IHEP/NumAMF2.jl.git) and [SeRA.jl](https://github.com/zhaoli-IHEP/SeRA2.jl.git).

Please clone the git repositories of these julia packages manually.

However, `NumAMF.jl` and `SeRA.jl` are not generated as the standard julia package. Please do it yourself, or check my forks as [NumAMF.jl](https://github.com/Fenyutanchan/NumAMF2.jl.git) and `SeRA.jl`[SeRA.jl](https://github.com/Fenyutanchan/SeRA2.jl.git). Notice that the packages of `NumAMF.jl` and `SeRA.jl` are shipped to version 2 but there is no modification to the `package-name.jl/src/pakcage-name.jl`, so you should ensure that the names of these packages are `NumAMF.jl` and `SeRA.jl`, respectively.

Finally, please using `Pkg.instantiate()` to make sure the environment installed correctly.
