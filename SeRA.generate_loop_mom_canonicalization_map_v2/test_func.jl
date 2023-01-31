using SeRA
using SymEngine
using Test


@testset "Test generate_loop_mom_canonicalization_map_v2" begin

  loop_den_list = [ Basic("Den(q1-q3,0,0)"), Basic("Den(q2+q3-k1,mt,0)") ]
  norm_dict = SeRA.generate_loop_mom_canonicalization_map_v2( 3, loop_den_list )
  new_loop_den_list = SeRA.normalize_loop_mom( subs.( loop_den_list, norm_dict... ) )

  new_loop_den_list_bench = [ Basic("Den(q1 + q3, 0, 0)"), Basic("Den(k1 + q2 + q3, mt, 0)") ]
  @test new_loop_den_list == new_loop_den_list_bench

end # @testset
