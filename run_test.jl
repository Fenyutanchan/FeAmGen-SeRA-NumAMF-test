for ii ∈ 1:Threads.nthreads()
    if isfile("$ii.log")
        rm("$ii.log")
    end
end

Threads.@threads for ii ∈ 1:100
    cmd_allow       =   @cmd "julia allow-disjoint-momenta/my-test.jl"
    cmd_disallow    =   @cmd "julia disallow-disjoint-momenta/my-test.jl"
    log             =   "$(Threads.threadid()).log"
    (run ∘ pipeline).(
        [cmd_allow, cmd_disallow];
        stdout=log, stderr=log, append=true
    )
end
