for ii ∈ 1:Threads.nthreads()
    rm("$ii.log")
end

Threads.@threads for ii ∈ 1:100
    cmd =   @cmd "julia my-test.jl"
    log =   "$(Threads.threadid()).log"
    (run ∘ pipeline)(cmd; stdout=log, append=true)
end
