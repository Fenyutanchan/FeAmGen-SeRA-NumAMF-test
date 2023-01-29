for ii ∈ 1:Threads.nthreads()
    if isfile("$ii.log")
        rm("$ii.log")
    end
end

Threads.@threads for ii ∈ 1:100
    cmd =   @cmd "julia my-test.jl"
    log =   "$(Threads.threadid()).log"
    (run ∘ pipeline)(cmd; stdout=log, stderr=log, append=true)
end
