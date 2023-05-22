a = []
lk = ReentrantLock() 
Threads.@threads for i = 1:10
    push!(a,Threads.threadid())
end



