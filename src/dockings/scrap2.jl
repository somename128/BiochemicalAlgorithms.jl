using Base.Threads

a = zeros(10)

@threads for i in eachindex(a)
    a[i] = threadid()
end

print(a)








