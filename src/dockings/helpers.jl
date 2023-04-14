function zero_small!(M, tol)
    for ι in eachindex(M)
        if abs(M[ι]) ≤ tol
            M[ι] = 0
        end
    end
    M
end