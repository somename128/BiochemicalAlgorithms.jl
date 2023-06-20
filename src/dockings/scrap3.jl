using DataFrames

# generating matrix for initalizing scoring table
R = Matrix3{Float32}([0 0 0; 0 0 0; 0 0 0]) 
R1 = Matrix3{Float32}([1 0 0; 0 1 0; 0 0 1]) 
# initialize scoring table
scoring_table = DataFrame(α=[zero(Float32)], β=[zero(Float32)], γ=[zero(Float32)], R=[R], score=[zero(Float32)])

for i in 1:10
    record = (α=Float32(1.0), β=Float32(1.0), γ=Float32(1.0), R=R1, score=Float32(i))
    push!(scoring_table,record)
end

sort!(scoring_table, [:score], rev=[true])
scoring_table[1:5,:]






