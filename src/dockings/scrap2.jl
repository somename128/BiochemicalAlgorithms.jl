using TypedTables

a = ones(2,2,2)
a[1,2,1] = 9

println(a)

t = Table(α=[], β=[], γ=[], R=[],score=[])

record = generate

