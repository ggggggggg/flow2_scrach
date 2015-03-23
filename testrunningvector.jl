using RunningVectors, Base.Test

r = RunningVector(Float64)
append!(r,[rand() for i=1:4])
@test length(r) == 4
freeuntil!(r,1)
@test length(r) == 4
r[2]=44
@test r[2] == 44
@test eltype(r)==Float64

rb = RunningBitVector()
append!(rb, randbool(10))
@test length(rb) == 10
freeuntil!(rb,1)
@test length(rb) == 10
rb[2] = true
@test rb[2] == true


rsb = RunningSumBitVector()
a = randbool(100)
a[100]=true
append!(rsb,a)
@test sum(a) == sum(rsb)
@test rsb[100]==true
freeuntil!(rsb,50)
@test sum(a) == sum(rsb)
rsb[100]=false
@test sum(rsb) == sum(a)-1
@test available(rsb) == 51:100
@test eltype(rsb)==Bool
@test rsb[100]==false


