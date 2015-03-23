type RunningVector{T<:Number}
	v::Vector{T}
	offset::Int64
end
RunningVector(v::Vector) = RunningVector(v,0)


Base.length(v::RunningVector) = length(v.v)
Base.endof(v::RunningVector) = length(v)
Base.getindex(v::RunningVector, r::Union(AbstractVector, Integer)) = v[r-v.offset]
Base.setindex!(v::RunningVector,x, r::Union(AbstractVector, Integer)) = v[r-v.offset] = x
Base.append!(v::RunningVector, x::Vector) = (append!(v.v, x);v)
Base.start(v::RunningVector) = Base.start(v.v)
Base.next(v::RunningVector, j) = Base.next(v.v,j)
Base.done(v::RunningVector, j) = Base.done(v.v,j)
firstind(v::RunningVector) = v.offset+1
function deluntil!(v::RunningVector, n::Integer)
	v.offset >= n && error("$v has no elements before $n left to delete")
	length(v)>= n || error("length of v is less than n=$n")
	deleteat!(v.v, 1:n-v.offset)
	v.offset = n
	v
end

function Base.show(io::IO, v::RunningVector) 
print(io, "$(length(v))-element $(typeof(v)):\n")
print(io, "offset: $(v.offset)\n")
print(io, v.v)
end

a = RunningVector([1:5])