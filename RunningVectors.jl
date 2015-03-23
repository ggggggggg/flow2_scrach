module RunningVectors
abstract AbstractRunningVector

type RunningVector{T} <: AbstractRunningVector
	d::Vector{T}
	offset::Int
end
type RunningBitVector <: AbstractRunningVector
	d::BitVector
	offset::Int

end
abstract AbstractRunningSumVector <: AbstractRunningVector
type RunningSumVector{T} <: AbstractRunningSumVector
	d::Vector{T}
	offset::Int
	sum::Int
end
type RunningSumBitVector <: AbstractRunningSumVector
	d::BitVector
	offset::Int
	sum::Int
end
RunningSumBitVector() = RunningSumBitVector(BitVector(0),0,0)
RunningVector{T}(::Type{T}) = RunningVector(Array(T,0),0)
RunningBitVector() = RunningBitVector(BitVector(0),0)


Base.getindex(r::AbstractRunningVector, inds) = getindex(r.d, inds-r.offset)
Base.setindex!(r::AbstractRunningVector, v, inds) = setindex!(r.d, v, inds-r.offset)
Base.length(r::AbstractRunningVector) = r.offset + length(r.d)
Base.endof(r::AbstractRunningVector) = length(r)
Base.size(r::AbstractRunningVector) = (length(r),)
Base.resize!(r::AbstractRunningVector,n) = resize!(r.d, n-r.offset)
Base.eltype(r::AbstractRunningVector) = eltype(r.d)
runningendof(x) = typemax(Int)
runningendof(r::AbstractRunningVector) = endof(AbstractRunningVector)
available(r::AbstractRunningVector) = r.offset+1:length(r)

function freeuntil!(r::AbstractRunningVector, i::Int) 
	r.offset >= i && return r
	i > length(r) && error("trying to free past existing data in, length(r) = $(length(r)), i=$i")
	todelete = i-r.offset
	deleteat!(r.d, 1:todelete)
	r.offset = i
	r
end

Base.push!(r::AbstractRunningVector, items...) = push!(r.d, items...)
Base.append!(r::AbstractRunningVector, items) = append!(r.d, items)
function Base.append!(r::AbstractRunningSumVector, items)
	r.sum+=sum(items)
	append!(r.d,items)
end
function Base.push!(r::AbstractRunningSumVector, items...)
	r.sum+=sum(items)
	push!(r.d,items...)
end
function Base.setindex!(r::AbstractRunningSumVector, v, inds)
	oldv = getindex(r.d, inds-r.offset)
	r.sum-=sum(oldv)
	r.sum+=sum(v)
	setindex!(r.d,v, inds-r.offset)
end
Base.sum(r::AbstractRunningSumVector) = r.sum


export RunningVector, RunningBitVector, RunningSumVector, RunningSumBitVector, freeuntil!, runningendof, available,
		AbstractRunningVector, AbstractRunningSumVector
end
