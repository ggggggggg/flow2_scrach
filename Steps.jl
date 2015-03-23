# module RunningSteps

using RunningVectors

type RunningStep
	func::Symbol
	ins::NTuple{Symbol}
	outs::NTuple{Symbol}
	is1to1::NTuple{Bool}
end
s = RunningStep(:+, (:c,:d), (:a,:b), (false, true))
# a,b=f(c,d) with a d 1to1 is represented as
# func = :f
# outs = (:a,:b)
# ins = (:c, :d)
# is1to1 = [false, true]
# outputs are always 1to1 

insdefined(s::RunningStep) = all(isdefined, s.ins)

function outsfins(e::Expr)
	outs = Symbol[]
	ins = Symbol[]
	ein = copy(e)
	if e.head == :tuple
		for i = 1:length(e.args)-1
			push!(outs, e.args[i])
		end
		e=e.args[end]
	end
	if e.head == :(=)
		push!(outs, e.args[1])
		e=e.args[2]
	else
	 	error("did you really pass an Expr like a,b=f(c,d,e) $e?")
	end
	if e.head == :call
		if length(e.args) >1
			ins = convert(Vector{Symbol},e.args[2:length(e.args)])
		end
	else
	 	error("did you really pass an Expr like a,b=f(c,d,e) $e?")
	end
	f = e.args[1]

	return outs, f, ins
end

function makeincrement(s::RunningStep)
	lastout = s.outs[end]
	ex = :( $(s.func)($(s.ins...)) )
	for j in 1:length(s.ins) # modify running inputs
		if s.is1to1[j]
			ex.args[j+1]=:($(s.ins[j])[i])
		end
	end
	# add output
	ex = :($lastout[i] = $ex)
	# add other outputs if neccesary
	for j in length(s.outs)-1:-1:1
		ex = Expr(:tuple, :( $(s.outs[j])[i] ),ex)
	end
	return ex

end
ex=makeincrement(s)

macro makeincrement(e)
	outs, f, ins = outsfins(e)
	s = RunningStep(f, outs, ins, tuple([rand(Bool) for i in ins]) )
	makeincrement(s)
end

b = RunningVector(Float64)
c = RunningVector(Float64)
resize!(c,1)
d = RunningVector(Float64)
resize!(d,1)
a=1
push!(b,1)
i=1

c[i],d[i] = (+(a,b[i]), a*b[i])
expr = quote c[i],d[i] = (+(a,b[i]), a*b[i]) end

function outputlengths(s::RunningStep)
	[length(o) for eval(:($o)) in s.outs]
end






# export #types/methods to export

#  #module body

# end