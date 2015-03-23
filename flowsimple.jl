using RunningVectors, Graphs, HDF5, JLD, GraphViz
import JLD: JldGroup, JldFile, JldDataset
include("mockpulses.jl")
include("summarize.jl")

const DONETHRU_MAX = typemax(Int)-1

type Histogram
	x::FloatRange{Float64}
	counts::Vector{Int}
	seen::Int # total events seen, including those not in counts
end
Histogram(x::Range) = Histogram(convert(FloatRange{Float64},x),zeros(Int, length(x)-1),0)
bin_edges(h::Histogram) = h.x
bin_centers(h::Histogram) = convert(FloatRange{Float64},range(first(h.x)*step(h.x)/2,step(h.x),length(h.x)-1))
addcounts!(h::Histogram, events::Array) = (h.counts+=hist(events,bin_edges(h))[2];h.seen+=length(events))
addseen!(h::Histogram,n)=h.seen+=n
counts(h::Histogram) = h.counts
misses(h::Histogram) = h.seen-sum(h.counts)
donethru(h::Histogram) = h.seen
inexists(h::Histogram)=true

type Calibration
	a::Float64
	b::Float64
	iscalibrated::Bool
end
donethru(c::Calibration) = c.iscalibrated?DONETHRU_MAX:0
inexists(c::Calibration) = c.iscalibrated

typealias Channel Dict{Symbol,Any}
const perpulse_symbols =  Set{Symbol}()
isperpulse(s::Symbol) = s in perpulse_symbols # TODO: memoize?

filename(c::Channel) = joinpath(pwd(), "flowsimple_test.jld")



function update_histogram(h::Histogram, selection, x)
	addcounts!(h,x[selection])
	addseen!(h,length(selection)-sum(selection))
end

mock_calibration(x) = Calibration(10,100,true)


abstract AbstractStep

type PerPulseStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
	outputs::Vector{Symbol}
end
PerPulseStep(func, inputs, outputs) = PerPulseStep(func, inputs, [is_])
type ToJLDStep <:AbstractStep
	inputs::Vector{Symbol}
end	
type HistogramStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
end
type ThresholdStep <: AbstractStep
	func::Function
	inputs::Vector{Symbol}
	outputs::Vector{Symbol}
	to_watch::Symbol
	to_watch_func::Function
	threshold::Int
	do_if_able::Bool
end
type MockPulsesStep <: AbstractStep
	pg::TwoExponentialPulseGenerator{Int}
	pulses_per_step::Int
	outputs::Vector{Symbol}
	last_rowstamp::Int
end

getfunction(s::AbstractStep) = s.func
graphlabel(s::AbstractStep) = string(getfunction(s))
inputs(s::AbstractStep) = s.inputs
outputs(s::AbstractStep) = s.outputs
inputs(s::AbstractStep, c::Channel) = [c[q] for q in inputs(s)]
outputs(s::AbstractStep, c::Channel) = [c[q] for q in outputs(s)]
inputs(s::AbstractStep, c::Channel, r::UnitRange{Int}) = 
[isperpulse(q)?c[q][r]:c[q] for q in s.inputs]
perpulse_inputs(s::AbstractStep) = inputs(s)[[isperpulse(q) for q in inputs(s)]]
perpulse_outputs(s::AbstractStep) = outputs(s)[[isperpulse(q) for q in outputs(s)]]
perpulse_inputs(s::AbstractStep, c::Channel) = [c[q] for q in perpulse_inputs(s)]
perpulse_outputs(s::AbstractStep, c::Channel) = [c[q] for q in perpulse_outputs(s)]
other_inputs(s::AbstractStep) = inputs(s)[[!isperpulse(q) for q in inputs(s)]]
other_outputs(s::AbstractStep) = outputs(s)[[!isperpulse(q) for q in outputs(s)]]
other_inputs(s::AbstractStep, c::Channel) = [c[q] for q in other_inputs(s)]
other_outputs(s::AbstractStep, c::Channel) = [c[q] for q in other_outputs(s)]
mindonethru(x) = minimum(map(donethru,x))
Base.range(s::AbstractStep, c::Channel) = 1+mindonethru(perpulse_outputs(s,c)):mindonethru(perpulse_inputs(s,c))
other_inputs_exist(s::AbstractStep, c::Channel) = all([haskey(c,q) for q in other_inputs(s)])
# placeholder versions of exists, probably would use Nullable types here
function dostep(s::PerPulseStep,c::Channel)
	f = getfunction(s)
	r = range(s,c)
	other_inputs_exist(s,c) || (return last(r):-1)
	length(r)==0 && (return r)
	fout = f(inputs(s,c,r)...) # it doesn't know the types here, tho f will, also I could write a macro so it knew types
	if length(outputs(s)) == 1 # fout will just be a value
		append!(outputs(s,c)[1],fout)
	elseif length(outputs(s)) > 1 # fout will be a tuple of values, one for each of outputs(s)
		@assert(length(fout) == length(outputs(s)))
		for (fo, co) in zip(fout, outputs(s,c))
			# fo is function out
			# co is channel output eg the thing store in the channel object
			append!(co,fo)
			assert(length(co)==last(r)) # append isn't the right choice if the lengths differ
			# I'm not sure, but there may be no good reason to support differing lengths
		end
	end
	# freeuntil needs information on other dependencies
	r
end

outputs(s::HistogramStep) = other_inputs(s) # for now I'm conflating the idea of perpulse and inputs vs inputs/outputs
perpulse_outputs(s::HistogramStep) = Symbol[]
other_outputs(s::HistogramStep) = outputs(s)
function inputs(s::HistogramStep, c::Channel, r::Range)
	out = Any[]
	for inp in inputs(s,c)
		push!(out, isa(inp,Histogram)?inp:inp[r])
	end
	out
end
function Base.range(s::HistogramStep, c::Channel)
	# for now I'm conflating the idea of per_pulse and histogram on inputs
	start = mindonethru(other_inputs(s,c))+1
	stop = mindonethru(perpulse_inputs(s,c))
	start:stop
end
function dostep(s::HistogramStep, c::Channel)
	f = getfunction(s)
	r = range(s,c)
	other_inputs_exist(s,c) || (return last(r):-1)
	fout = f(inputs(s,c,r)...) # it doesn't know the types here, tho f will, also I could write a macro so it knew types
	r
end

function dostep(s::ThresholdStep, c::Channel)
	s.do_if_able || (return false)
	other_inputs_exist(s,c) || (return false)
	n = s.to_watch_func(c[s.to_watch])
	if n >= s.threshold
		s.do_if_able = false
		f = getfunction(s)
		r=1:s.threshold
		fout = f(inputs(s,c,r)...)
		if length(outputs(s)) == 1 # fout will be a value
			c[outputs(s)[1]] =  fout
		elseif length(outputs(s)) > 1
			@assert(length(fout) == length(outputs(s)))
			for (fo, co) in zip(fout, outputs(s))
				c[fo]=co
			end
		end
		return true
	else
		return false
	end
end

graphlabel(s::MockPulsesStep) = "MockPulsesStep"
inputs(s::MockPulsesStep) = []
function dostep(s::MockPulsesStep, c::Channel)
	pulses = gettriggeredpulse(s.pg,s.pulses_per_step)
	rowstamps = getrowstamp(s.pg, s.pulses_per_step, s.last_rowstamp)
	s.last_rowstamp = rowstamps[end]
	append!(c[s.outputs[1]], pulses)
	append!(c[s.outputs[2]], rowstamps)
	s.pulses_per_step
end


# todisk
function g_require(parent::Union(JldFile,JldGroup), name::ASCIIString)
    exists(parent, name) ? (return parent[name]) : g_create(parent.plain, name)
end
# Create a new or update an existing dataset within an HDF5 object
# extends the dataset if required
function d_extend(parent::Union(JldFile,JldGroup), name::ASCIIString, value::Vector, range::UnitRange)
	d = d_require(parent, name, value)
	set_dims!(parent[name].plain, (maximum(range),))
	d[range] = value
	d
end
d_extend(parent::Union(JldFile,JldGroup), name::ASCIIString, value::Vector) = d_extend(parent, name, value, endof(parent[name])+1:endof(parent[name])+length(value))
d_update(parent::Union(JldFile,JldGroup), name::ASCIIString, value::Vector) = d_extend(parent, name, value, 1:endof(value))
function d_require(parent::Union(JldFile,JldGroup), name, value::Vector,chunksize = 10000)
	dims = ((1,), (-1,)) # create a minimum size dataset, zero isn't allowed
	exists(parent,name) ? parent[name] : d_create(parent.plain, name, eltype(value), dims, "chunk", (chunksize,))
end
function d_require(parent::Union(JldFile,JldGroup), name, elementype::Type ,chunksize = 10000)
	dims = ((1,), (-1,)) # create a minimum size dataset, zero isn't allowed
	exists(parent,name) ? parent[name] : d_create(parent.plain, name, elementype, dims, "chunk", (chunksize,))
end
# Create a new or update an existing attribute within an HDF5 object
# a_require will create the attribute if it doesn't exist, or assert that the existing attribute is equal to value
function a_require(parent::Union(JldFile,JldGroup),name::ASCIIString,value)
    if exists(attrs(parent), name)
    	a_read(parent, name) == value ? (return value) : error("new value $value != existing value $(a_read(parent,name)) for attr $parent[$name]")
	end
    attrs(parent)[name] = value	
end
# a_update will create or replace and existing attribute with value
function a_update(parent::Union(JldFile,JldGroup), name::ASCIIString, h::Histogram)
	a_update(parent, name*"_bin_centers",collect(bin_centers(h)))
	a_update(parent, name*"_counts",counts(h))
end
function a_update(parent::Union(JldFile,JldGroup),name::ASCIIString,value)
    if exists(attrs(parent.plain), name)
        old_value = a_read(parent, name)
    	old_value == value && (return value)
    	a_delete(parent.plain, name)
	end
    attrs(parent.plain)[name] = value
end
graphlabel(s::ToJLDStep) = "to JLD"
outputs(s::ToJLDStep) = []
function dostep(s::ToJLDStep,c::Channel)
	jld = jldopen(filename(c),"r+")
	for (sym,value) in zip(perpulse_inputs(s),perpulse_inputs(s,c))
		if string(sym) in names(jld)
			d=d_require(jld, string(sym), eltype(value))
			r = length(d)+1:length(value)
		else # account for the fact that the dataset is created with 1 element, not zero 
			d=d_require(jld, string(sym), eltype(value))
			r = 1:length(value)
		end
		if length(r)>0
			d_extend(jld, string(sym), value[r],r)
		end
	end
	for (sym, value) in  zip(other_inputs(s),other_inputs(s,c))
		info("writing $sym to $jld, probably too often")
		a_update(jld, string(sym), value)
	end
	close(jld)
end #dostep

default_vertex_attrs = AttributeDict()
default_vertex_attrs["fontsize"]=18.0
default_vertex_attrs["color"]="red"
default_vertex_attrs["fillcolor"]="lightgrey"
default_vertex_attrs["style"]="filled"
label(v)=v.label
vertex_labels(g) = [label(v) for v in vertices(g)]
function vertex(v,label)
	i=first(indexin([label], vertex_labels(g)))
	i==0 && error("no vertex with label $label in $v")
	vertices(g)[i]
end
add_exvertex!(g,name::Symbol, attrs_in) = add_exvertex!(g,string(name), attrs_in)
function add_exvertex!(g, name,attrs_in)
    verts = vertices(g)
    vert_labels = [label(v) for v in verts]
    if !(name in vert_labels)
        v = ExVertex(num_vertices(g)+1, name)
        attrs = attributes(v,g)
        attrs["label"]=name
        for (key,val) in default_vertex_attrs
        	attrs[key]=val
        end
        for (key,val) in attrs_in
        	attrs[key]=val
        end
        add_vertex!(g, v)
        v
    else
        first(verts[indexin([name],vert_labels)])
    end
end
add_perpulse_data!(g,name)=add_exvertex!(g,name,[("shape","doubleoctagon"),("type","perpulse_data")])
add_other_data!(g,name)=add_exvertex!(g,name,[("shape","octagon"),("type","other_data")])
add_func!(g,name)=add_exvertex!(g,name,[("shape","box"),("type","function")])
function Graphs.add_edge!(g,label1, label2)
    verts = vertices(g)
    vert_labels = [label(v) for v in verts]
    i1 = first(indexin([label1], vert_labels))
    i2 = first(indexin([label2], vert_labels))
    if i1 < 1 || i2 < 1
        error("error label $label1 or $label2 is not in $g")
    end
    add_edge!(g, verts[i1], verts[i2])
end
function add_step!(g, s)
    vf = add_func!(g,graphlabel(s))
    for p in perpulse_inputs(s)
        v=add_perpulse_data!(g,p)
        add_edge!(g,v,vf)
    end
    for o in other_inputs(s)
        v=add_other_data!(g,o)
        add_edge!(g,v,vf)
    end
    for p in perpulse_outputs(s)
        v=add_perpulse_data!(g,p)
        add_edge!(g,vf,v)
    end
    for o in other_outputs(s)
        v=add_other_data!(g,o)
        add_edge!(g,vf,v)
    end
    g
end
function graph(steps::Vector{AbstractStep})
	g=inclist(ExVertex, is_directed=true)
	for s in steps
    	add_step!(g,s)
	end
	g
end
function savegraph(fname,g)
dot = to_dot(g)
gviz = GraphViz.Graph(dot)
GraphViz.layout!(gviz,engine="dot")
	open("$fname.svg","w") do f 
	       GraphViz.writemime(f, MIME"image/svg+xml"(),gviz)
	end #do
	open("$fname.dot","w") do f
		write(f,dot)
	end
end


function mock_apply_calibration(filt_value, c::Calibration)
	filt_value.*c.a+c.b
end
function compute_whitenoise_filter(pulse, selection_good) 
	filter = mean(pulse[selection_good])
	normalization = (maximum(filter)-minimum(filter))./dot(filter, filter)
	filter*normalization
end
function filter1lag(pulse, filter)
	out = Array(Float64, length(pulse))
	for i=1:length(pulse)
		out[i] = dot(pulse[i], filter)
	end
	out
end

function selectfromcriteria(x...) # placeholder, kinda ugly to use and probalby a bit slow
	iseven(length(x)) || error("x must be indicator,criteria,indicator,criteria...")
	out = trues(length(x[1]))
	for i = 1:2:length(x)
		low, high = x[i+1]
		out &= low .< x[i] .< high
	end
	out
end

steps = AbstractStep[]
push!(steps, MockPulsesStep(TwoExponentialPulseGenerator{Int}(520, 100, 50, 200,1, 13.3, 100000,520,Normal(1000,1),30,0), 100, [:pulse,:rowstamp],0))
push!(steps, PerPulseStep(compute_summary, [:pulse, :pre_samples, :frame_time],
	[:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, :peak_index, :peak_value, :min_value]))
push!(steps, PerPulseStep(mock_apply_calibration, [:filt_value, :calibration], [:energy]) )
push!(steps, PerPulseStep(selectfromcriteria, [:pretrig_rms, :pretrig_rms_criteria, :peak_index, :peak_index_criteria, :postpeak_deriv, :postpeak_deriv_criteria], [:selection_good]))
#push!(steps, PerPulseStep(select_forstep, [:filt_value, :selection_criteria2],[true,false], [:selection_good],[true]))
push!(steps, ToJLDStep([:filt_value,:pretrig_rms, :energy]))
push!(steps, HistogramStep(update_histogram, [:filt_value_hist, :selection_good, :filt_value]))
push!(steps, HistogramStep(update_histogram, [:energy_hist, :selection_good, :energy]))
push!(steps, ThresholdStep(mock_calibration, [:filt_value_hist],[:calibration],:selection_good, sum, 100, true))
push!(steps, ThresholdStep(compute_whitenoise_filter, [:pulse, :selection_good], [:whitenoise_filter], :selection_good, sum, 100, true))
push!(steps, PerPulseStep(filter1lag, [:pulse, :whitenoise_filter], [:filt_value]))

push!(perpulse_symbols, :filt_value, :selection_good, :energy, :pulse, :rowstamp,
	:pretrig_mean, :pretrig_rms, :pulse_average, :pulse_rms, :rise_time, :postpeak_deriv, 
	:peak_index, :peak_value, :min_value, :selection_good)

g = graph(steps)

c=Dict{Symbol,Any}()
c[:pretrig_mean] = RunningVector(Float64)
c[:pretrig_rms] = RunningVector(Float64)
c[:pulse_average] = RunningVector(Float64)
c[:pulse_rms] = RunningVector(Float64)
c[:rise_time] = RunningVector(Float64)
c[:postpeak_deriv] = RunningVector(Float64)
c[:peak_index] = RunningVector(Uint16)
c[:peak_value] = RunningVector(Uint16)
c[:min_value] = RunningVector(Uint16)
c[:filt_value] = RunningVector(Float64)
c[:pretrig_rms] = RunningVector(Float64)
c[:selection_good] = RunningSumBitVector()
c[:energy] = RunningVector(Float64)
c[:filt_value_hist] = Histogram(0:1:20000)
c[:energy_hist] = Histogram(0:1:20000)
c[:calibration] = Calibration(0,0,false)
c[:pulse] = RunningVector(Vector{Int})
c[:rowstamp] = RunningVector(Int)
c[:pre_samples] = 100
c[:frame_time] = 1/100000
c[:peak_index_criteria] = (183,198)
c[:pretrig_rms_criteria] = (0.0,1.5)
c[:postpeak_deriv_criteria] = (0.0,10.0)

# figure out what we can free up
donethru(x::AbstractRunningVector) = length(x)
donethru(x::Vector{Float64}) = DONETHRU_MAX
earliest_needed_index(x) = donethru(x)+1
function earliest_needed_index(c::Channel, q::Symbol,p::Symbol)
	if q in names(c)
		return earliest_needed_index(c[q])
	elseif q==:to_disk
		return donethru_jld(c::Channel,q,p)+1
	end
end
function donethru_jld(c::Channel,q::Symbol,p::Symbol)
	jld = jldopen(filename(c),"r+")
	l = 0
	if string(p) in names(jld)
		l = length(jld[string(p)])
	end
	close(jld)
	l
end
function earliest_needed_index(parent::Symbol, c::Channel, g::AbstractGraph) 
	v = vertex(g,string(parent))
	children = visited_vertices(g,BreadthFirst(),v) # consider memoizing this
	shift!(children) # remove the first element, which is always v
	children_sym = [symbol(label(u)) for u in children]
	eni = [earliest_needed_index(c,q,parent) for q in children_sym]
	# println([(q,earliest_needed_index(c,q,parent)) for q in children_sym])
	filter!(x->x != nothing, eni)
	isempty(eni)?length(c[parent])+1:minimum(eni)
end

# savegraph("graph",g)
stepelapsed = Array(Float64, length(steps))
workdone = Array(Any, length(steps))
close(jldopen(filename(c),"w")) # wipe the test file
sumn = 0
for i = 1:5
	println("loop iteration $i")
	#do steps
	for (j,s) in enumerate(steps)
		tstart = time()
		workdone[j] = dostep(s,c)
		stepelapsed[i] = time()-tstart
	end

    # free unneeded data
	for q in perpulse_symbols
		d=c[q]
		# freeuntil!(d,min(earliest_needed_index(q,c,g)-1,length(d)))
	end
	println(workdone)
	println(stepelapsed)
end
jld = jldopen(filename(c), "r+")


