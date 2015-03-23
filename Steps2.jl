using RunningVectors

type RunningStep
	func::Function
	outs::NTuple
	ins::NTuple
	is1to1::Vector{Bool}
end
a=RunningVector(Float64)
b=RunningVector(Float64)
c=RunningVector(Float64)
d=4
f(c,d) = (c+d,c-d)
s = RunningStep(f, (a,b), (c,d), [true, false])


inputlengths1to1(s::RunningStep) = Int[length(o) for o in s.ins[s.is1to1]]
outputlengths1to1(s::RunningStep) = Int[length(o) for o in s.outs]

inputlengths1to1(s)
outputlengths1to1(s)

function range1to1(s)
	ol = outputlengths1to1(s)
	assert(all(ol.==ol[1])) # all output lengths should be the same
	il = inputlengths1to1(s)
	ol[1]+1:minimum(il)
end


getins(s,i) = [s.is1to1[j] ? s.ins[j][i] : s.ins[j] for j in 1:length(s.ins)]

function increment(n)
for j in 1:n
	push!(c,4)
end
r=range1to1(s)
println(r)
for i in r
	println("did $i")
	out = s.func(getins(s,i)...)
	for (k,o) in enumerate(out)
		push!(s.outs[k],o)
	end
end
end



function update_channel(ljhfile)
	getnewpulses(pulses!, ljhfile)

	#summarize
	oldsummarizelen = length(p_pretrig_mean)
	newsummarizelen = length(pulses)
	for j=oldsummarizelen+1:newsummarizelen
		summaryj = summarize(pulses[j])
		p_pretrig_mean[j] = summaryj[1]
		p_pulse_mean[j] = summaryj[2]
		p_rise_time[j] = summary[3]
	end

	# cuts
	oldcutlen = length(good)
	newcutlen = length(pulses)
	for j = oldcutlen+1:newcutlen
		good[j] = pt_low < p_pretrig_mean[j] < pt_high && pm_low < p_pulse_mean < pt_high
	end

	# calibrate
	if !alreadycalibrated()
		if npulses > num_pulses_for_calibration
			p_filt_value_calibration = calibrate(p_filt_value[good])
			lastneeded = inf # note taht we no longer have a hold on past p_filt_value entries
		end
	end

	# apply calibration
	if alreadycalibrated()
		oldenergylen = length(p_energy)
		newenergylen = length(p_filt_value)
		for j = oldenergylen+1:newenergylen
			p_energy[j] = applycal(p_filt_value_calibration, p_filt_value[j])
		end
	end

	# drift correction
	if alreadydriftcorrected()
		olddclen = length(p_filt_value_dc)
		newdclen = length(p_filt_value)
		for j in olddclen+1:newdclen
			applyt(dcinfo, p_filt_value[j])
		end
	end

	# histogram
	counts, edges = histogram(p_filt_value_dc[alreadyhistogrammed:1:end],edges)
	p_filt_value_histogram += counts 

	# freeup old stuff
	freeuntil(pulses, lastneeded(pulses))
	freeuntil(p_filt_value, lastneeded(p_filt_value))
end



end #update_channel
