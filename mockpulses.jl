using Distributions
abstract PulseGenerator

immutable TwoExponentialPulseGenerator{T<:Real} <: PulseGenerator
	record_length::Int
    pre_rise_points::Int
	rise_points::Float64
	fall_points::Float64
	noise_fwhm::Float64
	event_rate_hz::Float64
    samples_per_second::Float64
    min_points_since_last_pulse::Float64
    amplitude_distribution::Distribution
    numrows::Int
    rownumber::Int
    function TwoExponentialPulseGenerator(record_length,pre_rise_points,rise_points,fall_points,noise_fwhm,event_rate_hz,samples_per_second,min_points_since_last_pulse,amplitude_distribution,numrows,rownumber)
    	rise_points<fall_points || error("rise_points should be less than fall points")
    	new(record_length,pre_rise_points,rise_points,fall_points,noise_fwhm,event_rate_hz,samples_per_second,min_points_since_last_pulse,amplitude_distribution,numrows,rownumber)
    end
end

function two_exponential_pulses(
	points::Int, rise_points, fall_points, arrival_point::Tuple, amplitude::Tuple)
	pulse = zeros(Float64,points)
	@assert length(arrival_point) == length(amplitude)
	k=kfactor(1/rise_points,1/fall_points) 
	for (ap, amp) in zip(arrival_point, amplitude)
		@inbounds for i = 1:points
			t = i-ap
			pulse[i] += t>0.0?amp*k*(exp(-t/rise_points)-exp(-t/fall_points)):0.0
		end
	end
	pulse
end

# factor for calculating mathematical parameters from physical parameters and maintainging
# constant amplitude
# Otto-von-Guericke Univ. Magdeburg, Magdeburg, Germany
# IEEE Transactions on Electromagnetic Compatibility (Impact Factor: 1.35). 12/2010; DOI: 10.1109/TEMC.2010.2052621 
# Estimation of the Mathematical Parameters of Double-Exponential Pulses Using the Nelder–Mead Algorithm - ResearchGate. Available from: http://www.researchgate.net/publication/224153186_Estimation_of_the_Mathematical_Parameters_of_Double-Exponential_Pulses_Using_the_NelderMead_Algorithm [accessed Mar 20, 2015].kfactor(a,b) = 1/(exp(-a*(log(a/b)/(a-b)))-exp(-b*(log(a/b)/(a-b))))
kfactor(a,b) = 1/(exp(-a*(log(a/b)/(a-b)))-exp(-b*(log(a/b)/(a-b)))) # memoize?
const X = [0.842411, 1.660852, 8.887335, 0.624129]
const Y = [2.167370, 0.360349, 0.013289, 1.358701, 0.137057, 122.802861, 1.314850] 
alpha(t_rise, t_fwhm) = (X[1]/((t_fwhm/t_rise)^X[2]-X[3])^X[4])/t_rise
beta(t_rise, t_fwhm) = (Y[1]-Y[2]*exp(-t_fwhm*Y[3]/t_rise)-Y[4]*exp(-t_fwhm*Y[5]/t_rise)-Y[6]*exp(-t_fwhm*Y[7]/t_rise))/t_rise

white_noise(points::Int, noise_fwhm) = rand(points)*noise_fwhm/(2*sqrt(2*log(2))) # fwhm to std
function white_noise!(pulse, noise_fwhm) 
	@inbounds for i = 1:length(pulse) pulse[i]+=randn()*noise_fwhm end
	pulse
end
to_type_and_white_noise{T<:Integer}(::Type{T}, noise_fwhm, pulse) = convert(Vector{T}, round(white_noise!(pulse, noise_fwhm)))


function getcleanpulse{T}(pg::TwoExponentialPulseGenerator{T}, amp)
	to_type_and_white_noise(T, pg.noise_fwhm, 
		two_exponential_pulses(pg.record_length, pg.rise_points, pg.fall_points, (pg.pre_rise_points,), (amp,)))
end
function gettriggeredpulse{T}(pg::TwoExponentialPulseGenerator{T})
	last_arrival, next_arrival = rand(Exponential(pg.samples_per_second/pg.event_rate_hz),2)
	amplitudes = tuple(rand(pg.amplitude_distribution,3)...)
	arrival_points = (pg.pre_rise_points-pg.min_points_since_last_pulse-last_arrival
		, pg.pre_rise_points, pg.pre_rise_points+next_arrival)
	pulse = two_exponential_pulses(pg.record_length, pg.rise_points, pg.fall_points, arrival_points, amplitudes)
	to_type_and_white_noise(T, pg.noise_fwhm, pulse)
end
gettriggeredpulse(pg,n) = [gettriggeredpulse(pg) for i=1:n]
times_s(pg) = range(0,1/pg.samples_per_second, pg.record_length)

# Base.start(pg::TwoExponentialPulseGenerator) = 0.0
# function Base.next(pg::TwoExponentialPulseGenerator, last_arrival::Float64)
# 	timestamp = last_arrival+rand(Exponential(pg.samples_per_second/pg.event_rate_hz),2)
# 	pulse = gettriggeredpulse(pg)
# 	return ((timestamp,pulse), timestamp)
# end
# Base.done(pg::TwoExponentialPulseGenerator, last_arrival::Float64) = false

function getrowstamp(pg::TwoExponentialPulseGenerator,n::Int, last_rowstamp::Int=0)
	rowstamp_diffs = pg.numrows*rand(Exponential(pg.samples_per_second/pg.event_rate_hz),n)
	last_rowstamp = pg.numrows*div(last_rowstamp,pg.numrows)+pg.rownumber # make sure rowstamp reflect the correct row
	rowstamp_diffs[1]+=last_rowstamp
	cumsum(convert(Vector{Int}, round(rowstamp_diffs)))
end


