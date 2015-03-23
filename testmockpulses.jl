include("mockpulses.jl")

pg = TwoExponentialPulseGenerator{Int}(520, 100, 50, 200,1, 100, 100000,520,Normal(1000,1))
getcleanpulse(pg,1000) 

using PyPlot
plot(getcleanpulse(pg,1000))
plot(gettriggeredpulse(pg))
