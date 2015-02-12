module GainPatterns

# package code goes here
export normalize, normalize!, sampleGains, crosscorrelate, plotgains
import PGFPlots: PolarAxis, Plots, save

# Returns a normalized distribution
# This normalization was shown in Graefenstein 2009
function normalize(gains::Vector{Float64})
	newGains = Array(Float64, size(gains))
	return normalize!(gains)
end

# In place version of the above.
function normalize!(gains::Vector{Float64})
	m = mean(gains)
	s = std(gains)
	l = length(gains)
	for i = 1:l
		gains[i] = (gains[i] - m) / s
	end
	return gains
end

# Sample from a given distribution of gains
# Returns the gains sampled, and the angles at which these samples occurred
function sampleGains(gains, step::Int64, bearing::Int64)

	noise = 0.25*randn(360)

	# antenna_angles are angles antenna is oriented at
	# sampled_angles is angle antenna observes in its measurements
	antenna_angles = [0:step:359]
	sampled_angles = mod(antenna_angles - bearing, 360)
	sampled_gains = gains[sampled_angles + 1] + noise[sampled_angles + 1]
	normalize!(sampled_gains)
	return sampled_gains, antenna_angles
end

# Function carries out cross-correlation to determine the bearing
# angles is the vector of angles corresponding to the samples
#
# Returns vector containing shifts at every degree
#	ith index corresponds to (i-1) shift
# TODO: Make this more efficient
#
# The idea is that the reference is of length 360, 
#  but the sample need not be.
# Therefore, a vector of sampled angles is provided.
function crosscorrelate(ref::Vector{Float64}, sample::Vector{Float64}, angles)
	shifts = [0:359]
	cs = zeros(length(shifts))
	for shift in shifts
		cs[shift+1] = cc_helper(shift, ref, sample, angles)
	end
	return cs
end


# Helper function to perform cross-correlation
# This is the cross-correlation function used in Graefenstein 2009
# Inputs are:
#	rarkk
function cc_helper(shift, ref, sample, angles)
	i = 1
	c = 0
	for deg in angles
		c += sample[i] * ref[mod(deg - shift, 360) + 1]
		i += 1
	end
	return c
end


# Plots gains and angles
function plotgains(angles::Vector{Float64}, gains::Vector{Float64})

	# Make copies before we make some changes
	plot_gains = copy(gains)
	plot_angles = copy(angles)

	# First, we must handle negatives in gains
	labelgain = 0
	mingain = minimum(gains)
	if mingain < 0
		# If negative gains, add this to gains...
		# This does not mess with the gains passed in
		plot_gains -= mingain
		labelgain = mingain
	end

	# Last point must be same as first point to complete the plot
	# If not, it will be missing a section between last point and first
	if angles[1] != angles[end]
		push!(plot_angles, plot_angles[1])
		push!(plot_gains, plot_gains[1])
	end

	# Finally create the plot
	# TODO: Add yticklabel stuff to make it good
	p = Plots.Linear(plot_angles, plot_gains, mark="none")
	pa = PolarAxis(p)

	# Save to a pdf file...
	# TODO: give user option to specify this...
	save("temp.pdf", pa)
end

end # module
