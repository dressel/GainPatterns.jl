module GainPatterns

# package code goes here
export GainPattern, validgain, plot
export normalize, normalize!, sampleGains, crosscorrelate, plotgains
export PolarAxis, save, Axis
import PGFPlots: PolarAxis, Plots, save, Axis

# angles is the list of angles
# meangains is the list of gains
# samples is a vector of gain vectors.
#  The mean of each of these makes up mean gains
type GainPattern

	angles::Vector{Float64}
	meangains::Vector{Float64}
	samples::Vector{Vector{Float64}}

	# Default constructor
	# TODO: Create a real default... for all three to be specified
	GainPattern(angles::Vector{Float64}, gains::Vector{Float64}) = new(angles, gains)

	# Allow vectors of any reals to work
	GainPattern{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{T2}) = new(float(angles), float(gains))

	# Allow user to input matrix
	# Makes no assumption about the validity of all the inputs
	function GainPattern{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Matrix{T2})

		# Determine number of angles and samples we are looking at
		n_angles, n_samples = size(gains)

		# Create two of the main data structures we need (angles given to us)
		meangains = Array(Float64, n_angles)
		samples = Array(Vector{Float64}, n_angles)

		# Loop over all samples to determine
		for i = 1:n_angles
			samples[i] = Array(Float64, 0)
			rowavg = 0.0
			numvalid = 0
			for j = 1:n_samples
				tempgain = gains[i,j]
				if validgain(tempgain)
					rowavg += tempgain
					numvalid += 1
					push!(samples[i], tempgain)
				end
			end
			meangains[i] = rowavg / numvalid
		end

		# We now have all the data structures we need
		return new(angles, meangains, samples)
	end

	# Makes the assumption that all of these entries are valid
	function GainPattern{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{Vector{T2}})
		n_angles = length(angles)
		meangains = Array(Float64, n_angles)
		for i = 1:n_angles
			meangains[i] = mean(gains[i])
		end

		return new(angles, meangains, gains)
	end

	# Create gain pattern from csv file
	# Assumes first column is angles
	# This cuts out any null observations..
	function GainPattern(datafile::String)

		# Read the data
		alldata = readcsv(datafile)
		angles = float(alldata[:, 1])
		gains = float(alldata[:, 2:end])

		# Call the constructor we all know and love
		return GainPattern(angles, gains)
	end

end # type GainPattern


# Function determining if gain value is real
# The idea is that the user could over write this
function validgain(gain::Real)
	nullobs = 2147483647.
	if gain != nullobs
		return true
	else
		return false
	end
end

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


# Allows plotting with gain pattern
function plot(gp::GainPattern; ymin = 0.0)
	plotgains(gp.angles, gp.meangains)
end

# Plots gains and angles on a polar plot.
# Can handle negative values, which radial plots normally cannot.
# Saves the plot to temp.pdf
#
# ymin = minimum y (gain or radial) value. If it is positive, it is ignored.
#  If it is less than the minimum value of gains, also ignored.
# 
# TODO: Test the mingain stuff, including the default
function plotgains{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{T2}, ymin::Real)

	# Make copies before we make some changes
	plot_gains = copy(gains)
	plot_angles = copy(angles)

	# Set mingain to minimum(ymin, minimum(gains), 0)
	mingain = minimum(gains)
	mingain = (ymin < mingain ? ymin : mingain)
	mingain = (mingain > 0. ? 0. : mingain)

	# Shift the plotgains by mingain
	plot_gains -= mingain

	# Last point must be same as first point to complete the plot
	# If not, it will be missing a section between last point and first
	if angles[1] != angles[end]
		push!(plot_angles, plot_angles[1])
		push!(plot_gains, plot_gains[1])
	end

	# Finally create the plot
	p = Plots.Linear(plot_angles, plot_gains, mark="none", style="red,thick")
	pa = PolarAxis(p, yticklabel="{\\pgfmathparse{$mingain+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}")
end

# Allow gains to be plotted without specifying minimum gain
function plotgains{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{T2})
	plotgains(angles, gains, minimum(gains))
end

# Plots errors and stuff
# User passes in a matrix.
# Each row corresponds to the measurements for a specific angle.
function plotgains{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Matrix{T2}, ymin::Real)

	# Number of angles, number of samples per angle comes up often
	nangles = size(gains, 1)
	nsamples = size(gains, 2)

	# Find the means and make a copy of the angles
	meangains = reshape(mean(gains,2), nangles)
	plot_angles = copy(angles)

	# Create vector of linear plots
	# One for each angle plus one for all of them
	plot_array = Array(Plots.Linear, nangles + 1)

	# Once the means are found, determine absolute smallest gain
	# Set mingain to minimum(ymin, minimum(gains), 0)
	mingain = minimum(gains)
	mingain = (ymin < mingain ? ymin : mingain)
	mingain = (mingain > 0. ? 0. : mingain)

	# Shift all gains by this...
	meangains -= mingain

	# Last point must be same as first point to complete the plot
	# If not, it will be missing a section between last point and first
	if angles[1] != angles[end]
		push!(plot_angles, plot_angles[1])
		push!(meangains, meangains[1])
	end

	# Finally create the plot
	plot_array[1] = Plots.Linear(plot_angles, meangains, mark="none", style="red, thick, solid")
	p = plot_array[1]

	# Plot all the errors...
	for i = 1:length(angles)
		plot_array[i+1] = Plots.Linear(angles[i] * ones(nsamples), reshape(gains[i,:], nsamples)-mingain, mark="x", style="blue,smooth")
	end

	# Finally, create the polar axis
	pa = PolarAxis(plot_array, yticklabel="{\\pgfmathparse{$mingain+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}")
	#pa = PolarAxis(p, yticklabel="{\\pgfmathparse{$mingain+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}")
end

function plotgains{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Matrix{T2})
	plotgains(angles, gains, minimum(gains))
end

function plot_error{T<:Real}(gains::Vector{T}, angle::Real)
	p = Plots.Linear(angle*ones(length(gains)), gains)
end

end # module
