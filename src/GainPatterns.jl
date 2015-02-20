module GainPatterns

# package code goes here
export GainPattern, validgain, plot, rotate!
export normalize, normalize!, sampleGains, crosscorrelate
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
	# TODO: make use of the refine function
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

	# Create gain pattern from samples
	# First, we refine the samples (remove invalid entries, compute mean)
	# Then, just call the default constructor
	function GainPattern{T1<:Real,T2<:Real}(angles::Vector{T1}, samples::Vector{Vector{T2}})
		meangains, refined_samples = refine(samples)
		return new(angles, meangains, refined_samples)
	end

	# Create gain pattern from csv file
	# Assumes first column is angles
	# This cuts out any null observations..
	function GainPattern(datafile::String)

		# Read the data, creating the angles vector and an array of data
		alldata = readcsv(datafile)
		angles = float(alldata[:, 1])
		gains = alldata[:, 2:end]

		# Create a vector containing samples
		num_angles = length(angles)
		samples = Array(Vector{Float64}, num_angles)

		# Create a vector for each sample from the matrix of data
		# If some angles have fewer samples, we will have empty quotes, ""
		# We assume that these are after the actual samples
		for i = 1:num_angles
			tempvec = vec(gains[i,:])
			while tempvec[end] == ""
				pop!(tempvec)
			end
			samples[i] = tempvec
		end

		# Call the constructor 
		return GainPattern(angles, samples)
	end

end # type GainPattern


# Determes if gain value is real
# The idea is that the user could over write this
function validgain(gain::Real)
	nullobs = 2147483647.
	if gain != nullobs
		return true
	else
		return false
	end
end


# Eliminates any invalid gain values, and computes the mean
# Returns vector with means, and new vector with samples
function refine(samples::Vector{Vector{Float64}})

	# Calculate number of angles and create required data structures
	num_angles = length(samples)
	refined_samples = Array(Vector{Float64}, num_angles)
	meangains = Array(Float64, num_angles)

	# Loop through all sample vectors, 
	#  eliminating invalid gain values and calculating the mean
	for i = 1:num_angles
		refined_samples[i] = Array(Float64, 0)
		rowavg = 0.0
		numvalid = 0
		num_samples = length(samples[i])
		for j = 1:num_samples
			tempgain = samples[i][j]
			if validgain(tempgain)
				rowavg += tempgain
				numvalid += 1
				push!(refined_samples[i], tempgain)
			end
		end
		meangains[i] = rowavg / numvalid
	end

	return meangains, refined_samples
end


# Returns a normalized distribution
# This normalization was shown in Graefenstein 2009
# Creates a copy of the gains then calls the in-place version
function normalize{T<:Real}(gains::Vector{T})
	new_gains = float(gains)	# create copy of gains, while ensuring float
	normalize!(new_gains)
	return new_gains
end

# In place version of the above.
function normalize!(gains::Vector{Float64})

	# Calculate mean, std deviation, and vector length
	m = mean(gains)
	s = std(gains)
	l = length(gains)

	# Perform the normalization
	for i = 1:l
		gains[i] = (gains[i] - m) / s
	end
end

# TODO: Normalize function that takes in gain pattern
#  How do we normalize both the means and the samples?
#  Do we just use the means and std deviations of the meangains vector?

# Rotates the gain pattern by a specified number of degrees
# We must also call mod so the degrees stay in the proper range
function rotate!(gp::GainPattern, degrees::Real)
	for i = 1:length(gp.angles)
		gp.angles[i] = mod(gp.angles[i] + degrees, 360)
	end
end


# Sample from a given distribution of gains
# Returns the gains sampled, and the angles at which these samples occurred
# TODO: Look into this, see if it is still working
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
function plot(gp::GainPattern; ymin = 0.0, showsamples::Bool=false)

	n_angles = length(gp.angles)
	mingain = minimum(gp.meangains)

	# TODO: check that samples exists
	if showsamples
		# Plot showing the samples...
		# First we need to determine the mean gain...

		# Once the means are found, determine absolute smallest gain
		# Set mingain to minimum(ymin, minimum(gains), 0)
		# Loop through all samples...
		# TODO: Really, mingain should be called ymin. Clean that up
		for i = 1:n_angles
			tempmin = minimum(gp.samples[i])
			mingain = (tempmin < mingain ? tempmin : mingain)
		end
		mingain = (ymin < mingain ? ymin : mingain)
		mingain = (mingain > 0. ? 0. : mingain)

		gain_plot = plotgains(gp.angles, gp.meangains, mingain)
		plot_array = plotsamples(gp.angles, gp.samples, mingain)
		push!(plot_array, gain_plot)
		pa = PolarAxis(plot_array, yticklabel="{\\pgfmathparse{$mingain+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}")

	else
		# Determine the min gain...
		mingain = (ymin < mingain ? ymin : mingain)
		mingain = (mingain > 0. ? 0. : mingain)

		gain_plot = plotgains(gp.angles, gp.meangains, mingain)
		pa = PolarAxis(gain_plot, yticklabel="{\\pgfmathparse{$mingain+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}")
	end

	return pa

end

# Creates Plots.Linear object (from PGFPlots) with gains vs angles
# Can handle negative values, which radial plots normally cannot.
#
# ymin = minimum y (gain or radial) value. If it is positive, it is ignored.
#  If it is less than the minimum value of gains, also ignored.
function plotgains{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{T2}, ymin::Real)

	# Make copies before we make some changes
	# Remember to shift gains by the minimum value to be plotted
	plot_gains = copy(gains) - ymin
	plot_angles = copy(angles)

	# Last point must be same as first point to complete the plot
	# If not, it will be missing a section between last point and first
	if angles[1] != angles[end]
		push!(plot_angles, plot_angles[1])
		push!(plot_gains, plot_gains[1])
	end

	# Create a linear plot and return it
	return Plots.Linear(plot_angles, plot_gains, mark="none", style="red,thick")
end

# Returns an array of plots
function plotsamples{T1<:Real,T2<:Real}(angles::Vector{T1}, samples::Vector{Vector{T2}}, ymin::Real)

	# Create an array of plots
	num_angles = length(angles)
	plot_array = Array(Plots.Linear, num_angles)

	# Create a linear plot for each sample
	for i = 1:num_angles
		plot_array[i] = Plots.Linear(angles[i]*ones(length(samples[i])), samples[i]-ymin, mark="x", style="blue,smooth")
	end

	# Return the array of plots
	return plot_array
end


end # module
