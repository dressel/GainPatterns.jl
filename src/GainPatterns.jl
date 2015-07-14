module GainPatterns

# package code goes here
export GainPattern, validgain, plot, rotate!
export normalize, normalize!, sampleGains, bearing_cc, addsamples!, csv
export PolarAxis, save, Axis
import PGFPlots: PolarAxis, Plots, save, Axis
import StatsBase.sample
export sample

_nullobs = 2147483647.

# angles is the list of angles
# meangains is the list of gains
# samples is a vector of gain vectors.
#  The mean of each of these makes up mean gains
type GainPattern

	angles::Vector{Float64}
	meangains::Vector{Float64}
	samples::Vector{Vector{Float64}}


	# Constructor where user passes in angles and mean gains
	# The samples array uses the single gain at each angle as that gain
	# TODO: check that these gains are valid
	# TODO: just don't insert that entry if invalid
	function GainPattern(angles::Vector{Float64}, gains::Vector{Float64})

		# Error if the angles and gains vectors are of different lengths
		if length(angles) != length(gains)
			error("Number of angles must match number of gains")
		end

		# Create samples vectors from gains
		num_gains = length(gains)
		samples = Array(Vector{Float64}, num_gains)
		for i = 1:num_gains
			if validgain(gains[i])
				samples[i] = [gains[i]]
			else
				samples[i] = []
				gains[i] = NaN
			end
		end
		return new(angles, gains, samples)
	end

	# Constructor taking in only angles and mean gains
	# Allows vectors of type real, and converts to float arrays
	function GainPattern{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{T2})
		GainPattern(float(angles), float(gains))
	end

	# Allow user to input matrix
	# Makes no assumption about the validity of all the inputs
	function GainPattern{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Matrix{T2})

		# Determine number of angles and samples we are looking at
		num_angles, num_samples = size(gains)

		# Turn gains matrix into proper samples vector, then refine and return
		samples = Array(Vector{Float64}, num_angles)
		for i = 1:num_angles
			samples[i] =  vec(gains[i,:])
		end
		meangains, refined_samples = refine(samples)
		return new(angles, meangains, refined_samples)
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

		return GainPattern(angles, samples)
	end

end # type GainPattern



# Determes if gain value is valid
# In my research, I receive a value of 2147483647 if the gain is invalid
# If this is the value of the gain, I know it is not valid
# The user can override this to meet their own data needs
function validgain(gain::Real)
	#nullobs = 2147483647.
	gain == _nullobs ? false : true
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
			push!(refined_samples[i], tempgain)
			if validgain(tempgain)
				rowavg += tempgain
				numvalid += 1
			end
		end
		meangains[i] = rowavg / numvalid
	end

	return meangains, refined_samples
end


# Subtracts gp2 from gp1
# Only works if they both have the same number of angles
# Assumes that the angles are the same
# TODO: Should I check that the angles are the same?
# If the samples are all of equal length, these are subtracted
# Otherwise, only the means are subtracted,
#  and a new samples array is initialized (with only the mean as a sample)
function -(gp1::GainPattern, gp2::GainPattern)
	num_angles1 = length(gp1.angles)
	num_angles2 = length(gp2.angles)

	if num_angles1 != num_angles2
		error("Number of angles must be the same!")
	end

	# Loop through samples, ensuring they are the same length
	same_sample_length = true
	num_angles = num_angles1
	gains = Array(Float64, num_angles)
	for i = 1:num_angles
		if length(gp1.samples[i]) != length(gp2.samples[i])
			same_sample_length = false
		end
		gains[i] = gp1.meangains[i] + gp2.meangains[i]
	end

	# If the samples are all the same length, add them
	# Otherwise, just create a new gain pattern with the calculated means
	samples = Array(Vector{Float64}, num_angles)
	if same_sample_length
		for i = 1:num_angles
			#samples[i] = gp1.samples[i] + gp2.samples[i]
			samples[i] = Array(Float64, 0)
			for j = 1:length(gp1.samples[1])
				if !validgain(gp1.samples[i][j]) || !validgain(gp2.samples[i][j])
					push!(samples[i], _nullobs)
				else
					push!(samples[i], gp1.samples[i][j] - gp2.samples[i][j])
				end
			end
		end
		gp_new = GainPattern(gp1.angles, samples)
	else
		gp_new = GainPattern(gp1.angles, gains)
	end

	return gp_new
end


# Addition function for gain patterns.
# Identical to subtraction, except for + sign in a few spots instead of -
# Leave invalid observations alone
function +(gp1::GainPattern, gp2::GainPattern)
	num_angles1 = length(gp1.angles)
	num_angles2 = length(gp2.angles)

	if num_angles1 != num_angles2
		error("Number of angles must be the same!")
	end

	# Loop through samples, ensuring they are the same length
	same_sample_length = true
	num_angles = num_angles1
	gains = Array(Float64, num_angles)
	for i = 1:num_angles
		if length(gp1.samples[i]) != length(gp2.samples[i])
			same_sample_length = false
		end
		gains[i] = gp1.meangains[i] + gp2.meangains[i]
	end

	# If the samples are all the same length, add them
	# Otherwise, just create a new gain pattern with the calculated means
	samples = Array(Vector{Float64}, num_angles)
	if same_sample_length
		for i = 1:num_angles
			#samples[i] = gp1.samples[i] + gp2.samples[i]
			samples[i] = Array(Float64, 0)
			for j = 1:length(gp1.samples[1])
				if !validgain(gp1.samples[i][j]) || !validgain(gp2.samples[i][j])
					push!(samples[i], _nullobs)
				else
					push!(samples[i], gp1.samples[i][j] + gp2.samples[i][j])
				end
			end
		end
		gp_new = GainPattern(gp1.angles, samples)
	else
		gp_new = GainPattern(gp1.angles, gains)
	end

	return gp_new
end


# Adding a constant term to a gain pattern
# Leave any invalid observations alone
function +(gp::GainPattern, c::Real)
	gp_new = deepcopy(gp)
	num_angles = length(gp.angles)
	for i = 1:num_angles
		gp_new.meangains[i] += c
		for j = 1:length(gp_new.samples[i])
			if validgain(gp_new.samples[i][j])
				gp_new.samples[i][j] += c
			end
		end
	end

	return gp_new
end

+(c::Real, gp::GainPattern) = +(gp, c)
-(gp::GainPattern, c::Real) = +(gp, -c)

# Appends the samples of gp2 to the samples vector of gp1
# Recalculates the mean gains of gp1
function addsamples!(gp1::GainPattern, gp2::GainPattern)
	@assert gp1.angles == gp2.angles
	num_angles = length(gp1.angles)
	for i = 1:num_angles
		for j = 1:length(gp2.samples[i])
			push!(gp1.samples[i], gp2.samples[i][j])
		end
	end
	meangains, refined_samples = refine(gp1.samples)
	gp1.meangains = meangains
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
# TODO: Normalize function that takes in gain pattern
#  How do we normalize both the means and the samples?
#  Do we just use the means and std deviations of the meangains vector?
function normalize!(gains::Vector{Float64})
	m = mean(gains)
	s = std(gains)
	len = length(gains)

	for i = 1:len
		gains[i] = (gains[i] - m) / s
	end
end


# Calculates mean and standard deviation of meangains
# Applies the normalization (x-m) / s
# Also applies normalization to each of the samples, but only valid ones
function normalize!(gp::GainPattern)
	m = mean(gp.meangains)
	s = std(gp.meangains)
	num_angles = length(gp.angles)

	for i = 1:num_angles
		gp.meangains[i] = (gp.meangains[i]-m)/s
		for j = 1:length(gp.samples[i])
			if validgain(gp.samples[i][j])
				gp.samples[i][j] = (gp.samples[i][j] - m) / s
			end
		end
	end
end


# Rotates the gain pattern by a specified number of degrees
# We must also call mod so the degrees stay in the proper range
function rotate!(gp::GainPattern, degrees::Real)
	for i = 1:length(gp.angles)
		gp.angles[i] = mod(gp.angles[i] + degrees, 360)
	end
end


###########################################################################
# SAMPLING
###########################################################################
# Samples at a specified angle
# TODO: allow for angles not specified in pattern. Requires interpolation
function sample(gp::GainPattern, angle::Real)
	i = findfirst(gp.angles, angle)
	i == 0 ? error("Currently, you must sample from angle for which the gain pattern has samples") : nothing
	sample(gp.samples[i])
end

# Samples at a vector of specified angles
# These angles must all exist in the gain pattern's angle vector
# TODO: maybe allow for angles not specified in pattern?
function sample{T<:Real}(gp::GainPattern, angles::Vector{T})
	num_angles = length(angles)
	gains = Array(Float64, num_angles)
	for i = 1:num_angles
		idx = findfirst(gp.angles, angles[i])
		idx == 0 ? error("Currently, you must sample from angle for which the gain pattern has samples") : nothing
		gains[i] = sample(gp.samples[idx])
	end
	return gains
end

# Sample a gain for every angle in gp.
function sample(gp::GainPattern)
	sample(gp, gp.angles)
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
# TODO: handle null observations
# TODO: give this a thorough scrubbing. It is important to me that it works
function bearing_cc(angles, gains, ref_angles, ref_gains)
	cval = 0.0
	max_cval = -100000
	max_i = 0.0

	for i = 0:359
		cval = cross_correlate(angles, gains, ref_angles, ref_gains, i)
		if (cval > max_cval)
			max_cval = cval
			max_i = i
		end
	end
	return 360.0 - max_i
end


# Helper function to perform cross-correlation
# This is the cross-correlation function used in Graefenstein 2009
# Inputs are:
#	rarkk
function cross_correlate(angles, gains, ref_angles, ref_gains, shift)
	c = 0.0
	len = length(angles)
	for i = 1:len
		c += gains[i] * interp_gain(angles[i], shift, ref_angles, ref_gains)
	end
	return c
end

function interp_gain(angle, shift, ref_angles, ref_gains)
	# Add the shift to the desired angle, be sure to mod it
	angle = mod(angle+shift, 360)

	# Find the two ref_angles this is between
	ref_len = length(ref_angles)

	# final i value will be between i, i+1
	# find the i value that works...
	# this is index of lower close value
	i = 1
	while (i < ref_len) && (angle > ref_angles[i+1])
		i += 1
	end

	#if i == ref_len, between last and first
	if i == ref_len
		val = (360.0 - angle) * ref_gains[i]
		val += (angle - ref_angles[i]) * ref_gains[1]
		val = val / (360.0 - ref_angles[i])
	else
		val = (ref_angles[i+1] - angle) * ref_gains[i]
		val += (angle - ref_angles[i]) * ref_gains[i+1]
		val = val / (ref_angles[i+1] - ref_angles[i])
	end
	return val
end


###########################################################################
# PLOTTING
###########################################################################
# Returns the PolarAxis object containing a plot of the gain pattern
function plot(gp::GainPattern; ymin::Real=0.0, ymax=nothing, showsamples::Bool=false, lastleg::Bool=true, style=nothing, degrees::Bool=false)
	plot([gp], ymin=ymin, ymax=ymax, showsamples=showsamples, lastleg=lastleg, styles=[style],degrees=degrees)
end

function plot(gp_array::Vector{GainPattern}; ymin::Real=0.0, ymax=nothing, showsamples::Bool=false, lastleg::Bool=true, legendentries=nothing, colors=nothing, styles=nothing, degrees::Bool=false)

	# Create an array with length of angles for each gain pattern
	# Create an array with minimum mean gain for each pattern
	# Initialize the minimum_gain to the lowest of these minimum gains
	num_gp = length(gp_array)
	mingain_array = zeros(num_gp)
	nangles_array = zeros(Int, num_gp)
	for i = 1:num_gp
		nangles_array[i] = length(gp_array[i].angles)
		mingain_array[i] = minimum(gp_array[i].meangains)
	end
	mingain = minimum(mingain_array)

	# Currently, only plot samples for one gain pattern
	# TODO: Don't allow legend entries for showing samples
	#  Or do, but it has to be done well...
	if showsamples && (num_gp == 1)
		gp = gp_array[1]
		for i = 1:nangles_array[1]
			tempmin = minimum(gp.samples[i])
			mingain = (tempmin < mingain ? tempmin : mingain)
		end

		ymin = min(ymin, mingain, 0.0)
		if typeof(ymax) != Nothing
			ymax -= ymin
			emsg = "ymax is smaller than smallest gain, 0, and specified ymin"
			ymax < ymin ? error(emsg) : nothing
		end

		gain_plot = plotgains(gp.angles, gp.meangains, ymin, lastleg)
		plot_array = plotsamples(gp.angles, gp.samples, ymin)
		push!(plot_array, gain_plot)
		pa = PolarAxis(plot_array, ymax=ymax, yticklabel="{\\pgfmathparse{$ymin+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}")
	else
		ymin = min(ymin, mingain, 0.0)
		if typeof(ymax) != Nothing
			ymax -= ymin
			emsg = "ymax is smaller than smallest gain, 0, and specified ymin"
			ymax < ymin ? error(emsg) : nothing
		end

		# legendentries, styles must be indexable
		if legendentries == nothing
			legendentries = Array(Nothing, num_gp)
		end
		if styles == nothing
			styles = Array(Nothing, num_gp)
		end

		# Do we include the degrees
		if degrees
			xtl = "{\$\\pgfmathprintnumber{\\tick}^{\\circ}\$}"
		else
			xtl = nothing
		end

		plot_array = Array(Plots.Linear, num_gp)
		for i = 1:num_gp
			gp = gp_array[i]
			plot_array[i] = plotgains(gp.angles, gp.meangains, ymin, lastleg, legendentries[i], styles[i])
		end
		#xticklabel=$\pgfmathprintnumber{\tick}^\circ$
		pa = PolarAxis(plot_array, ymax=ymax, yticklabel="{\\pgfmathparse{$ymin+\\tick} \\pgfmathprintnumber{\\pgfmathresult}}", xticklabel=xtl)

	end

	# Return the polar axis object
	return pa

end

# Creates Plots.Linear object (from PGFPlots) with gains vs angles
# Can handle negative values, which radial plots normally cannot.
#
# ymin = minimum y (gain or radial) value. If it is positive, it is ignored.
#  If it is less than the minimum value of gains, also ignored.
function plotgains{T1<:Real,T2<:Real}(angles::Vector{T1}, gains::Vector{T2}, ymin::Real, lastleg::Bool, legendentry, style)

	# Make copies before we make some changes
	# Remember to shift gains by the minimum value to be plotted
	plot_gains = copy(gains) - ymin
	plot_angles = copy(angles)

	# Last point must be same as first point to complete the plot
	# If not, it will be missing a section between last point and first
	if (angles[1] != angles[end]) && lastleg
		push!(plot_angles, plot_angles[1])
		push!(plot_gains, plot_gains[1])
	end

	# Create a linear plot and return it
	#return Plots.Linear(plot_angles, plot_gains, mark="none", style="red,thick")
	return Plots.Linear(plot_angles, plot_gains, mark="none", legendentry=legendentry, style=style)
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

# Writes a gain pattern to a file
# TODO: check that the file exists
# TODO: Check that they add csv file extension if not, add it
# Currently only prints out meangains, but could theoretically do samples too
function csv(gp::GainPattern, filename::String="temp.csv"; samples=true)
	num_angles = length(gp.angles)
	outfile = open(filename, "w")

	if samples
		for i = 1:num_angles
			num_samples = length(gp.samples[i])
			write(outfile, "$(gp.angles[i])")
			for j = 1:num_samples
				write(outfile, ",$(gp.samples[i][j])")
			end
			write(outfile, "\n")
		end
	else
		for i = 1:num_angles
			write(outfile, "$(gp.angles[i]),$(gp.meangains[i])\n")
		end
	end

	close(outfile)
end

# Finds minimum and maximum gains in gain pattern's samples, in one sweep
# Checks if each gain is valid
# I don't know how I feel about that...
function Base.extrema(gp::GainPattern)
	num_angles = length(gp.angles)
	min_sample = Inf
	max_sample = -Inf
	for i = 1:num_angles
		num_samples = length(gp.samples[i])
		for j = 1:num_samples
			temp = gp.samples[i][j]
			if validgain(temp)
				max_sample = max(temp, max_sample)
				min_sample = min(temp, min_sample)
			end
		end
	end
	return min_sample, max_sample
end

# Returns the maximum sample in the gain pattern's samples
function Base.maximum(gp::GainPattern)
	extrema(gp)[2]
end

# Returns the minimum sample in the gain pattern's samples
function Base.minimum(gp::GainPattern)
	extrema(gp)[1]
end



end # module
