module GainPatterns

# package code goes here
import StatsBase.sample
export GainPattern, validgain, rotate!
export angular_error, angular_error_rel
export bearing_ls, bearing_cc, bearing_mle, bearing_sls, bearing_max
export bearing_ccn
export sample, rotate!, csv, normalize, normalize!, addsamples!, sampleGains

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
	# WHY DO I MAKE IT A NAN??
	# Maybe it isn't totally useless info?
	function GainPattern(angles::Vector{Float64}, gains::Vector{Float64})

		# Error if the angles and gains vectors are of different lengths
		if length(angles) != length(gains)
			error("Number of angles must match number of gains")
		end

		# Create samples vectors from gains
		num_gains = length(gains)
		#samples = Array(Vector{Float64}, num_gains)
		samples = Array(Vector{Float64}, 0)
		new_angles = Array(Float64, 0)
		new_gains = Array(Float64, 0)
		for i = 1:num_gains
			if validgain(gains[i])
				push!(samples, [gains[i]])
				push!(new_angles, angles[i])
				push!(new_gains, gains[i])
				#samples[i] = [gains[i]]
			else
				#samples[i] = []
				#gains[i] = NaN
			end
		end
		return new(new_angles, new_gains, samples)
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
	function GainPattern(datafile::AbstractString)

		# Read the data, creating the angles vector and an array of data
		alldata = readcsv(datafile)
		angles = float(alldata[:, 1])
		gains = alldata[:, 2:end]

		if size(gains, 2) == 1
			return GainPattern(angles, vec(gains))
		else
			return GainPattern(angles, gains)
		end

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


# Normalizes a vector by subtracting the mean and dividing by std dev
# Pretends invalid gains aren't even there
function normalize!(gains::Vector{Float64})

	temp_gains = copy(gains)
	len = length(temp_gains)
	for j = len:-1:1
		if !validgain(temp_gains[j])
			deleteat!(temp_gains, j)
		end
	end

	m = mean(temp_gains)
	s = std(temp_gains)

	for i = 1:len
		if validgain(gains[i])
			gains[i] = (gains[i] - m) / s
		end
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
# ANGULAR ERROR
# Error returned is always a mangitude (positive quantity)
# Assumes both errors are between zero and 360
###########################################################################
function angular_error(angle1::Real, angle2::Real)
	angular_error(float(angle1), float(angle2))
end
function angular_error(angle1::Float64, angle2::Float64)
	err = angle2 - angle1
	if err > 180
		err = 360 - err
	elseif err < -180
		err += 360
	end
	if err < 0
		err *= -1
	end
	return err
end
function angular_error{T<:Real}(angles1::Vector{T}, angles2::Vector{T})
	angular_error(float(angles1), float(angles2))
end
function angular_error(angles1::Vector{Float64}, angles2::Vector{Float64})
	@assert length(angles1) == length(angles2)
	err_len = length(angles1)
	err_arr = zeros(err_len)
	for i = 1:err_len
		err_arr[i] = angular_error(angles1[i], angles2[i])
	end
	return err_arr
end
function angular_error{T<:Real}(angles1::Vector{T}, angle2::Real)
	err_len = length(angles1)
	angles2 = angle2 * ones(err_len)
	return angular_error(float(angles1), angles2)
end
function angular_error{T<:Real}(angle1::Real, angles2::Vector{T})
	return angular_error(angles2, angle1)
end

# OLD! Inaccurate!
function angular_error_relold{T1<:Real, T2<:Real}(v1::Vector{T1}, v2::Vector{T2})
	err_arr = v1 - v2
	temp = v1 - (360 + v2)
	indy = abs(temp) .< abs(err_arr)
	err_arr[indy] = temp[indy]
	return err_arr
end

# Finds the shortest distance between points, but preserves a sign
# say v2[i] = 30, v1[i] = 15 => ans[i] = 15
# say v2[i] = 15, v1[i] = 30 => ans[i] = -15
function angular_error_rel{T1<:Real, T2<:Real}(v1::Vector{T1}, v2::Vector{T2})
	eps_val = 1e-6
	ref = angular_error(v1, v2)
	len = length(v1)
	for i = 1:len
		if abs(v2[i] - v1[i] - ref[i]) < eps_val
			# we know it is positive
			ref[i] = ref[i]
		elseif abs(360 + v2[i] - v1[1] - ref[i]) < eps_val
			ref[i] = ref[i]
		else
			ref[i] = -ref[i]
		end

	end
	return ref
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



# Assumes gp only has one sample (looks at meangains)
# Doesn't handle null obs... does it have to?
function bearing_max(gp::GainPattern)
	num_angles = length(gp.angles)
	temp_max = -Inf
	max_bearing = 0
	for i = 1:num_angles
		if validgain(gp.meangains[i])
			if gp.meangains[i] > temp_max
				temp_max = gp.meangains[i]
				max_bearing = gp.angles[i]
			end
		end
	end
	return max_bearing
end

###########################################################################
# Cross-correlation
###########################################################################
# Function carries out cross-correlation to determine the bearing
# angles is the vector of angles corresponding to the samples
#
# Returns vector containing shifts at every degree
#	ith index corresponds to (i-1) shift
#
# The idea is that the reference is of length 360, 
#  but the sample need not be.
# Therefore, a vector of sampled angles is provided.
# TODO: give this a thorough scrubbing. It is important to me that it works
function bearing_cc(angles, gains, ref_angles, ref_gains)
	cval = 0.0
	max_cval = -Inf
	max_i = 0.0

	norm_gains = copy(gains)
	normalize!(norm_gains)

	for i = 0:359
		cval = cross_correlate(angles, norm_gains, ref_angles, ref_gains, i)
		if (cval > max_cval)
			max_cval = cval
			max_i = i
		end
	end
	return 360.0 - max_i
end
# returns max correlation coefficient
function bearing_ccn(angles, gains, ref_angles, ref_gains)
	cval = 0.0
	max_cval = -Inf
	max_i = 0.0

	norm_gains = copy(gains)
	normalize!(norm_gains)

	for i = 0:359
		cval = cross_correlate(angles, norm_gains, ref_angles, ref_gains, i)
		if (cval > max_cval)
			max_cval = cval
			max_i = i
		end
	end
	return max_cval
end

function bearing_cc(angles, gains, ref_file::AbstractString)
	ref_gp = GainPattern(ref_file)
	return bearing_cc(angles, gains, ref_gp.angles, ref_gp.meangains)
end

# Returns array of bearing values
# Go through all samples in test_file
# Assumes all samples are the same length in test_gp
function bearing_cc(test_file::AbstractString, ref_file::AbstractString)
	test_gp = GainPattern(test_file)
	ref_gp = GainPattern(ref_file)

	num_samples = length(test_gp.angles)
	num_exp = length(test_gp.samples[1])

	aoa_array = zeros(num_exp)
	angles = zeros(num_samples)
	gains = zeros(num_samples)
	for i = 1:num_exp
		for j = 1:num_samples
			# Loop through all samples in this experiment...
			angles[j] = test_gp.angles[j]
			gains[j] = test_gp.samples[j][i]
		end
		aoa_array[i] = bearing_cc(angles, gains, ref_gp.angles, ref_gp.meangains)
	end
	return aoa_array
end

function bearing_ccn(test_file::AbstractString, ref_file::AbstractString)
	test_gp = GainPattern(test_file)
	ref_gp = GainPattern(ref_file)

	num_samples = length(test_gp.angles)
	num_exp = length(test_gp.samples[1])

	aoa_array = zeros(num_exp)
	angles = zeros(num_samples)
	gains = zeros(num_samples)
	for i = 1:num_exp
		for j = 1:num_samples
			# Loop through all samples in this experiment...
			angles[j] = test_gp.angles[j]
			gains[j] = test_gp.samples[j][i]
		end
		aoa_array[i] = bearing_ccn(angles, gains, ref_gp.angles, ref_gp.meangains)
	end
	return aoa_array
end


# Helper function to perform cross-correlation
# This is the cross-correlation function used in Graefenstein 2009
# Handle null observations
function cross_correlate(angles, gains, ref_angles, ref_gains, shift)
	c = 0.0
	len = length(angles)
	for i = 1:len
		if validgain(gains[i])
			c += gains[i] * interp_linear(angles[i], shift, ref_angles, ref_gains)
		end
	end
	return c
end

###########################################################################
# Least-squares
###########################################################################
function bearing_ls(angles, gains, ref_angles, ref_gains)
	@assert length(angles) == length(gains)
	num_gains = length(gains)
	best_err = Inf
	best_shift = 0
	for shift = 0:359
		err = 0.0
		for j = 1:num_gains
			if validgain(gains[j])
				ref_gain = interp_linear(angles[j], shift, ref_angles, ref_gains)
				err += (gains[j] - ref_gain)^2
			end
		end
		if err < best_err
			best_err = err
			best_shift = shift
		end
	end
	return 360 - best_shift
end

# Uses least-square method
# samples is how many samples you want to use (maybe you dont want whole rotation)
# sense is direction of rotation: -1 for negative sense, 1 for positive sense
function bearing_ls(dir_file::AbstractString, omni_file::AbstractString, ref_file::AbstractString; samples::Int64=0, sense::Int64=1)
	test_gp = GainPattern(dir_file) - GainPattern(omni_file) + 5
	ref_gp = GainPattern(ref_file)

	num_samples = length(test_gp.angles)
	num_exp = length(test_gp.samples[1])

	# Allowing user to specify count
	if (samples > 0) && (samples < num_samples)
		num_samples = samples
	end

	aoa_array = zeros(num_exp)
	angles = zeros(num_samples)
	gains = zeros(num_samples)
	for i = 1:num_exp
		if sense > 0
			for j = 1:num_samples
				# Loop through all samples in this experiment...
				angles[j] = test_gp.angles[j]
				gains[j] = test_gp.samples[j][i]
			end
		elseif sense == -1
			# Backwards loop, starting from last measurement
			for j = 1:num_samples
				angles[j] = test_gp.angles[end - j + 1]
				gains[j] = test_gp.samples[end - j + 1][i]
			end
		elseif sense == -2
			# Backwards loop, starting from origin
			for j = 2:num_samples
				angles[j] = test_gp.angles[end - j + 2]
				gains[j] = test_gp.samples[end - j + 2][i]
			end
			angles[1] = test_gp.angles[1]
			gains[1] = test_gp.samples[1][i]
		end
		aoa_array[i] = bearing_ls(angles, gains, ref_gp.angles, ref_gp.meangains)
	end
	return aoa_array
end


###########################################################################
# Scaled Least-squares
###########################################################################
function bearing_sls(angles, gains, ref_angles, ref_gains)
	@assert length(angles) == length(gains)
	num_gains = length(gains)
	best_err = Inf
	best_shift = 0
	best_scale = 0
	temp_gains = copy(gains)
	for scale_factor = 0:100
		temp_gains += 1
		for shift = 0:359
			err = 0.0
			for j = 1:num_gains
				if validgain(gains[j])
					ref_gain = interp_linear(angles[j], shift, ref_angles, ref_gains)
					err += (temp_gains[j] - ref_gain)^2
				end
			end
			if err < best_err
				best_err = err
				best_shift = shift
				best_scale = scale_factor
			end
		end
	end
	return 360 - best_shift
end

function bearing_sls(dir_file::AbstractString, ref_file::AbstractString; samples::Int64=0, sense::Int64=1)
	test_gp = GainPattern(dir_file)
	ref_gp = GainPattern(ref_file)

	num_samples = length(test_gp.angles)
	num_exp = length(test_gp.samples[1])

	# Allowing user to specify count
	if (samples > 0) && (samples < num_samples)
		num_samples = samples
	end

	aoa_array = zeros(num_exp)
	angles = zeros(num_samples)
	gains = zeros(num_samples)
	for i = 1:num_exp
		if sense > 0
			for j = 1:num_samples
				# Loop through all samples in this experiment...
				angles[j] = test_gp.angles[j]
				gains[j] = test_gp.samples[j][i]
			end
		elseif sense == -1
			# backwards loop, starting from last measurement
			for j = 1:num_samples
				angles[j] = test_gp.angles[end - j + 1]
				gains[j] = test_gp.samples[end - j + 1][i]
			end
		elseif sense == -2
			# backwards loop, starting from origin
			for j = 2:num_samples
				angles[j] = test_gp.angles[end - j + 2]
				gains[j] = test_gp.samples[end - j + 2][i]
			end
			angles[1] = test_gp.angles[1]
			gains[1] = test_gp.samples[1][i]
		end
		aoa_array[i] = bearing_sls(angles, gains, ref_gp.angles, ref_gp.meangains)
	end
	return aoa_array
end


###########################################################################
# MLE
###########################################################################
# TODO: ensure that you are within bounds for obs_matrix indices
function bearing_mle(angles, gains, obs_matrix, min_gain)
	@assert length(angles) == length(gains)
	null_index = size(obs_matrix, 2)
	num_gains = length(gains)
	best_prob = 0.0
	best_shift = -1
	for shift = 0:359
		prob = 1.0
		for j = 1:num_gains
			# Converting theta_rel to index
			theta_rel = angles[j] - shift
			if theta_rel < 0
				theta_rel += 360
			end
			theta_rel = mod(int(theta_rel/10.0), 36) + 1
			if validgain(gains[j])
				# Convert gain to index
				index = max(gains[j] - min_gain + 1, 1)
			else
				index = null_index
			end
			prob *= max(obs_matrix[theta_rel, index], 1e-6)
		end
		if prob > best_prob
			best_prob = prob
			best_shift = shift
		end
	end
	return best_shift
end


###########################################################################
# INTERPOLATION
# Assumes ref_angles always given in 0-360 order
# First ref_angle will always be zero, but last will be less than 360
###########################################################################
# Finds the value (from ref_gains) at angle + shift
# Looks like linear interpolation to me
function interp_linear(angle, shift, ref_angles, ref_gains)
	# Add the shift to the desired angle, be sure to mod it
	angle = mod(angle+shift, 360)

	# Find the two ref_angles this is between
	# final i value will be between i, i+1
	# this is index of lower close value
	i = 1
	ref_len = length(ref_angles)
	while (i < ref_len) && (angle > ref_angles[i+1])
		i += 1
	end

	if i == ref_len
		# Between last and first ref_angles
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

function interp_nearest(angle, shift, ref_angles, ref_gains)

	# Determine angle and ref_angle indices it falls between
	angle = mod(angle+shift, 360)
	i = 1
	ref_len = length(ref_angles)
	while (i < ref_len) && (angle > ref_angles[i+1])
		i += 1
	end

	if i == ref_len
		# between last ref_angle and first ref_angle
		d1 = abs(ref_angles[i] - angle)
		d2 = abs(360 - angle)
		if d2 > d1
			val = ref_gains[i]
		else
			val = ref_gains[1]
		end
	else
		# Determine which you are closer too...
		d1 = abs(ref_angles[i] - angle)
		d2 = abs(ref_angles[i+1] - angle)
		if d2 > d2
			val = ref_gains[i]
		else
			val = ref_gains[i+1]
		end
	end
	return val
end



# Writes a gain pattern to a file
# TODO: check that the file exists
# TODO: Check that they add csv file extension if not, add it
# Currently only prints out meangains, but could theoretically do samples too
function csv(gp::GainPattern, filename::AbstractString="temp.csv"; samples=true)
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


###########################################################################
# EXTREMA
# I don't remember what I used this for, but it appears unused
###########################################################################
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
