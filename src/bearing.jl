######################################################################
# bearing.jl
# Bearing calculation methods:
# bearing_cc - 
# bearing_half
######################################################################

###########################################################################
# Cross-correlation
#
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


#################################################################
# half-way point
function bearing_half(gp_file::AbstractString, half_val::Real=3)
	gp = GainPattern(gp_file)
	return bearing_half(gp)
end
function bearing_half(gp::GainPattern, half_val::Real=3)

	# Determine index and value of maximum gain
	max_idx = indmax(gp.meangains)
	max_gain = gp.meangains[max_idx]
	max_angle = gp.angles[max_idx]
	desired_gain = max_gain - half_val

	# Determine how many angles on a side
	num_angles = length(gp.angles)
	side_angles = div(num_angles, 2)

	# Find the -3 db point
	best_sub = Inf
	best_add = Inf
	sub_angle = 0.0
	add_angle = 0.0
	for i = 1:side_angles
		sub_idx = max_idx - i
		if sub_idx <= 0
			sub_idx += num_angles
		end
		add_idx = max_idx + i
		if add_idx > num_angles
			add_idx -= num_angles
		end
		if abs(gp.meangains[sub_idx] - desired_gain) < best_sub
			sub_angle = gp.angles[sub_idx]
			best_sub = abs(gp.meangains[sub_idx] - desired_gain)
		end
		if abs(gp.meangains[add_idx] - desired_gain) < best_add
			add_angle = gp.angles[add_idx]
			best_add = abs(gp.meangains[add_idx] - desired_gain)
		end
	end

	# return the half way point between the two angles
	#println("sub_angle = ", sub_angle)
	#println("add_angle = ", add_angle)
	ret_angle = add_angle + angular_error_rel(add_angle, sub_angle) / 2.
	if ret_angle < 0.1
		ret_angle += 360.0
	end
	return ret_angle
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
