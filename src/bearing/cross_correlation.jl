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
