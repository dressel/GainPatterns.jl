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
