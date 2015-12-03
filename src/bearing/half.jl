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
