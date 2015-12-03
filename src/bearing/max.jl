#################################################################
# maximum bearing method 
#
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
bearing_max(filename::AbstractString) = bearing_max(GainPattern(filename))
