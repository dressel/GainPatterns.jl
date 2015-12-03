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
