function [output] = ecalUtilde(what_vec, mu_vec, sigma_mat, gamma_c)
	output=(mu_vec')*what_vec-(gamma_c/2)*((what_vec')*(sigma_mat*what_vec));
end

