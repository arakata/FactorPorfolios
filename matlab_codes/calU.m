function [output] = calU(w_vec, mu_vec, sigma_mat, gamma_c)
	output=(mu_vec')*w_vec-(gamma_c/2)*((w_vec')*(sigma_mat*w_vec));
end

