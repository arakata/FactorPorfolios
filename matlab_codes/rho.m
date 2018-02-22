function [output] = rho(ws_vec, wshat_vec, mu_vec,sigma_mat,gamma_c)
	output=calU(ws_vec, mu_vec, sigma_mat, gamma_c)-calUtilde(wshat_vec, mu_vec, sigma_mat, gamma_c);
end

