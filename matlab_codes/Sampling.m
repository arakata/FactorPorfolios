UsingZeroInvestment

[NTIME NASSET] = size(ret_mat);

muhat_vec 		= mean(ret_mat)';
sigmahat_mat  	= mycov1(ret_mat);
sigmatilde_mat  = mycov2(ret_mat);
rhat_mat		= inv(sigmahat_mat)-(1/((one_vec')*(sigmahat_mat\one_vec)))*(sigmahat_mat\one_vec)*((sigmahat_mat\one_vec)');
rtilde_mat		= inv(sigmatilde_mat)-(1/((one_vec')*(sigmatilde_mat\one_vec)))*(sigmatilde_mat\one_vec)*((sigmatilde_mat\one_vec)');

% desicion rules 

wqhat_vec = (sigmahat_mat\one_vec)/((one_vec')*(sigmahat_mat\one_vec))+(1/gamma_c)*(rhat_mat*muhat_vec);
wghat_vec = (sigmahat_mat\one_vec)/((one_vec')*(sigmahat_mat\one_vec));
whhat_vec = (1/gamma_c)*(rhat_mat*muhat_vec);

wqtilde_vec = (sigmatilde_mat\one_vec)/((one_vec')*(sigmatilde_mat\one_vec))+(1/gamma_c)*(rtilde_mat*muhat_vec);
wgtilde_vec = (sigmatilde_mat\one_vec)/((one_vec')*(sigmatilde_mat\one_vec));
whtilde_vec = (1/gamma_c)*(rtilde_mat*muhat_vec);

% out of sample performance
% a)
eUtildeofwqhatGIVENsigma 	= mug_c -(gamma_c/2)*sigma2g_c+(psi2_c^2)/(2*gamma_c)-(NASSET-1)/(2*gamma_c*NASSET);
rhoofwqwqhatGIVENsigma 		= (NASSET-1)/(2*gamma_c*NTIME);
% b)
eUtildeofwqhatGIVENmu		= mug_c-(gamma_c/2)*((NTIME-2)/(MTIME-NASSET-1))*sigma2g_c+(1-K1)*(psi2_c/(2*gamma_c));