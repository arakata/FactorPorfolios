%%%%%%%%%%%%%
% question 1: 
%%%%%%%%%%%%%

% is it possible to do better than wt_vec? by selecting another mixture?

% Assume: Parameters Knwon

conststar_c	=(1/gamma_c)*(mug_c/sigma2g_c);
sutilstar_c =calU(conststar_c*wg_vec+wh_vec, mu_vec, sigma_mat, gamma_c);

m_const 	=50;
sutil_vec 	=ones(m_const,1);
for i=1:m_const
	const_c 		=i/15;
	sutil_vec(i) 	=calU(const_c*wg_vec+wh_vec, mu_vec, sigma_mat, gamma_c);
end

plot([sutil_vec, sutilstar_c*ones(m_const,1)])

% no. it is not possible. conststar_c=1.6 (approx) is the best option.   

T_c 	=200;
nsim_c 	=20000;

% Assume: Parameters Unknown

	% Set 1 of Returns 
	sret1_mat 	= mvnrnd(mu_vec,sigma_mat,T_c);
	emu1_vec 	= mean(sret1_mat)';
	esigma1_mat = cov(sret1_mat,1);
	x1_vec 		= esigma1_mat\one_vec;
	esigma2g1_c = ((one_vec'*x1_vec)^(-1));
	er1_mat 	= inv(esigma1_mat)-esigma2g1_c*(x1_vec*x1_vec');
	z1_vec 		= er1_mat*emu1_vec;
	ewg1_vec 	= esigma2g1_c*x1_vec;
	emug1_c 	= emu1_vec'*ewg1_vec;
	y1_vec 		= esigma1_mat\emu1_vec;

    % Investment Rules 
	iesigma1_mat= inv(esigma1_mat);

    ewt1_vec 	= (1/gamma_c)*iesigma1_mat*emu1_vec;
	ewg1_vec	= (1/((one_vec')*(esigma1_mat\one_vec)))*(esigma1_mat\one_vec);
	ewh1_vec 	= (1/gamma_c)*er1_mat*emu1_vec;
	ewq1_vec	= (1/(one_vec'*(esigma1_mat\one_vec)))*(esigma1_mat\one_vec)+(1/gamma_c)*er1_mat*emu1_vec;
    
	conststar1_c = (1/gamma_c)*(emug1_c/esigma2g1_c);
	sutilstar1_c = calU(conststar_c*ewg1_vec+ewh1_vec, emu1_vec, esigma1_mat, gamma_c);

const_vec	= linspace(0.1,5,m_const);
sutil2_mat 	= ones(m_const,nsim_c);
for i=1:nsim_c  

	% Set 2 of Returns (independent of Set 1)
	sret2_mat 	= mvnrnd(mu_vec,sigma_mat,T_c);
    emu2_vec 	= mean(sret2_mat)';
    esigma2_mat = cov(sret2_mat,1);
    x2_vec 		= esigma2_mat\one_vec;
    esigma2g2_c	= ((one_vec'*x2_vec)^(-1));
    er2_mat 	= inv(esigma2_mat)-esigma2g2_c*(x2_vec*x2_vec');
    Z2_vec 		= er2_mat*emu2_vec;
    ewg2_vec 	= esigma2g2_c*x2_vec;
    emug2_c 	= emu2_vec'*ewg2_vec;
    
    % Investment Rules 
	iesigma2_mat= inv(esigma2_mat);

    ewt2_vec 	= (1/gamma_c)*iesigma2_mat*emu2_vec;
	ewg2_vec	= (1/((one_vec')*(esigma2_mat\one_vec)))*(esigma2_mat\one_vec);
	ewh2_vec 	= (1/gamma_c)*er2_mat*emu2_vec;
	ewq2_vec	= (1/(one_vec'*(esigma2_mat\one_vec)))*(esigma2_mat\one_vec)+(1/gamma_c)*er2_mat*emu2_vec;

	for j=1:m_const 
		sutil2_mat(j,i) = calU(const_vec(j)*ewg2_vec+ewh2_vec, emu1_vec, esigma1_mat, gamma_c);
	end

end

[sutil2val,sutil2index] = max(sutil2_mat);
boxplot(sutil2val-sutilstar1_c)

% It is not possible to improve the performance just by chance. 
