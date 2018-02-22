	% Set 1 of Returns 
	%sret1_mat 	= mvnrnd(mu_vec,sigma_mat,T_c);
	%emu1_vec 	= mean(sret1_mat)';
	%esigma1_mat = mycov1(sret1_mat);
	%x1_vec 		= esigma1_mat\one_vec;
	%esigma2g1_c = ((one_vec'*x1_vec)^(-1));
	%er1_mat 	= inv(esigma1_mat)-esigma2g1_c*(x1_vec*x1_vec');
	%z1_vec 		= er1_mat*emu1_vec;
	%ewg1_vec 	= esigma2g1_c*x1_vec;
	%emug1_c 	= emu1_vec'*ewg1_vec;
	%y1_vec 		= esigma1_mat\emu1_vec;

    % Investment Rules 
	%iesigma1_mat= inv(esigma1_mat);

    %ewt1_vec 	= (1/GAMMA)*iesigma1_mat*emu1_vec;
	%ewg1_vec	= (1/((one_vec')*(esigma1_mat\one_vec)))*(esigma1_mat\one_vec);
	%ewh1_vec 	= (1/GAMMA)*er1_mat*emu1_vec;
	%ewq1_vec	= (1/(one_vec'*(esigma1_mat\one_vec)))*(esigma1_mat\one_vec)+(1/GAMMA)*er1_mat*emu1_vec;
    
	%conststar1_c = (1/GAMMA)*(emug1_c/esigma2g1_c);
	%sutilstar1_c = calU(CONSTs*ewg1_vec+ewh1_vec, emu1_vec, esigma1_mat, GAMMA);

