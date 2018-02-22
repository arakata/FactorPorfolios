% cd \Users\lf.rosalesm\Desktop\matlab_codes

UsingZeroInvestment

%%%%%%%%%%%%%
% question 1: 
%%%%%%%%%%%%%

% is it possible to do better than wt_vec? by selecting another mixture?

% Assume: Parameters Knwon

CONSTs		= (1/GAMMA)*(MUg/SIGMA2g);
UTILs 		= calU(CONSTs*wg_vec+wh_vec, mu_vec, sigma_mat, GAMMA);

const_vec	= linspace(0.1,5,50);
util_vec 	= 0*const_vec;

for i=1:length(const_vec)
	util_vec(i) 	= calU(const_vec(i)*wg_vec+wh_vec, mu_vec, sigma_mat, GAMMA);
end

plot([util_vec', UTILs*ones(length(const_vec),1)])

% no. it is not possible. CONSTs=1.6 (approx) is the best option.   

T 	= 40;
NS 	= 20000;

% Assume: Parameters Unknown

const_vec	= linspace(0.1,5,50);
util_mat 	= ones(50,NS);
for i=1:NS  

	% Set of Returns
	retsim_mat 		= mvnrnd(mu_vec,sigma_mat,T);

   	muhat_vec 		= mean(retsim_mat)';
   	mutilde_vec 	= mean(retsim_mat)'; % note that muhat_vec and mutilde_vec are identical.
    sigmahat_mat 	= cov(retsim_mat,1);
    sigmahat_mat 	= mycov1(retsim_mat);
    sigmatilde_mat 	= mycov2(retsim_mat);
    
    % with hat
    x_vec 		 	= sigmahat_mat\one_vec;
    SIGMAHAT2g		= ((one_vec'*x_vec)^(-1));
    rhat_mat 		= inv(sigmahat_mat)-SIGMAHAT2g*(x_vec*x_vec');
    z_vec 			= rhat_mat*muhat_vec;
    wghat_vec 		= SIGMAHAT2g*x_vec;
    MUHATg 			= muhat_vec'*wghat_vec;
	isigmahat_mat	= inv(sigmahat_mat);

    wthat_vec 	= (1/GAMMA)*isigmahat_mat*muhat_vec;
	wghat_vec	= (1/((one_vec')*(sigmahat_mat\one_vec)))*(sigmahat_mat\one_vec);
	whhat_vec 	= (1/GAMMA)*rhat_mat*muhat_vec;
	wqhat_vec	= (1/(one_vec'*(sigmahat_mat\one_vec)))*(sigmahat_mat\one_vec)+(1/GAMMA)*rhat_mat*muhat_vec;
    
    %with tilde
    x_vec 				= sigmatilde_mat\one_vec;
    SIGMATILDE2g		= ((one_vec'*x_vec)^(-1));
    rtilde_mat 			= inv(sigmatilde_mat)-SIGMATILDE2g*(x_vec*x_vec');
    z_vec 				= rtilde_mat*mutilde_vec;
    wgtilde_vec 		= SIGMATILDE2g*x_vec;
    MUTILDEg 			= mutilde_vec'*wgtilde_vec;
	isigmatilde_mat		= inv(sigmatilde_mat);

    wttilde_vec 	= (1/GAMMA)*isigmatilde_mat*mutilde_vec;
	wgtilde_vec		= (1/((one_vec')*(sigmatilde_mat\one_vec)))*(sigmatilde_mat\one_vec);
	whtilde_vec 	= (1/GAMMA)*rtilde_mat*mutilde_vec;
	wqtilde_vec		= (1/(one_vec'*(sigmatilde_mat\one_vec)))*(sigmatilde_mat\one_vec)+(1/GAMMA)*rtilde_mat*mutilde_vec;

    % Investment Rules 

	for j=1:length(const_vec)
		utilhat_mat(j,i) = calU(const_vec(j)*wghat_vec+whhat_vec, mu_vec, sigma_mat, GAMMA);
		utiltilde_mat(j,i) = calU(const_vec(j)*wgtilde_vec+whtilde_vec, mu_vec, sigma_mat, GAMMA);		
	end

end

[UTILH,UTILHX] = max(utilhat_mat);
[UTILT,UTILTX] = max(utiltilde_mat);

x = [UTILH'-UTILs';UTILT'-UTILs'];
g = [ones(size(UTILT'-UTILs')); 2*ones(size(UTILT'-UTILs'))];
boxplot(x,g)

% In general it is not possible to improve the performance just by chance (see T=200). 
% If T decreases considerably, say T=40, then it is possible to find such alternative. 