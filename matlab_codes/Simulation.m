T 	= 100;
NS 	= 20000;

% Assume: Parameters Unknown

const_vec	= linspace(0.1,5,50);
util_mat 	= ones(50,NS);
tfr1_vec    = ones(1,NS);
tfr2_vec    = ones(1,NS);
tfrr_vec    = ones(1,NS);

for i = 1:NS  

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

    THETA2HAT       = (muhat_vec')*(sigmahat_mat\muhat_vec);
    TAU2HAT         = (((muhat_vec')*(sigmahat_mat\one_vec))^2)/((one_vec')*(sigmahat_mat\one_vec));
    PSI2HAT         = (muhat_vec')*(rhat_mat*muhat_vec);

    % Three fund rules

    Ds                      = (1/GAMMA)*(MUHATg/SIGMAHAT2g)*((M-N-1)/(M-2));
    ECALUTILDE_WQ3IHATCsDs  = (TAU2HAT/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2HAT/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2HAT/(PSI2HAT+(N-1)/M));

    Css                         = Cs;
    Dss                         = (1/GAMMA)*MUHATg*(((M-N-1)*(M-N-4))/(M*(M-2)));

    ECALUTILDE_WQ3IIHATCssDss   = (TAU2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2HAT/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2HAT)/(PSI2HAT+(N-1)/M)); 

    % Three fund rule revisited

    CKZ3s                           = (((M-N-1)*(M-N-4))/(M*(M-2)))*((PSI2HAT)/(PSI2HAT+N/M));
    DKZ3s                           = (1/GAMMA)*(((M-N-1)*(M-N-4))/(M*(M-2)))*((N/M)/(PSI2HAT+N/M))*MUHATg;

    ECALUTILDE_WKZ3HATCKZ3sDKZ3s    = (THETA2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2HAT+(THETA2HAT/PSI2HAT)*(N/M)));

    tfr1hat_vec(1,i) = ECALUTILDE_WQ3IHATCsDs;
    tfr2hat_vec(1,i) = ECALUTILDE_WQ3IIHATCssDss;
    tfrrhat_vec(1,i) = ECALUTILDE_WKZ3HATCKZ3sDKZ3s;

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

    THETA2TILDE       = (mutilde_vec')*(sigmatilde_mat\mutilde_vec);
    TAU2TILDE         = (((mutilde_vec')*(sigmatilde_mat\one_vec))^2)/((one_vec')*(sigmatilde_mat\one_vec));
    PSI2TILDE         = (mutilde_vec')*(rtilde_mat*mutilde_vec);

    % Three fund rules

    Ds                      = (1/GAMMA)*(MUTILDEg/SIGMATILDE2g)*((M-N-1)/(M-2));
    ECALUTILDE_WQ3IHATCsDs  = (TAU2TILDE/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2TILDE/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2TILDE/(PSI2TILDE+(N-1)/M));

    Css                         = Cs;
    Dss                         = (1/GAMMA)*MUTILDEg*(((M-N-1)*(M-N-4))/(M*(M-2)));

    ECALUTILDE_WQ3IIHATCssDss   = (TAU2TILDE/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2TILDE/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2TILDE)/(PSI2TILDE+(N-1)/M)); 

    % Three fund rule revisited

    CKZ3s                           = (((M-N-1)*(M-N-4))/(M*(M-2)))*((PSI2TILDE)/(PSI2TILDE+N/M));
    DKZ3s                           = (1/GAMMA)*(((M-N-1)*(M-N-4))/(M*(M-2)))*((N/M)/(PSI2TILDE+N/M))*MUTILDEg;

    ECALUTILDE_WKZ3HATCKZ3sDKZ3s    = (THETA2TILDE/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2TILDE+(THETA2TILDE/PSI2TILDE)*(N/M)));

    tfr1tilde_vec(1,i) = ECALUTILDE_WQ3IHATCsDs;
    tfr2tilde_vec(1,i) = ECALUTILDE_WQ3IIHATCssDss;
    tfrrtilde_vec(1,i) = ECALUTILDE_WKZ3HATCKZ3sDKZ3s;


    % Investment Rules 
	for j = 1:length(const_vec)
		utilhat_mat(j,i) = calU(const_vec(j)*wghat_vec+whhat_vec, mu_vec, sigma_mat, GAMMA);
		utiltilde_mat(j,i) = calU(const_vec(j)*wgtilde_vec+whtilde_vec, mu_vec, sigma_mat, GAMMA);		
	end

end

[UTILH,UTILHX] = max(utilhat_mat);
[UTILT,UTILTX] = max(utiltilde_mat);

x = [UTILH'-UTILs';UTILT'-UTILs'];
g = [ones(size(UTILT'-UTILs')); 2*ones(size(UTILT'-UTILs'))];

x1 = [tfr1hat_vec';tfr2hat_vec';tfrrhat_vec'];
g1 = [ones(size(tfr1hat_vec')); 2*ones(size(tfr2hat_vec')); 3*ones(size(tfrrhat_vec'))];
boxplot(x1,g1)

x2 = [tfr1tilde_vec';tfr2tilde_vec';tfrrtilde_vec'];
g2 = [ones(size(tfr1tilde_vec')); 2*ones(size(tfr2tilde_vec')); 3*ones(size(tfrrtilde_vec'))];
boxplot(x2,g2)

% In general it is not possible to improve the performance just by chance (see T=200). 
% If T decreases considerably, say T=40, then it is possible to find such alternative. 

subplot(1,2,1)
boxplot(x1,g2,'Labels',{'TFRI','TFRII','TFRR'})
ylim([min([x1',x2']),max([x1',x2'])])
title('HAT Estimator')

subplot(1,2,2)
boxplot(x2,g2,'Labels',{'TFRI','TFRII','TFRR'})
ylim([min([x1',x2']),max([x1',x2'])])
title('Tilde Estimator')

print('FillPageFigure','-dpdf')