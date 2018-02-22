T 	= 200;
NS 	= 20000;

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

%%%%%%%%%%%
% corollary 
%%%%%%%%%%%

wghat_vec'*sigma_mat*whhat_vec	

% should be near zero, and it is not really if T is low.

%%%%%%%%%%%%%%%%
% Two fund rules
%%%%%%%%%%%%%%%%

% it is possible to replicate wqhat_vec by making a mixing of wghat and whhat (with a weight).
% note the linear relation between wqhat_vec-wghat_vec and whhat_vec in

plot(wqhat_vec-wghat_vec,whhat_vec)

% optimal constant C=C*=Cs

Cs 					= ((M-N)*(M-N-3)/(M*(M-2)))*(PSI2HAT/(PSI2HAT+(N-1)/M));

ECALUTILDE_WQHATCs	= MUHATg-(GAMMA/2)*((M-2)/(M-N-1))*SIGMAHAT2g + (PSI2HAT/(2*GAMMA))*((M-N)*(M-N-3)/((M-2)*(M-N-1)))*(PSI2HAT/(PSI2HAT+(N-1)/M));


% two fund rules of Kan and Zhou

CKZ2s					= ((M-N-1)*(M-N-4)/(M*(M-2)))*(THETA2HAT/(THETA2HAT+N/M));
ECALUTILDE_WKZ2HATCKZ2s	= (THETA2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(THETA2HAT/(THETA2HAT+N/M));

%%%%%%%%%%%%%%%%%%%
% Three fund rules
%%%%%%%%%%%%%%%%%%%

Ds 						= (1/GAMMA)*(MUHATg/SIGMAHAT2g)*((M-N-1)/(M-2));
ECALUTILDE_WQ3IHATCsDs	= (TAU2HAT/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2HAT/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2HAT/(PSI2HAT+(N-1)/M));

Css 						= Cs;
Dss 						= (1/GAMMA)*MUHATg*(((M-N-1)*(M-N-4))/(M*(M-2)));

ECALUTILDE_WQ3IIHATCssDss 	= (TAU2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2HAT/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2HAT)/(PSI2HAT+(N-1)/M)); 

% Three fund rule revisited

CKZ3s 							= (((M-N-1)*(M-N-4))/(M*(M-2)))*((PSI2HAT)/(PSI2HAT+N/M));
DKZ3s 							= (1/GAMMA)*(((M-N-1)*(M-N-4))/(M*(M-2)))*((N/M)/(PSI2HAT+N/M))*MUHATg;

ECALUTILDE_WKZ3HATCKZ3sDKZ3s	= (THETA2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2HAT+(THETA2HAT/PSI2HAT)*(N/M)));

