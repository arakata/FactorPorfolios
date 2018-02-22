
% some constants
K1 		= (M/(M-N-1))*(2-(M*(M-2))/((M-N)*(M-N-3)));
K2 		= ((N-1)*M*(M-2))/((M-N)*(M-N-1)*(M-N-3));

K1TILDE = 2-((M-2)*(M-N-1))/((M-N)*(M-N-3));
K2TILDE = ((N-1)*(M-2)*(M-N-1))/(M*(M-N)*(M-N-3));

%%%%%%%%%%%%%%%%
% proposition 1
%%%%%%%%%%%%%%%%

% a
ECALUTILDE_WQHAT_GIVENSIGMA 	= MUg -(GAMMA/2)*SIGMA2g+PSI2/(2*GAMMA)-(N-1)/(2*GAMMA*M);
RHO_WQHAT_GIVENSIGMA			= (N-1)/(2*GAMMA*M);

% b
ECALUTILDE_WQHAT_GIVENMU		= MUg-(GAMMA/2)*((M-2)/(M-N-1))*SIGMA2g+K1*(PSI2/(2*GAMMA));
RHO_WQHAT_GIVENMU				= (GAMMA/2)*((N-1)/(M-N-1))*SIGMA2g+(1-K1)*(PSI2/(2*GAMMA));

% c
ECALUTILDE_WQHAT 				= MUg-(GAMMA/2)*((M-2)/(M-N-1))*SIGMA2g+K1*(PSI2/(2*GAMMA))-K2*(1/(2*GAMMA));
RHO_WQHAT  						= (GAMMA/2)*((N-1)/(M-N-1))*SIGMA2g+(1-K1)*(PSI2/(2*GAMMA))+K2*(1/(2*GAMMA));

% proposition 2
ECALUTILDE_WQTILDE 				= MUg-(GAMMA/2)*((M-2)/(M-N-1))*SIGMA2g+K1TILDE*(PSI2/(2*GAMMA))-K2TILDE*(1/(2*GAMMA));
RHO_WQTILDE  					= (GAMMA/2)*((N-1)/(M-N-1))*SIGMA2g+(1-K1TILDE)*(PSI2/(2*GAMMA))+K2TILDE*(1/(2*GAMMA));

% corollary 

wghat_vec'*sigma_mat*whhat_vec	

% should be near zero, and it is not really if T is low.

%%%%%%%%%%%%%%%%
% Two fund rules
%%%%%%%%%%%%%%%%

% it is possible to replicate wqhat_vec by making a mixing of wghat and whhat (with a weight).
% note the linear relation between wqhat_vec-wghat_vec and whhat_vec in

plot(wqhat_vec-wghat_vec,whhat_vec)

% optimal constant C=C*=Cs

Cs 					= ((M-N)*(M-N-3)/(M*(M-2)))*(PSI2/(PSI2+(N-1)/M));

ECALUTILDE_WQHATCs	= MUg-(GAMMA/2)*((M-2)/(M-N-1))*SIGMA2g + (PSI2/(2*GAMMA))*((M-N)*(M-N-3)/((M-2)*(M-N-1)))*(PSI2/(PSI2+(N-1)/M));


% two fund rules of Kan and Zhou

CKZ2s					= ((M-N-1)*(M-N-4)/(M*(M-2)))*(THETA2/(THETA2+N/M));
ECALUTILDE_WKZ2HATCKZ2s	= (THETA2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(THETA2/(THETA2+N/M));

%%%%%%%%%%%%%%%%%%%
% Three fund rules
%%%%%%%%%%%%%%%%%%%

Ds 						= (1/GAMMA)*(MUg/SIGMA2g)*((M-N-1)/(M-2));
ECALUTILDE_WQ3IHATCsDs	= (TAU2/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2/(PSI2+(N-1)/M));

Css 						= Cs;
Dss 						= (1/GAMMA)*MUg*(((M-N-1)*(M-N-4))/(M*(M-2)));

ECALUTILDE_WQ3IIHATCssDss 	= (TAU2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2)/(PSI2+(N-1)/M)); 

% Three fund rule revisited

CKZ3s 							= (((M-N-1)*(M-N-4))/(M*(M-2)))*((PSI2)/(PSI2+N/M));
DKZ3s 							= (1/GAMMA)*(((M-N-1)*(M-N-4))/(M*(M-2)))*((N/M)/(PSI2+N/M))*MUg;

ECALUTILDE_WKZ3HATCKZ3sDKZ3s	= (THETA2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2+(THETA2/PSI2)*(N/M)));

