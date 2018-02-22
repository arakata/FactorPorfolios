%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USING A ZERO-INVESTMENT PORTFOLIO TO REDUCE ESTIMATION 
%RISK IN MEAN-VARIANCE PORTFOLIO OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DRAFT: 	LUIS CARLOS CHAVEZ-BEDOYA
% DATE:  	JULY 24TH, 2017
% CODE:  	JAN 2018
% EDITED:	FRANCISCO ROSALES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREAMBLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
format long
clc

% input
% file:       	datakanzhou.xlsx
% details     	pages 5 and 6
% ret_mat:		matrix with returns
%rf_vec:		vector with risk free asset info    

ret_mat 	= xlsread('datakanzhou.xlsx',5, 'B2:Z1021');
rf_vec 		= xlsread('datakanzhou.xlsx',6, 'B2:B1021');

% constants
[M,N] = size(ret_mat); 		% number of possible assets
GAMMA = 3;		% risk attitude

for k = 1:N
    rf_mat(:,k) = (1/100)*rf_vec;
end
ret_mat = (1/100)*ret_mat-rf_mat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASIC STATS & INSTRUMENTAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stats
% we assume these are the real values of the process

mu_vec		= mean(ret_mat)';
sigma_mat	= cov(ret_mat,1);

save mu_vec
save sigma_mat

% auxiliary

isigma_mat	= inv(sigma_mat);
one_vec		= ones(N,1);
r_mat 		= isigma_mat-(1/((one_vec')*(sigma_mat\one_vec)))*(sigma_mat\one_vec)*((sigma_mat\one_vec)');


% theoretical results

THETA2 		= (mu_vec')*(sigma_mat\mu_vec);
TAU2 		= (((mu_vec')*(sigma_mat\one_vec))^2)/((one_vec')*(sigma_mat\one_vec));
PSI2 		= (mu_vec')*(r_mat*mu_vec);

wt_vec 	= (1/GAMMA)*isigma_mat*mu_vec;
wg_vec	= (1/((one_vec')*(sigma_mat\one_vec)))*(sigma_mat\one_vec);
wh_vec 	= (1/GAMMA)*r_mat*mu_vec;
wq_vec	= (1/(one_vec'*(sigma_mat\one_vec)))*(sigma_mat\one_vec)+(1/GAMMA)*r_mat*mu_vec;

MUg		= (wg_vec')*(mu_vec);
SIGMA2g	= (wg_vec')*(sigma_mat*wg_vec);

% Q: is it possible to do better than wt_vec? by selecting another mixture?

CONSTs		= (1/GAMMA)*(MUg/SIGMA2g);
UTILs 		= calU(CONSTs*wg_vec+wh_vec, mu_vec, sigma_mat, GAMMA);

const_vec	= linspace(0.1,5,50);
util_vec 	= 0*const_vec;

for i=1:length(const_vec)
	util_vec(i) 	= calU(const_vec(i)*wg_vec+wh_vec, mu_vec, sigma_mat, GAMMA);
end

plot([util_vec', UTILs*ones(length(const_vec),1)])

% ANS: no. it is not possible. CONSTs=1.6 (approx) is the best option.   

M=1000;

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
