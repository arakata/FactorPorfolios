% data kan and zhou
clear all
format long
clc
% N=25 size portfolios
ret=xlsread('datakanzhou.xlsx',5, 'B2:Z1021');
% N=10 portfolios
%ret=xlsread('datakanzhou.xlsx',7, 'B2:K1021');
% N=32 portfolios
%ret=xlsread('datakanzhou.xlsx',8, 'B2:AG643');
% reading the risk-free rate
rf=xlsread('datakanzhou.xlsx',6, 'B2:B1021');
% fixing the size
N=25;
for k=1:N
    rf_mat(:,k)=(1/100)*rf;
end
ret=(1/100)*ret-rf_mat;
mu=mean(ret)';
save mu
Sigma=cov(ret,1);
save Sigma
invSigma=inv(Sigma);
onevec=ones(N,1);
sigma2_g=((onevec'*(Sigma\onevec))^(-1))
R=invSigma-sigma2_g*((Sigma\onevec)*(Sigma\onevec)');
theta=sqrt(mu'*(Sigma\mu))
wg=sigma2_g*(Sigma\onevec); %min variance portfolio
mu_g=mu'*wg
psi=sqrt((mu-mu_g*onevec)'*(Sigma\(mu-mu_g*onevec)))
psi=sqrt(mu'*R*mu)