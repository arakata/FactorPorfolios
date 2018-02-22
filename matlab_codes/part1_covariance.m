% code to estimate out of sample
clear all
clc
% we read the true mean vector and covariance matrix
load mu.mat
load Sigma.mat
% we define certain variables
Sigmah=sqrtm(Sigma);
%invSigmah=inv(Sigmah);
N=10;
onevec=ones(N,1);
% sharpe ratio tangent% number of observations
theta2=mu'*(Sigma\mu);
% sharpe ratio hedge
sigma2_g=((onevec'*(Sigma\onevec))^(-1));
wg=sigma2_g*(Sigma\onevec);
mu_g=mu'*wg;
psi2=((mu-mu_g*onevec)'*(Sigma\(mu-mu_g*onevec)));
T=300;
nu=sqrt(sigma2_g)*(Sigmah\onevec);
eta=(1/sqrt(psi2))*(Sigmah\(mu-mu_g*onevec));
% optimal performance
%Eteo=(theta2-psi2)*T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4))
% simulation starts
% number of simulations
n_sim=1;
V=0;
B=0;
W=0;
    for i=1:n_sim
        ret_sim=mvnrnd(mu,Sigma,T);
        mue=mean(ret_sim)';
        Sigmae=cov(ret_sim,1);
        X=Sigmae\onevec;
        Y=Sigmae\mue;
        YY=Sigmae\mu;
        sigma2_ge=((onevec'*X)^(-1));
        Re=inv(Sigmae)-sigma2_ge*(X*X');
        Z=Re*mue;
        wge=sigma2_ge*X;
        mue_g=mue'*wge;
        Wpinv=Sigmah*(inv(Sigmae))*Sigmah;
        P=(sqrt(psi2)/sqrt(sigma2_g))*(nu'*Wpinv*Wpinv*eta);
        V=mu_g*X'*Sigma*X;
        %V=V+mu_g*X'*Sigma*X*(onevec'*Y)/(onevec'*X);
        B=X'*Sigma*YY;
    end
P
B-V
X'*Sigma*X
(1/sigma2_g)*nu'*Wpinv*Wpinv*nu



