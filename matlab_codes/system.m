% code to estimate out of sample
clear all
clc
% we read the true mean vector and covariance matrix
load mu.mat
load Sigma.mat
%invSigmah=inv(Sigmah);
N=10;
T=150
onevec=ones(N,1);
evec=[1; zeros(N-2,1)];
% sharpe ratio tangent% number of observations
theta2=mu'*(Sigma\mu);
% sharpe ratio hedge
sigma2_g=((onevec'*(Sigma\onevec))^(-1));
wg=sigma2_g*(Sigma\onevec);
mu_g=mu'*wg;
psi2=((mu-mu_g*onevec)'*(Sigma\(mu-mu_g*onevec)));
const1=onevec'*(Sigma\mu);
% matrix input
a11=T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4))*(1/sigma2_g);
a12=T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4))*const1;
a13=0;
a22=T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4))*theta2+N*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4));
a23=T*T*(T-2)/((T-N)*(T-N-2)*(T-N-3))*psi2+(N-1)*T*(T-2)/((T-N)*(T-N-2)*(T-N-3));
a33=T*T*(T-2)/((T-N)*(T-N-1)*(T-N-3))*psi2+(N-1)*T*(T-2)/((T-N)*(T-N-1)*(T-N-3));
a14=0*(1/sigma2_g);
a24=0*const1+1;
a34=0;
a41=0;
a42=1;
a43=0;
a44=0;
b1=T*const1/(T-N-2);
b2=T*theta2/(T-N-2);
b3=T*psi2/(T-N-1);
b4=((T-N-1)*(T-N-4)/(T*(T-2)))*(psi2/(psi2+(N/T)));
A=[a11 a12 a13;
    a12 a22 a23;
    a13 a23 a33]
A2=[a11 a12 a13 a14;
    a12 a22 a23 a24;
    a13 a23 a33 a34;
    a41 a42 a43 a44]
b=[b1;b2;b3]
b2=[b;b4]
ct=(T-N)*(T-N-3)/(T*(T-2))*(psi2/(psi2+(N-1)/T))
dt=(T-N-1)*(T-N-4)*mu_g/(T*(T-2))
ckzt=((T-N-1)*(T-N-4)/(T*(T-2)))*(psi2/(psi2+(N/T)))
dkz=((T-N-1)*(T-N-4)/(T*(T-2)))*((N/T)/(psi2+(N/T)))*mu_g
et=(T-N)*(T-N-3)/(T*(T-2))*(psi2/(psi2+(N-1)/T))*(1-(T-N-1)*(T-N-1)*(T-N-4)*(psi2+(N-1)/T)/((T-N)*(T-N-2)*(T-N-3)*(psi2+(N)/T)))
et=ct-((T-N-1)*(T-N-1)*(T-N-4)/(T*(T-2)*(T-N-2)))*(psi2/(psi2+(N/T)));
True=[dkz;ckzt;et]
