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
evec=[1; zeros(N-2,1)];
% sharpe ratio tangent% number of observations
theta2=mu'*(Sigma\mu);
% sharpe ratio hedge
sigma2_g=((onevec'*(Sigma\onevec))^(-1));
wg=sigma2_g*(Sigma\onevec);
mu_g=mu'*wg;
psi2=((mu-mu_g*onevec)'*(Sigma\(mu-mu_g*onevec)));
T=500;
% optimal performance
n_sim=5000;
E1=0;
E2=0;
E3=0;
    for i=1:n_sim
        % simulation normal method
        ret_y=mvnrnd(mu,Sigma,T);
        mue=mean(ret_y)';
        Sigmae=cov(ret_y,1);
        X=Sigmae\onevec;
        Y=Sigmae\mue;
        YY=Sigmae\mu;
        V1=onevec'*X;
        mue_g=onevec'*Y/V1;
        E3=E3+mue_g*X'*Sigma*Y;
        Z1=X'*Sigma*YY;
        B1=onevec'*YY;
        E1=E1+Z1*B1/V1;
        E2=E2+X'*Sigma*inv(Sigmae)*Sigma*X/(V1);
    end
    Esperadoform1=E1/n_sim
    Teo1=(T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4)))*(theta2-psi2*(T-N-1)*(T-N-4)/((T-N)*(T-N-3)))
    Esperadoform2=E2/n_sim
    %Teo2=(T*(T-2)*(N-1)/((T-N)*(T-N-1)*(T-N-3)))*(T/(T-N-2))+(T*T/((T-N-2)*(T-N-4)))*(((T+N-3)/(T-N-1))+((N-1)*(N+1)/((T-N-1)*(T-N-3))))
    Teo2=(T*T*(T-2)*(N-1)/((T-N)*(T-N-1)*(T-N-3)*(T-N-2)))+(T*T*(T-2)*(T-4)/((T-N-1)*(T-N-2)*(T-N-3)*(T-N-4)))
    Esperadofinal=E3/n_sim
    Teofinal=Teo1+(1/T)*Teo2

