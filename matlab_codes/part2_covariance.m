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
T=50;
nu=sqrt(sigma2_g)*(Sigmah\onevec);
eta=(1/sqrt(psi2))*(Sigmah\(mu-mu_g*onevec));
% optimal performance
%Eteo=(theta2-psi2)*T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4))
% simulation starts
% number of simulations
n_sim=200000;
E=0;
F=0;
E2=0;
    for i=1:n_sim
        % simulation orthogonal method
        ret_sim=mvnrnd(zeros(N-1,1),eye(N-1),T);
        z=mvnrnd(zeros(N-1,1),eye(N-1)/(T))';
        u=chi2rnd(T-N)/T;
        A=cov(ret_sim,1);
        Ah=sqrtm(A);
        M1=(evec'*(Ah\z));
        %F1=((1/(u*u))*(M1)+(1/u)*(evec'*((A*Ah)\z))+(1/(u*u))*(M1)*(z'*(A\z)));
        %F1=M1;
        F1=1;
        %F1=(T*T/((T-N-2)*(T-N-4)))*(M1*M1)+(T/(T-N-2))*(evec'*((A*Ah)\z))*M1+(T*T/((T-N-2)*(T-N-4)))*(M1*M1)*(z'*(A\z));
        %F1=(1/(u*u))*(M1*M1)+(1/(u))*(evec'*((A*Ah)\z))*M1+(1/(u*u))*(M1*M1)*(z'*(A\z));
        F11=mu_g*((1/(u*u))*(M1)+(1/(u))*(evec'*((A*Ah)\z))+(1/(u*u))*(M1)*(z'*(A\z)));
        %F2=sigma2_g*u;
        %F3=(mu_g/(u*sigma2_g))+(sqrt(psi2)/sqrt(sigma2_g))*(evec'*(A\z))*(1/u);
        %F=F+(u*sigma2_g)*F1*F3;
        F=F+0*sqrt(psi2)*sqrt(sigma2_g)*F1+F11;
        % simulation normal method
        ret_y=mvnrnd(mu,Sigma,T);
        mue=mean(ret_y)';
        Sigmae=cov(ret_y,1);
        X=Sigmae\onevec;
        YY=Sigmae\mu;
        Wpinv=Sigmah*(inv(Sigmae))*Sigmah;
        P=(nu'*Wpinv*Wpinv*eta);
        V=onevec'*X;
        B=onevec'*YY;
        %Q=(sqrt(sigma2_g)*(1/sqrt(psi2))*X'*Sigma*(Sigmae\(mu-mu_g*onevec)))^2;
        E=E+P*B/V;
        %E2=E2+Q;
        %F=F+(YY'*Sigma*X)*(B/V);
    end
    Esperadohoja=F/n_sim
    Esperadoform=E/n_sim
    %Teo=T*(T-2)/((T-N)*(T-N-1)*(T-N-3))
    %Teo=((T-N)*(N-1)-2*(N-2)+2*(T-2))/((T-N)*(T-N-1)*(T-N-3))
    %Teo2=(N+1)/((T-N-1)*(T-N-3))
    A1=T*T/((T-N-2)*(T-N-4));
    Teo1=sqrt(psi2)*sqrt(sigma2_g)*((A1*(1/(T-N-1)))+(T/(T-N-2))*(T*(T-2)/((T-N)*(T-N-1)*(T-N-3)))+A1*((N+1)/((T-N-1)*(T-N-3))))
    Teo2=sqrt(psi2)*sqrt(sigma2_g)*(2*T*T*(T-2)/((T-N)*(T-N-1)*(T-N-3)*(T-N-4)))




