clear all
clc
clear all
clc
T=9;
N=3;
mu=[0.1;0.2;0.05];
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
C1=inv(Sigma);
C2=((T-1)/(T-N-2))*inv(Sigma);
ER=((T-1)/(T-N-1))*(C1-C1*ones(3,1)*ones(1,3)*C1/(ones(1,3)*C1*ones(3,1)));
Teo=-mu'*ER*mu+mu'*C2*mu
theta=mu'*C1*mu;
psi=theta-((mu'*C1*ones(3,1))^2)/(ones(1,3)*C1*ones(3,1));
Teo2=((T-1)/(T-N-2))*(theta-psi+psi/(T-N-1))
SSim=0;
for k=1:100000
    X=mvnrnd(zeros(3,1),Sigma,T);
    C=inv(cov(X));
    mus=mean(X);
    w_mins=C*ones(3,1)/(ones(1,3)*C*ones(3,1));
    mu_gs=mu'*w_mins;
    val_sim=mu_gs*ones(1,3)*C*mu;
    SSim=SSim+val_sim;
end
Savg=(1/k)*SSim

