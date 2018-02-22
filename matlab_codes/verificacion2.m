clear all
clc
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
mu=[0.1;0.2;0.05];
C1=inv(Sigma);
w_min=C1*ones(3,1)/(ones(1,3)*C1*ones(3,1));
mu_g=mu'*w_min;
Teo=mu_g*mu'*C1*ones(3,1);
SSim=zeros(1,1);
for k=1:100000
    X=mvnrnd(zeros(3,1),Sigma,8);
    mus=mean(X);
    C=inv(cov(X));
    w_mins=C1*ones(3,1)/(ones(1,3)*C1*ones(3,1));
    mu_gs=mus*w_mins;
    val_sim=mu_gs*mus*C*ones(3,1);
    SSim=SSim+val_sim;
end
SSim=(1/k)*(SSim)
Teo=Teo
