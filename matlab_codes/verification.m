clear all
clc
clear all
clc
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
C1=inv(Sigma);
Teo=C1-C1*ones(3,1)*ones(1,3)*C1/(ones(1,3)*C1*ones(3,1));
SSim=zeros(3,3);
for k=1:100000
    X=mvnrnd(zeros(3,1),Sigma,8);
    C=inv(cov(X));
    A=C-C*ones(3,1)*ones(1,3)*C/(ones(1,3)*C*ones(3,1));
    SSim=SSim+A;
end
Savg=(1/k)*SSim
Steo=((8-1)/(8-3-1))*Teo
