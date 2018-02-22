clear all
clc
clear all
clc
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
A=[1 1 1];
b=[1];
C1=inv(Sigma);
A1=inv(A*C1*A');
wa=C1*A'*A1*b
w_acum=zeros(3,1);
for k=1:10000
    X=mvnrnd(zeros(3,1),Sigma,10);
    C1sim=inv(cov(X,1));
    A1sim=inv(A*C1sim*A');
    was=C1sim*A'*A1sim*b;
    w_acum=w_acum+was;
end
wavg=(1/k)*w_acum
%Steo=((8-1)/(8-3-1))*Teo
