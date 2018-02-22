clear all
clc
clear all
clc
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
A=[1 1 1;0.5 -0.5 0];
b=[1;0];
C1=inv(Sigma);
A1=inv(A*C1*A');
F1=C1-C1*A'*A1*A*C1;
trazaF1=trace(A'*A1*A*C1)
wa=C1*A'*A1*b
K=40000
w_acum=zeros(K,3);
T=5
for k=1:K
    X=mvnrnd(zeros(3,1),Sigma,T);
    C1sim=inv(cov(X,1));
    A1sim=inv(A*C1sim*A');
    was=C1sim*A'*A1sim*b;
    w_acum(k,:)=was';
end
S_sim=cov(w_acum)
S_teo=(1/(T-3+2-2))*((C1-C1*A'*(A1)*A*C1)*(b'*A1*b))
%Steo=((8-1)/(8-3-1))*Teo
