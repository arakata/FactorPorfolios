clear all
clc
clear all
clc
mu=[0.01;0.02;0.03]
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
e=[1;1;1];
A=[1 1 1; 2 -1 3];
invSigma=inv(Sigma)
C=invSigma*A'*inv(A*invSigma*A')
c1=invSigma*A(1,:)'/(A(1,:)*invSigma*A(1,:)')
