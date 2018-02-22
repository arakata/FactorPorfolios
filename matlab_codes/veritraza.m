clear all
clc
clear all
clc
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
e=[1;0;0];
A=[1 1 1;0.5 -0.5 0];
%b=[1;0];
C1=inv(Sigma);
A1=inv(A*C1*A');
F1=C1-C1*A'*A1*A*C1;
trazaF1=trace(F1*e*e'*F1*Sigma)
F1(2,2)
