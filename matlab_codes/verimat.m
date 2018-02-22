clear all
clc
clear all
clc
mu=[0.01;0.02;0.03]
Sigma=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
e=[1;1;1];
C1=inv(Sigma)
wg=(1/(e'*C1*e))*(C1*e)
mug=mu'*wg;
sigma_g=sqrt(wg'*Sigma*wg);
SR2g=(mug*mug)/(sigma_g^2)
R=C1-(1/(e'*C1*e))*(C1*e*e'*C1)
alpha=mu-mug*e
Prueba=R-(1/(100+alpha'*R*alpha))*R*alpha*alpha'*R
