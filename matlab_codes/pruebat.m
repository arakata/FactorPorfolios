clear all
clc
T=200
N=50;
Teo1=N*T*(T-2)*(T-N-3)/((T-N-1)*(T-N-2)*(T-N-4)*(T-N-3));
Teo2=-(N-1)*T*(T-2)/((T-N)*(T-N-1)*(T-N-2)*(T-N-3));
Teo3=-T*(T-2)*(T-4)/((T-N-1)*(T-N-2)*(T-N-3)*(T-N-4));
Teo=Teo1+Teo2+Teo3
Teover=T*(T-2)*(N-1)/((T-N)*(T-N-2)*(T-N-3))