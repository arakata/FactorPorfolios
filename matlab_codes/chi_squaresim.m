clear all
clc
N=8;
T=600;
nsim=40000;
Xsum=0;
Xsum2=0;
for i=1:nsim
    u=chi2rnd(T-N)/T;
    X=1/u;
    X2=1/(u*u);
    Xsum=X+Xsum;
    Xsum2=X2+Xsum2;
end
Xbar1=Xsum/nsim
Xteo=T/(T-N-2)
Xbar2=Xsum2/nsim
Xteo2=T*T/((T-N-2)*(T-N-4))