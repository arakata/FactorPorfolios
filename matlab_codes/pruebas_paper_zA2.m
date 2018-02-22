clear all 
clc
N=20;
T=200;
nsim=20000;
Xsum=0;
Ysum=0;
for i=1:nsim
    ret_sim=mvnrnd(zeros(N-1,1),eye(N-1),T);
    ret_y=mvnrnd(zeros(N,1),eye(N),T);
    z=mvnrnd(zeros(N-1,1),eye(N-1)/(T))';
    u=chi2rnd(T-N)/T;
    A=cov(ret_sim,1);
    Sigmae=cov(ret_y,1);
    X=(1/(u*u))*((1+z'*(A\z))+(z'*(A\z))+(z'*(A\z))^2);
    X=X+(1/u)*(z'*((A*A)\z));
    f=Sigmae\ones(N,1);
    invSigmae=inv(Sigmae);
    Y=ones(1,N)*invSigmae*invSigmae*invSigmae*ones(N,1);
    Y=(1/(f'*ones(N,1)))*Y;
    Xsum=X+Xsum;
    Ysum=Y+Ysum;
end
Xbar=Xsum/nsim
Ybar=Ysum/nsim
Xteo=(T*(T-2)*(N-1)/((T-N)*(T-N-1)*(T-N-3)))*(T/(T-N-2));
Xteo=Xteo+(T*T/((T-N-2)*(T-N-4)))*(((T+N-3)/(T-N-1))+((N-1)*(N+1)/((T-N-1)*(T-N-3))))
Xteo2=(1/((T-N)*(T-N-1)*(T-N-2)*(T-N-3)*(T-N-4)))*T*T*(T-2)*(T*T-5*T-N*N+N+4)