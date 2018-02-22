% code to estimate out of sample
clear all
clc
% we read the true mean vector and covariance matrix
load mu.mat
load Sigma.mat
% we define certain variables
N=25;
onevec=ones(N,1);
gamma=3;
% sharpe ratio tangent% number of observations
theta2=mu'*(Sigma\mu);
% sharpe ratio hedge
sigma2_g=((onevec'*(Sigma\onevec))^(-1));
wg=sigma2_g*(Sigma\onevec);
mu_g=mu'*wg;
psi2=((mu-mu_g*onevec)'*(Sigma\(mu-mu_g*onevec)));
T=200;
% optimal kz3 performance
EUkz3=100*(theta2/(2*gamma))*((T-N-1)*(T-N-4)/((T-2)*(T-N-2)))*(1-(N/T)/(theta2+(theta2/(psi2))*(N/T)));
%cteo=((T-N-1)*(T-N-4)/(T*(T-2)))*(psi2/(psi2+(N/T)));
%EUkz31t=100*((cteo*theta2/gamma)*(T/(T-N-2))-(cteo*cteo*0.5/gamma)*(theta2+(N/T))*((T*T*(T-2)/((T-N-1)*(T-N-2)*(T-N-4)))));
EUkz31=100*(1/(2*gamma))*((T-N-1)*(T-N-4)/((T-2)*(T-N-2)))*(psi2/(psi2+((N)/T)))*(2*theta2-((psi2)/(psi2+(N/T)))*(theta2+(N/T)));
dteo=(1/gamma)*((T-N-1)*(T-N-4)/(T*(T-2)))*((N/T)/(psi2+(N/T)))*mu_g;
EUkz32=100*((T/(T-N-2))*dteo*(onevec'*(Sigma\mu))-gamma*0.5*dteo*dteo*T*T*(T-2)*(onevec'*(Sigma\onevec))/((T-N-1)*(T-N-2)*(T-N-4)));
EUkz3teo=[EUkz31 EUkz32 (EUkz3-EUkz31-EUkz32) EUkz3]
% optimal 3fund performance
EUq3=100*(psi2/(2*gamma))*((T-N)*(T-N-3)/((T-2)*(T-N-1)))*(psi2/(psi2+((N-1)/T)))+100*((theta2-psi2)/(2*gamma))*((T-N-1)*(T-N-4)/((T-2)*(T-N-2)));
EUq31=100*(psi2/(2*gamma))*((T-N)*(T-N-3)/((T-2)*(T-N-1)))*(psi2/(psi2+((N-1)/T)));
EUq32=100*((theta2-psi2)/(2*gamma))*((T-N-1)*(T-N-4)/((T-2)*(T-N-2)));
EUq3teo=[EUq31 EUq32 (EUq3-EUq31-EUq32) EUq3]
% optimal 2 fund full invested
EUq2=100*(psi2/(2*gamma))*((T-N)*(T-N-3)/((T-2)*(T-N-1)))*(psi2/(psi2+((N-1)/T)))+100*(mu_g-0.5*gamma*((T-2)/(T-N-1))*sigma2_g)
% simulation starts
% number of simulations
n_sim=2000;
Ukz3_acum=0;
Uq3_acum=0;
Uq4_acum=0;
    for i=1:n_sim
        ret_sim=mvnrnd(mu,Sigma,T);
        mue=mean(ret_sim)';
        Sigmae=cov(ret_sim,1);
        X=Sigmae\onevec;
        Y=Sigmae\mue;
        sigma2_ge=((onevec'*X)^(-1));
        Re=inv(Sigmae)-sigma2_ge*(X*X');
        Z=Re*mue;
        wge=sigma2_ge*X;
        mue_g=mue'*wge;
        % estimate for theta2
        thetae2=mue'*Y;
        betatheta2=betainc((thetae2/(1+thetae2)),(N)/2,(T-N)/2)*beta((N)/2,(T-N)/2);
        thetaa2=(((T-N-2)*thetae2-(N))/T)+(2/(T*betatheta2))*((thetae2)^((N)/2))*((1+thetae2)^(-(T-2)/2));
        % estimate for psi2 
        psie2=(mue-mue_g*onevec)'*(Sigmae\(mue-mue_g*onevec));
        betapsi2=betainc((psie2/(1+psie2)),(N-1)/2,(T-N+1)/2)*beta((N-1)/2,(T-N+1)/2);
        psia2=(((T-N-1)*psie2-(N-1))/T)+(2/(T*betapsi2))*((psie2)^((N-1)/2))*((1+psie2)^(-(T-2)/2));
        % portfolio rules
        whkz3=((T-N-1)*(T-N-4)/(T*(T-2)*gamma))*((psia2/((N/T)+psia2))*Y+((N/T)/(psia2+(N/T)))*mue_g*X);
        whq3=((T-N)*(T-N-3)/(T*(T-2)*gamma))*((psia2/(((N-1)/T)+psia2)))*Z+((T-N-1)*(T-N-4)/(T*(T-2)*gamma))*mue_g*X;
        %whq4=whkz3+(((T-N)*(T-N-3)/(T*(T-2)*gamma))*(psia2/(((N-1)/T)+psia2))-((T-N-1)*(T-N-1)*(T-N-4)/(T*(T-2)*(T-N-2)*gamma))*(psia2/((N/T)+psia2)))*Z;
        %whq3=((T-N)*(T-N
        
        -3)/(T*(T-2)*gamma))*(psia2/(((N-1)/T)+psia2))*Z+((T-N-1)*(T-N-2)/(T*(T-2)*gamma))*(thetaa2-psia2)*(1/mue_g)*wge;
        Ukz3_sim=whkz3'*mu-(gamma/2)*(whkz3'*Sigma*whkz3);
        Uq3_sim=whq3'*mu-(gamma/2)*(whq3'*Sigma*whq3);
        %Uq4_sim=whq4'*mu-(gamma/2)*(whq4'*Sigma*whq4);
        % accumulation
        Ukz3_acum=Ukz3_acum+Ukz3_sim*100;       
        Uq3_acum=Uq3_acum+Uq3_sim*100;
        %Uq4_acum=Uq4_acum+Uq4_sim*100;
end
Ukz3_est=Ukz3_acum/n_sim;
Uq3_est=Uq3_acum/n_sim;
%Uq4_est=Uq4_acum/n_sim;
%res=[Ukz3_est Uq3_est Uq4_est]
res=[Ukz3_est Uq3_est]

