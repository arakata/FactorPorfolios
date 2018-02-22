clear all
clc
% we introduce the parameters for Kan and Zhou (2007)
gamma=3;
theta=0.344;
psi=0.267;
tau=sqrt(theta^2-psi^2);
mu_g=0.00889;
sigma_g=sqrt(mu_g*mu_g/(theta^2-psi^2));
N=25;
EU_opt=theta^2/(2*gamma)*100;
EUg_opt=100*(mu_g-0.5*gamma*sigma_g^2+(1/(2*gamma))*psi^2);
%EUmvg=100*(mu_g-0.5*gamma*sigma_g^2);
z=1;
for T=60:10:480
    f1=1-(N/T)/(theta^2+(theta^2/(psi^2))*(N/T));
    f2=((T-N-1)*(T-N-4)/((T-2)*(T-N-2)));
    EU_opt3=100*((theta^2/(2*gamma))*f1*f2);
    EUg_opt3=100*((tau^2/(2*gamma))*(((T-N-1)*(T-N-4))/((T-2)*(T-N-2)))+(psi^2/(2*gamma))*(psi^2/(psi^2+((N-1)/T)))*((T-N)*(T-N-3)/((T-2)*(T-N-1))));
    data(z,:)=[T EU_opt EUg_opt EU_opt3 EUg_opt3 ];
    z=z+1;
end
createfigure_three(data(:,1), data(:,2:end))



    