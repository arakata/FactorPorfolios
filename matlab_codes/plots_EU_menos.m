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
EUmvg=100*(mu_g-0.5*gamma*sigma_g^2);
z=1;
for T=60:10:480
    %k1=((T-1)/(T-N-2))*(2-(T-1)*(T-2)/((T-N-1)*(T-N-4)));
    %EU_plug=100*(k1*theta^2/(2*gamma)-(1/(2*gamma*T))*(N*(T-1)*(T-1)*(T-2))/((T-N-1)*(T-N-2)*(T-N-4)));
    EU_opt2f=100*(theta^2/(2*gamma))*(theta^2/(theta^2+(N/T)))*((T-N-1)*(T-N-4)/((T-2)*(T-N-2)));
    %k1_g=((T-1)/(T-N-1))*(2-(T-1)*(T-2)/((T-N)*(T-N-3)));
    %k2_g=((N-1)*(T-1)*(T-1)*(T-2))/(T*(T-N)*(T-N-1)*(T-N-3));
    %EUg_plug=100*(mu_g-0.5*gamma*(1+(N-1)/(T-N-1))*sigma_g^2+k1_g*psi^2/(2*gamma)-(1/(2*gamma))*k2_g);
    %EUg_optg=100*((tau*tau*0.5/gamma)*(T-N-1)/(T-2)+(k2_g*psi^2/(2*gamma)-(1/(2*gamma*T))*);
    EUg_opth=100*(mu_g-0.5*gamma*(1+(N-1)/(T-N-1))*sigma_g^2+(psi^2/(2*gamma))*(psi^2/(theta^2+((N-1)/T)))*((T-N)*(T-N-3)/((T-2)*(T-N-1))));
    data(z,:)=[T EU_opt EUg_opt EU_opt2f EUg_opth EUmvg];
    z=z+1;
end
createfigure_menos(data(:,1), data(:,2:end))



    