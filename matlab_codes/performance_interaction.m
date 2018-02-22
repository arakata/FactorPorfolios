clear all
clc
gamma=3;
theta=0.341799088677402;
psi=0.273168755272346;
%psi=0.4;
tau=sqrt(theta^2-psi^2);
mu_g=0.008455877392722;
sigma_g=sqrt(mu_g*mu_g/(theta^2-psi^2));
N=25;
EUg_opt=mu_g-0.5*gamma*sigma_g^2+(1/(2*gamma))*psi^2
z=1;
for T=60:0.5:480
    k1_g=(T/(T-N-1))*(2-T*(T-2)/((T-N)*(T-N-3)));
    k2_g=((N-1)*T*(T-2))/((T-N)*(T-N-1)*(T-N-3));
    EUg_u=mu_g-0.5*gamma*sigma_g^2+(1/(2*gamma))*(psi^2)-(1/(2*gamma))*((N-1)/T);
    EUg_s=mu_g-0.5*gamma*(1+(N-1)/(T-N-1))*sigma_g^2+k1_g*psi^2/(2*gamma);
    EUg=mu_g-0.5*gamma*(1+(N-1)/(T-N-1))*sigma_g^2+k1_g*(psi^2)/(2*gamma)-(1/(2*gamma))*k2_g;
    
    Lmvg=mu_g-0.5*gamma*((T-2)/(T-N-1))*sigma_g^2;
    Lhdg=(k1_g)*(psi^2)/(2*gamma)-(1/(2*gamma))*k2_g;
    Lu=(1-EUg_u/EUg_opt)*100;
    Lg=(1-EUg_s/EUg_opt)*100;
    Lt=(1-EUg/EUg_opt)*100;
    Li=Lt-Lu-Lg;
    %data2(z,:)=[T (1-(Lmvg/(mu_g-0.5*gamma*sigma_g^2)))*100 (1-(Lhdg/((1/(2*gamma))*(psi^2))))*100];
    data(z,:)=[T Lu Lg Li];
    z=z+1;
end
figure(2)
area(data(:,1),data(:,2:end))
axis([60 480 0 2000])
h=legend('$L{\mu}$ ($\Sigma$ known)','$L_{\Sigma}$ ($\mu$ known)','Interaction')
set(h,'Interpreter','latex')
% Create xlabel
xlabel('T (number of periods)','FontSize',16.5);
% Create ylabel
ylabel('% Loss','FontSize',16.5);
% Create title
title('Percentage loss of expected out-of-sample performance Q','FontWeight','bold',...
    'FontSize',18);
hold on
grid on
%figure(3)
%grid on
%hold on
%area(data2)



    