clear all
clc
gamma=3;
theta=0.341799088677402;
%gamma=3;
%theta=0.4;
N=25;
EUg_opt=(1/(2*gamma))*(theta^2)
z=1;
for T=60:1:480
    k1_g=(T/(T-N-2))*(2-(T*(T-2)/((T-N-1)*(T-N-4))));
    k2_g=(T*T*(T-2))/((T-N-1)*(T-N-2)*(T-N-4));
    EUg_u=(1/(2*gamma))*theta^2-(1/(2*gamma))*(N)/T;
    EUg_s=k1_g*theta^2/(2*gamma);
    EUg=k1_g*theta^2/(2*gamma)-(N/(2*gamma*T))*k2_g;
    Lu=(1-EUg_u/EUg_opt)*100;
    Lg=(1-EUg_s/EUg_opt)*100;
    Lt=(1-EUg/EUg_opt)*100;
    Li=Lt-Lu-Lg;
    data(z,:)=[T Lu Lg Li];
    z=z+1;
end
figure(1)
area(data(:,1),data(:,2:end))
axis([60 480 0 2000])
h=legend('$\hat{\mu}$ ($\Sigma$ known)','$\hat{\Sigma}$ ($\mu$ known)','$\hat{\mu}$ and $\hat{\Sigma}$ - Interaction')
set(h,'Interpreter','latex')
% Create xlabel
xlabel('T (number of periods)','FontSize',16.5);
% Create ylabel
ylabel('% Loss','FontSize',16.5);
% Create title
title('Percentage loss of expected out-of-sample performance T','FontWeight','bold',...
    'FontSize',18);
hold on
grid on