clear all
clc
clear all
clc
gamma=12;
S=[0.02778 0.00387 0.00021;
    0.00387 0.01112 -0.00020;
    0.00021 -0.00020 0.0015];
mu=[0.03; 0.01; 0.01];
Sinv=inv(S);
R=Sinv-Sinv*ones(3,1)*ones(1,3)*Sinv/(ones(1,3)*Sinv*ones(3,1));
wg=(1/(ones(3,1)'*Sinv*ones(3,1)))*Sinv*ones(3,1);
wg1=Sinv*ones(3,1);
wt=(1/(ones(3,1)'*Sinv*mu))*Sinv*mu;
wt1=Sinv*mu;
wh=R*mu;
wopt=(1/gamma)*Sinv*mu;
T=1;
Q=1;
for t=0:0.005:0.02
    for v=-3:0.05:1.5
        wp=(1/gamma)*(t*wg1+v*wt1);
        mu_p=wp'*mu;
        sig_p=sqrt(wp'*S*wp);
        data2(Q,:)=[sig_p^2 mu_p];
        Q=Q+1;
    end
end
for alpha=-2:0.1:3
    we=alpha*wg+(1-alpha)*wt;
    ws=wg+alpha*0.15*wh;
    ws1=(1/gamma)*(mu'*wg/(wg'*S*wg))*wg+alpha*0.15*wh;
    mu_e=we'*mu;
    sig_e=sqrt(we'*S*we);
    mu_s=ws'*mu;
    sig_s=sqrt(ws'*S*ws);
    mu_s1=ws1'*mu;
    sig_s1=sqrt(ws1'*S*ws1);
    data(T,:)=[sig_e^2 mu_e sig_s^2 mu_s sig_s1^2 mu_s1];
    T=T+1;
end
data
plot(data(:,1),data(:,2))
hold on
plot(data2(:,1),data2(:,2),'k*')
hold on
plot(data(:,3),data(:,4),'r.')
hold on
plot(data(:,5),data(:,6),'c.')
hold on
grid on
plot(sqrt(wg'*S*wg)^2,wg'*mu,'*')
plot(sqrt(wt'*S*wt)^2,wt'*mu,'o')
plot(sqrt(wh'*S*wh)^2,wh'*mu,'sq')
plot(sqrt(wopt'*S*wopt)^2,wopt'*mu,'r*')
plot(0,0,'^')








