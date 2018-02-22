KK=0;
for T = [50 80 120 200 500 1000]
KK = KK+1;
    display(T)
    NS 	= 20000;

    % Assume: Parameters Unknown

%    tfr1_vec    = ones(NS,6);
%    tfr2_vec    = ones(NS,6);
%    tfrr_vec    = ones(NS,6);

    for i = 1:NS  
       % theoretical parameters
       TheoreticalModule

	   % Set of Returns
	   retsim_mat 		= mvnrnd(mu_vec,sigma_mat,T);

       	muhat_vec 		= mean(retsim_mat)';
       	mutilde_vec 	= mean(retsim_mat)'; % note that muhat_vec and mutilde_vec are identical.
        sigmahat_mat 	= cov(retsim_mat,1);
        sigmahat_mat 	= mycov1(retsim_mat);
        sigmatilde_mat 	= mycov2(retsim_mat);
    
        % with hat
        x_vec 		 	= sigmahat_mat\one_vec;
        SIGMAHAT2g		= ((one_vec'*x_vec)^(-1));
        rhat_mat 		= inv(sigmahat_mat)-SIGMAHAT2g*(x_vec*x_vec');
        z_vec 			= rhat_mat*muhat_vec;
        wghat_vec 		= SIGMAHAT2g*x_vec;
        MUHATg 			= muhat_vec'*wghat_vec;
    	isigmahat_mat	= inv(sigmahat_mat);

        wthat_vec 	= (1/GAMMA)*isigmahat_mat*muhat_vec;
    	wghat_vec	= (1/((one_vec')*(sigmahat_mat\one_vec)))*(sigmahat_mat\one_vec);
    	whhat_vec 	= (1/GAMMA)*rhat_mat*muhat_vec;
    	wqhat_vec	= (1/(one_vec'*(sigmahat_mat\one_vec)))*(sigmahat_mat\one_vec)+(1/GAMMA)*rhat_mat*muhat_vec;

        THETA2HAT       = (muhat_vec')*(sigmahat_mat\muhat_vec);
        TAU2HAT         = (((muhat_vec')*(sigmahat_mat\one_vec))^2)/((one_vec')*(sigmahat_mat\one_vec));
        PSI2HAT         = (muhat_vec')*(rhat_mat*muhat_vec);

        % Three fund rules
        Cs                      = ((M-N)*(M-N-3)/(M*(M-2)))*(PSI2HAT/(PSI2HAT+(N-1)/M));

        Ds                      = (1/GAMMA)*(MUHATg/SIGMAHAT2g)*((M-N-1)/(M-2));
        ECALUTILDE_WQ3IHATCsDs  = (TAU2HAT/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2HAT/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2HAT/(PSI2HAT+(N-1)/M));
        ECALUTILDE_WQ3ICsDs     = (TAU2/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2/(PSI2+(N-1)/M));

        Css                         = Cs;
        Dss                         = (1/GAMMA)*MUHATg*(((M-N-1)*(M-N-4))/(M*(M-2)));

        ECALUTILDE_WQ3IIHATCssDss   = (TAU2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2HAT/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2HAT)/(PSI2HAT+(N-1)/M)); 
        ECALUTILDE_WQ3IICssDss      = (TAU2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2)/(PSI2+(N-1)/M)); 

        % Three fund rule revisited

        CKZ3s                           = (((M-N-1)*(M-N-4))/(M*(M-2)))*((PSI2HAT)/(PSI2HAT+N/M));
        DKZ3s                           = (1/GAMMA)*(((M-N-1)*(M-N-4))/(M*(M-2)))*((N/M)/(PSI2HAT+N/M))*MUHATg;

        ECALUTILDE_WKZ3HATCKZ3sDKZ3s    = (THETA2HAT/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2HAT+(THETA2HAT/PSI2HAT)*(N/M)));
        ECALUTILDE_WKZ3CKZ3sDKZ3s       = (THETA2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2+(THETA2/PSI2)*(N/M)));

        tfr1hat_vec(i,KK) = ECALUTILDE_WQ3IHATCsDs;
        tfr2hat_vec(i,KK) = ECALUTILDE_WQ3IIHATCssDss;
        tfrrhat_vec(i,KK) = ECALUTILDE_WKZ3HATCKZ3sDKZ3s;

        tfr1_vec(i,KK) = ECALUTILDE_WQ3ICsDs;
        tfr2_vec(i,KK) = ECALUTILDE_WQ3IICssDss;
        tfrr_vec(i,KK) = ECALUTILDE_WKZ3CKZ3sDKZ3s;

        %with tilde
        x_vec 				= sigmatilde_mat\one_vec;
        SIGMATILDE2g		= ((one_vec'*x_vec)^(-1));
        rtilde_mat 			= inv(sigmatilde_mat)-SIGMATILDE2g*(x_vec*x_vec');
        z_vec 				= rtilde_mat*mutilde_vec;
        wgtilde_vec 		= SIGMATILDE2g*x_vec;
        MUTILDEg 			= mutilde_vec'*wgtilde_vec;
    	isigmatilde_mat		= inv(sigmatilde_mat);

        wttilde_vec 	= (1/GAMMA)*isigmatilde_mat*mutilde_vec;
    	wgtilde_vec		= (1/((one_vec')*(sigmatilde_mat\one_vec)))*(sigmatilde_mat\one_vec);
    	whtilde_vec 	= (1/GAMMA)*rtilde_mat*mutilde_vec;
    	wqtilde_vec		= (1/(one_vec'*(sigmatilde_mat\one_vec)))*(sigmatilde_mat\one_vec)+(1/GAMMA)*rtilde_mat*mutilde_vec;

        THETA2TILDE       = (mutilde_vec')*(sigmatilde_mat\mutilde_vec);
        TAU2TILDE         = (((mutilde_vec')*(sigmatilde_mat\one_vec))^2)/((one_vec')*(sigmatilde_mat\one_vec));
        PSI2TILDE         = (mutilde_vec')*(rtilde_mat*mutilde_vec);

        % Three fund rules
        Cs                      = ((M-N)*(M-N-3)/(M*(M-2)))*(PSI2TILDE/(PSI2TILDE+(N-1)/M));
        Ds                      = (1/GAMMA)*(MUTILDEg/SIGMATILDE2g)*((M-N-1)/(M-2));
        ECALUTILDE_WQ3IHATCsDs  = (TAU2TILDE/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2TILDE/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2TILDE/(PSI2TILDE+(N-1)/M));
        ECALUTILDE_WQ3ICsDs     = (TAU2/(2*GAMMA))*((M-N-1)/(M-2))+(PSI2/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*(PSI2/(PSI2+(N-1)/M));

        Css                         = Cs;
        Dss                         = (1/GAMMA)*MUTILDEg*(((M-N-1)*(M-N-4))/(M*(M-2)));

        ECALUTILDE_WQ3IIHATCssDss   = (TAU2TILDE/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2TILDE/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2TILDE)/(PSI2TILDE+(N-1)/M)); 
        ECALUTILDE_WQ3IICssDss      = (TAU2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))+(PSI2/(2*GAMMA))*(((M-N)*(M-N-3))/((M-2)*(M-N-1)))*((PSI2)/(PSI2+(N-1)/M)); 

        % Three fund rule revisited

        CKZ3s                           = (((M-N-1)*(M-N-4))/(M*(M-2)))*((PSI2TILDE)/(PSI2TILDE+N/M));
        DKZ3s                           = (1/GAMMA)*(((M-N-1)*(M-N-4))/(M*(M-2)))*((N/M)/(PSI2TILDE+N/M))*MUTILDEg;

        ECALUTILDE_WKZ3HATCKZ3sDKZ3s    = (THETA2TILDE/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2TILDE+(THETA2TILDE/PSI2TILDE)*(N/M)));
        ECALUTILDE_WKZ3CKZ3sDKZ3s       = (THETA2/(2*GAMMA))*(((M-N-1)*(M-N-4))/((M-2)*(M-N-2)))*(1-(N/M)/(THETA2+(THETA2/PSI2)*(N/M)));

        tfr1tilde_vec(i,KK) = ECALUTILDE_WQ3IHATCsDs;
        tfr2tilde_vec(i,KK) = ECALUTILDE_WQ3IIHATCssDss;
        tfrrtilde_vec(i,KK) = ECALUTILDE_WKZ3HATCKZ3sDKZ3s;

        tfr1_vec(i,KK) = ECALUTILDE_WQ3ICsDs;
        tfr2_vec(i,KK) = ECALUTILDE_WQ3IICssDss;
        tfrr_vec(i,KK) = ECALUTILDE_WKZ3CKZ3sDKZ3s;

    end
end

% save simulations

save tfr1_vec;
save tfr2_vec;
save tfrr_vec;

save tfr1hat_vec;
save tfr2hat_vec;
save tfrrhat_vec;

save tfr1tilde_vec;
save tfr2tilde_vec;
save tfrrtilde_vec;

% plots
subplot(3,1,1)
boxplot(...
    [tfr1hat_vec(:,1),tfr1tilde_vec(:,1),tfr1_vec(:,1),...
    tfr1hat_vec(:,2),tfr1tilde_vec(:,2),tfr1_vec(:,2),...
    tfr1hat_vec(:,3),tfr1tilde_vec(:,3),tfr1_vec(:,3),...
    tfr1hat_vec(:,4),tfr1tilde_vec(:,4),tfr1_vec(:,4),...
    tfr1hat_vec(:,5),tfr1tilde_vec(:,5),tfr1_vec(:,5),...        
    tfr1hat_vec(:,6),tfr1tilde_vec(:,6),tfr1_vec(:,6)],...
    'Labels',...
    {'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True'})
set(gca,'FontSize',10,'XTickLabelRotation',90)
title('Three Fund Rule I')
line([3.5 3.5], get(gca, 'ylim'));
line([6.5 6.5], get(gca, 'ylim'));
line([9.5 9.5], get(gca, 'ylim'));
line([12.5 12.5], get(gca, 'ylim'));
line([15.5 15.5], get(gca, 'ylim'));
text(1.5,0.15,'T=50')
text(4.5,0.15,'T=80')
text(7.5,0.15,'T=120')
text(10.5,0.15,'T=200')
text(13.5,0.15,'T=500')
text(16.5,0.15,'T=1000')

subplot(3,1,2)
boxplot(...
    [tfr2hat_vec(:,1),tfr2tilde_vec(:,1),tfr2_vec(:,1),...
    tfr2hat_vec(:,2),tfr2tilde_vec(:,2),tfr2_vec(:,2),...
    tfr2hat_vec(:,3),tfr2tilde_vec(:,3),tfr2_vec(:,3),...
    tfr2hat_vec(:,4),tfr2tilde_vec(:,4),tfr2_vec(:,4),...
    tfr2hat_vec(:,5),tfr2tilde_vec(:,5),tfr2_vec(:,5),...
    tfr2hat_vec(:,6),tfr2tilde_vec(:,6),tfr2_vec(:,6)],...
    'Labels',...
    {'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True'})
set(gca,'FontSize',10,'XTickLabelRotation',90)
title('Three Fund Rule II')
line([3.5 3.5], get(gca, 'ylim'));
line([6.5 6.5], get(gca, 'ylim'));
line([9.5 9.5], get(gca, 'ylim'));
line([12.5 12.5], get(gca, 'ylim'));
line([15.5 15.5], get(gca, 'ylim'));
text(1.5,0.15,'T=50')
text(4.5,0.15,'T=80')
text(7.5,0.15,'T=120')
text(10.5,0.15,'T=200')
text(13.5,0.15,'T=500')
text(16.5,0.15,'T=1000')

subplot(3,1,3)
boxplot(...
    [tfrrhat_vec(:,1),tfrrtilde_vec(:,1),tfrr_vec(:,1),...
    tfrrhat_vec(:,2),tfrrtilde_vec(:,2),tfrr_vec(:,2),...
    tfrrhat_vec(:,3),tfrrtilde_vec(:,3),tfrr_vec(:,3),...
    tfrrhat_vec(:,4),tfrrtilde_vec(:,4),tfrr_vec(:,4),...
    tfrrhat_vec(:,5),tfrrtilde_vec(:,5),tfrr_vec(:,5),...        
    tfrrhat_vec(:,6),tfrrtilde_vec(:,6),tfrr_vec(:,6)],...
    'Labels',...
    {'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True',...
    'Hat', 'Tilde', 'True'})
set(gca,'FontSize',10,'XTickLabelRotation',90)
title('Three Fund Rule Revisited')
line([3.5 3.5], get(gca, 'ylim'));
line([6.5 6.5], get(gca, 'ylim'));
line([9.5 9.5], get(gca, 'ylim'));
line([12.5 12.5], get(gca, 'ylim'));
line([15.5 15.5], get(gca, 'ylim'));
text(1.5,0.15,'T=50')
text(4.5,0.15,'T=80')
text(7.5,0.15,'T=120')
text(10.5,0.15,'T=200')
text(13.5,0.15,'T=500')
text(16.5,0.15,'T=1000')
print('FillPageFigure','-dpdf','-fillpage')
