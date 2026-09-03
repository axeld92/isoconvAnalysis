 data1 = St01;
 data2 = St06;
 data3 = St11;
 data4 = St16;

 %Deconvolve signals


[out_1_05,out_2_05,out_3_05,c_05] = deconvolve(data1);
[out_1_10,out_2_10,out_3_10,c_10] = deconvolve(data2);
%[out_1_15,out_2_15,out_3_15,c_15] = deconvolve(data3);
[out_1_20,out_2_20,out_3_20,c_20] = deconvolve(data4);


figure
hold on
plot(out_1_05(:,1),out_1_05(:,3))
plot(out_1_10(:,1),out_1_10(:,3))
%plot(out_1_15(:,1),out_1_15(:,3))
plot(out_1_20(:,1),out_1_20(:,3))
hold off
grid
legend('5','10','20')
title('hemicellulose')
% 

figure
hold on
plot(out_2_05(:,1),out_2_05(:,3))
plot(out_2_10(:,1),out_2_10(:,3))
%plot(out_2_15(:,1),out_2_15(:,3))
plot(out_2_20(:,1),out_2_20(:,3))
hold off
grid
legend('5','10','20')
title('cellulose')


figure
hold on
plot(out_3_05(:,1),out_3_05(:,3))
plot(out_3_10(:,1),out_3_10(:,3))
%plot(out_3_15(:,1),out_3_15(:,3))
plot(out_3_20(:,1),out_3_20(:,3))
hold off
grid
legend('5','10','20')
title('lignin')
%%
% % Calculate activation energy of individual steps
% 
aa = 0.2:0.01:0.8;

EA_v1 = Starink3(aa,out_1_05,out_1_10,out_1_20);
ea1_max = max(EA_v1(:,1));
ea1_min = min(EA_v1(:,1));
E_1 = mean(EA_v1(:,1));
variation_1 = 100*(ea1_max-ea1_min)/E_1
% 
% 
figure
hold on
errorbar(aa,EA_v1(:,1),EA_v1(:,2),'s-')
hold off
grid
title('1')
% % 
% % 
% 
EA_v2 = Starink3(aa,out_2_05,out_2_10,out_2_20);
ea2_max = max(EA_v2(:,1));
ea2_min = min(EA_v2(:,1));
E_2 = mean(EA_v2(:,1));
variation_2 = 100*(ea2_max-ea2_min)/E_2

% 
figure
hold on
errorbar(aa,EA_v2(:,1),EA_v2(:,2),'s-')
hold off
grid
title('2')

EA_v3 = Vyazovkinfor3(aa,out_3_05,out_3_10,out_3_20);
ea3_max = max(EA_v3(:,1));
ea3_min = min(EA_v3(:,1));
E_3 = mean(EA_v3(:,1));
variation_3 = 100*(ea3_max-ea3_min)/E_3

% 
figure
hold on
errorbar(aa,EA_v3(:,1),EA_v3(:,2),'s-')
hold off
grid
title('3')

%%

% % Calculate pre-exponential factor using compensation effect

A_1 = mean([compeffects(E_1,out_1_05),compeffects(E_1,out_1_10),compeffects(E_1,out_1_20)]);
A_2 = mean([compeffects(E_2,out_2_05),compeffects(E_2,out_2_10),compeffects(E_2,out_2_20)]);
A_3 = mean([compeffects(E_3,out_3_05),compeffects(E_3,out_3_10),compeffects(E_3,out_3_20)]);

% % Determine model

nmp_1 = sblin(A_1,E_1,out_1_10);
nmp_2 = sblin(A_2,E_2,out_2_10);
nmp_3 = sblin(A_3,E_3,out_3_10);

% % 


%%
% % Full kinetic model

R = 8.314E-3;

c = mean([c_05';c_10';c_20']);

dadt = @(t,aT)    [ A_1*exp(-E_1/R/aT(4))*nmp_1(1)*aT(1)^nmp_1(2)*(1-aT(1))^nmp_1(3) ;          %1
                    A_2*exp(-E_2/R/aT(4))*nmp_2(1)*aT(2)^nmp_2(2)*(1-aT(2))^nmp_2(3) ;          %2
                    A_3*exp(-E_3/R/aT(4))*nmp_3(1)*aT(3)^nmp_3(2)*(1-aT(3))^nmp_3(3) ;          %3
                    10 ];
                    
                    
                   
%M = eye(6);M(6,6) = 0;
            
% DeltaT = 780-130;
% time = 0:DeltaT*60/5;
% Temp = 130+273.15 + 5*time;
% dadtt = dadt()                
 

St06c = cleandata(St06);

%options = odeset('Mass',M);               
 T0 = St06c(1,2);
 Tf = St06c(end,2);


 [t,aT] = ode15s(dadt,[T0 Tf],[1E-3 1E-3 1E-3 400]);
%  aT(aT~=real(aT)) = 0;

 
a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3);



 
 
 figure
 hold on
 plot(aT(:,4)-9.7,a_calc,'-')
 plot(St06c(:,1),St06c(:,4),'--')
 hold off
 

%%


dadT =@(Temp,cA_1,cA_2,cA_3,E_1,E_2,E_3,n_1,n_2,n_3,m_1,m_2,m_3) [cA_1*exp(-E_1/R./Temp)*a_1.^n_1.*(1-a_1).^m_1 +...
                                                                cA_2*exp(-E_2/R./Temp)*a_2.^n_2.*(1-a_2).^m_2 +...
                                                                cA_3*exp(-E_3/R./Temp)*a_3.^n_3.*(1-a_3).^m_3];
                                                            
funct = @(p) dadT(Temp,p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10),p(11),p(12));


Temp = St06c(:,1);
aexp = St06c(:,4);




