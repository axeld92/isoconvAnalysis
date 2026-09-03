 data1 = St02;
 data2 = St07;
 data3 = St12;
 data4 = St17;

 %Deconvolve signals


[out_1_05,out_2_05,out_3_05,out_4_05,c_05] = deconvolve4(data1);
[out_1_10,out_2_10,out_3_10,out_4_10,c_10] = deconvolve4(data2);
[out_1_15,out_2_15,out_3_15,out_4_15,c_15] = deconvolve4(data3);
[out_1_20,out_2_20,out_3_20,out_4_20,c_20] = deconvolve4(data4);




%%
% % Calculate activation energy of individual steps
% 
aa = 0.1:0.1:0.9;

EA_v1 = Vyazovkin(aa,out_1_05,out_1_10,out_1_15,out_1_20);
ea1_max = max(EA_v1(:,1));
ea1_min = min(EA_v1(:,1));
E_1 = median(EA_v1(:,1));
variation_1 = 100*(ea1_max-ea1_min)/E_1
% 
% 
% figure
% hold on
% errorbar(aa,EA_v1(:,1),EA_v1(:,2),'s-')
% hold off
% grid
% title('1')
% xlabel('\alpha')
% ylabel('E_{\alpha}')
% % 
% % 
% %
EA_v2 = Vyazovkin(aa,out_2_05,out_2_10,out_2_15,out_2_20);
ea2_max = max(EA_v2(:,1));
ea2_min = min(EA_v2(:,1));
E_2 = median(EA_v2(:,1));
variation_2 = 100*(ea2_max-ea2_min)/E_2

% 
% figure
% hold on
% errorbar(aa,EA_v2(:,1),EA_v2(:,2),'s-')
% hold off
% grid
% title('2')
% xlabel('\alpha')
% ylabel('E_{\alpha}')

EA_v3 = Vyazovkin(aa,out_3_05,out_3_10,out_3_15,out_3_20);
ea3_max = max(EA_v3(:,1));
ea3_min = min(EA_v3(:,1));
E_3 = median(EA_v3(:,1));
variation_3 = 100*(ea3_max-ea3_min)/E_3

% figure
% hold on
% errorbar(aa,EA_v3(:,1),EA_v3(:,2),'s-')
% hold off
% grid
% title('3')
% xlabel('\alpha')
% ylabel('E_{\alpha}')


EA_v4 = Vyazovkin(aa,out_4_05,out_4_10,out_4_15,out_4_20);
ea4_max = max(EA_v4(:,1));
ea4_min = min(EA_v4(:,1));
E_4 = median(EA_v4(:,1));
variation_4 = 100*(ea4_max-ea4_min)/E_4


% 
figure
hold on
errorbar(aa,EA_v1(:,1),EA_v1(:,2),'o-')
errorbar(aa,EA_v2(:,1),EA_v2(:,2),'s-')
errorbar(aa,EA_v3(:,1),EA_v3(:,2),'^-')
errorbar(aa,EA_v4(:,1),EA_v4(:,2),'v-')
hold off
grid

xlabel('\alpha')
ylabel('E_{\alpha}')
legend('1','2','3','4')


%%

% % Calculate pre-exponential factor using compensation effect

As_1 = [compeffects(E_1,out_1_05),compeffects(E_1,out_1_10),compeffects(E_1,out_1_15),compeffects(E_1,out_1_20)];
As_2 = [compeffects(E_2,out_2_05),compeffects(E_2,out_2_10),compeffects(E_2,out_2_15),compeffects(E_2,out_2_20)];
As_3 = [compeffects(E_3,out_3_05),compeffects(E_3,out_3_10),compeffects(E_3,out_3_15),compeffects(E_3,out_3_20)];
As_4 = [compeffects(E_4,out_4_05),compeffects(E_4,out_4_10),compeffects(E_4,out_4_15),compeffects(E_4,out_4_20)];



A_1 = real(mean(As_1));
A_2 = real(mean(As_2));
A_3 = real(mean(As_3));
A_4 = real(mean(As_4));

%% Determine model

nmp_1 = [sblin(A_1,E_1,out_1_05)' ; sblin(A_1,E_1,out_1_10)' ; sblin(A_1,E_1,out_1_15)' ; sblin(A_1,E_1,out_1_20)']
nmp_2 = [sblin(A_2,E_2,out_2_05)' ; sblin(A_2,E_2,out_2_10)' ; sblin(A_2,E_2,out_2_15)' ; sblin(A_2,E_2,out_2_20)']
nmp_3 = [sblin(A_3,E_3,out_3_05)' ; sblin(A_3,E_3,out_3_10)' ; sblin(A_3,E_3,out_3_15)' ; sblin(A_3,E_3,out_3_20)']
nmp_4 = [sblin(A_4,E_4,out_4_05)' ; sblin(A_4,E_4,out_4_10)' ; sblin(A_4,E_4,out_4_15)' ; sblin(A_4,E_4,out_4_20)']

nmp_1 = mean(nmp_1);
nmp_2 = mean(nmp_2);
nmp_3 = mean(nmp_3);
nmp_4 = mean(nmp_4);

% % 


%%
% % Full kinetic model

R = 8.314E-3;

c = mean([c_05';c_10';c_15';c_20']);

dadt = @(t,aT)    [ A_1*exp(-E_1/R/aT(5))*nmp_1(1)*aT(1)^nmp_1(2)*(1-aT(1))^nmp_1(3) ;          %1
                    A_2*exp(-E_2/R/aT(5))*nmp_2(1)*aT(2)^nmp_2(2)*(1-aT(2))^nmp_2(3) ;          %2
                    A_3*exp(-E_3/R/aT(5))*nmp_3(1)*aT(3)^nmp_3(2)*(1-aT(3))^nmp_3(3) ;          %3
                    A_4*exp(-E_4/R/aT(5))*nmp_4(1)*aT(4)^nmp_4(2)*(1-aT(4))^nmp_4(3) ;          %4
                    10 ];
                    
 St06c = cleandata(St06);
 
 
 
 
t0 = St06c(1,2);
tf = St06c(end,2);
T0 = St06c(1,1);


 [t,aT] = ode15s(dadt,[0 tf],[1E-3 1E-3 1E-3 1E-3 T0]);
%  aT(aT~=real(aT)) = 0;
aT = real(aT);
 
a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3) + c(4)*aT(:,4);

a_fr = interp1(St06c(:,1),St06c(:,4),aT(:,5));
da_fr = interp1(St06c(:,1),St06c(:,3),aT(:,5));
resid = a_fr - a_calc;

figure 
subplot(5,1,1:4);
hold on
plot(aT(:,5),a_calc,'-')
plot(St06c(:,1),St06c(:,4),'--')
hold off
legend('Modelo','Datos')
grid
ylabel('\alpha')
subplot(5,1,5);
plot(aT(:,5),resid)
xlabel('Temperatura [K]')
grid 
 %set(gca, 'xtick', [] );
 %p = get(gca,'Position');
 %p_diff = p(4)*0.1;
 %p(4) = p(4) + p_diff
 %p(2) = p(2) - p_diff
