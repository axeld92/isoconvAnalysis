 data1 = St05;
 data2 = St010;
 data3 = St15;
 data4 = St20;
%  data1 = cleandata(St01);
%  data2 = cleandata(St06);
%  data3 = cleandata(St11);
%  data4 = cleandata(St16);

% figure
% hold on
%  plot(data1(:,1),data1(:,3))
%  plot(data2(:,1),data2(:,3))
%  plot(data3(:,1),data3(:,3))
%  plot(data4(:,1),data4(:,3))
%  hold off
%  
%  %%
 %Deconvolve signals


[out_1_05,out_2_05,out_3_05,c_05] = deconvolve(data1);
[out_1_10,out_2_10,out_3_10,c_10] = deconvolve(data2);
[out_1_15,out_2_15,out_3_15,c_15] = deconvolve(data3);
[out_1_20,out_2_20,out_3_20,c_20] = deconvolve(data4);


%%
% 
figure
hold on
plot(out_1_05(:,1),out_1_05(:,3))
plot(out_1_10(:,1),out_1_10(:,3))
plot(out_1_15(:,1),out_1_15(:,3))
plot(out_1_20(:,1),out_1_20(:,3))
hold off
grid
legend('5','10','15','20')
title('1')
% 

figure
hold on
plot(out_2_05(:,1),out_2_05(:,3))
plot(out_2_10(:,1),out_2_10(:,3))
plot(out_2_15(:,1),out_2_15(:,3))
plot(out_2_20(:,1),out_2_20(:,3))
hold off
grid
legend('5','10','15','20')
title('2')

figure
hold on
plot(out_3_05(:,1),out_3_05(:,3))
plot(out_3_10(:,1),out_3_10(:,3))
plot(out_3_15(:,1),out_3_15(:,3))
plot(out_3_20(:,1),out_3_20(:,3))
hold off
grid
legend('5','10','15','20')
title('3')
%%
% % Calculate activation energy of individual steps
% 
aa = 0.1:0.01:0.9;

EA_v1 = Vyazovkin(aa,out_1_05,out_1_10,out_1_15,out_1_20);
ea1_max = max(EA_v1(:,1));
ea1_min = min(EA_v1(:,1));
E_1 = median(EA_v1(:,1));
variation_1 = 100*(ea1_max-ea1_min)/E_1
% 
 
% figure
% hold on
% errorbar(aa,EA_v1(:,1),EA_v1(:,2),'s-')
%hold off
%grid
%title('1')
%xlabel('\alpha')
%ylabel('E_{\alpha}')



EA_v2 = Vyazovkin(aa,out_2_05,out_2_10,out_2_15,out_2_20);
ea2_max = max(EA_v2(:,1));
ea2_min = min(EA_v2(:,1));
E_2 = median(EA_v2(:,1));
variation_2 = 100*(ea2_max-ea2_min)/E_2

% 
%figure
%hold on
% errorbar(aa,EA_v2(:,1),EA_v2(:,2),'s-')
%hold off
%grid
%title('2')
%xlabel('\alpha')
%ylabel('E_{\alpha}')

EA_v3 = Vyazovkin(aa,out_3_05,out_3_10,out_3_15,out_3_20);
ea3_max = max(EA_v3(:,1));
ea3_min = min(EA_v3(:,1));
E_3 = median(EA_v3(:,1));
variation_3 = 100*(ea3_max-ea3_min)/E_3


 data1c = cleandata(data1);
 data2c = cleandata(data2);
 data3c = cleandata(data3);
 data4c = cleandata(data4);




EA_t = Vyazovkin(aa,data1c,data2c,data3c,data4c);



% 
%figure
%hold on
figure
hold on
errorbar(aa,EA_v1(:,1),EA_v1(:,2),'o-')
errorbar(aa,EA_v2(:,1),EA_v2(:,2),'s-')
errorbar(aa,EA_v3(:,1),EA_v3(:,2),'^-')
errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
hold off
grid
%title('3')
xlabel('\alpha')
ylabel('E_{\alpha}')
legend('1','2','3','Total')
%%

% % Calculate pre-exponential factor using compensation effect

As_1 = [compeffects(E_1,out_1_05),compeffects(E_1,out_1_10),compeffects(E_1,out_1_15),compeffects(E_1,out_1_20)];
As_2 = [compeffects(E_2,out_2_05),compeffects(E_2,out_2_10),compeffects(E_2,out_2_15),compeffects(E_2,out_2_20)];
As_3 = [compeffects(E_3,out_3_05),compeffects(E_3,out_3_10),compeffects(E_3,out_3_15),compeffects(E_3,out_3_20)];


A_1 = real(mean(As_1));
A_2 = real(mean(As_2));
A_3 = real(mean(As_3));

%%

% % Determine model

nmp_1 = [sblin(A_1,E_1,out_1_05)' ; sblin(A_1,E_1,out_1_10)' ; sblin(A_1,E_1,out_1_15)' ; sblin(A_1,E_1,out_1_20)']
nmp_2 = [sblin(A_2,E_2,out_2_05)' ; sblin(A_2,E_2,out_2_10)' ; sblin(A_2,E_2,out_2_15)' ; sblin(A_2,E_2,out_2_20)']
nmp_3 = [sblin(A_3,E_3,out_3_05)' ; sblin(A_3,E_3,out_3_10)' ; sblin(A_3,E_3,out_3_15)' ; sblin(A_3,E_3,out_3_20)']


nmp_1 = real(mean(nmp_1));
nmp_2 = real(mean(nmp_2));
nmp_3 = real(mean(nmp_3));




% % 


%%
% % Full kinetic model

R = 8.314E-3;

c = mean([c_05';c_10';c_15';c_20']);

dadt = @(t,aT)    [ A_1*exp(-E_1/R/aT(4))*nmp_1(1)*aT(1)^nmp_1(2)*(1-aT(1))^nmp_1(3) ;          %1
                    A_2*exp(-E_2/R/aT(4))*nmp_2(1)*aT(2)^nmp_2(2)*(1-aT(2))^nmp_2(3) ;          %2
                    A_3*exp(-E_3/R/aT(4))*nmp_3(1)*aT(3)^nmp_3(2)*(1-aT(3))^nmp_3(3) ;          %3
                    10 ];
                    
St06c = cleandata(St06);

t0 = St06c(1,2);
tf = St06c(end,2);
T0 = St06c(1,1);

 options = odeset('RelTol',1e-10,'AbsTol',1e-10) ;
 [t,aT] = ode23s(dadt,[0 tf],[1E-3 1E-3 1E-3 T0],options);
%  aT(aT~=real(aT)) = 0;
aT = real(aT);
 
a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3);

a_fr = interp1(St06c(:,1),St06c(:,4),aT(:,4));
da_fr = interp1(St06c(:,1),St06c(:,3),aT(:,4));
resid = a_fr - a_calc;
%%
figure 
subplot(5,1,1:4);
hold on
plot(aT(:,4),a_calc,'-')
plot(St06c(:,1),St06c(:,4),'--')
hold off
legend('Modelo','Datos')
grid
ylabel('\alpha')
subplot(5,1,5);
plot(aT(:,4),resid)
xlabel('Temperatura [K]')
grid 
 %set(gca, 'xtick', [] );
 %p = get(gca,'Position');
 %p_diff = p(4)*0.1;
 %p(4) = p(4) + p_diff
 %p(2) = p(2) - p_diff


