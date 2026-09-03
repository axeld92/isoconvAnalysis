 data1 = St03;
 data2 = St08;
 data3 = St13;
 data4 = St18;
%%
 %Deconvolve signals
 [out,C] = deconvolution4(data1,data2,data3,data4);

 
  %% Plot results
data1c = cleandata(data4);
hfig = figure; 
hold on
plot(data1c(:,1),data1c(:,3),'LineWidth',1.5)
plot(out{4}{3}(:,1),out{4}{1}(:,3)+out{4}{2}(:,3)+out{4}{3}(:,3) + out{4}{4}(:,3),'LineWidth',1.5)
plot(out{4}{1}(:,1),out{4}{1}(:,3),'--','LineWidth',1.5)
plot(out{4}{4}(:,1),out{4}{4}(:,3),'--','LineWidth',1.5)
plot(out{4}{3}(:,1),out{4}{3}(:,3),'--','LineWidth',1.5)
plot(out{4}{2}(:,1),out{4}{2}(:,3),'--','LineWidth',1.5)


hold off
limy = [-1  1100*max(data1c(:,3))];
limy = limy*1E-3;
ylim(limy)
xlabel('Temperatura [K]')
ylabel('d\alpha/dt [min^{-1}]')
legend('Experimental','Ajustado','Pico 1','Pico 2','Pico 3','Pico 4')


picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(hfig,'grafico11','-dpdf','-r0')
 
%%
for i = 1:4
for j = 1:4

    out{i}{j}(:,3) = out{i}{j}(:,3)*C(i,j);
end
end


%% 
%%% plot conversion curves
n = 4;
for j = 1:4
figure
hold on
    for i = 1:n
        plot(out{i}{j}(:,1),out{i}{j}(:,4))
    end
hold off
grid
legend('5','10','15','20')
end

%%% plot conversion rate curves

for j = 1:4
figure
hold on
    for i = 1:n
        plot(out{i}{j}(:,1),out{i}{j}(:,3))
    end
hold off
grid
legend('5','10','15','20')
end

% % grid
% % legend('5','10','15','20')
% % title('1')
% % 
%%
% % Calculate activation energy of individual steps
% 
aa = 0.1:0.05:0.9;
EA_v = {};
ea_max = {};

for j = 1:4
EA_v{j} = Vyazovkinv(aa,out{1}{j},out{2}{j},out{3}{j},out{4}{j});
ea_max{j} = max(EA_v{j}(:,1));
ea_min{j} = min(EA_v{j}(:,1));
E{j} = mean(EA_v{j}(:,1));
variation{j} = 100*(ea_max{j}-ea_min{j})/E{j};
end% 
%  
% figure
% hold on
% plot(aa,EA_v1(:,1),'s-')
% hold off
% grid
% title('1')
% xlabel('\alpha')
% ylabel('E_{\alpha}')


% 
% EA_v2 = Vyazovkinv(aa,out{1}{2},out{2}{2},out{3}{2},out{4}{2});
% ea2_max = max(EA_v2(:,1));
% ea2_min = min(EA_v2(:,1));
% E_2 = median(EA_v2(:,1));
% variation_2 = 100*(ea2_max-ea2_min)/E_2
% 
% % 
% figure
% hold on
% %errorbar(aa,EA_v2(:,1),EA_v2(:,2),'s-')
% plot(aa,EA_v2(:,1),'s-')
% hold off
% grid
% title('2')
% xlabel('\alpha')
% ylabel('E_{\alpha}')
% %%
% EA_v3 = Vyazovkin(aa,out{1}{3},out{2}{3},out{3}{3},out{4}{3});
% ea3_max = max(EA_v3(:,1));
% ea3_min = min(EA_v3(:,1));
% E_3 = median(EA_v3(:,1));
% variation_3 = 100*(ea3_max-ea3_min)/E_3
% 
%%
 data1c = cleandata(data1);
 data2c = cleandata(data2);
 data3c = cleandata(data3);
 data4c = cleandata(data4);

EA_t = Vyazovkinv(aa,data1c,data2c,data3c,data4c);
% figure
% hold on
% yyaxis left
% errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
% yyaxis right
% plot(data1c(:,4),data1c(:,3))
% hold off

%figure
%hold on
%%
hfig = figure;
hold on
errorbar(aa,EA_v{1}(:,1),EA_v{1}(:,2),'o-','LineWidth',1.5)
errorbar(aa,EA_v{4}(:,1),EA_v{4}(:,2),'s-','LineWidth',1.5)


errorbar(aa,EA_v{3}(:,1),EA_v{3}(:,2),'^-','LineWidth',1.5)
errorbar(aa,EA_v{2}(:,1),EA_v{2}(:,2),'v-','LineWidth',1.5)
%errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
hold off
%grid
%title('3')
xlabel('\alpha')
ylabel('E_{\alpha}')
%legend('1','2','3','4')%,'Total')

legend('1','2','3','4','Location','northwest');%,'Total')

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(hfig,'grafico11','-dpdf','-r0')
%%
figure
hold on
plot(out{1}{1}(:,3))
plot(out{1}{2}(:,3))
plot(out{1}{3}(:,3))
plot(out{1}{4}(:,3))
hold off;legend('1','2','3','4')

%%

% % Calculate pre-exponential factor using compensation effect
for j = 1:4
As{j} = [compeffects(E{j},out{1}{j}),compeffects(E{j},out{2}{j}),compeffects(E{j},out{3}{j}),compeffects(E{j},out{4}{j})];
end
% As_2 = [compeffects(E_2,out_2_05),compeffects(E_2,out_2_10),compeffects(E_2,out_2_15),compeffects(E_2,out_2_20)];
% As_3 = [compeffects(E_3,out_3_05),compeffects(E_3,out_3_10),compeffects(E_3,out_3_15),compeffects(E_3,out_3_20)];


A_1 = real(mean(As{1}));
A_2 = real(mean(As{2}));
A_3 = real(mean(As{3}));
A_4 = real(mean(As{4})); 
%close all
%%

% % Determine model

nmp_1 = [sbrlin(A_1,E{1},out{1}{1})' ; sbrlin(A_1,E{1},out{2}{1})' ; sbrlin(A_1,E{1},out{3}{1})' ; sbrlin(A_1,E{1},out{4}{1})'];
nmp_2 = [sbrlin(A_2,E{2},out{1}{2})' ; sbrlin(A_2,E{2},out{2}{2})' ; sbrlin(A_2,E{2},out{3}{2})' ; sbrlin(A_2,E{2},out{4}{2})'];
nmp_3 = [sbrlin(A_3,E{3},out{1}{3})' ; sbrlin(A_3,E{3},out{2}{3})' ; sbrlin(A_3,E{3},out{3}{3})' ; sbrlin(A_3,E{3},out{4}{3})'];
nmp_4 = [sbrlin(A_4,E{4},out{1}{4})' ; sbrlin(A_4,E{4},out{2}{4})' ; sbrlin(A_4,E{4},out{3}{4})' ; sbrlin(A_4,E{4},out{4}{4})'];

nmp_1 = real(mean(nmp_1));
nmp_2 = real(mean(nmp_2));
nmp_3 = real(mean(nmp_3));
nmp_4 = real(mean(nmp_4));


close all
% % 


%%
% % Full kinetic model

R = 8.314E-3;

c = mean([C]);

R = 8.314E-3;

c = mean([C]);

datareal = cleandata(St15);
% 
t0 = datareal(1,2);
tf = datareal(end,2);
T0 = datareal(1,1);
Tf = datareal(end,1);
b = (Tf-T0)/(tf-t0);
% 


dadt = @(t,aT)    [ A_1*exp(-E{1}/R/aT(5))*real(aT(1)^nmp_1(1)*(1-aT(1))^nmp_1(2)) ;          %1
                    A_2*exp(-E{2}/R/aT(5))*real(aT(2)^nmp_2(1)*(1-aT(2))^nmp_2(2)) ;          %2
                    A_3*exp(-E{3}/R/aT(5))*real(aT(3)^nmp_3(1)*(1-aT(3))^nmp_3(2)) ;          %3
                    A_4*exp(-E{4}/R/aT(5))*real(aT(4)^nmp_4(1)*(1-aT(4))^nmp_4(2)) ;          %3
                    b ];
                    


% dadt = @(t,aT)    [ A_1*exp(-E{1}/R/aT(4))*aT(1)^nmp_1(1)*(1-aT(1))^nmp_1(2)*(-log(1-aT(1)))^nmp_1(3) ;          %1
%                     A_2*exp(-E{2}/R/aT(4))*aT(2)^nmp_2(1)*(1-aT(2))^nmp_2(2)*(-log(1-aT(2)))^nmp_1(3) ;          %2
%                     A_3*exp(-E{3}/R/aT(4))*aT(3)^nmp_3(1)*(1-aT(3))^nmp_3(2)*(-log(1-aT(3)))^nmp_1(3) ;          %3
%                     10 ];

% dadt = @(t,aT)    [ A_1*exp(-E{1}/R/aT(4))*real(nmp_1(1)*aT(1)^nmp_1(2)*(1-aT(1))^nmp_1(3)*(-log(1-aT(1)))^nmp_1(4)) ;          %1
%                     A_2*exp(-E{2}/R/aT(4))*real(nmp_2(1)*aT(2)^nmp_2(2)*(1-aT(2))^nmp_2(3)*(-log(1-aT(2)))^nmp_1(4)) ;          %2
%                     A_3*exp(-E{3}/R/aT(4))*real(nmp_3(1)*aT(3)^nmp_3(2)*(1-aT(3))^nmp_3(3)*(-log(1-aT(3)))^nmp_1(4)) ;          %3
%                     10 ];



 options = odeset('RelTol',1e-8,'AbsTol',1e-8) ;
 [t,aT] = ode15s(dadt,[0 tf],[1E-3 1E-3 1E-3 1E-3 T0],options);
%  aT(aT~=real(aT)) = 0;
aT = real(aT);
 
a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3) + c(4)*aT(:,4);
da_calc = diff_alpha(aT(:,4),a_calc);



figure
hold on
plot(t,a_calc,'--','LineWidth',2)
plot(datareal(:,2),datareal(:,4),':','LineWidth',2)
hold off
legend('Simulación','Datos experimentales')



a_fr = interp1(datareal(:,1),datareal(:,4),aT(:,4));
da_fr = interp1(datareal(:,1),datareal(:,3),aT(:,4));
resid = a_fr - a_calc;
%%

figure
hold on
plot(t,da_calc)
plot(datareal(:,2),datareal(:,3))
hold off
ylim([0,0.012])


%%
figure 
subplot(5,1,1:4);
hold on
plot(aT(:,4),a_calc,'-')
plot(datareal(:,1),datareal(:,4),'--')
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

%%

a1 = aT(:,1);
a2 = aT(:,2);
a3 = aT(:,3);
T = aT(:,4);

a = interp1(datareal(:,1),datareal(:,4),T);
% 
% M = [ones(size(T)), -1/R./T,  log(a1),  log(1-a1),   ones(size(T)), -1/R./T,  log(a2),  log(1-a2),...
% ones(size(T)), -1/R./T,  log(a3),  log(1-a3)];
% M(M==-inf) = -36.736800569677100;

% %%
% 
% 
% MM = M'*M;
% Ma = M'*a;
% 
% y = MM\Ma;



objfun =@(p) sum((a - kineticmodel(T,a1,a2,a3,p)).^2);

p0 = [c(1)*A_1,E{1},nmp_1(2),nmp_1(3),c(2)*A_2,E{2},nmp_2(2),nmp_2(3),c(3)*A_3,E{3},nmp_3(2),nmp_3(3)];

idx = p0<0;
lb = zeros(size(p0));
lb(idx) = p0(idx) + p0(idx)*0.5;
lb(~idx) = p0(~idx) - p0(~idx)*0.5;

ub = zeros(size(p0));
ub(idx) = p0(idx) - p0(idx)*0.5;
ub(~idx) = p0(~idx) + p0(~idx)*0.5;


opts = optimoptions('particleswarm','SwarmSize',1000);
 [fitparam , fval] = particleswarm(objfun,12,lb,ub,opts);
 checklb = fitparam == ub
 checkub = fitparam == lb



%%

calcrate = kineticmodel(T,a1,a2,a3,fitparam);
figure
plot(T,a,T,calcrate)









%%
function rate = kineticmodel(T,a1,a2,a3,params)
cA1 = params(1);
E1 = params(2);
n1 = params(3);
m1 = params(4);
cA2 = params(5);
E2 = params(6);
n2 = params(7);
m2 = params(8);
cA3 = params(9);
E3 = params(10);
n3 = params(11);
m3 = params(12);
R = 8.314E-3;

drate = cA1.*exp(-E1/R./T).*a1.^n1.*(1-a1).^m1 + cA2.*exp(-E2/R./T).*a2.^n2.*(1-a2).^m2 + cA3.*exp(-E3/R./T).*a3.^n3.*(1-a3).^m3;
rate = cumtrapz(T,drate);
end
