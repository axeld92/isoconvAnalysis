 data1 = St01;
 data2 = St06;
 data3 = St11;
 data4 = St16;
%%
 %Deconvolve signals
 [out,C] = deconvolution(data1,data2,data3,data4);

 %% Plot results
data1c = cleandata(data3);
hfig = figure; 
hold on
plot(data1c(:,1),data1c(:,3),'LineWidth',1.5)
plot(out{3}{1}(:,1),out{3}{1}(:,3)+out{3}{2}(:,3)+out{3}{3}(:,3),'LineWidth',1.5)
plot(out{3}{1}(:,1),out{3}{1}(:,3),'--','LineWidth',1.5)
plot(out{3}{2}(:,1),out{3}{2}(:,3),'--','LineWidth',1.5)
plot(out{3}{3}(:,1),out{3}{3}(:,3),'--','LineWidth',1.5)

hold off
limy = [-0.5  11];
limy = limy*1E-3;
ylim(limy)
xlabel('Temperature (K)')
ylabel('d\alpha/dt (1/min)')
legend('Exp.','Fit','Peak 1','Peak 2','Peak 3')


picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional

set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(hfig,'grafico11','-dpdf','-r0')
%% 
data1c = cleandata(data4);
ii = 4;
hfig = figure; 
hold on
plot(data1c(:,1),data1c(:,3),'LineWidth',1.5)
plot(out{ii}{3}(:,1),out{ii}{1}(:,3)+out{ii}{2}(:,3)+out{ii}{3}(:,3),'LineWidth',1.5)
plot(out{ii}{3}(:,1),out{ii}{3}(:,3),'--','LineWidth',1.5)
plot(out{ii}{1}(:,1),out{ii}{1}(:,3),'--','LineWidth',1.5)
plot(out{ii}{2}(:,1),out{ii}{2}(:,3),'--','LineWidth',1.5)

hold off
limy = [-0.5  10.2];
limy = limy*1E-3;
ylim(limy)
xlabel('Temperature (K)')
ylabel('d\alpha/dt (1/min)')
legend('Exp.','Fit','Peak 1','Peak 2','Peak 3')


picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


%%

datar = {data1,data2,data3,data4};
for ii = 1:4;

dataexp = cleandata(datar{ii});
dataexp = dataexp(:,3);
maxdata = max(dataexp);
dataexp = dataexp/maxdata;
datacalc = out{ii}{1}(:,3)+out{ii}{2}(:,3)+out{ii}{3}(:,3);
datacalc = datacalc/maxdata;
ECM(ii) = sum((dataexp - datacalc).^2)/(length(dataexp)-12);
Fitperc(ii) = (1-sqrt(ECM(ii))/max(dataexp))*100;
end
ECM
Fitperc




%%
% 

n = 4;
for j = 1:3
figure
hold on
    for i = 1:n
        plot(out{i}{j}(:,1),out{i}{j}(:,4)/C(i,j))
    end
hold off
grid
legend('5','10','15','20')
end

for j = 1:3
figure
hold on
    for i = 1:n
        plot(out{i}{j}(:,1),out{i}{j}(:,3)/C(i,j))
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
aa = 0.1:0.1:0.9;
EA_v = {};
ea_max = {};
for j = 1:3
EA_v{j} = Vyazovkinv(aa,out{1}{j},out{2}{j},out{3}{j},out{4}{j});
ea_max{j} = max(EA_v{j}(:,1));
ea_min{j} = min(EA_v{j}(:,1));
weights{j} =1./(EA_v{j}(:,3)-EA_v{j}(:,2)).^2;
E{j} = sum(weights{j}.*EA_v{j}(:,1))/sum(weights{j});
EAic_lo{j} = E{j} - sum(weights{j}.*EA_v{j}(:,2))/sum(weights{j});
EAic_up{j} = -E{j} + sum(weights{j}.*EA_v{j}(:,3))/sum(weights{j});
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
%%
%figure
%hold on
hfig = figure;
hold on
err_lo = EA_v{2}(:,1)-EA_v{2}(:,2);
err_up = -EA_v{2}(:,1)+EA_v{2}(:,3);
errorbar(aa,EA_v{3}(:,1),EA_v{2}(:,1)-EA_v{2}(:,2),-EA_v{2}(:,1)+EA_v{2}(:,3),'o-','LineWidth',1.5)
errorbar(aa,EA_v{1}(:,1),EA_v{1}(:,1)-EA_v{1}(:,2),-EA_v{1}(:,1)+EA_v{1}(:,3),'s-','LineWidth',1.5)
errorbar(aa,EA_v{2}(:,1),EA_v{3}(:,1)-EA_v{3}(:,2),-EA_v{3}(:,1)+EA_v{3}(:,3),'^-','LineWidth',1.5)
% errorbar(aa,EA_v{2}(:,1),EA_v{2}(:,1)-EA_v{2}(:,2),-EA_v{2}(:,1)+EA_v{2}(:,3),'o-','LineWidth',1.5)
% errorbar(aa,EA_v{1}(:,1),EA_v{1}(:,2),EA_v{1}(:,3),'s-','LineWidth',1.5)
% errorbar(aa,EA_v{3}(:,1),EA_v{2}(:,2),EA_v{3}(:,3),'^-','LineWidth',1.5)

%errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
hold off
%grid
%title('3')
xlabel('\alpha')
ylabel('E_{\alpha} (kJ/mol)')
legend('1','2','3','Location','best');%,'Total')

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%%

% % Calculate pre-exponential factor using compensation effect
for j = 1:3
As{j} = [compeffects(E{j},out{1}{j}),compeffects(E{j},out{2}{j}),compeffects(E{j},out{3}{j}),compeffects(E{j},out{4}{j})];
end
%%
% As_2 = [compeffects(E_2,out_2_05),compeffects(E_2,out_2_10),compeffects(E_2,out_2_15),compeffects(E_2,out_2_20)];
% As_3 = [compeffects(E_3,out_3_05),compeffects(E_3,out_3_10),compeffects(E_3,out_3_15),compeffects(E_3,out_3_20)];

A_1 = real(mean(As{1}))
A_2 = real(mean(As{2}))
A_3 = real(mean(As{3}))

lnA_1 = log(A_1)
lnA_2 = log(A_2)
lnA_3 = log(A_3)


%%
[damax, Tindex] = max(out{1}{1}(:,3));
Tps = out{1}{1}(:,1);
T_p = Tps(Tindex);
R = 8.314;

DeltaH_1 = (E{1}*1000 - R*T_p)/1000
DeltaG_1 = (E{1}*1000 + R*T_p*log(1.381E-23*T_p/6.626E-34/A_1))/1000
DeltaS_1 = ((DeltaH_1 - DeltaG_1)/T_p)





%%

% % Determine model  REVISAR QUE se corresponden las energías de activación
% y los factores preexponenciales y n y m. Idea, no importa porque todos se
% juntan al final en la suma error posiblemente en C

nmp_1 = [sbrlin(A_1,E{1},out{1}{1})' ; sbrlin(A_1,E{1},out{2}{1})' ; sbrlin(A_1,E{1},out{3}{1})' ; sbrlin(A_1,E{1},out{4}{1})'];
nmp_2 = [sbrlin(A_2,E{2},out{1}{2})' ; sbrlin(A_2,E{2},out{2}{2})' ; sbrlin(A_2,E{2},out{3}{2})' ; sbrlin(A_2,E{2},out{4}{2})'];
nmp_3 = [sbrlin(A_3,E{3},out{1}{3})' ; sbrlin(A_3,E{3},out{2}{3})' ; sbrlin(A_3,E{3},out{3}{3})' ; sbrlin(A_3,E{3},out{4}{3})'];


% figure
% plot(nmp_1);
% figure
% plot(nmp_2);
% figure
% plot(nmp_3);
dev_nmp_1 = real(std(nmp_1))
nmp_1 = real(mean(nmp_1));
dev_nmp_2 = real(std(nmp_2))
nmp_2 = real(mean(nmp_2));
dev_nmp_3 = real(std(nmp_3))
nmp_3 = real(mean(nmp_3));

%close all

%% determine model with log

nmp_1 = [sbrlin(A_1,E{1},out{1}{1})' ; sbrlin(A_1,E{1},out{2}{1})' ; sbrlin(A_1,E{1},out{3}{1})' ; sbrlin(A_1,E{1},out{4}{1})'];
nmp_2 = [sbrlin(A_2,E{2},out{1}{2})' ; sbrlin(A_2,E{2},out{2}{2})' ; sbrlin(A_2,E{2},out{3}{2})' ; sbrlin(A_2,E{2},out{4}{2})'];
nmp_3 = [sbrlin(A_3,E{3},out{1}{3})' ; sbrlin(A_3,E{3},out{2}{3})' ; sbrlin(A_3,E{3},out{3}{3})' ; sbrlin(A_3,E{3},out{4}{3})'];


nmp_1 = real(mean(nmp_1));
nmp_2 = real(mean(nmp_2));
nmp_3 = real(mean(nmp_3));

%%
% % Full kinetic model

R = 8.314E-3;

c = mean([C]);

datareal = cleandata(St16);
% 
t0 = datareal(1,2);
tf = datareal(end,2);
T0 = datareal(1,1);
Tf = datareal(end,1);
b = (Tf-T0)/(tf-t0);




dadt = @(t,aT)    [ A_1*exp(-E{1}/R/aT(4))*real(aT(1)^nmp_1(1)*(1-aT(1))^nmp_1(2)) ;          %1
                    A_2*exp(-E{2}/R/aT(4))*real(aT(2)^nmp_2(1)*(1-aT(2))^nmp_2(2)) ;          %2
                    A_3*exp(-E{3}/R/aT(4))*real(aT(3)^nmp_3(1)*(1-aT(3))^nmp_3(2)) ;          %3
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
 [t,aT] = ode15s(dadt,[0 tf],[1E-3 1E-3 1E-3 T0],options);
%  aT(aT~=real(aT)) = 0;
aT = real(aT);

figure
plot(t,aT(:,4),datareal(:,2),datareal(:,1))
legend('mod','dat')
 
a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3);
da_calc = diff_alpha(aT(:,4),a_calc);

a_fr = interp1(datareal(:,1),datareal(:,4),aT(:,4));
a_fr(isnan(a_fr))=1;
da_fr = interp1(datareal(:,1),datareal(:,3),aT(:,4));
resid = a_fr - a_calc;

MSE_model = sum(resid.^2)/(length(resid))
%%

tplot = datareal(:,2);
daplot = datareal(:,3);
aplot = datareal(:,4);
nn = length(tplot);
idx = 1:100:nn;
tplot = tplot(idx);
daplot = daplot(idx);
aplot = aplot(idx);
figure

hold on
plot(t,a_calc,'LineWidth',1)
plot(tplot,aplot,'--','LineWidth',1)
hold off
%ylim([0,0.012])
%% peak by peak


R = 8.314E-3;
A = {A_1,A_2,A_3};
c = mean([C]);
for k = 1:3
    
datar = out{1}{k};

T = datar(:,1);
t = datar(:,2);
da = datar(:,3);
a = datar(:,4);

t0 = t(1);
tend = t(end);
T0 = T(1);
Tend = T(end);
b = (Tend-T0)/(tend-t0);



%plot(T,a)

dadT =@(t,aT)  [A{k}.*exp(-E{k}/R./aT(2))*real(aT(1).^nmp_1(1).*(1-aT(1)).^nmp_1(2));b] ;




options = odeset('RelTol',1e-8,'AbsTol',1e-8) ;
[tt,aT] = ode15s(dadT,[0 tend],[1E-3 T0],options);
dadTT = diff_alpha(aT(:,2),aT(:,1));
figure
hold on
plot(t,da)
plot(tt,dadTT)
hold off

end


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



objfun =@(p) sum((a - ratelaw(T,a1,a2,a3,p)).^2);

p0 = [c(1)*A_1,E{1},nmp_1(1),nmp_1(2),c(2)*A_2,E{2},nmp_2(1),nmp_2(2),c(3)*A_3,E{3},nmp_3(1),nmp_3(2)];

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
 ac = kineticmodel(T,a1,a2,a3,fitparam);
 ac = ac/max(ac);
 
 
 plot(T,a,T,ac)
 
 
%%

calcrate = kineticmodel(T,a1,a2,a3,fitparam);
figure
plot(T,a,T,calcrate)

%%

function rate = ratelaw(t,T,a,params)
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

drate=@(t,X)    [cA1.*exp(-E1/R./T).*X(1).^n1.*(1-X(1)).^m1;
                cA2.*exp(-E2/R./T).*X(2).^n2.*(1-X(2)).^m2;
                cA3.*exp(-E3/R./T).*X(3).^n3.*(1-X(3)).^m3;
                15;]
 options = odeset('RelTol',1e-8,'AbsTol',1e-8) ;            
 [t,X] = ode15s(drate,[0 tf],[1E-3 1E-3 1E-3 T0],options);
 rate = X(:,1) + X(:,2)+ X(:,3) 
 
end







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

%drate = cA1.*exp(-E1/R./T).*a1.^n1.*(1-a1).^m1 + cA2.*exp(-E2/R./T).*a2.^n2.*(1-a2).^m2 + cA3.*exp(-E3/R./T).*a3.^n3.*(1-a3).^m3;
%rate = cumtrapz(T,drate);
rate = cA1.*exp(-E1/R./T).*a1.^n1.*(1-a1).^m1 + cA2.*exp(-E2/R./T).*a2.^n2.*(1-a2).^m2 + cA3.*exp(-E3/R./T).*a3.^n3.*(1-a3).^m3;

end
