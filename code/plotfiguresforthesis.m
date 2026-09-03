% Producir figuras
%close all
data1 = cleandata(St01);
data2 = cleandata(St05);
d = 0.01;
hfig = figure;
yyaxis left
plot(data1(:,1)+273,data1(:,4),'LineWidth',1.5);
d = 0.01;
ylim([min(data1(:,4))-d,max(data1(:,4))+d])
ylabel('Conversion, \alpha')
yyaxis right
plot(data1(:,1)+273,data1(:,3),'--','LineWidth',1.5);
ylabel('Rate of conversion, d\alpha/dt')
d = 0.0001;
ylim([min(data1(:,3))-d,max(data1(:,3))+d])
xlabel('Temperatura [K]')
legend('TGA','DTG','Location','east')
xlim([min(data1(:,1)+273),max(data1(:,1)+273)])

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional

set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

%%
%set(gca,'fontname','times')

hfig = figure;
yyaxis left
plot(data2(:,1)+273,data2(:,4),'LineWidth',1.5);
%ylim([-70,0])
ylabel('Pérdida de masa [%]')
yyaxis right
plot(data2(:,1)+273,data2(:,3),'--','LineWidth',1.5);
ylabel('Velocidad de pérdida de masa [%/min]')
%ylim([-0.7,0.1])
xlabel('Temperatura [K]')
legend('TGA','DTG','Location','east')
%grid on
%set(gca,'fontname','times')
picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%%








%%

n_1 = 1;
E_1 = 150;
A_1 = 10^10;

R = 8.314E-3;

dXdt =@(t,X) [  A_1*exp(-E_1/R/X(2))*(1-X(1)).^n_1;
                10];
%+ 0.5*A_2*exp(-E_2/R/X(2))*(1-X(1))^0.5
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,x] = ode15s(dXdt,[0,25],[0.001,600],options);
dx1 =  A_1*exp(-E_1/R./x(:,2)).*(1-x(:,1));

hfig = figure;
yyaxis left
plot(x(:,2),x(:,1),'LineWidth',1.5);
ylim([-0.01,1.01])
ylabel('\alpha')
yyaxis right
plot(x(:,2),dx1,'--','LineWidth',1.5);
ylabel('d\alpha/dt')
ylim([-0.001,max(dx1)+0.001])
xlabel('Temperature')
legend('TGA','DTG','Location','east')
%grid on
%set(gca,'fontname','times')
picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional

set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%%







%%

hfig = figure;  % save the figure handle in a variable
% t = 0:0.02:10; x = t.*sin(2*pi*t)+ 2*rand(1,length(t)); % data
% plot(t,x,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
% xlabel('time $t$ (s)')
% ylabel('$\Omega$ (V)')
%figure
plot(data1(:,1),data1(:,4))
ylabel('Pérdida de masa [%]')
yyaxis right
plot(data1(:,1),data1(:,3))
ylabel('Velocidad de pérdida de masa')
xlabel('Temperatura [°C]')
fname = 'myfigure';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')