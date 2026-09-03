function plotmassloss(data1,data2,data3,data4)

data1 = cleandatam(data1);
data2 = cleandatam(data2);
data3 = cleandatam(data3);
data4 = cleandatam(data4);

hfig = figure;
colororder({'k','k'})
hold on
plot(data1(:,1)+273.15,data1(:,4),'color','#0072BD','LineStyle','-','LineWidth',1.5)
plot(data2(:,1)+273.15,data2(:,4),'color','#D95319','LineStyle','-','LineWidth',1.5)
plot(data3(:,1)+273.15,data3(:,4),'color','#EDB120','LineStyle','-','LineWidth',1.5)
plot(data4(:,1)+273.15,data4(:,4),'color','#7E2F8E','LineStyle','-','LineWidth',1.5)
%ylim([-80,0]);
ylabel('Mass loss (%)')
ylim([-70,1])
yyaxis right


plot(data1(:,1)+273.15,data1(:,3),'color','#0072BD','LineStyle','--','LineWidth',1.5)
plot(data2(:,1)+273.15,data2(:,3),'color','#D95319','LineStyle','--','LineWidth',1.5)
plot(data3(:,1)+273.15,data3(:,3),'color','#EDB120','LineStyle','--','LineWidth',1.5)
plot(data4(:,1)+273.15,data4(:,3),'color','#7E2F8E','LineStyle','--','LineWidth',1.5)
ylim([-0.70,0.05])
%plot(data1(:,1),zeros(size(data1(:,3))),'-')
legend('5 K/min','10 K/min','15 K/min','20 K/min','5 K/min','10 K/min','15 K/min','20 K/min','Location','east')
%legend('-','TG','--','DTG','Location','east')
ylabel('Mass loss rate (%/min)')
%grid
%title('Curvas termogravimétricas para % neumáticos')
xlabel('Temperature (K)')

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % 
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
