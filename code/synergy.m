function synergy(data_1,data_2,data_3,x)

data_1 = cleandata(data_1);
data_2 = cleandata(data_2);
data_3 = cleandata(data_3);

Temp_1 = data_1(:,1);
Temp_2 = data_2(:,1);
Temp_3 = data_3(:,1);

a_1 = data_1(:,4);
a_2 = data_2(:,4);
a_3 = data_3(:,4);

dadT_1 = data_1(:,3);
dadT_2 = data_2(:,3);
dadT_3 = data_3(:,3);




aa_1 = interp1(Temp_1,a_1,Temp_3);
aa_2 = interp1(Temp_2,a_2,Temp_3);
aa_3 = interp1(Temp_3,a_3,Temp_3);
daadT_1 = interp1(Temp_1,dadT_1,Temp_3);
daadT_2 = interp1(Temp_2,dadT_2,Temp_3);
daadT_3 = interp1(Temp_3,dadT_3,Temp_3);


aa_s = (1-x)*aa_1 + x*aa_2;
daa_s = (1-x)*daadT_1 + x*daadT_2;

dif = aa_3-aa_s;
ddif = (daadT_3 - daa_s);%./daadT_3;%./daa_s;

hfig = figure ;
subplot(5,1,1:4);
hold on
plot(Temp_3,daadT_3,'LineWidth',1.5);
plot(Temp_3,daa_s,'--','LineWidth',1.5);
hold off
legend('Real','Teórico')
%grid
ylabel('d\alpha/dT')
subplot(5,1,5);
plot(Temp_3,ddif,'LineWidth',1.5)
%ylim([-0.5 0.5])
xlabel('Temperatura [K]')
%grid 
picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% 
% figure 
% subplot(5,1,1:4);
% hold on
% plot(Temp_3,aa_3);
% plot(Temp_3,aa_s,'--');
% hold off
% legend('Real','Teórico')
% grid
% ylabel('d\alpha/dT')
% subplot(5,1,5);
% plot(Temp_3,dif)
% %ylim([-0.5 0.5])
% xlabel('Temperatura [K]')
% grid 







%plot(Temp_3,a_1(Temp_3),Temp_3,a_2(Temp_3),Temp_3,a_3(Temp_3),Temp_2,aa_s)
% figure
% hold on
% plot(Temp_1,a_1)
% plot(Temp_2,a_2)
% plot(Temp_3,a_3)
% plot(Temp_3,aa_s)
% hold off
% legend('Carozo','Neumáticos', 'Real','Teórico')
% grid


% difference =100*(a_3./aa_s - 1) ;
% difference = smoothdata(difference,'sgolay');
%  figure
%  hold on
% plot(a_3,difference)
%  xlim([0.05 0.95])
%  hold off
% grid
% 
% ddifference = 100*(dadT_3-daa_s)./daa_s;
% ddifference = smoothdata(ddifference,'sgolay');
% figure
% hold on
% plot(a_3,ddifference)
% xlim([0.05 0.95])
% hold off
% grid
% 
% % plot(Temp_3,dadT_3(Temp_3))
% % xlim([500 950])
