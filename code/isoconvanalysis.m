%Script general para el análisis isoconversional
function isoconvanalysis(data_1,data_2,data_3,data_4)

%limpiar datos
data_1 = cleandata(data_1);
data_2 = cleandata(data_2);
data_3 = cleandata(data_3);
data_4 = cleandata(data_4);

%energia de activacion segun vyazovkin
alphas = 0.15:0.05:0.9;
EA_s = Starink(alphas,data_1,data_2,data_3,data_4);
EA_v = Vyazovkin(alphas,data_1,data_2,data_3,data_4);

hfig = figure;
hold on
errorbar(alphas,EA_s(:,1),EA_s(:,2),'-o','LineWidth',1.5)
errorbar(alphas,EA_v(:,1),EA_v(:,2),'-^','LineWidth',1.5)
% plot(alphas,EA_s(:,1))
% plot(alphas,EA_v(:,1))
ylim([0 800])

%grid;
%title('Energía de activación para  % neumáticos')
ylabel('E_{\alpha} [kJ/mol]')
xlabel('\alpha')
legend('Starink','Vyazovkin')
picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional



 