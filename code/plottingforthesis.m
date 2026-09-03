%%

p1 = polyfit(x1,y1,1);
p5 = polyfit(x5,y5,1);
p9 = polyfit(x9,y9,1);

xx1 = [min(x1) max(x1)];
xx5 = [min(x5) max(x5)];
xx9 = [min(x9) max(x9)];

yy1 = polyval(p1,xx1);
yy5 = polyval(p5,xx5);
yy9 = polyval(p9,xx9);

E1 = E{1};
E5 = E{2};
E9 = E{3};

lnA1 = polyval(p1,E1);
lnA5 = polyval(p5,E5);
lnA9 = polyval(p9,E9);




hfig = figure;
hold on
plot(x9,y9,'color','#0072BD','Marker','o','LineWidth',1.5,'LineStyle','none')
plot(x1,y1,'color','#D95319','Marker','o','LineWidth',1.5,'LineStyle','none')
plot(x5,y5,'color','#EDB120','Marker','o','LineWidth',1.5,'LineStyle','none')

plot(xx9,yy9,'color','#0072BD','LineStyle','-','LineWidth',1.5)
plot(xx1,yy1,'color','#D95319','LineStyle','-','LineWidth',1.5)
plot(xx5,yy5,'color','#EDB120','LineStyle','-','LineWidth',1.5)

plot(E9,lnA9,'MarkerFaceColor','#0072BD','Marker','s','LineWidth',1.5,'LineStyle','none','MarkerSize',10,'MarkerEdgeColor','k')
plot(E1,lnA1,'MarkerFaceColor','#D95319','Marker','s','LineWidth',1.5,'LineStyle','none','MarkerSize',10,'MarkerEdgeColor','k')
plot(E5,lnA5,'MarkerFaceColor','#EDB120','Marker','s','LineWidth',1.5,'LineStyle','none','MarkerSize',10,'MarkerEdgeColor','k')

hold off
xlabel('E_j (kJ/mol)')
ylabel('ln A_j')
legend('Peak 1','Peak 2','Peak 3','Location','southeast')
% xlim([-100 1000])
% legend('Peak 1 - Theoretical pairs','Peak 2 - Theoretical pairs','Peak 3 - Theoretical pairs','Peak 1 - Compensation effect line',...
%     'Peak 2 - Compensation effect line','Peak 3 - Compensation effect line','Peak 1 - Isoconversional pair','Peak 2 - Isoconversional pair'...
%     ,'Peak 3 - Isoconversional pair','Location','southeast')
picturewidth = 20; % set this parameter and keep it forever
% set(gcf,'position',[0, 0, 10, 5]); 
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional

%%

hfig = figure;
hold on
plot(a9,lfa9,'Color','#0072BD','LineStyle','--','LineWidth',1.5)
plot(a9,fitlfa9,'Color','#0072BD','LineStyle','-','LineWidth',1.5)
plot(a1,lfa1,'Color','#D95319','LineStyle','--','LineWidth',1.5)
plot(a1,fitlfa1,'Color','#D95319','LineStyle','-','LineWidth',1.5)
plot(a5,lfa5,'Color','#EDB120','LineStyle','--','LineWidth',1.5)
plot(a5,fitlfa5,'Color','#EDB120','LineStyle','-','LineWidth',1.5)
hold off
legend('Peak 1 - Exp.','Peak 1 - SB model','Peak 2 - Exp.','Peak 2 - SB model','Peak 3 - Exp.','Peak 3 - SB model','Location','southwest')
xlabel('\alpha')
ylabel('ln[f(\alpha)]')
picturewidth = 20; % set this parameter and keep it forever
% set(gcf,'position',[0, 0, 10, 5]); 
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional

