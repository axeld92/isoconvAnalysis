function plotconversion(data1,data2,data3,data4)

data1 = cleandata(data1);
data2 = cleandata(data2);
data3 = cleandata(data3);
data4 = cleandata(data4);

figure
hold on
plot(data1(:,2),data1(:,4))
plot(data2(:,2),data2(:,4))
plot(data3(:,2),data3(:,4))
plot(data4(:,2),data4(:,4))
legend('1','2','3','4')

yyaxis right

% 
% plot(data1(:,1),data1(:,3))
% plot(data2(:,1),data2(:,3))
% plot(data3(:,1),data3(:,3))
% plot(data4(:,1),data4(:,3))
% plot(data1(:,1),zeros(size(data1(:,3))),'-')
% legend('5K/min','10K/min','15k/min','20K/min','5K/min','10K/min','15k/min','20K/min','zero')
% grid
% title('Curvas de conversión para % neumáticos')