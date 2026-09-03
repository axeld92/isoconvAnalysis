function plotmasslossi(data1,data2,data3)

data1 = cleandatami(data1);
data2 = cleandatami(data2);
data3 = cleandatami(data3);
% data4 = cleandatam(data4);

figure
hold on
plot(data1(:,2),data1(:,4))
plot(data2(:,2),data2(:,4))
plot(data3(:,2),data3(:,4))
% plot(data4(:,1),data4(:,4))


yyaxis right


plot(data1(:,2),data1(:,3))
plot(data2(:,2),data2(:,3),':')
plot(data3(:,2),data3(:,3),'-.')
% plot(data4(:,),data4(:,3),'--')
% plot(data1(:,2),zeros(size(data1(:,3))),'-')
legend('25','50','75','25','50','75')

grid
title('curvas termogravimetricas isotermicas ')