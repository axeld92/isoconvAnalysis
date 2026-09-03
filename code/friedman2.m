function EA = friedman2(alphas,data_1,data_2,data_3,data_4)

%Get temperatures
Temp_1 = data_1(:,1)+273.15;
Temp_2 = data_2(:,1)+273.15;
Temp_3 = data_3(:,1)+273.15;
Temp_4 = data_4(:,1)+273.15;

%Get times
time_1 = data_1(:,2)-data_1(1,2);
time_2 = data_2(:,2)-data_2(1,2);
time_3 = data_3(:,2)-data_3(1,2);
time_4 = data_4(:,2)-data_4(1,2);

%get massloss

massloss_1 = data_1(:,4);
massloss_2 = data_2(:,4);
massloss_3 = data_3(:,4);
massloss_4 = data_4(:,4);

%calculate conversion
alpha_1 = (massloss_1(1)-massloss_1)/(massloss_1(1)-massloss_1(end));
alpha_2 = (massloss_2(1)-massloss_2)/(massloss_2(1)-massloss_2(end));
alpha_3 = (massloss_3(1)-massloss_3)/(massloss_3(1)-massloss_3(end));
alpha_4 = (massloss_4(1)-massloss_4)/(massloss_4(1)-massloss_4(end));

%get rid of data outside of alpha = 0.05 to alpha = 0.095

ind = getind(alpha_1);
ind = logical(ind);
alpha_1 = alpha_1(ind);
Temp_1 = Temp_1(ind);
time_1 = time_1(ind);

ind = getind(alpha_2);
ind = logical(ind);
alpha_2 = alpha_2(ind);
Temp_2 = Temp_2(ind);
time_2 = time_2(ind);

ind = getind(alpha_3);
ind = logical(ind);
alpha_3 = alpha_3(ind);
Temp_3 = Temp_3(ind);
time_3 = time_3(ind);

ind = getind(alpha_4);
ind = logical(ind);
alpha_4 = alpha_4(ind);
Temp_4 = Temp_4(ind);
time_4 = time_4(ind);

beta = [(Temp_1(end)-Temp_1(1))/(time_1(end)-time_1(1)) , 
        (Temp_2(end)-Temp_2(1))/(time_2(end)-time_2(1)) , 
        (Temp_3(end)-Temp_3(1))/(time_3(end)-time_3(1)) , 
        (Temp_4(end)-Temp_4(1))/(time_4(end)-time_4(1))];

dadT_1 = diff_alpha(Temp_1,alpha_1);
dadT_2 = diff_alpha(Temp_2,alpha_2);
dadT_3 = diff_alpha(Temp_3,alpha_3);
dadT_4 = diff_alpha(Temp_4,alpha_4);

%figure
% hold on
% plot(Temp_1,dadT1(Temp_1))
% plot(Temp_2,dadT1(Temp_2))
% plot(Temp_3,dadT1(Temp_3))
% plot(Temp_4,dadT1(Temp_4))
% hold off
% grid
% title('d\alpha/dT vs T')

R = 8.314;
%alphas = 0.2:0.01:0.8;
figure
hold on
for i = 1:length(alphas)
    aa = alphas(i);

TT = [T_1(aa) T_2(aa) T_3(aa) T_4(aa)];
dadTT = [dadT1(TT(1)) dadT2(TT(2)) dadT3(TT(3)) dadT4(TT(4))];

x = [1/TT(1) ; 1/TT(2) ; 1/TT(3) ;1/TT(4)];
    y = [log(beta(1)*dadTT(1)); log(beta(2)*dadTT(2)) ; log(beta(3)*dadTT(3)) ; log(beta(4)*dadTT(4))];
    [p,S]=polyfit(x,y,1);
    ci = polyparci(p,S);
    y1 = polyval(p,x);
    plot(x,y,'o')
    plot(x,y1)
    Ee(i) = -p(1)*8.314/1000;
    CI(i) = abs(-ci(2,1)*8.314/1000-Ee(i));
end
EA = [Ee' CI'];
hold off
grid
title('Friedman')
xlabel('T_{\alpha,i}^{-1}')
ylabel('ln(\beta_i (d\alpha/dT))')



