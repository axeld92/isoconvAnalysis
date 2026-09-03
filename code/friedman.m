function EA = friedman(alphas,data_1,data_2,data_3,data_4)

lower = 0.01;
upper = 0.99;


ind = (data_1(:,4)>lower).*(data_1(:,4)<upper);ind = logical(ind);data_1 = data_1(ind,:);
ind = (data_2(:,4)>lower).*(data_2(:,4)<upper);ind = logical(ind);data_2 = data_2(ind,:);
ind = (data_3(:,4)>lower).*(data_3(:,4)<upper);ind = logical(ind);data_3 = data_3(ind,:);
ind = (data_4(:,4)>lower).*(data_4(:,4)<upper);ind = logical(ind);data_4 = data_4(ind,:);




%Get temperatures
Temp_1 = data_1(:,1);
Temp_2 = data_2(:,1);
Temp_3 = data_3(:,1);
Temp_4 = data_4(:,1);

%Get times
time_1 = data_1(:,2)-data_1(1,2);
time_2 = data_2(:,2)-data_2(1,2);
time_3 = data_3(:,2)-data_3(1,2);
time_4 = data_4(:,2)-data_4(1,2);


beta = [(Temp_1(end)-Temp_1(1)/(time_1(end)-time_1(1))) ...
    (Temp_2(end)-Temp_2(1)/(time_2(end)-time_2(1))) ... 
    (Temp_3(end)-Temp_3(1)/(time_3(end)-time_3(1))) ... 
    (Temp_4(end)-Temp_4(1)/(time_4(end)-time_4(1)))];

%get massloss

a_1 = data_1(:,4);
a_2 = data_2(:,4);
a_3 = data_3(:,4);
a_4 = data_4(:,4);

%calculate conversion
da_1 = data_1(:,3);
da_2 = data_2(:,3);
da_3 = data_3(:,3);
da_4 = data_4(:,3);


R = 8.314;
%alphas = 0.2:0.01:0.8;
figure
hold on
for i = 1:length(alphas)
    aa = alphas(i);

TT = [interp1(a_1,Temp_1,aa) interp1(a_2,Temp_2,aa) interp1(a_3,Temp_3,aa) interp1(a_4,Temp_4,aa)];
dadTT = [interp1(a_1,da_1,aa) interp1(a_2,da_2,aa) interp1(a_3,da_3,aa) interp1(a_4,Temp_4,aa)];

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



