data_1 = out_3_05;
data_2 = out_3_10;
data_3 = out_3_15;
data_4 = out_3_20;

A = A_3;

Ea = ea3_avg;


T_1 = data_1(:,1);
t_1 = data_1(:,2);
da_1 = data_1(:,3);
a_1 = data_1(:,4);
b_1 = (T_1(end)-T_1(1))/(t_1(end)-t_1(1));

T_2 = data_2(:,1);
t_2 = data_2(:,2);
da_2 = data_2(:,3);
a_2 = data_2(:,4);
b_2 = (T_2(end)-T_2(1))/(t_2(end)-t_2(1));

T_3 = data_3(:,1);
t_3 = data_3(:,2);
da_3 = data_3(:,3);
a_3 = data_3(:,4);
b_3 = (T_3(end)-T_3(1))/(t_3(end)-t_3(1));

T_4 = data_4(:,1);
t_4 = data_4(:,2);
da_4 = data_4(:,3);
a_4 = data_4(:,4);
b_4 = (T_4(end)-T_4(1))/(t_4(end)-t_4(1));

R = 8.314;

% Calculate model from data
fa_1 = da_1./(A*exp(-Ea./(R*T_1)));
fa_2 = da_2./(A*exp(-Ea./(R*T_2)));
fa_3 = da_3./(A*exp(-Ea./(R*T_3)));
fa_4 = da_4./(A*exp(-Ea./(R*T_4)));
%a = 0:0.01:1;
%a=a';

%fa = a.^34.*(1-a).^27;
figure
hold on
plot(a_1,fa_1)
plot(a_2,fa_2)
plot(a_3,fa_3)
plot(a_4,fa_4)
hold off
legend('5','10','15','20')
