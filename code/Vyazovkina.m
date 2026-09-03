function EA = Vyazovkina(alphas,data_1,data_2,data_3,data_4)

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

%get massloss

massloss_1 = data_1(:,4);
massloss_2 = data_2(:,4);
massloss_3 = data_3(:,4);
massloss_4 = data_4(:,4);

%calculate conversion
alpha_1 = (massloss_1-massloss_1(1))/(massloss_1(end)-massloss_1(1));
alpha_2 = (massloss_2-massloss_2(1))/(massloss_2(end)-massloss_2(1));
alpha_3 = (massloss_3-massloss_3(1))/(massloss_3(end)-massloss_3(1));
alpha_4 = (massloss_4-massloss_4(1))/(massloss_4(end)-massloss_4(1));

%clean data
[time_1,Temp_1,alpha_1] = removeerror(time_1,Temp_1,alpha_1);
[time_2,Temp_2,alpha_2] = removeerror(time_2,Temp_2,alpha_2);
[time_3,Temp_3,alpha_3] = removeerror(time_3,Temp_3,alpha_3);
[time_4,Temp_4,alpha_4] = removeerror(time_4,Temp_4,alpha_4);


%plot conversions vs time

figure
plot(Temp_1,alpha_1)
hold on
plot(Temp_2,alpha_2)
plot(Temp_3,alpha_3)
plot(Temp_4,alpha_4)
hold off
grid
title('Conversion vs temperature')
legend('1','2','3','4')
xlabel('temperature')
ylabel('conversion')

m=length(alphas);
delta_a = 1/(m+1);

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



%create cfit objects for temperature and time as a function of conversion

T_1 = Tfit(alpha_1,Temp_1);
T_2 = Tfit(alpha_2,Temp_2);
T_3 = Tfit(alpha_3,Temp_3);
T_4 = Tfit(alpha_4,Temp_4);

t_1 = timefit(alpha_1,time_1);
t_2 = timefit(alpha_2,time_2);
t_3 = timefit(alpha_3,time_3);
t_4 = timefit(alpha_4,time_4);
% 
% figure
% plot(T_1(alpha_1),alpha_1)
% hold on
% %plot(Temp_1,alpha_1)
% plot(T_2(alpha_2),alpha_2)
% plot(T_3(alpha_3),alpha_3)
% plot(T_4(alpha_4),alpha_4)
% hold off
% grid
% title('Conversion vs temperature')
% legend('1','2','3','4')
% xlabel('temperature')
% ylabel('conversion')


%alphas = 0.2:0.01:0.8;

%create matrix for times and temperatures

%delta_a = 0.01;
% m=length(alphas);
% delta_a = 1/(m+1);
z=10;
for i = 1:length(alphas)

aa = alphas(i);
tt = [t_1(aa-delta_a:delta_a/z:aa),t_2(aa-delta_a:delta_a/z:aa),t_3(aa-delta_a:delta_a/z:aa),t_4(aa-delta_a:delta_a/z:aa)];
TT = [T_1(aa-delta_a:delta_a/z:aa),T_2(aa-delta_a:delta_a/z:aa),T_3(aa-delta_a:delta_a/z:aa),T_4(aa-delta_a:delta_a/z:aa)];

%perform vyazovkin calculation

E = 1000;

phi =@(e) phiofe2(tt,TT,e);
  options = optimset('MaxFunEvals',100000000,'TolFun',1e-15,'TolX',1e-15);
[Ea,fval] = fminsearch(phi,E,options);

errorp = abs(fval-12)*100/12;
%errorp = fval;
Eea(i,:) = [Ea, errorp];
end
EA = Eea;