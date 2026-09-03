function EA = Vyazovkinfor3(alphas,data_1,data_2,data_3)

lower = 0.05;
upper = 0.95;


ind = (data_1(:,4)>lower).*(data_1(:,4)<upper);ind = logical(ind);data_1 = data_1(ind,:);
ind = (data_2(:,4)>lower).*(data_2(:,4)<upper);ind = logical(ind);data_2 = data_2(ind,:);
ind = (data_3(:,4)>lower).*(data_3(:,4)<upper);ind = logical(ind);data_3 = data_3(ind,:);


%Get temperatures
Temp_1 = data_1(:,1);
Temp_2 = data_2(:,1);
Temp_3 = data_3(:,1);


%Get times
time_1 = data_1(:,2);
time_2 = data_2(:,2);
time_3 = data_3(:,2);

%get massloss

% massloss_1 = data_1(:,4);
% massloss_2 = data_2(:,4);
% massloss_3 = data_3(:,4);
% massloss_4 = data_4(:,4);

alpha_1 = data_1(:,4);
alpha_2 = data_2(:,4);
alpha_3 = data_3(:,4);


%create cfit objects for temperature and time as a function of conversion

T_1 = Tfit(alpha_1,Temp_1);
T_2 = Tfit(alpha_2,Temp_2);
T_3 = Tfit(alpha_3,Temp_3);


t_1 = timefit(alpha_1,time_1);
t_2 = timefit(alpha_2,time_2);
t_3 = timefit(alpha_3,time_3);

figure
plot(T_1(alpha_1),alpha_1)
hold on
plot(T_2(alpha_2),alpha_2)
plot(T_3(alpha_3),alpha_3)
hold off
grid
title('Conversion vs temperature')
legend('1','2','3')
xlabel('temperature')
ylabel('conversion')




%alphas = 0.2:0.01:0.8;

%create matrix for times

m=length(alphas);
delta_a = 0.02;
z=10;
Eea = zeros(length(alphas),2);

for i = 1:length(alphas)

aa = alphas(i);
tt = [t_1(aa-delta_a:delta_a/z:aa),t_2(aa-delta_a:delta_a/z:aa),t_3(aa-delta_a:delta_a/z:aa)];
TT = [T_1(aa-delta_a:delta_a/z:aa),T_2(aa-delta_a:delta_a/z:aa),T_3(aa-delta_a:delta_a/z:aa)];

%perform vyazovkin calculation

E = 200;

phi =@(e) abs(phiofe3(tt,TT,e)-12);
%  options = optimset('MaxFunEvals',10000000,'TolFun',1e-15,'TolX',1e-15);
[Ea,fval] = fminsearch(phi,E);

%[Ea,fval] = particleswarm(phi,1);


errorp = fval*100/12;
%errorp = fval;
Eea(i,:) = [Ea, errorp];
end
EA = Eea;