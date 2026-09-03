clear
clc
n_1 = 1;
E_1 = 200;
A_1 = 10^20;

% A_2 = 10^20;
% E_2 = 200*1.1;
% n_2 = 1
 R = 8.314E-3;

dXdt =@(t,X) [  A_1*exp(-E_1/R/X(5))*(1-X(1)).^n_1;
                A_1*exp(-E_1/R/X(6))*(1-X(2)).^n_1;
                A_1*exp(-E_1/R/X(7))*(1-X(3)).^n_1;
                A_1*exp(-E_1/R/X(8))*(1-X(4)).^n_1;
                5;10;15;20];
%+ 0.5*A_2*exp(-E_2/R/X(2))*(1-X(1))^0.5
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

[t,x] = ode15s(dXdt,[0,30],[0.00001,0.00001,0.00001,0.00001, 400, 400, 400, 400],options);

dx1 =  A_1*exp(-E_1/R./x(:,5)).*(1-x(:,1))/5;
dx2 =  A_1*exp(-E_1/R./x(:,6)).*(1-x(:,2))/10;
dx3 =  A_1*exp(-E_1/R./x(:,7)).*(1-x(:,3))/15;
dx4 =  A_1*exp(-E_1/R./x(:,8)).*(1-x(:,4))/20;

figure
%yyaxis left
hold on
plot(x(:,5),x(:,1))
plot(x(:,6),x(:,2))
plot(x(:,7),x(:,3))
plot(x(:,8),x(:,4))
hold off
ylabel('\alpha')
legend('5','10','15','20')

% yyaxis right
figure
hold on
plot(x(:,5),dx1(:,1),'--')
plot(x(:,6),dx2(:,1),'--')
plot(x(:,7),dx3(:,1),'--')
plot(x(:,8),dx4(:,1),'--')
ylabel('d\alpha /dt')


data1 = [x(:,5) t dx1 x(:,1)];
data2 = [x(:,6) t dx1 x(:,2)];
data3 = [x(:,7) t dx1 x(:,3)];
data4 = [x(:,8) t dx1 x(:,4)];

xlabel('Temperatura [K]')
grid on



aa = 0.1:0.1:0.9;
EAS = Vyazovkinv(aa,data1,data2,data3,data4);
figure
errorbar(aa,EAS(:,1),EAS(:,2))
EA = mean(EAS(:,1))
AA = [compeffects(EA,data1) compeffects(EA,data2) compeffects(EA,data3) compeffects(EA,data4)];
AAm = mean(AA)
nm = [sbrlin(AA(1),EA,data1)' ; sbrlin(AA(2),EA,data2)' ; sbrlin(AA(3),EA,data3)' ; sbrlin(AA(4),EA,data4)'];
nmm = mean(nm);


