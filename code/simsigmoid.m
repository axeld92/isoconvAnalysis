clear
clc
n_1 = 1;
E_1 = 200;
A_1 = 10^20;

A_2 = 10^20;
E_2 = 200*1.1;
n_2 = 1;
R = 8.314E-3;

dXdt =@(t,X) [  A_1*exp(-E_1/R/X(3))*X(1)*(1-X(1)).^n_1;
                A_2*exp(-E_2/R/X(3))*(1-X(2)).^n_2;
                10];
%+ 0.5*A_2*exp(-E_2/R/X(2))*(1-X(1))^0.5
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

[t,x] = ode15s(dXdt,[0,50],[0.00001,0.00001,400],options);

dx1 =  A_1*exp(-E_1/R./x(:,3)).*(1-x(:,1));
dx2 =  A_2*exp(-E_2/R./x(:,3)).*(1-x(:,2));

figure
%yyaxis left
hold on
plot(x(:,3),x(:,1))
plot(x(:,3),x(:,2))
hold off
ylabel('\alpha')

% yyaxis right
% 
% plot(x(:,3),dx1(:,1),'--')
% plot(x(:,3),dx2(:,1),'--')
% ylabel('d\alpha /dt')
% 
% legend('\alpha','d\alpha /dt')
xlabel('Temperatura [K]')
grid on
