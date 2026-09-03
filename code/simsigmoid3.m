clear
clc
n_1 = 0;
n_2 = 1;
C = 1;
E_1 = 200;
A_1 = 10^15;

% A_2 = 10^20;
% E_2 = 200*1.1;
% n_2 = 1
 R = 8.314E-3;
b = [5 10 15 20];
T0 = 500;
dXdt =@(t,X) [  A_1*exp(-E_1/R/X(5))*C.*X(1)^n_1.*(1-X(1)).^n_2;
                A_1*exp(-E_1/R/X(6))*C.*X(2)^n_1.*(1-X(2)).^n_2;
                A_1*exp(-E_1/R/X(7))*C.*X(3)^n_1.*(1-X(3)).^n_2;
                A_1*exp(-E_1/R/X(8))*C.*X(4)^n_1.*(1-X(4)).^n_2;
                b(1);b(2);b(3);b(4)];
%+ 0.5*A_2*exp(-E_2/R/X(2))*(1-X(1))^0.5
options = odeset('RelTol',1e-20,'AbsTol',1e-20);

[t,x] = ode15s(dXdt,[0 50],[0.0000001,0.0000001,0.0000001,0.0000001,T0,T0,T0,T0],options);

%global dx1 dx2 dx3 dx4 a1 a2 a3 a4 T1 T2 T3 T4
% dx1 =  A_1*exp(-E_1/R./x(:,5)).*(1-x(:,1));%/5;
% dx2 =  A_1*exp(-E_1/R./x(:,6)).*(1-x(:,2));%/10;
% dx3 =  A_1*exp(-E_1/R./x(:,7)).*(1-x(:,3));%/15;
% dx4 =  A_1*exp(-E_1/R./x(:,8)).*(1-x(:,4));%/20;

dx1 =  A_1*exp(-E_1/R./x(:,5)).*(1-x(:,1))/5;
dx2 =  A_1*exp(-E_1/R./x(:,6)).*(1-x(:,2))/10;
dx3 =  A_1*exp(-E_1/R./x(:,7)).*(1-x(:,3))/15;
dx4 =  A_1*exp(-E_1/R./x(:,8)).*(1-x(:,4))/20;


a1 = x(:,1);
a2 = x(:,2);
a3 = x(:,3);
a4 = x(:,4);

T1 = x(:,5);
T2 = x(:,6);
T3 = x(:,7);
T4 = x(:,8);

% 
% yy = [log(dx1)-log(a1.^n_1.*(1-a1).^n_2);log(dx2)-log(a2.^n_1.*(1-a2).^n_2)];
% rT = [1./T1 ; 1./T2];
%  plot(rT,yy,'.')

 %%

figure
%yyaxis left
hold on
plot(t,x(:,1))
plot(t,x(:,2))
plot(t,x(:,3))
plot(t,x(:,4))
hold off
ylabel('\alpha')
legend('5','10','15','20')

%%

% figure
% %yyaxis left
% hold on
% plot(t,x(:,1))
% plot(t,x(:,2))
% plot(t,x(:,3))
% plot(t,x(:,4))
% hold off
% ylabel('\alpha')
% legend('5','10','15','20')

% % yyaxis right
figure
hold on
plot(1./x(:,5),log(dx1(:,1)),'.')
plot(1./x(:,6),log(dx2(:,1)),'.')
plot(1./x(:,7),log(dx3(:,1)),'.')
plot(1./x(:,8),log(dx4(:,1)),'.')
ylabel('d\alpha /dt')
xlabel('Temperatura [K]')
grid on
hold off

% figure
% hold on
% plot(t,dx1(:,1),'--')
% plot(t,dx2(:,1),'--')
% plot(t,dx3(:,1),'--')
% plot(t,dx4(:,1),'--')
% ylabel('d\alpha /dt')
% xlabel('Temperatura [K]')
% grid on
% hold off

%%

data1 = [T1 t dx1 a1];
data2 = [T2 t dx1 a2];
data3 = [T3 t dx1 a3];
data4 = [T4 t dx1 a4];



lower = 0.2;
upper = 0.8;

     idx = (data1(:,4)>lower).*(data1(:,4)<upper);
     idx = logical(idx);
     data1 = data1(idx,:);

     idx = (data2(:,4)>lower).*(data2(:,4)<upper);
     idx = logical(idx);
     data2 = data2(idx,:);
     
     idx = (data3(:,4)>lower).*(data3(:,4)<upper);
     idx = logical(idx);
     data3 = data3(idx,:);
     
     idx = (data4(:,4)>lower).*(data4(:,4)<upper);
     idx = logical(idx);
     data4 = data4(idx,:);

     
%%
  

data1 = out{1}{2};
data2 = out{2}{2};
data3 = out{3}{2};
data4 = out{4}{2};


lower = 0.1;
upper = 0.9;

     idx = (data1(:,4)>lower).*(data1(:,4)<upper);
     idx = logical(idx);
     data1 = data1(idx,:);

     idx = (data2(:,4)>lower).*(data2(:,4)<upper);
     idx = logical(idx);
     data2 = data2(idx,:);
     
     idx = (data3(:,4)>lower).*(data3(:,4)<upper);
     idx = logical(idx);
     data3 = data3(idx,:);
     
     idx = (data4(:,4)>lower).*(data4(:,4)<upper);
     idx = logical(idx);
     data4 = data4(idx,:);

%%

dx1 = data1(:,3);
dx2 = data2(:,3);
dx3 = data3(:,3);
dx4 = data4(:,3);

a1 = data1(:,4);
a2 = data2(:,4);
a3 = data3(:,4);
a4 = data4(:,4);

T1 = data1(:,1);
T2 = data2(:,1);
T3 = data3(:,1);
T4 = data4(:,1);




repT = [1./T1 ; 1./T2 ; 1./T3 ; 1./T4];
yy1 = logfun(0,1,data1);
yy2 = logfun(0,1,data2);
yy3 = logfun(0,1,data3);
yy4 = logfun(0,1,data4);

pp = polyfit(1./T4,yy4,1);
E = -pp(1)*8.314E-3

%%
figure
hold on
plot(1./T1,yy1,'.')
plot(1./T2,yy2,'.')
plot(1./T3,yy3,'.')
plot(1./T4,yy4,'.')
hold off
 
%%
%prueba = logfun(1,0,data1,data2,data3,data4);


%figure1
%plot(repT,prueba,'.')

objfun =@(mn) -rfromreg(mn(1),mn(2),repT,data1,data2,data3,data4);

% optoptions = optimset('TolFun',1E-8,'TolX',1E-8,'Display','iter','PlotFcns',@optimplotfval);
% [mno,ropt] = fminsearch(objfun,[0,1],optoptions);
optionsps = optimoptions('patternsearch','TolMesh',1e-10,'TolFun',1e-10,'Display','iter','PlotFcn',@psplotbestf);
[mno,ropt] = patternsearch(objfun,[10,10],[],[],[],[],[0,0],[],[],optionsps);
%       x = patternsearch(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
 prueba = logfun(mno(1),mno(2),data1,data2,data3,data4);
 
 figure
 plot(repT,prueba,'.')

%rsqrd = objfunct(0.5,1)

function l = logfun(m,n,varargin)
num = length(varargin);
data = varargin;
l = [];
for i=1:num
T{i} = data{i}(:,1);
da{i} = data{i}(:,3);
a{i} = data{i}(:,4);
l = [l ; log(da{i})-log(a{i}.^m.*(1-a{i}).^n)];
end

end

function r = rfromreg(m,n,repT,varargin)
data = varargin;
y_i = logfun(m,n,data{1},data{2},data{3},data{4});
y_bar = mean(y_i);
[p,S] = polyfit(repT,y_i,1);
y_fit = polyval(p,repT);
Rsqrd = 1 - (S.normr/norm(y_i - mean(y_i)))^2;
%r = 1-sum((y_i-y_fit).^2)/sum((y_i-y_bar).^2);
r = Rsqrd
E = -p(1)*8.314E-3
lncA = p(2)
end




