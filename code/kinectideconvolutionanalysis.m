data = St16;
data = cleandata(data);

% 
t = data(:,2);
T = data(:,1);
a = data(:,4);
b = (Tf-T0)/(tf-t0);

objfun =@(p) sum((a - ratelaw(t,T,p)).^2);

p0 = [c(1)*A_1,E{1},nmp_1(1),nmp_1(2),c(2)*A_2,E{2},nmp_2(1),nmp_2(2),c(3)*A_3,E{3},nmp_3(1),nmp_3(2)];

% idx = p0<0;
% lb = zeros(size(p0));
% lb(idx) = p0(idx) + p0(idx)*0.5;
% lb(~idx) = p0(~idx) - p0(~idx)*0.5;
% 
% ub = zeros(size(p0));
% ub(idx) = p0(idx) - p0(idx)*0.5;
% ub(~idx) = p0(~idx) + p0(~idx)*0.5;
% 
% options = optimset('PlotFcns',@optimplotfval);
% %opts = optimoptions('particleswarm','SwarmSize',1000);
%  %[fitparam , fval] = particleswarm(objfun,12,lb,ub,opts);
%  [fitparam , fval] = fminsearch(objfun,p0,options);
%  checklb = fitparam == ub
%  checkub = fitparam == lb

respuesta = ratelaw(t,T,p0);




function rate = ratelaw(t,T,params)
cA1 = params(1);
E1 = params(2);
n1 = params(3);
m1 = params(4);
cA2 = params(5);
E2 = params(6);
n2 = params(7);
m2 = params(8);
cA3 = params(9);
E3 = params(10);
n3 = params(11);
m3 = params(12);
R = 8.314E-3;
t0 = t(1);
tf = t(end);
T0 = T(1);
Tf = T(end);
b = (Tf-T0)/(tf-t0);

drate =@(t,X)    [cA1.*exp(-E1/R./X(4)).*X(1).^n1.*(1-X(1)).^m1;
                cA2.*exp(-E2/R./X(4)).*X(2).^n2.*(1-X(2)).^m2;
                cA3.*exp(-E3/R./X(4)).*X(3).^n3.*(1-X(3)).^m3;
                b];
 options = odeset('RelTol',1e-8,'AbsTol',1e-8) ;            
 [t,X] = ode15s(drate,[0 tf],[1E-3 1E-3 1E-3 T0]);
 rate = X(:,1) + X(:,2)+ X(:,3) ;
 
end