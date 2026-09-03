
function mnpc = sb(A,Ea,data)

ind = (data(:,4)>0.1).*(data(:,4)<0.9);
ind = logical(ind);
data = data(ind,:);

T = data(:,1);
t = data(:,2);
da = data(:,3);
a = data(:,4);

%plot(T,a)


b = (T(end)-T(1))/(t(end)-t(1));

R = 8.314/1000;

% Calculate model from data
fa = b*da./(A*exp(-Ea./(R*T)));

%define sestak berggren model

%sesberg =@(alpha,mnp) alpha.^mnp(1).*(1-alpha).^mnp(2) ;

%mnp0 = [5;4];


obj_fun =@(mnp) sum((fa - sesberg(a,mnp)).^2) ;

lb = [-100,-100,-100];
ub = [100,100,100];

rng default
opts = optimoptions('particleswarm','SwarmSize',1000,'FunctionTolerance',1E-10,'MaxStallIterations',50,'MaxIterations',10000);

[mnpc , fval] = particleswarm(obj_fun,3,lb,ub,opts);

mnpc
fval


% options = optimset('MaxFunEvals',100000000,'TolFun',1e-20,'TolX',1e-15,'MaxIter',1000000,'PlotFcns',@optimplotfval);
% mnpc = fminsearch(obj_fun,mnp0,options);

figure
hold on
plot(a,fa)
% plot(a,sesberg(a,mnp0))
plot(a,(sesberg(a,mnpc)),'--')
legend('fa','sb0','sb')


    function  s = sesberg(alpha,mnp)
s = alpha.^mnp(1).*(1-alpha).^mnp(2).*(-log(1-alpha)).^mnp(3) ;
end
     
%     function s = sesberg(alpha,mnp)
%         s = mnp(1)*log(alpha) + mnp(2)*log(1-alpha) - mnp(3)*log(log(1-alpha));
%          s(s~=real(s))=0;
%     end
end
