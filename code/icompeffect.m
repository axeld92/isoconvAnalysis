function A = icompeffect(Ea,data)

ind = (data(:,4)>0.05).*(data(:,4)<0.9);
ind = logical(ind);
data = data(ind,:);


T = data(:,1);
t = data(:,2);
da = data(:,3);
a = data(:,4);

b = (T(end)-T(1))/(t(end)-t(1));
%b = 1;

a_1 = a.^(1/4); % powerlaw 1
a_2 = a.^(1/3); % powerlaw 2
a_3 = a.^(1/2); % powerlaw 3
a_4 = a.^(3/2); % powerlaw 4
a_5 = a.^2 ; %1D diffusion
a_6 = -log(1-a); % Mampel 1st order
a_7 = (-log(1-a)).^(1/4); % Avrami-Erofeev 1
a_8 = (-log(1-a)).^(1/3); % Avrami-Erofeev 2
a_9 = (-log(1-a)).^(1/2); % Avrami-Erofeev 3
a_10 = (1-(1-a).^(1/3)).^2; %3D diffusion
a_11 = 1-(1-a).^(1/3); %contracting sphere
a_12 = 1-(1-a).^(1/2); %Contracting cylinder

% 
% figure
% hold on
% plot(a,fa_1,a,fa_2,a,fa_3,a,fa_4,a,fa_5)
% legend('1','2','3','4','5')


R = 8.314/1000;
EA0 = [10^15,Ea];




options = optimset('MaxFunEvals',100000000,'MaxIter',1E5,'PlotFcn',@optimplotfval);
    
    ga_1 =@(EA) cumtrapz(a_1,EA(1)*exp(-EA(2)/R./T)/b./da);
    ga_2 =@(EA) cumtrapz(a_2,EA(1)*exp(-EA(2)/R./T)/b./da);
    ga_3 =@(EA) cumtrapz(a_3,EA(1)*exp(-EA(2)/R./T)/b./da);
    ga_4 =@(EA) cumtrapz(a_4,EA(1)*exp(-EA(2)/R./T)/b./da);
    ga_5 =@(EA) cumtrapz(a_5,EA(1)*exp(-EA(2)/R./T)/b./da);
    ga_6 =@(EA) cumtrapz(a_6,EA(1)*exp(-EA(2)/R./T)/b./da);
    
    
    difference1 =@(EA) (a-ga_1(EA)).^2;
    difference2 =@(EA) (a-ga_2(EA)).^2;
    difference3 =@(EA) (a-ga_3(EA)).^2;
    difference4 =@(EA) (a-ga_4(EA)).^2;
    difference5 =@(EA) (a-ga_5(EA)).^2;
    difference6 =@(EA) (a-ga_6(EA)).^2;
    
    obj_fun1 =@(EA) sum(difference1(EA));
    obj_fun2 =@(EA) sum(difference2(EA));
    obj_fun3 =@(EA) sum(difference3(EA));
    obj_fun4 =@(EA) sum(difference4(EA));
    obj_fun5 =@(EA) sum(difference5(EA));
    obj_fun6 =@(EA) sum(difference6(EA));
    
    EAA1 = fminsearch(obj_fun1,EA0,options);
    EAA2 = fminsearch(obj_fun2,EA0,options);
    EAA3 = fminsearch(obj_fun3,EA0,options);
    EAA4 = fminsearch(obj_fun4,EA0,options);
    EAA5 = fminsearch(obj_fun5,EA0,options);
    EAA6 = fminsearch(obj_fun6,EA0,options);
    
    figure
    hold on
    plot(a,da)
    plot(a,ga_1(EAA1))
    plot(a,ga_2(EAA2))
    plot(a,ga_3(EAA3))
    plot(a,ga_4(EAA4))
    plot(a,ga_5(EAA5))
    plot(a,ga_6(EAA6))
    hold off
    legend('data','1','2','3','4','5','6')

    
    
    
    
    
    
   x =      [EAA1(2) EAA2(2) EAA3(2) EAA4(2) EAA5(2) EAA6(2)]
   y = log( [EAA1(1) EAA2(1) EAA3(1) EAA4(1) EAA5(1) EAA6(1)])
   
   
   
   p=polyfit(x,y,1);
   logA = polyval(p,Ea);
   figure
   plot(x,y,'o',Ea,logA,'*')
   A = exp(logA);
%    
%    figure
%    hold on
%    plot(T,da);
%    plot(T,dT_1c([A,Ea]));
%    plot(T,dT_2c(EAA6));
%    hold off









