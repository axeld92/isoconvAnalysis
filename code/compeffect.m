function A = compeffect(Ea,data)

ind = (data(:,4)>0.05).*(data(:,4)<0.9);
ind = logical(ind);
data = data(ind,:);


T = data(:,1);
t = data(:,2);
da = data(:,3);
a = data(:,4);
da = log(da);
b = (T(end)-T(1))/(t(end)-t(1));
%b = 1;




fa_1 = 4*a.^(3/4); %PL

fa_2 = 3*a.^(2/3); %PL

fa_3 = 2*a.^(1/2); %PL

fa_4 = 2/3*a.^(-1/2); %PL

fa_5 = (a.^(-1))/2; %1D diffusion

fa_6 = 1-a; %1st order

fa_7 = 4*(1-a).*(-log(1-a)).^(3/4); %A4

fa_8 = 3*(1-a).*(-log(1-a)).^(2/3); %A3

fa_9 = 2*(1-a).*(-log(1-a)).^(1/2); %A2

fa_10 = 3/2 * (1-a).^(2/3).*(1-(1-a).^(1/3)).^(-1) ; %3D diffusion

fa_11 = 3*(1-a).^(2/3); % Contracting sphere

fa_12 = 2*(1-a).^(1/2); % contracting cylinder

fa_13 = (-log(1-a)).^(-1); % 2D diffusion

fa_7(fa_7~=real(fa_7)) = 0;
fa_8(fa_8~=real(fa_8)) = 0;
fa_9(fa_9~=real(fa_9)) = 0;
fa_13(fa_13~=real(fa_13)) = 0;




% figure
% hold on
% plot(a,fa_1,a,fa_2,a,fa_3,a,fa_4,a,fa_5,a,fa_6,a,fa_7,a,fa_8,a,fa_9,a,fa_10,a,fa_11,a,fa_12,a,fa_13)
% legend('1','2','3','4','5','6','7','8','9','10','11','12','13')


R = 8.314E-3;
EA0 = [10^15,Ea];





    
    dT_1c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_1;
    dT_2c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_2;
    dT_3c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_3;
    dT_4c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_4;
    dT_5c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_5;
    dT_6c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_6;
    dT_7c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_7;
    dT_8c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_8;
    dT_9c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_9;
    dT_10c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_10;
    dT_11c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_11;
    dT_12c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_12;
    dT_13c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_13;


    dT_1c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_1/b);
    dT_2c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_2/b);
    dT_3c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_3/b);
    dT_4c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_4/b);
    dT_5c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_5/b);
    dT_6c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_6/b);
    dT_7c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_7/b);
    dT_8c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_8/b);
    dT_9c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_9/b);
    dT_10c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_10/b);
    dT_11c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_11/b);
    dT_12c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_12/b);
    dT_13c =@(EA) EA(1) - EA(2)./(R*T) + log(fa_13/b);
    
    
    
    difference1 =@(EA) (da-dT_1c(EA)).^2;
    difference2 =@(EA) (da-dT_2c(EA)).^2;
    difference3 =@(EA) (da-dT_3c(EA)).^2;
    difference4 =@(EA) (da-dT_4c(EA)).^2;
    difference5 =@(EA) (da-dT_5c(EA)).^2;
    difference6 =@(EA) (da-dT_6c(EA)).^2;
    difference7 =@(EA) (da-dT_7c(EA)).^2;
    difference8 =@(EA) (da-dT_8c(EA)).^2;
    difference9 =@(EA) (da-dT_9c(EA)).^2;
    difference10 =@(EA) (da-dT_10c(EA)).^2;
    difference11 =@(EA) (da-dT_11c(EA)).^2;
    difference12 =@(EA) (da-dT_12c(EA)).^2;
    difference13 =@(EA) (da-dT_13c(EA)).^2;
    
    obj_fun1 =@(EA) sum(difference1(EA));
    obj_fun2 =@(EA) sum(difference2(EA));
    obj_fun3 =@(EA) sum(difference3(EA));
    obj_fun4 =@(EA) sum(difference4(EA));
    obj_fun5 =@(EA) sum(difference5(EA));
    obj_fun6 =@(EA) sum(difference6(EA));
    obj_fun7 =@(EA) sum(difference7(EA));
    obj_fun8 =@(EA) sum(difference8(EA));
    obj_fun9 =@(EA) sum(difference9(EA));
    obj_fun10 =@(EA) sum(difference10(EA));
    obj_fun11 =@(EA) sum(difference11(EA));
    obj_fun12 =@(EA) sum(difference12(EA));
    obj_fun13 =@(EA) sum(difference13(EA));
    
 
%     options = optimset('MaxFunEvals',100000000,'MaxIter',1E5);
%     
%     EAA1 = fminsearch(obj_fun1,EA0,options);
%     EAA2 = fminsearch(obj_fun2,EA0,options);
%     EAA3 = fminsearch(obj_fun3,EA0,options);
%     EAA4 = fminsearch(obj_fun4,EA0,options);
%     EAA5 = fminsearch(obj_fun5,EA0,options);
%     EAA6 = fminsearch(obj_fun6,EA0,options);
%     EAA7 = fminsearch(obj_fun7,EA0,options);
%     EAA8 = fminsearch(obj_fun8,EA0,options);
%     EAA9 = fminsearch(obj_fun9,EA0,options);
%     EAA10 = fminsearch(obj_fun10,EA0,options);
%     EAA11 = fminsearch(obj_fun11,EA0,options);
%     EAA12 = fminsearch(obj_fun12,EA0,options);
%     EAA13 = fminsearch(obj_fun13,EA0,options);
%options = optimoptions('particleswarm','FunctionTolerance',1E-40,'MaxGenerations',10000,'PopulationSize',200);

options = optimoptions('particleswarm','SwarmSize',200,'FunctionTolerance',1E-20);


% rng default
    EAA1 = particleswarm(obj_fun1,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA2 = particleswarm(obj_fun2,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA3 = particleswarm(obj_fun3,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA4 = particleswarm(obj_fun4,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA5 = particleswarm(obj_fun5,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA6 = particleswarm(obj_fun6,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA7 = particleswarm(obj_fun7,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA8 = particleswarm(obj_fun8,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA9 = particleswarm(obj_fun9,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA10 = particleswarm(obj_fun10,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA11 = particleswarm(obj_fun11,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA12 = particleswarm(obj_fun12,2,[],[],options);%,[],[],[],[],[],[],[],options);
    EAA13 = particleswarm(obj_fun13,2,[],[],options);%,[],[],[],[],[],[],[],options);

    
    
    
%     
%     figure
%     hold on
%     plot(a,da)
%     plot(a,dT_1c(EAA1))
%     plot(a,dT_2c(EAA2))
%     plot(a,dT_3c(EAA3))
%     plot(a,dT_4c(EAA4))
%     plot(a,dT_5c(EAA5))
%     plot(a,dT_6c(EAA6))
%     hold off
%     legend('data','1','2','3','4','5','6')






    
    
    
    
    
    
   x =      [EAA1(2) EAA2(2) EAA3(2) EAA4(2) EAA5(2) EAA6(2) EAA7(2) EAA8(2) EAA9(2) EAA10(2) EAA11(2) EAA12(2) EAA13(2)];
   y = [EAA1(1) EAA2(1) EAA3(1) EAA4(1) EAA5(1) EAA6(1) EAA7(1) EAA8(1) EAA9(1) EAA10(1) EAA11(1) EAA12(1) EAA13(1)];
   
   
   
    p=polyfit(x,y,1);
    logA = polyval(p,Ea);
    
    xx = [x Ea];
    lx = [min(xx) max(xx)];
    ly = polyval(p,lx);
   figure
   hold on
   plot(x,y,'o')
   plot(Ea,logA,'*')
   plot(lx,ly,'-')
   hold off
xlabel('E_{\alpha}')
ylabel('ln(A_{\alpha})')
legend('Puntos de los modelos','Punto correcto',"l\'inea de m\'inimos cuadrados")
   A = exp(logA);

   
   yy = abs(y-logA); 
   
   [minloga,indx] = min(yy);
   closest = indx
   
   
   
   %    
%    figure
%    hold on
%    plot(T,da);
%    plot(T,dT_1c([A,Ea]));
%    plot(T,dT_2c(EAA6));
%    hold off









