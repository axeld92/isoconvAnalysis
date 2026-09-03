function p = comeffect(T,a,b)

%T_1 = Tfit(T5_0,a5_0);
dT = smoothdadT2(T,a);
%dT(1:87) = 0;

    fa_1 = (1-a);
   fa_2 =2*(1-a).*(-log(1-a)).^(1/2);
    fa_3 = 3*(1-a).*(-log(1-a)).^(2/3);
    fa_4 = 4*(1-a).*(-log(1-a)).^(3/4);
    fa_2(fa_2~=real(fa_2)) = 0;
      fa_3(fa_3~=real(fa_3)) = 0;
        fa_4(fa_4~=real(fa_4)) = 0;
    
%     fa_2 =  (1-a).^2;%2*(1-a5_0).*(-log(1-a5_0)).^(1/2);
%     fa_3 =  (1-a).^3;%3*(1-a5_0).*(-log(1-a5_0)).^(2/3);
%     fa_4 =  (1-a).^4;%4*(1-a5_0).*(-log(1-a5_0)).^(3/4);
%     
    R = 8.314/1000;
    EA0 = [3,20];
    
    options = optimset('MaxFunEvals',10000);
    
    dT_1c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_1;
    dT_2c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_2;
    dT_3c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_3;
    dT_4c =@(EA) (1/b)*EA(1)*exp(-EA(2)./(R*T)).*fa_4;
    
    difference1 =@(EA) (dT-dT_1c(EA)).^2;
    difference2 =@(EA) (dT-dT_2c(EA)).^2;
    difference3 =@(EA) (dT-dT_3c(EA)).^2;
    difference4 =@(EA) (dT-dT_4c(EA)).^2;
    
    obj_fun1 =@(EA) sum(difference1(EA));
    obj_fun2 =@(EA) sum(difference2(EA));
    obj_fun3 =@(EA) sum(difference3(EA));
    obj_fun4 =@(EA) sum(difference4(EA));
    
    EAA1 = fminsearch(obj_fun1,EA0);
    EAA2 = fminsearch(obj_fun2,EA0);
    EAA3 = fminsearch(obj_fun3,EA0);
    EAA4 = fminsearch(obj_fun4,EA0);
    
   x = [EAA1(2) EAA2(2) EAA3(2) EAA4(2)];
   y = log([EAA1(1) EAA2(1) EAA3(1) EAA4(1)]);
%    figure
%    plot(x,y)
   
   p=polyfit(x,y,1);
end
  