function p = comeffect2(data_1)

%Get temperatures
Temp_1 = data_1(:,1)+273.15;
% Temp_2 = data_2(:,1)+273.15;
% Temp_3 = data_3(:,1)+273.15;
% Temp_4 = data_4(:,1)+273.15;

%Get times
time_1 = data_1(:,2)-data_1(1,2);
% time_2 = data_2(:,2)-data_2(1,2);
% time_3 = data_3(:,2)-data_3(1,2);
% time_4 = data_4(:,2)-data_4(1,2);

%get massloss

massloss_1 = data_1(:,4);
% massloss_2 = data_2(:,4);
% massloss_3 = data_3(:,4);
% massloss_4 = data_4(:,4);

%calculate conversion
alpha_1 = (massloss_1(1)-massloss_1)/(massloss_1(1)-massloss_1(end));
% alpha_2 = (massloss_2(1)-massloss_2)/(massloss_2(1)-massloss_2(end));
% alpha_3 = (massloss_3(1)-massloss_3)/(massloss_3(1)-massloss_3(end));
% alpha_4 = (massloss_4(1)-massloss_4)/(massloss_4(1)-massloss_4(end));

%get rid of data outside of alpha = 0.05 to alpha = 0.095

ind = getind(alpha_1);
ind = logical(ind);
alpha_1 = alpha_1(ind);
Temp_1 = Temp_1(ind);
time_1 = time_1(ind);

% ind = getind(alpha_2);
% ind = logical(ind);
% alpha_2 = alpha_2(ind);
% Temp_2 = Temp_2(ind);
% time_2 = time_2(ind);
% 
% ind = getind(alpha_3);
% ind = logical(ind);
% alpha_3 = alpha_3(ind);
% Temp_3 = Temp_3(ind);
% time_3 = time_3(ind);
% 
% ind = getind(alpha_4);
% ind = logical(ind);
% alpha_4 = alpha_4(ind);
% Temp_4 = Temp_4(ind);
% time_4 = time_4(ind);
% 
 beta = (Temp_1(end)-Temp_1(1))/(time_1(end)-time_1(1)); 
%         (Temp_2(end)-Temp_2(1))/(time_2(end)-time_2(1)) , 
%         (Temp_3(end)-Temp_3(1))/(time_3(end)-time_3(1)) , 
%         (Temp_4(end)-Temp_4(1))/(time_4(end)-time_4(1))];

%create cfit objects for temperature and time as a function of conversion

T_1 = Tfit(Temp_1,alpha_1);
% T_2 = Tfit(alpha_2,Temp_2);
% T_3 = Tfit(alpha_3,Temp_3);
% T_4 = Tfit(alpha_4,Temp_4);

%create cfit objects for dalpha as a function of temperature

dadT1 = smoothdadT2(Temp_1,alpha_1);
% dadT2 = smoothdadT2(Temp_2,alpha_2);
% dadT3 = smoothdadT2(Temp_3,alpha_3);
% dadT4 = smoothdadT2(Temp_4,alpha_4);
 a = alpha_1;

fa_1 = (1-a);
   %fa_2 =2*(1-a).*(-log(1-a)).^(1/2);
   fa_2 = (1-a).^2; 
   fa_3 = 3*(1-a).*(-log(1-a)).^(2/3);
    fa_4 = 4*(1-a).*(-log(1-a)).^(3/4);
    %fa_2(fa_2~=real(fa_2)) = 0;
      fa_3(fa_3~=real(fa_3)) = 0;
        fa_4(fa_4~=real(fa_4)) = 0;
    
%     fa_2 =  (1-a).^2;%2*(1-a5_0).*(-log(1-a5_0)).^(1/2);
%     fa_3 =  (1-a).^3;%3*(1-a5_0).*(-log(1-a5_0)).^(2/3);
%     fa_4 =  (1-a).^4;%4*(1-a5_0).*(-log(1-a5_0)).^(3/4);
%     
    R = 8.314/1000;
    EA0 = [10^8,250];
    b = 5;
    dT = dadT1(Temp_1);
    T = Temp_1;
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









