function [output_1,output_2,output_3] = fs_out(data)
%data = St01; %selecciona datos de una corrida en particular
Temp = data(:,1)+273.15; %pasa temperatura a Kelvin
massloss = data(:,4); %Selecciona datos de pérdida de masa
time = data(:,2);
alpha_raw = (massloss(1)-massloss)/(massloss(1)-massloss(end)); % calcula conversión
alpha = smoothdata(alpha_raw,'sgolay',2);
dT = -data(:,3);
%alpha = alpha_raw;

% [time,Temp,alpha] = removeerror(time,Temp,alpha);
 
 % ind = getind(alpha); %indices para eliminar conversiones muy bajas y muy altas
% ind = logical(ind); 
% alpha = alpha(ind);
% Temp = Temp(ind);
% time = time(ind);

% dadT = [0 ; diff(alpha)./diff(Temp)]; %calcular derivada de la conversion
% dT = smoothdata(dadT,'sgolay',100); %vector de datos de derivada de conversión

% RdT = dadT;

% height_1 = parameters_1(1);
% skew_1  = parameters_1(2);
% position_1 = parameters_1(3);
% width_1


parameters0 =   [0.003,-0.15,250+273.15,50,...
                0.006,-0.15,320+273.15,30,...
                0.001,-0.15,390+273.15,200];

% % 
% lb = [0,-0.5,530,0,         0,-0.5,590,0,             0,-0.5,550,0];
% ub = [0.012,0.6,580,250,    0.012,0.9,640,250,       0.0025,0.6,720,360];

lb = [0,-.33,273,50,...
      0,-.33,563,0,...
      0,-.29,603,160];
  
ub = [2,0.6,560,100,...
      2,0.25,653,60,...
      2,2,720,1000];

% lb = [];
% ub = [];
  
 
obj_fun =@(parameters) sum( (dT - fs_mixture(Temp,parameters(1:4),parameters(5:8),parameters(6:12) ) ).^2);
rng default % For reproducibility
%opts = 
% gs = GlobalSearch('FunctionTolerance',1e-10, 'BasinRadiusFactor', 0.2,'XTolerance',1e-8,'Display','iter')
% problem = createOptimProblem('fmincon','x0',parameters0,'objective',obj_fun,'lb',lb,'ub',ub);
% fsfit = run(gs,problem)
 options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',10000000000,...
    'FunctionTolerance',1E-15,'MaxIterations',10000,'StepTolerance',1E-20);
% [fsfit,error] = lsqcurvefit(@(parameters0,Temp)fs_mixture2(parameters0,Temp),parameters0,Temp,dT,lb,ub,options)
[fsfit,fval] = lsqcurvefit(@(parameters0,Temp)fs_mixture2(parameters0,Temp),parameters0,Temp,dT,lb,ub,options)

%fval



out = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
out_1 = fs_function(Temp,fsfit(1:4));
out_2 = fs_function(Temp,fsfit(5:8));
out_3 = fs_function(Temp,fsfit(9:12));
% 
% hemicellulose = fsfit(1:4);
% cellulose = fsfit(5:8);
% ligning = fsfit(9:12)



output_1 = [Temp, time, out_1, zeros(size(time))];
output_2 = [Temp, time, out_2, zeros(size(time))];
output_3 = [Temp, time, out_3, zeros(size(time))];

for i=2:length(Temp)
    output_1(i,4) = trapz(Temp(1:i),out_1(1:i));
end
   
for i=2:length(Temp)
    output_2(i,4) = trapz(Temp(1:i),out_2(1:i));
end

for i=2:length(Temp)
    output_3(i,4) = trapz(Temp(1:i),out_3(1:i));
end

precent_hemic = output_1(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4))
precent_cellulose = output_2(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4))
precent_lignin = output_3(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4))

output_1(:,4) = (output_1(:,4) - output_1(1,4))/(output_1(end,4)-output_1(1,4));
output_2(:,4) = (output_2(:,4) - output_2(1,4))/(output_2(end,4)-output_2(1,4));
output_3(:,4) = (output_3(:,4) - output_3(1,4))/(output_3(end,4)-output_3(1,4));



% 
figure
hold on
plot(Temp,out)
plot(Temp,dT,':')
plot(Temp,out_1,'--')
plot(Temp,out_2,'--')
plot(Temp,out_3,'--');
hold off
grid
legend('Calc','real','out_1','out_2','out_3')

% 
% figure
% hold on
% %plot(Temp,out)
% %plot(Temp,dT,':')
% plot(Temp,output_1(:,4),'--')
% plot(Temp,output_2(:,4),'--')
% plot(Temp,output_3(:,4),'--');
% hold off
% grid
% legend('hemi','cell','lign')


%output = [out_1 , out_2 , out_3];
