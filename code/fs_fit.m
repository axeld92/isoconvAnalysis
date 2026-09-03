function [Temp_out,fs]  = fs_fit(data)

Temp = data(:,1)+273.15; %pasa temperatura a Kelvin
massloss = data(:,4); %Selecciona datos de pérdida de masa
alpha = (massloss(1)-massloss)/(massloss(1)-massloss(end)); % calcula conversión
ind = getind(alpha); %indices para eliminar conversiones muy bajas y muy altas
ind = logical(ind); 
alpha = alpha(ind);
Temp = Temp(ind);

dadT = diff_alpha(Temp,alpha); %crear objeto cfit para la derivada de la conversion
figure
plot(Temp,dadT)
dT = smoothdata(dadT,'sgolay',500); %vector de datos de derivada de conversión

parameters0 = [0.003,-0.15,250+273.15,50,0.006,-0.15,320+273.15,30,0.001,-0.15,670,200];

lb = [0,-10,0,0,0,-10,0,0,0,-10,0,0];
ub = [0.015,1000,1000,1000,0.015,1000,1000,1000,0.005,1000,1000,1000];

options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',10000000000,...
    'FunctionTolerance',1E-10,'MaxIterations',10000,'StepTolerance',1E-15);
[fsfit,error] = lsqcurvefit(@(parameters0,Temp)fs_mixture2(parameters0,Temp),parameters0,Temp,dT,lb,ub,options);


out = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
out_1 = fs_function(Temp,fsfit(1:4));
out_2 = fs_function(Temp,fsfit(5:8));
out_3 = fs_function(Temp,fsfit(9:12));

fs = [out out_1 out_2 out_3];
Temp_out = Temp;