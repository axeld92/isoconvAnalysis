data = St02; %selecciona datos de una corrida en particular
Temp = data(:,1)+273.15; %pasa temperatura a Kelvin
massloss = data(:,4); %Selecciona datos de pérdida de masa
alpha = (massloss(1)-massloss)/(massloss(1)-massloss(end)); % calcula conversión
ind = getind(alpha); %indices para eliminar conversiones muy bajas y muy altas
ind = logical(ind); 
alpha = alpha(ind);
Temp = Temp(ind);
% % 
% % dadT = smoothdadT2(Temp,alpha); %crear objeto cfit para la derivada de la conversion
% % dT = dadT(Temp); %vector de datos de derivada de conversión
% % 
% % %fsmix =@(parameters) fs_mixture(Temp,parameters(1:4),parameters(5:8),parameters(6:12));
% % 
% % %obj_fun =@(parameters) abs(sum((dT -  fsmix(parameters)).^2));
% % 
% % %Función Objetivo
% % obj_fun =@(parameters) sum( ( 1000*dT - 1000*fs_mixture(Temp,parameters(1:4),parameters(5:8),parameters(6:12) ) ).^2);
% % 
% % % height_1 = parameters_1(1);
% % % skew_1  = parameters_1(2);
% % % position_1 = parameters_1(3);
% % % width_1
% % 
% % 
% % %parameters0 = [1,-1,500,40,1,-1,600,31,1,-1,700,100];
% % parameters0 = [0.006,-0.01,530,70,0.012,-0.01,611,20,0.0015,-0.01,550,300];
% % 
% % 
% % 
% % % 
% % options = optimset('MaxFunEvals',1000000,'TolFun',1e-9,'TolX',1e-9,'PlotFcns',@optimplotfval);
% % fsfit = fminsearch(obj_fun,parameters0,options);
% % %out = fsmix(fsfit);
% % % out = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
% % % out_1 = fs_function(Temp,fsfit(1:4));
% % % out_2 = fs_function(Temp,fsfit(5:8));
% % % out_3 = fs_function(Temp,fsfit(9:12));
% % 
% % 
% % % figure
% % % 
% % % 
% % % plot(Temp,out,Temp,dT,Temp,out_1,Temp,out_2,Temp,out_3);
% % % grid
% % % legend('Calc','real','out_1','out_2','out_3')
% % 
% % %initial guess
% % %parameters0 = [0.3,-0.1,550,64,0.003,-0.3,620,24,0.01,-0.006,600,100];
% % 
% % lowerBound = [0,-10,0,0,0,-10,0,0,0,-10,0,0];
% % upperBound = [0.015,1000,1000,1000,0.015,1000,1000,1000,0.005,1000,1000,1000];
% % 
% % % rng default % For reproducibility
% % % gs = GlobalSearch;
% % % problem = createOptimProblem('fmincon','x0',parameters0,'objective',obj_fun,'lb',lowerBound,'ub',upperBound);
% % % fsfit = run(gs,problem)
% % 
% % out = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
% % out_1 = fs_function(Temp,fsfit(1:4));
% % out_2 = fs_function(Temp,fsfit(5:8));
% % out_3 = fs_function(Temp,fsfit(9:12));
% % 
% % figure
% % 
% % 
% % plot(Temp,out,Temp,dT,Temp,out_1,Temp,out_2,Temp,out_3);
% % grid
% % legend('Calc','real','out_1','out_2','out_3')
% % 

dadT = diff_alpha(Temp,alpha); %crear objeto cfit para la derivada de la conversion
dT = smoothdata(dadT,'sgolay',1000); %vector de datos de derivada de conversión

%fsmix =@(parameters) fs_mixture(Temp,parameters(1:4),parameters(5:8),parameters(6:12));

%obj_fun =@(parameters) abs(sum((dT -  fsmix(parameters)).^2));

%Función Objetivo
%obj_fun =@(parameters) sum( (dT - fs_mixture(Temp,parameters(1:4),parameters(5:8),parameters(6:12) ) ).^2);

% height_1 = parameters_1(1);
% skew_1  = parameters_1(2);
% position_1 = parameters_1(3);
% width_1


%parameters0 = [1,-1,500,40,1,-1,600,31,1,-1,700,100];
parameters0 = [0.004,-0.1,550.15,50,    0.006,-0.1,620,30,      0.0001,-0.1,700,250];


lb = [0,-0.5,540,0,         0,-0.5,590,0,             0,-0.5,550,0];
ub = [0.012,0.5,580,250,    0.012,0.5,640,250,       0.012,0.5,710,350];



% 
% options = optimset('MaxFunEvals',1000000,'TolFun',1e-9,'TolX',1e-9,'PlotFcns',@optimplotfval);
% fsfit = fminsearch(obj_fun,parameters0,options);
% % out = fsmix(fsfit);
% out = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
% out_1 = fs_function(Temp,fsfit(1:4));
% out_2 = fs_function(Temp,fsfit(5:8));
% out_3 = fs_function(Temp,fsfit(9:12));


% figure
% 
% 
% plot(Temp,out,Temp,dT,Temp,out_1,Temp,out_2,Temp,out_3);
% grid
% legend('Calc','real','out_1','out_2','out_3')

%initial guess
%parameters0 = [0.3,-0.1,550,64,0.003,-0.3,620,24,0.01,-0.006,600,100];
% 

% 
% %rng default % For reproducibility
% gs = GlobalSearch;
% problem = createOptimProblem('fmincon','x0',parameters0,'objective',obj_fun,'lb',lowerBound,'ub',upperBound);
% fsfit = run(gs,problem);

options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',10000000000,...
    'FunctionTolerance',1E-10,'MaxIterations',10000,'StepTolerance',1E-15);
[fsfit,error] = lsqcurvefit(@(parameters0,Temp)fs_mixture2(parameters0,Temp),parameters0,Temp,dT,lb,ub,options)


out = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
out_1 = fs_function(Temp,fsfit(1:4));
out_2 = fs_function(Temp,fsfit(5:8));
out_3 = fs_function(Temp,fsfit(9:12));


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


