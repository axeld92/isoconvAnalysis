function [output_1,output_2,output_3] = fs_out2(data)
%data = St01; %selecciona datos de una corrida en particular
Temp = data(:,1); %pasa temperatura a Kelvin
massloss = data(:,4); %Selecciona datos de pérdida de masa
time = data(:,2);
alpha_raw = (massloss(1)-massloss)/(massloss(1)-massloss(end)); % calcula conversión
alpha = smoothdata(alpha_raw,'sgolay',2);
dT = data(:,3);

% height_1 = parameters_1(1);
% skew_1  = parameters_1(2);
% position_1 = parameters_1(3);
% width_1

% 
% parameters0 =   [0.03,0.05,570,50,...
%                 0.06,-0.05,605,30,...
%                 0,-0.05,659,200];

%parameters0 = [0.0363 0.0682 560.0000 55.1830 0.0597 -0.2543 632.3738 35.0086 0.0118 0.9939 609.8258 160.0000];
parameters0 = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];

%parameters0 = [0.0543    0.4773  559.5178   56.8689    0.0643    0.0008  633.3581   43.2677    0.0406    0.5249  686.6261  173.0819];


% lb = [0,-0.5,530,0,         0,-0.5,590,0,             0,-0.5,550,0];
% ub = [0.012,0.6,580,250,    0.012,0.9,640,250,       0.0025,0.6,720,360];

lb = [0,    -.33,   273,    50,...
      0,    -.33,   563,    0,...
      0,    -.29,   603,    100];
  
ub = [0.1,  0.1,    590,    100,...
      0.1,  0.25,   653,    70,...
      0.1,  2,      720,    400];

% lb = [];
% ub = [];
  
 
obj_fun =@(parameters) (sum((fs_mixture(Temp,parameters(1:4),parameters(5:8),parameters(6:12))-dT).^2))/(length(Temp)-12);

% obj_fun =@(parameters)    sum(([0.1*(fs_mixture(Temp(1:102),parameters(1:4),parameters(5:8),parameters(6:12))-dT(1:102));
%                            0.8*(fs_mixture(Temp(102:300),parameters(1:4),parameters(5:8),parameters(6:12))-dT(102:300));
%                            0.7*(fs_mixture(Temp(300:end),parameters(1:4),parameters(5:8),parameters(6:12))-dT(300:end))]).^2);

% rng default % For reproducibility
% opts = optimoptions(@fmincon,'Algorithm','sqp')
%   gs = GlobalSearch('StartPointsToRun','all','FunctionTolerance',1e-10 ,'BasinRadiusFactor', 0.2,...
%       'XTolerance',1e-8,'NumStageOnePoints',1000)%,'PlotFcn',@gsplotbestf,,'NumTrialPoints',10000)
%   problem = createOptimProblem('fmincon','x0',parameters0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
%   fsfit = run(gs,problem)
%   
% problem = createOptimProblem('fmincon','x0',parameters0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
%   ms = MultiStart;
%  [fsfit,fval] = run(ms,problem,200)
  
  
%  opts = optimoptions('ga','MaxGenerations',1000,'FunctionTolerance',1E-10,'PopulationSize',1000,'UseParallel',true);
% [fsfit,fval] = ga(obj_fun,12,[],[],[],[],lb,ub,[],opts)

% opts = optimoptions('particleswarm','SwarmSize',1000);
% [fsfit , fval] = particleswarm(obj_fun,12,lb,ub,opts)


% opts = optimoptions('patternsearch','MaxIterations',1E10,'MeshTolerance',1E-15,'StepTolerance',1E-15,...
%     'MaxFunctionEvaluations',1E20,'MeshExpansionFactor',2,'PollMethod','GPSPositiveBasisNp1');
% [fsfit,fval] = patternsearch(obj_fun,parameters0,[],[],[],[],lb,ub,opts)

% opts = optimoptions('surrogateopt','MinSurrogatePoints',10000)
% [fsfit,fval] = surrogateopt(obj_fun,lb,ub)


 options = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',10000000000,...
    'FunctionTolerance',1E-15,'MaxIterations',10000,'StepTolerance',1E-20,'OptimalityTolerance',1E-15);
[fsfit,fval] = lsqcurvefit(@(parameters0,Temp)fs_mixture2(parameters0,Temp),parameters0,Temp,dT,lb,ub,options)



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
