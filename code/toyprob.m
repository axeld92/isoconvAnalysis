% %Parameters to create data
% %pa = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];
% %pa = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];
% pa = [0.006487 0.0682 547.8 55.1830 0.01162 -0.2543 611.6 35.0086 0.0012 0.9939 649 160];
% Temp = (400:1000)';
% 
% %create data
% fs_data = fsm(Temp,pa);
% fs_int = cumtrapz(Temp,fs_data);
% 
% plot(Temp,fs_int)
% 
% 





%%
data = St01;
data = cleandata(data);

fs_data = data(:,3);
Temp = data(:,1);
time = data(:,2);
fs_int = data(:,4);


Tempr = Temp;
Temp = (Tempr-Tempr(1))/(Tempr(end)-Tempr(1));
fs_datar = fs_data;
hmax = max(fs_data);
fs_data = fs_datar/hmax;



lb = [0.01,    -0.5,   0.01,    0.01,...
      0.01,    -0.5,   0.01,    0.01,...
      0.01,    -0.5,   0.01,    0.01];
  
ub = [0.99,  0.5,   0.99,    0.99,...
      0.99,  0.5,   0.99,    0.99,...
      0.99,  0.5,   0.99,    0.99];

  
%Set initial guess

p0 = (lb + ub)/2;


%Set objective function
obj_fun =@(parameters) (sum((fsm(Temp,parameters)-fs_data).^2))/(length(Temp)-12);


%Set optmization problem
 rng default % For reproducibility
 opts = optimoptions(@fmincon,'Algorithm','sqp');
 gs = GlobalSearch;
 problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
% 
%fit model to data
[fsfit,fval] = run(gs,problem)

Temp = Tempr;
fs_data = fs_datar;
ind = [3 7 11];
fsfit(ind) = fsfit(ind)*(Temp(end)-Temp(1)) + Temp(1); 
ind = [4 8 12];
fsfit(ind)  = fsfit(ind)*(Temp(end)-Temp(1));
ind = [1 5 9];
fsfit(ind) = fsfit(ind)*hmax;

% problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
%   ms = MultiStart;
%  [fsfit,fval] = run(ms,problem,200)

% % % PATTERN SEARCH
% rng default
% opts = optimoptions('patternsearch','MaxIterations',1E10,'MeshTolerance',1E-15,'StepTolerance',1E-15,...
%    'MaxFunctionEvaluations',1E20,'MeshExpansionFactor',2);
% 
% [fsfit,fval] = patternsearch(obj_fun,p0,[],[],[],[],lb,ub,opts)


plot(Temp,fs_int,Temp,cumtrapz(Temp,fsm(Temp,fsfit)))







%%

% options = optimoptions('fmincon','Display','iter','Algorithm','sqp','OptimalityTolerance',1E-15,...
%     'StepTolerance',1E-20,'MaxFunctionEvaluations',10000);
% [fsfit , fval] = fmincon(obj_fun,fsfit,[],[],[],[],lb,ub,[],options);


figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fsm(Temp,fsfit))
plot(Temp,fs(Temp,fsfit(1:4)))
plot(Temp,fs(Temp,fsfit(5:8)))
plot(Temp,fs(Temp,fsfit(9:12)))
plot(Temp,fs(Temp,pa(1:4)),'--')
plot(Temp,fs(Temp,pa(5:8)),'--')
plot(Temp,fs(Temp,pa(9:12)),'--')
hold off
legend('real','calc','1','2','3','1c','2c','3c')
%%
figure
hold on
plot(Temp,fs_int)
plot(Temp,intfsm(Temp,fsfit),':')
hold off

% 
% figure
% plot(Temp,fsm(Temp,fsfit)-fs_data)





%Define Fraser Suzuki Function
function y = fs(T,p)

 in = 2 * p(2) * ((T-p(3)) / p(4));
 out = - log(2) / p(2)^2;
 
 y = p(1) * exp(out .* (log(1+in)).^2);
 y(y~=real(y))=0;
end

%Define Fraser suzuki mixture

function y = fsm(T,pm)

y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12));

end

function y = intfsm(T,pm)

y = cumtrapz(T,fsm(T,pm));

end


% % Set optmization problem using globalsearch
% rng default % For reproducibility
%  opts = optimoptions(@fmincon,'OptimalityTolerance',1E-8,...
%       'StepTolerance',1E-10,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
% gs = GlobalSearch;
% problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
% 
% % % fit model to data
% tic;
% [fsfit,fval] = run(gs,problem)
% disp('GS')
% toc
% 
% % % MultiStart
% rng default
% problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
% ms = MultiStart;
% tic
% [fsfit,fval] = run(ms,problem,200)
% disp('MS')
% toc

% % CURVE FIT (LOCAL)
%  options = optimoptions(@lsqcurvefit,'Algorithm','trust-region-reflective','MaxFunctionEvaluations',10000000000,...
%     'FunctionTolerance',1E-15,'MaxIterations',10000,'StepTolerance',1E-20,'OptimalityTolerance',1E-15);
% [fsfit,fval] = lsqcurvefit(@(p0,Temp)fs_mixture2(p0,Temp),p0,Temp,fs_data,lb,ub,options)

% % PATTERN SEARCH
% rng default
% opts = optimoptions('patternsearch','MaxIterations',1E10,'MeshTolerance',1E-15,'StepTolerance',1E-15,...
%     'MaxFunctionEvaluations',1E20,'MeshExpansionFactor',2);
% tic;
% [fsfit,fval] = patternsearch(obj_fun,p0,[],[],[],[],lb,ub,opts)
% disp('PS')
% toc
% % SURROGATE OPTIMIZATION
% rng default
% opts = optimoptions('surrogateopt','MinSurrogatePoints',10000)
% tic;
% [fsfit,fval] = surrogateopt(obj_fun,lb,ub)
% disp('SO')
% toc
% 
% % GENETIC ALGORITHM
% rng default  
% opts = optimoptions('ga','MaxGenerations',1000,'FunctionTolerance',1E-10,'PopulationSize',1000);
% tic;
%   [fsfit,fval] = ga(obj_fun,12,[],[],[],[],lb,ub,[],opts)
% disp('GA')
% toc
% 
% % % SIMULATED ANNEALING
% rng default
% options = optimoptions('simulannealbnd','FunctionTolerance',1E-10);
% tic;
% [fsfit,fval] = simulannealbnd(obj_fun,p0,lb,ub,options)
% disp('SA')
% toc