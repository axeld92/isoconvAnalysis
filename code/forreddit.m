%Parameters to create data
pa = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];

% Create vector for Temperatures
Temp = (400:1000)';

%create data
fs_data = fsm(Temp,pa);

%Set lower and upper bounds
lb = [0,    -.33,   273,    50,...
      0,    -.33,   563,    0,...
      0,    -.29,   603,    100];
  
ub = [0.1,  0.1,    590,    100,...
      0.1,  0.25,   653,    70,...
      0.1,  2,      720,    400];
  
%Set initial guess
%p0 = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];
p0 = (lb + ub)/2

%Set objective function
obj_fun =@(parameters) (sum((fsm(Temp,parameters)-fs_data).^2))/(length(Temp)-12);

%Set optmization problem
rng default % For reproducibility
opts = optimoptions(@fmincon,'OptimalityTolerance',1E-10,...
     'StepTolerance',1E-10,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);

%fit model to data
[fsfit,fval] = run(gs,problem)

%Plot results
figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fsm(Temp,fsfit))
hold off
legend('real','calc')

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