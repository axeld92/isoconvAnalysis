%function deconvolve(data)

 data = cleandata(St05);

fs_data = data(:,3);
Temp = data(:,1);
time = data(:,2);

Tempr = Temp;
Temp = (Tempr-Tempr(1))/(Tempr(end)-Tempr(1));

%find peaks
% 
% [peaks , locp, w, prom ]= findpeaks(fs_data,'MinPeakProminence',0.0005)
%  Tp = Temp(locp);
fs_datar = fs_data;
hmax = max(fs_data);
fs_data = fs_datar/hmax;

%set lower and upper bounds
% lb = [0.01,    -1,   0.01,    0.01,...
%       0.01,    -1,   0.01,    0.01,...
%       0.01,    -1,   0.01,    0.01];
%   
% ub = [0.99,  -0.0001,   0.99,    0.99,...
%       0.99,  -0.0001,   0.99,    0.99,...
%       0.99,  -0.0001,   0.99,    0.99];
lb = [0.01,    -10,   0.01,    0.01,...
      0.01,    -10,   0.01,    0.01,...
      0.01,    -10,   0.01,    0.01];
  
ub = [0.99,  10,   0.99,    0.99,...
      0.99,  10,   0.99,    0.99,...
      0.99,  10,   0.99,    0.99];


%Set initial guess

 p0 = (lb + ub)/2;

%Set objective function
obj_fun =@(parameters) (sum((fsm(Temp,parameters)-fs_data).^2))/(length(Temp)-12);

% % PARTICLE SWARM
% rng default;
% options = optimoptions('fmincon','OptimalityTolerance',1E-15,...
%     'StepTolerance',1E-30,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
% opts = optimoptions('particleswarm','SwarmSize',1000,'HybridFcn',{@fmincon,options});
% [fsfit , fval] = particleswarm(obj_fun,12,lb,ub,opts)

rng default % For reproducibility
 opts = optimoptions(@fmincon,'OptimalityTolerance',1E-5,'StepTolerance',1E-5,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
[fsfit,fval] = run(gs,problem);
%%
funf =@(p)  fsm(Temp,p) - fs_data;

[x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(funf,fsfit)


%%

Covaria = resnorm*pinv(full(jacobian'*jacobian));

sj = sqrt(diag(Covaria));

ttinv = tinv(1-0.025,length(Temp)-12);

ICs = sj*ttinv;

% Plot errorbars
errorbar(1:12,x,ICs,'o','LineWidth',1.5)

%%
Corrij = 


%%

figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fsm(Temp,fsfit))
plot(Temp,fs(Temp,fsfit(1:4)))
plot(Temp,fs(Temp,fsfit(5:8)))
plot(Temp,fs(Temp,fsfit(9:12)))
hold off
ylim([-0.001 0.012])
legend('Real','Calculado','1','2','3')
xlabel('Temperatura [K]')
ylabel('d\alpha/dT')

% Calculate outputs




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


    function y = fsmi(T,pm)
        y = cumtrapz(T,fsm(T,pm));
    end

function y = fsm4(T,pm)

y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12)) + fs(T,pm(13:16));

end

