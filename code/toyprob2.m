%Parameters to create data
pa = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];
%pa = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];

Temp = (400:1000)';

%create data
fs_data = fsm(Temp,pa);

%find peaks
% [peaks , locp, w, prom ]= findpeaks(fs_data);
%  Tp = Temp(locp);
hmax = max(fs_data);


%set lower and upper bounds
lb = [0,    -.5,   500,    50,...
      0,    -.5,   500,    0,...
      0,    -.5,   500,    100];
  
ub = [hmax,  -0.00001,   800,    100,...
      hmax,  -0.00001,   800,    70,...
      hmax,  -0.00001,   800,    400];



%set lower and upper bounds
% lb = [prom(1),    -.33,   Tp(1)-5,    50,...
%       prom(1),    -.33,   Tp(2)-5,    0,...
%       0,    -.29,   603,    100];
%   
% ub = [peaks(1),  0.1,    Tp(1)+5,    100,...
%       peaks(2),  0.25,   Tp(2)+5,    70,...
%       0.1,  2,      720,    400];

lb = [0 0 0 0 0 0];
ub = [100 100 100 100 100 100];

%Set initial guess

p0 = (lb + ub)/2;


%p0 = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];
%p0 = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];




%Set objective function
obj_fun =@(parameters) (sum((fsl(Temp,parameters)-fs_data).^2));

% % PARTICLE SWARM
rng default;
options = optimoptions('fmincon','OptimalityTolerance',1E-15,...
    'StepTolerance',1E-30,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
opts = optimoptions('particleswarm','SwarmSize',1000,'HybridFcn',{@fmincon,options});

%tic;
[fsfit , fval] = particleswarm(obj_fun,6,lb,ub,opts)







% p0
% pa
% fsfit
% pdifference = abs((fsfit-pa)./pa)
% fval

figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fsl(Temp,fsfit))
plot(Temp,lognpdf(Temp,fsfit(1),fsfit(2)))
plot(Temp,lognpdf(Temp,fsfit(3),fsfit(4)))
plot(Temp,lognpdf(Temp,fsfit(5),fsfit(6)))
hold off
legend('real','calc','1','2','3')%,'1c','2c','3c')
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

function y = fsl(T,p)
y = lognpdf(T,p(1),p(2)) + lognpdf(T,p(3),p(4)) + lognpdf(T,p(5),p(6));
end

%Define Fraser suzuki mixture

function y = fsm(T,pm)

y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12));

end
