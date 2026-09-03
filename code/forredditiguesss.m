%Parameters to create data
pa = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];

%independent variable vector
Temp = (400:1000)';



% %create data
fs_data = fsm(Temp,pa);

plot(Temp,fs_data)

%%
% 
% data = cleandata(St01);
% 
% fs_data = data(:,3);
% Temp = data(:,1);

%set lower and upper bounds
lb = [0,    -.5,   0,    0,...
      0,    -.5,   0,    0,...
      0,    -.5,   0,    0];
  
ub = [0.1,  0.5,    800,    200,...
      0.1,  0.5,    800,    200,...
      0.1,  0.5,    800,    400];
  
%Set initial guess
p0 = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];


%Set objective function
obj_fun =@(parameters) sum((fsm(Temp,parameters)-fs_data).^2);

%Set optmization problem

% PATTERN SEARCH
rng default
[fsfit,fval] = patternsearch(obj_fun,p0,[],[],[],[],lb,ub)

% % Set optmization problem
% rng default % For reproducibility
% opts = optimoptions(@fmincon,'Algorithm','sqp');
% gs = GlobalSearch;
% problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);

% % fit model to data
% [fsfit,fval] = run(gs,problem)

difference = pa - fsfit


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