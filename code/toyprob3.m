%Parameters to create data using multiple
pa = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];
%pa = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];

Temp = (475:1000)';


%create data
fs_data = fsm(Temp,pa);

plot(Temp,fs_data)

Tempr = Temp;
Temp = (Tempr-Tempr(1))/(Tempr(end)-Tempr(1));

%%


lb = [0.1     0   0.2   0     0.2   0     0 0 0 ];
ub = [0.2   0.5   0.4   0.5   0.5   0.5   1 1 1];

%Set initial guess
% p0 = [560   55      632     35 660 160 0.5 0.5 0.5];
p0 = (lb+ub)/2;






%%
%Set objective function
obj_fun =@(parameters) (sum((fmgauss(Temp,parameters)-fs_data).^2));

 % Pattern search
rng default
opts = optimoptions('patternsearch','MaxIterations',1E10,'MeshTolerance',1E-15,'StepTolerance',1E-15,...
    'MaxFunctionEvaluations',1E20,'MeshExpansionFactor',2);
tic;
[fsfit,fval] = patternsearch(obj_fun,p0,[],[],[],[],lb,ub,opts)

% optss = optimset('TolFun',1E-20);%'PlotFcns',@optimplotfval);
% [fsfit,fval] = fminsearch(obj_fun,p0,optss)

figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fmgauss(Temp,fsfit))
plot(Temp,fsfit(7)*normpdf(Temp,fsfit(1),fsfit(2)))
plot(Temp,fsfit(8)*normpdf(Temp,fsfit(3),fsfit(4)))
plot(Temp,fsfit(9)*normpdf(Temp,fsfit(5),fsfit(6)))
hold off
legend('real','calc','1','2','3')%,'1c','2c','3c')



function y = fs(T,p)
%Fraser Suzuki function
 in = 2 * p(2) * ((T-p(3)) / p(4));
 out = - log(2) / p(2)^2;
 
 y = p(1) * exp(out .* (log(1+in)).^2);
 y(y~=real(y))=0;
end

function y = fmgauss(T,p)
%Multiple gaussian function

y_1 = normpdf(T,p(1),p(2));%y_1 = y_1/max(y_1);
y_2 = normpdf(T,p(3),p(4));%y_2 = y_2/max(y_2);
y_3 = normpdf(T,p(5),p(6));%y_3 = y_3/max(y_3);

y = p(7)*y_1 + p(8)*y_2 + p(9)*y_3;
end

%Define Fraser suzuki mixture

function y = fsm(T,pm)
%multiple fraser suzuki function
y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12));

end
