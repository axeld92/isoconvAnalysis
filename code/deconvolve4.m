function [output_1,output_2,output_3,output_4,c] = deconvolve4(data)

data = cleandata(data);

fs_data = data(:,3);
Temp = data(:,1);
time = data(:,2);

Tempr = Temp;
Temp = (Tempr-Tempr(1))/(Tempr(end)-Tempr(1));

%find peaks
% 
% [peaks , locp, w, prom ]= findpeaks(fs_data,'MinPeakProminence',0.0005)
%  Tp = Temp(locp);


hmax = max(fs_data);

lb = [0,    -0.5,   0.01,    0.01,...
      0,    -0.5,   0.01,    0.01,...
      0,    -0.5,   0.37,    0.01,...
      0,    -0.5,   0.43,    0.5];
  
ub = [hmax,  -.001,   0.99,    0.99,...
      hmax,  -.001,   0.99,    0.99,...
      hmax,  -.001,   0.99,    0.99,...
      hmax,  -.001,   0.99,    0.99];

%Set initial guess

 p0 = (lb + ub)/2;
 p0(13:16) = [0.0005    0.2494    0.4375    0.5756];

%Set objective function
obj_fun =@(parameters) (sum((fsm4(Temp,parameters)-fs_data).^2));

% % % PARTICLE SWARM
% rng default;
% options = optimoptions('fmincon','OptimalityTolerance',1E-15,...
%     'StepTolerance',1E-30,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
% opts = optimoptions('particleswarm','SwarmSize',1000,'UseParallel',true,'HybridFcn',{@fmincon,options});
% [fsfit , fval] = particleswarm(obj_fun,16,lb,ub,opts);

rng default % For reproducibility
 opts = optimoptions(@fmincon,'OptimalityTolerance',1E-8,'StepTolerance',1E-10,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
[fsfit,fval] = run(gs,problem);



A = [   fsfit(1:4);
        fsfit(5:8);
        fsfit(9:12);
        fsfit(13:16)];
    
AA = sortrows(A,3);

%fsfit = [AA(1,:) AA(2,:) AA(4,:) AA(3,:)];
fsfit = [AA(1,:) AA(2,:) AA(3,:) AA(4,:)];
Temp = Tempr;

ind = [3 7 11 15];
fsfit(ind) = fsfit(ind)*(Temp(end)-Temp(1)) + Temp(1); 
ind = [4 8 12 16];
fsfit(ind)  = fsfit(ind)*(Temp(end)-Temp(1));


fsfit
fval

% %  Plot results
figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fsm4(Temp,fsfit))
plot(Temp,fs(Temp,fsfit(1:4)))
plot(Temp,fs(Temp,fsfit(5:8)))
plot(Temp,fs(Temp,fsfit(9:12)))
plot(Temp,fs(Temp,fsfit(13:16)))
hold off
ylim([-0.001 0.012])
legend('Real','Calculado','1','2','3')
xlabel('Temperatura [K]')
ylabel('d\alpha/dT')


out_1 = fs_function(Temp,fsfit(1:4));
out_2 = fs_function(Temp,fsfit(5:8));
out_3 = fs_function(Temp,fsfit(9:12));
out_4 = fs_function(Temp,fsfit(13:16));

output_1 = [Temp, time, out_1, zeros(size(time))];
output_2 = [Temp, time, out_2, zeros(size(time))];
output_3 = [Temp, time, out_3, zeros(size(time))];
output_4 = [Temp, time, out_4, zeros(size(time))];

output_1(:,4) = cumtrapz(Temp,output_1(:,3));
output_2(:,4) = cumtrapz(Temp,output_2(:,3));
output_3(:,4) = cumtrapz(Temp,output_3(:,3));
output_4(:,4) = cumtrapz(Temp,output_4(:,3));

c = zeros(4,1);
c(1) = output_1(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));
c(2) = output_2(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));
c(3) = output_3(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));
c(4) = output_3(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));

output_1(:,4) = (output_1(:,4) - output_1(1,4))/(output_1(end,4)-output_1(1,4));
output_2(:,4) = (output_2(:,4) - output_2(1,4))/(output_2(end,4)-output_2(1,4));
output_3(:,4) = (output_3(:,4) - output_3(1,4))/(output_3(end,4)-output_3(1,4));
output_4(:,4) = (output_4(:,4) - output_4(1,4))/(output_4(end,4)-output_4(1,4));

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

function y = fsm4(T,pm)

y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12)) + fs(T,pm(13:16));

end

end
