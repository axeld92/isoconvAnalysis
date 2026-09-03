function [output_1,output_2,output_3,c] = deconvolveg(data)

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
% fs_datar = fs_data;
% hmax = max(fs_data);
% fs_data = fs_datar/hmax;
% figure


%set lower and upper bounds
lb =    [0  0   0   0   0   0   0   0   0];
ub =    [1  1   1   1   1   1   1   1   1];
%Set initial guess

 p0 = (lb + ub)/2;

%Set objective function
obj_fun =@(parameters) (sum((fmgauss(Temp,parameters)-fs_data).^2));
% % PARTICLE SWARM
% rng default;
% options = optimoptions('fmincon','OptimalityTolerance',1E-15,...
%     'StepTolerance',1E-30,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
% opts = optimoptions('particleswarm','SwarmSize',1000,'HybridFcn',{@fmincon,options});
% [fsfit , fval] = particleswarm(obj_fun,12,lb,ub,opts)

rng default % For reproducibility
 opts = optimoptions(@fmincon,'OptimalityTolerance',1E-15,'StepTolerance',1E-10,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
[fsfit,fval] = run(gs,problem);

A = [   fsfit(1:3);
        fsfit(4:6);
        fsfit(7:9)];
%         
%     
 AA = sortrows(A,1);
% 
 fsfit = [AA(2,:) AA(1,:) AA(3,:)];

Temp = Tempr;
% fs_data = fs_datar;
ind = [2 5 8];
fsfit(ind) = fsfit(ind)*(Temp(end)-Temp(1)) + Temp(1); 
ind = [3 6 9];
fsfit(ind)  = fsfit(ind)*(Temp(end)-Temp(1));
% ind = [1 4 7];
% fsfit(ind) = fsfit(ind)*hmax;
% 
fsfit
fval


% %  Plot results
figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fmgauss(Temp,fsfit))
plot(Temp,fgauss(Temp,fsfit(1:3)))
plot(Temp,fgauss(Temp,fsfit(4:6)))
plot(Temp,fgauss(Temp,fsfit(7:9)))
hold off
legend('real','calc','1','2','3')%,'1c','2c','3c')

% Calculate outputs

%outp = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
out_1 = fgauss(Temp,fsfit(1:3));
out_2 = fgauss(Temp,fsfit(4:6));
out_3 = fgauss(Temp,fsfit(7:9));

output_1 = [Temp, time, out_1, zeros(size(time))];
output_2 = [Temp, time, out_2, zeros(size(time))];
output_3 = [Temp, time, out_3, zeros(size(time))];

output_1(:,4) = cumtrapz(Temp,output_1(:,3));
output_2(:,4) = cumtrapz(Temp,output_2(:,3));
output_3(:,4) = cumtrapz(Temp,output_3(:,3));

c = zeros(3,1);
c(1) = output_1(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));
c(2) = output_2(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));
c(3) = output_3(end,4)/(output_1(end,4) + output_2(end,4) + output_3(end,4));

output_1(:,4) = (output_1(:,4) - output_1(1,4))/(output_1(end,4)-output_1(1,4));
output_2(:,4) = (output_2(:,4) - output_2(1,4))/(output_2(end,4)-output_2(1,4));
output_3(:,4) = (output_3(:,4) - output_3(1,4))/(output_3(end,4)-output_3(1,4));




%Define Fraser Suzuki Function
function y = fs(T,p)
%Fraser Suzuki function
 in = 2 * p(2) * ((T-p(3)) / p(4));
 out = - log(2) / p(2)^2;
 
 y = p(1) * exp(out .* (log(1+in)).^2);
 y(y~=real(y))=0;
end

function y = fmgauss(T,p)
%Multiple gaussian function

y = fgauss(T,p(1:3)) + fgauss(T,p(4:6)) + fgauss(T,p(7:9));
end

    function y = fgauss(T,p)

       y = normpdf(T,p(2),p(3));
       y = p(1)*y/max(y);
    end
%Define Fraser suzuki mixture

function y = fsm(T,pm)
%multiple fraser suzuki function
y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12));

end

end