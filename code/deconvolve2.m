function [output_1,output_2,output_3,c] = deconvolve2(data)

data = cleandata(data);

fs_data = data(:,3);
Temp = data(:,1);
time = data(:,2);

% plot(Temp,fs_data)
% %find peaks
% 
% [peaks , locp, w, prom ]= findpeaks(fs_data,'MinPeakProminence',0.0005)
%  Tp = Temp(locp);


 % % h s p w
%set lower and upper bounds
lb = [500   0       500     0       600     50 0 0 0 ];
ub = [650   100     700     100     800     300 1 1 1];

%Set initial guess
% p0 = [560   55      632     35 660 160 0.5 0.5 0.5];
p0 = (lb+ub)/2

% p0 = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];
%  p0 = [0.0036   -0.1939  Tp(1)   50.0000    0.0076   -0.3300  Tp(2)   58.9930    0.0014    1.0584  630.5029  165.8270]; 

obj_fun =@(parameters) (sum((fmgauss(Temp,parameters)-fs_data).^2));

% PATTERN SEARCH
rng default
opts = optimoptions('patternsearch','MaxIterations',1E10,'MeshTolerance',1E-15,'StepTolerance',1E-15,...
    'MaxFunctionEvaluations',1E20,'MeshExpansionFactor',2);
[fsfit,fval] = patternsearch(obj_fun,p0,[],[],[],[],lb,ub,opts)


fsfit

% options = optimoptions('fmincon','OptimalityTolerance',1E-15,...
%     'StepTolerance',1E-30,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
% [fsfit , fval] = fmincon(obj_fun,fsfit,[],[],[],[],lb,ub,[],options);


% %  Plot results
figure
hold on
plot(Temp,fs_data,'-')
plot(Temp,fmgauss(Temp,fsfit))
plot(Temp,fsfit(7)*normpdf(Temp,fsfit(1),fsfit(2)))
plot(Temp,fsfit(8)*normpdf(Temp,fsfit(3),fsfit(4)))
plot(Temp,fsfit(9)*normpdf(Temp,fsfit(5),fsfit(6)))
hold off
legend('real','calc','1','2','3')%,'1c','2c','3c')

% Calculate outputs

%outp = fs_mixture(Temp,fsfit(1:4),fsfit(5:8),fsfit(9:12));
out_1 = fs_function(Temp,fsfit(1:4));
out_2 = fs_function(Temp,fsfit(5:8));
out_3 = fs_function(Temp,fsfit(9:12));

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

end