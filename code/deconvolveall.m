function [output_1,output_2,output_3,c] = deconvolveall(data1,data2,data3,data4)

% clean data
data1 = cleandata(data1);
data2 = cleandata(data2);
data3 = cleandata(data3);
data4 = cleandata(data4);

%extract data
da1 = data1(:,3);
T1 = data1(:,1);
Tr1 = T1;
T1 = (Tr1-Tr1(1))/(Tr1(end)-Tr1(1));
dar1 = da1;
dam1 = max(dar1);
da1 = dar1/dam1;

da2 = data2(:,3);
T2 = data2(:,1);
Tr2 = T2;
T2 = (Tr2-Tr2(1))/(Tr2(end)-Tr2(1));
dar2 = da2;
dam2 = max(dar2);
da2 = dar2/dam2;

da3 = data3(:,3);
T3 = data3(:,1);
Tr3 = T3;
T3 = (Tr3-Tr3(1))/(Tr3(end)-Tr3(1));
dar3 = da3;
dam3 = max(dar3);
da3 = dar3/dam3;

da1 = data1(:,3);
T4 = data4(:,1);
Tr4 = T4;
T4 = (Tr4-Tr4(1))/(Tr4(end)-Tr4(1));
dar4 = da4;
dam4 = max(dar4);
da4 = dar4/dam4;

lb = [0.01,    -10,   0.01,    0.01,...
      0.01,    -10,   0.01,    0.01,...
      0.01,    -10,   0.01,    0.01];
  
ub = [0.99,  10,   0.99,    0.99,...
      0.99,  10,   0.99,    0.99,...
      0.99,  10,   0.99,    0.99];

p0 = (lb + ub)/2;

%Set objective function
obj_fun =@(parameters) (sum((fsm(Temp,parameters)-fs_data).^2));%/(length(Temp)-12);

rng default % For reproducibility
opts = optimoptions(@fmincon,'OptimalityTolerance',1E-8,'StepTolerance',1E-10,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
gs = GlobalSearch;
problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
[fsfit1,fval1] = run(gs,problem);

A = [   fsfit1(1:4);
        fsfit1(5:8);
        fsfit1(9:12)];
        
    
AA = sortrows(A,1);

fsfit1 = [AA(2,:) AA(3,:) AA(1,:)];
lb(:,3) = fsfit1(:,3);
ub(:,1) = fsfit1(:,1);




T1 = Tr1;
da1 = da1;
ind = [3 7 11];
fsfit1(ind) = fsfit1(ind)*(T1(end)-T(1)) + T(1); 
ind = [4 8 12];
fsfit1(ind)  = fsfit1(ind)*(Temp1(end)-Temp1(1));
ind = [1 5 9];
fsfit1(ind) = fsfit1(ind)*dam1;

fsfit1
fval1
MSE1 = sum((fsm(Temp,fsfit)-fs_data).^2)/(length(fs_data)-12)


























% %  Plot results
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