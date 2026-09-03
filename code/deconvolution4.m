function [output,c] = deconvolution4(varargin)

n = nargin;

%Extract data
das = [];
for i = 1:n
    data{i} = varargin{i};
    data{i} = cleandata(data{i});
    T{i} = data{i}(:,1);
    t{i} = data{i}(:,2);
    a{i} = data{i}(:,4);
    da{i} = data{i}(:,3);
    das = [das ; da{i}];
    %dar{i} = da{i};
    %da{i} = da{i}/max(da{i});
    Tr{i} = T{i};
    T{i} = (T{i}-T{i}(1))/(T{i}(end)-T{i}(1));
end

dmax = max(das(:));
dar = da;
for i = 1:n
da{i} = dar{i}/dmax;
end


%set lower and upper bounds
% lb1 =    [0  0   0   0   0   0   0   0   0 0 0 0];
% ub1 =    [1  1   1   1   1   1   1   1   1 1 1 1];
% %Set initial guess
% 
% % p01 = [0.5 0.24 0.5 0.5 0.34 0.5 0.5 0.457 0.5 0.5 0.55 0.5];
% 
%   p01 = [0.35 0.227 0.5 0.62 0.345 0.5 0.14 0.504 0.5 0.045 0.72 0.5];


%Set lb and ub
% % 
lb = [0.001,    -0.2,   0.,    0.001,...
      0.001,    -0.2,   0.301,    0.001,...
      0.001,    -0.2,   0.4,    0.05,...
  0.001,    -0.2,   0.3,    0.05];
  
ub = [0.999,  0.2,   0.3,    0.999,...
      0.999,  0.2,   0.38,    0.999,...
      0.999,  0.2,   0.6,    0.5,...
  0.999,  0.2,   0.4,    0.5];


%   lb = [0.01,    -10,   0.01,    0.01,...
%       0.01,    -10,   0.01,    0.01,...
%       0.01,    -10,   0.01,    0.01,...
%       0.01,    -10,   0.01,    0.01];
%   
% ub = [0.99,  10,   0.99,    0.99,...
%       0.99,  10,   0.99,    0.99,...
%       0.99,  10,   0.99,    0.99,...
%       0.99,  10,   0.99,    0.99];
%   
  



% lb = [0.01,    -1,   0.01,    0.13,...
%       0.01,    -1,   0.01,    0.05,...
%       0.01,    -0.2,   0.01,    0.43];
%   
% ub = [0.99,  0.000001,   0.99,    0.18,...
%       0.99,  1,   0.99,    0.082,...
%       0.99,  0.000001,   0.99,    0.6];
rng default 
opts = optimoptions(@fmincon,'OptimalityTolerance',1E-8,'StepTolerance',1E-10,'MaxFunctionEvaluations',10000,'MaxIterations',1000);
gs = GlobalSearch; 
fits = {};
fvals = {};
p0 = (lb + ub)/2;

% % 50:50
%p0 = [0.355 0.01 0.26 0.12 ...
%     0.883 0.01 0.393 0.1 ...
%     0.244 0.01 0.529 0.12 ...
%     0.06 0.01 0.7 0.7];

% %75:25
% p0 = [0.4 0.1 0.22 0.2...
%     0.79 0.1 0.346 0.2 ...
%     0.21 0.1 0.48 0.2 ...
%     0.05 0.1 0.605 0.7];
% 
% 
% idx = [1 3 4 5 7 8 9 11 12];
% 
% lb(idx) = p0(idx) - p0(idx)*0.1;
% ub(idx) = p0(idx) + p0(idx)*0.15;


for i=1:n
   
%        
% %     obj_fun1 =@(parameters) (sum((fmgauss(T{i},parameters)-da{i}).^2));
% %     problem = createOptimProblem('fmincon','x0',p01,'objective',obj_fun1,'lb',lb1,'ub',ub1,'options',opts);
% %     [fsg,fvalg] = run(gs,problem);
% %     fsg
% %     
% %     idx1 = [1 2 4 5 7 8 10 11];
% %     idx2 = [1 3 5 7 9 11 13 15];
% %     lb(idx2) = fsg(idx1) - fsg(idx1)*0.1;
% %     ub(idx2) = fsg(idx1) + fsg(idx1)*0.1;
% % %     p0(idx2) = fsg(idx1);
% % %     p01 = fsg;
% % %     lb1(idx1) = p01(idx1) - p01(idx1)*0.1;
% %     ub1(idx1) = p01(idx1) + p01(idx1)*0.1;
    
       
    obj_fun =@(parameters) (sum((fsm4(T{i},parameters)-da{i}).^2));% + sum((fsmi(T{i},parameters)-a{i}).^2));
    rng default % For reproducibility
    problem = createOptimProblem('fmincon','x0',p0,'objective',obj_fun,'lb',lb,'ub',ub,'options',opts);
    
    [fsfit,fval] = run(gs,problem);
    fits{i} = fsfit;
    fvals{i} = fval;
    
       MSE = fval/(length(T{i})-12)
    
    p0 = fsfit;
    idx = [1,3,4,5,7,8,9,11,12,13,15,16];
    lb(idx) = p0(idx) - p0(idx)*0.05;
    ub(idx) = p0(idx) + p0(idx)*0.05;
    ub([1,5,9,13]) = fsfit([1,5,9,13]);
    lb([3,7,11,15]) = fsfit([3,7,11,15]);
    lb([4,8,12,16]) = fsfit([4,8,12,16]);
  
end

T = Tr;
da = dar;


for i = 1:n
    fitss{i}{1} = fits{i}(1:4);
    fitss{i}{2} = fits{i}(5:8);
    fitss{i}{3} = fits{i}(9:12);
    fitss{i}{4} = fits{i}(13:16);
end


for i = 1:n
    for j = 1:4
        fitss{i}{j}(1) = dmax*fitss{i}{j}(1);
        
        fitss{i}{j}(3) = (T{i}(end) - T{i}(1))*fitss{i}{j}(3) + T{i}(1);
        fitss{i}{j}(4) = (T{i}(end) - T{i}(1))*fitss{i}{j}(4);
       
        daa{i}{j} = fs(T{i},fitss{i}{j});
        aa{i}{j} = cumtrapz(T{i},daa{i}{j});
            end
  for j = 1:4
         c(i,j) = aa{i}{j}(end)/(aa{i}{1}(end) + aa{i}{2}(end) + aa{i}{3}(end) + aa{i}{4}(end));
     end
%   for j = 1:4
%       aa{i}{j} = aa{i}{j}/c(i,j);
%          daa{i}{j} = daa{i}{j}/c(i,j);
%   end
    
end
    
    for i = 1:n
    for j = 1:4
    texto = ['tasa ',num2str(i),' pico ',num2str(j)];
    disp(texto)
        
    fitss{i}{j}    
    end
    end
    
    
   
for i = 1:n
output{i} = {[T{i},t{i},daa{i}{1},aa{i}{1}],[T{i},t{i},daa{i}{2},aa{i}{2}],[T{i},t{i},daa{i}{3},aa{i}{3}],[T{i},t{i},daa{i}{4},aa{i}{4}]};
end





for i=1:n
figure
hold on
plot(T{i},dar{i},'LineWidth',1)
plot(T{i},fsm4(T{i},[fitss{i}{1},fitss{i}{2},fitss{i}{3},fitss{i}{4}]),'LineWidth',1)
plot(T{i},daa{i}{1}*c(i,1),'LineWidth',1)
plot(T{i},daa{i}{2}*c(i,2),'LineWidth',1)
plot(T{i},daa{i}{3}*c(i,3),'LineWidth',1)
plot(T{i},daa{i}{4}*c(i,4),'LineWidth',1)
hold off
grid
%ylim([-0.002 0.012])
legend('Real','Ajustado','Pico 1','Pico 2','Pico 3','Pico 4')
ylabel('d\alpha/dt [1/min]')
xlabel('Temperatura [K]')
end

% for i=1:n
% figure
% hold on
% plot(T{i},da{i},'LineWidth',1)
% plot(T{i},fsm4(T{i},[fitss{i}{1},fitss{i}{2},fitss{i}{3},fitss{i}{4}]),'LineWidth',1)
% plot(T{i},daa{i}{1},'LineWidth',1)
% plot(T{i},daa{i}{2},'LineWidth',1)
% plot(T{i},daa{i}{3},'LineWidth',1)
% plot(T{i},daa{i}{4},'LineWidth',1)
% hold off
% grid
% ylim([-0.002 0.012])
% legend('Real','Ajustado','Pico 1','Pico 2','Pico 3','Pico 4')
% ylabel('d\alpha/dt [1/min]')
% xlabel('Temperatura [K]')
% end

% for i=1:n
% figure
% hold on
% plot(T{i},a{i})
% plot(T{i},fsmi(T{i},[fitss{i}{1},fitss{i}{2},fitss{i}{3}]))
% hold off
% %ylim([-0.002 0.012])
% legend('real','all')
% end



% for i=1:n
% figure
% hold on
% plot(T{i},a{i})
% %plot(T{i},cumtrapz(fsm(T{i},[fitss{i}{1},fitss{i}{2},fitss{i}{3}])))
% plot(T{i},aa{i}{1} + aa{i}{2} + aa{i}{3})
% hold off
% 
% legend('real','fs')
% end



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


 function y = fsmi(T,pm)
     
        y = cumtrapz(T,fsm(T,pm));
        
 end


function y = fmgauss(T,p)
%Multiple gaussian function

y = fgauss(T,p(1:3)) + fgauss(T,p(4:6)) + fgauss(T,p(7:9) + fgauss(T,p(10:12)));
end

    function y = fgauss(T,p)

       y = normpdf(T,p(2),p(3));
       y = p(1)*y/max(y);
    end



end
    
    




% ddadT = diff_alpha(T{i},da{i});
% ddadT = lowpass(ddadT,0.000000001,'Steepness',0.9999999);
% [mind,idx] = min(ddadT);
% Tmin = T{i}(idx:end);
% hmin = da{i}(idx:end);
% idx = hmin<0.11;
% hmin = hmin(idx);
% Tmin = Tmin(idx);
% datail = hmin;
% Ttail = Tmin;
% hmin = hmin(1);
% Tmin = Tmin(1);
% 
% obj_tail =@(ptail) sum((datail - fs(Ttail,ptail)).^2);
%  ptail0 = [hmin,-0.1,Tmin,Tmin];
% lbtail = [hmin-0.1,-0.2,Tmin-0.1,Tmin-0.1];
% ubtail = [hmin+0.1,0.2,Tmin+0.1,Tmin+0.1];
% ptailfit = fmincon(obj_tail,ptail0,[],[],[],[],lbtail,ubtail)
% figure
% plot(Ttail,datail,Ttail,fs(Ttail,ptailfit))
% end    
%     
%     
% for i=1:n    
    %     obj_fun1 =@(parameters) (sum((fmgauss(T{i},parameters)-da{i}).^2));
%     problem = createOptimProblem('fmincon','x0',p01,'objective',obj_fun1,'lb',lb1,'ub',ub1,'options',opts);
%     [fsg,fvalg] = run(gs,problem);
%     fsg
%     p0 = [fsg(1) 0.01 fsg(2) fsg(3)*3 fsg(4) 0.01 fsg(5)*3 fsg(6) fsg(7) 0.01 fsg(8) fsg(9)*3];
%     lb = p0 - 0.2*p0;
%     ub = p0 + 0.2*p0;
%    
% lb([9,11,12]) = ptailfit([1,3,4]) - ptailfit([1,3,4])*0.01;
% ub([9,11,12])= ptailfit([1,3,4]) + ptailfit([1,3,4])*0.01;
% p0([9,11,12]) = ptailfit([1,3,4]);
    





