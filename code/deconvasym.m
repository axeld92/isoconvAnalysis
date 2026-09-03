function [output,c] = deconvasym(varargin)

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
    da_todos = [] ;
    da_todos = [da_todos ; da{i}];
    dar{i} = da{i};
    %da{i} = da{i}/max(da{i});
    Tr{i} = T{i};
    T_todos = [];
    T_todos = [T_todos ; T{i}];
end

T_min = min(T_todos);
T_max = max(T_todos);
da_max = max(da_todos);

for i = 1:n
T{i} = (T{i}-T_min)/(T_max - T_min);
da{i} = da{i}/da_max;
end

% %%
% 
% T = 400:800;
% Tp = 600;
% theta = 10;
% w1 = 0;
% w2 = 50;
% w3 = 50*2;
% 
% 
% func = asymfun(T,[Tp,theta,w1,w2,w3]);
% plot(T,func)
% grid

%%

TT = T{1};
a_exp = a{1};
da_exp = da{1};

obj_fun =@(pa) sum((da_exp - asymfun3(TT,pa(1:5),pa(6:10),pa(11:15))).^2);
pa0 = 0.5*ones(1,15);
lb = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
ub = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];



%options = optimset('PlotFcns',@optimplotfval);
%pa_opt = fmincon(obj_fun,pa0,[],[],[],[],lb,ub);
pa_opt = surrogateopt(obj_fun,lb,ub)

pa_opt
da_calc = asymfun3(TT,pa_opt(1:5),pa_opt(6:10),pa_opt(11:15));
a_calc = cumtrapz(da_calc);
a_calc = a_calc/max(a_calc);

da_calc_1 = asymfun(TT,pa_opt(1:5));
da_calc_2 = asymfun(TT,pa_opt(6:10));
da_calc_3 = asymfun(TT,pa_opt(11:15));

figure
hold on
plot(TT,da_exp,'--')
plot(TT,da_calc);
plot(TT,da_calc_1);
plot(TT,da_calc_2);
plot(TT,da_calc_3);
hold off

figure
hold on
plot(TT,a_exp,'--')
plot(TT,a_calc)
hold off




%%

function f = asymfun(T,p)
Tp = p(1);
theta = p(2);
w1 = p(3);
w2 = p(4);
w3 = p(5);

f = theta./(1 + exp(-(T-Tp + w1/2)/w2)).*(1-1./(1 + exp(-(T-Tp + w1/2)/w3))) ;
end

function f = asymfun3(T,p1,p2,p3)

f = asymfun(T,p1) + asymfun(T,p2) + asymfun(T,p3);

end

end

