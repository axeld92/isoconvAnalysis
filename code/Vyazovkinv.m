function EA = Vyazovkinv(alphas,varargin)

lower = 0.01;
upper = 0.99;

n = length(varargin);
data = varargin;
 for i = 1:n
     idx = (data{i}(:,4)>lower).*(data{i}(:,4)<upper);
     idx = logical(idx);
     data{i} = data{i}(idx,:);
 end
Temp = {};
time = {};
alpha = {};
 
 
for i = 1:n
    Temp{i} = data{i}(:,1);
    time{i} = data{i}(:,2);
    alpha{i} = data{i}(:,4);
end

T = {};
t = {};
a = {};

for i = 1:n
    T{i} = Tfit(alpha{i},Temp{i});
    t{i} = timefit(alpha{i},time{i});

end




% figure
% hold on
% plot(alpha_1,Temp_1)
% plot(Temp_2,alpha_2)
% plot(Temp_3,alpha_3)
% plot(Temp_4,alpha_4)
% hold off

%create cfit objects for temperature and time as a function of conversion

%  T_1 = Tfit(alpha_1,Temp_1);
% %T_1 = griddedInterpolant(alpha_1,Temp_1);
% T_2 = Tfit(alpha_2,Temp_2);
% T_3 = Tfit(alpha_3,Temp_3);
% T_4 = Tfit(alpha_4,Temp_4);
% 
% t_1 = timefit(alpha_1,time_1);
% t_2 = timefit(alpha_2,time_2);
% t_3 = timefit(alpha_3,time_3);
% t_4 = timefit(alpha_4,time_4);
% 
% figure
% plot(T_1(alpha_1),alpha_1)
% hold on
% %plot(Temp_1,alpha_1)
% plot(T_2(alpha_2),alpha_2)
% plot(T_3(alpha_3),alpha_3)
% plot(T_4(alpha_4),alpha_4)
% hold off
% grid
% title('Conversion vs temperature')
% legend('1','2','3','4')
% xlabel('temperature')
% ylabel('conversion')
% 

%alphas = 0.2:0.01:0.8;

%create matrix for times and temperatures

%m=length(alphas);
delta_a = 0.02;
z=10;
Eea = zeros(length(alphas),3);

for i = 1:length(alphas)

aa = alphas(i);

for j = 1:n
    tt(:,j) = t{j}(aa-delta_a:delta_a/z:aa);
    TT(:,j) = T{j}(aa-delta_a:delta_a/z:aa);
end

E = 250;

phi =@(e) abs(phiofe2(tt,TT,e)-12);
  %options = optimset('MaxFunEvals',100000000,'TolFun',1e-15,'TolX',1e-15,'MaxIter',1000000,'PlotFcns',@optimplotfval);
 [Ea,fval] = fminsearch(phi,E);%options);
%[Ea,fval] = particleswarm(phi,1);
J_integrals = Jintegs(tt,TT,Ea);

siglevel = 0.05;
n = length(J_integrals);  % number of experiments

% Calculate minimum variance activation energy (Emin)
S2 = 0;
for ii = 1:n
    for jj = 1:n
        if ii ~= jj
            S2 = S2 + (J_integrals(ii)/J_integrals(jj) - 1)^2;
        end
    end
end
S2 = S2 / (n * (n-1));

F = finv(1-siglevel, n-1, n-1);
S2_lo = S2 * ((n-1) * F) / (n * (n-1) - F);
S2_up = S2 * ((n-1) * F) / (n * (n-1) + F);

E_lo = Ea - sqrt(S2_lo);
E_up = Ea + sqrt(S2_up);

errorp = [E_lo, E_up];
%errorp = fval;

Eea(i,:) = [Ea, errorp];
end
EA = Eea;
