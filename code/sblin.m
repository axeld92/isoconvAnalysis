function mnpc = sblin(A,Ea,data)

% data = out_1_05;
% A = A_1;
% Ea = ea1_avg;


ind = (data(:,4)>0.1).*(data(:,4)<0.9);
ind = logical(ind);
data = data(ind,:);

T = data(:,1);
t = data(:,2);
da = data(:,3);
a = data(:,4);

%plot(T,a)


b = (T(end)-T(1))/(t(end)-t(1));

R = 8.314/1000;

% Calculate model from data
fa = b*da./(A*exp(-Ea./(R*T)));

lfa = log(fa);

loga = log(a);
log1a = log(1-a);
%loglog1a = log(-log(1-a));

M = [ones(size(loga)) loga log1a];% loglog1a];
MM = M'*M;
yy = M'*lfa;

mnpc = MM\yy;

fitlfa = mnpc(1) + mnpc(2)*loga + mnpc(3)*log1a;% + mnpc(3)*loglog1a;

mnpc(1) = exp(mnpc(1));

figure
plot(a,lfa,'-',a,fitlfa,'--')
legend('Datos','Ajustado')
xlabel('\alpha')
ylabel('log[f(\alpha)]')



