function [mnpc,ci] = sbrlin(A,Ea,data,colorsrs)

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
loglog1a = log(-log(1-a));

%M = [ones(size(loga)) loga log1a loglog1a];
%M = [ones(size(loga)) loga log1a];
%M = [loga log1a loglog1a];
M = [loga log1a];
lmdl = fitlm(M,lfa,'Intercept',false)

MM = M'*M;
yy = M'*lfa;

mnpc = MM\yy;


SSE = lfa'*lfa - mnpc'*M'*lfa;
[nn,pp] = size(M);

sigma2 = SSE/(nn-pp)

imvMM = inv(MM);

varcoef = zeros(nn,1);
for i=1:nn
    varcoef(i) = sqrt(sigma2*(1+M(i,:)*imvMM*M(i,:)'));
    ci(i) = tinv(0.05/2,nn-pp)*varcoef(i);
end
ci = ci';




%fitlfa = mnpc(1) + mnpc(2)*loga + mnpc(3)*log1a + mnpc(4)*loglog1a;
%fitlfa = mnpc(1) + mnpc(2)*loga + mnpc(3)*log1a;
%fitlfa = mnpc(1)*loga + mnpc(2)*log1a + mnpc(3)*loglog1a;
fitlfa = mnpc(1)*loga + mnpc(2)*log1a;

%mnpc(1) = exp(mnpc(1));
%colorsrs = [0.9290 0.6940 0.1250]; 
%hfig = figure;
%hold on
% fill([a ; flipud(a)],[fitlfa+ci ; flipud(fitlfa-ci)],colorsrs,'FaceAlpha',0.1,'EdgeColor','none')
% plot(a,lfa,'--','LineWidth',1.5,'Color',	colorsrs)
% plot(a,fitlfa,'-','LineWidth',1.5,'Color',	colorsrs)
% drawnow
% hola = 1;
%plot(a,fitlfa+ci,'-.','LineWidth',1.5,'Color','#D95319')
%plot(a,fitlfa-ci,'-.','LineWidth',1.5,'Color','#D95319')
%hold off
%legend('Peak 1 - 95% CI','Peak 1 - Exp.','Peak 1 - SB model','Peak 2 - 95% CI','Peak 2 - Exp.','Peak 2 - SB model','Peak 3 - 95% CI','Peak 3 - Exp.','Peak 3 - SB model','Location','best')
%xlabel('\alpha')
%ylabel('log[f(\alpha)]')
%set(findall(hfig,'-property','FontSize'),'FontSize',13)
%set(findall(hfig,'-property','Box'),'Box','off') % optional


