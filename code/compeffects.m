function [A,CI_A] = compeffects(Ea,data)

ind = (data(:,4)>0.05).*(data(:,4)<0.9);
ind = logical(ind);
data = data(ind,:);


T = data(:,1);
t = data(:,2);
da = data(:,3);
a = data(:,4);
da = log(da);
b = (T(end)-T(1))/(t(end)-t(1));
%b = 1;


names = {'PL4','PL3','PL2','PL2/3','D1','D2','D3','GB','A4','A3','A2','A3/2','R3','R2','RS','F1','F2','F3','F4'};
% power laws
fa_1 = 4*a.^(3/4); %PL

fa_2 = 3*a.^(2/3); %PL

fa_3 = 2*a.^(1/2); %PL

fa_4 = 2/3*a.^(-1/2); %PL

%diffusion control

fa_5 = (a.^(-1))/2; %1D diffusion

fa_6 = (-log(1-a)).^(-1); % 2D diffusion

fa_7 = 3/2 * (1-a).^(2/3).*(1-(1-a).^(1/3)).^(-1) ; %3D diffusion

fa_8 = (3/2*(1-a).^(1/3)).^(-1) ; %ginstling - brounshtein

%nucleation

fa_9 = 4*(1-a).*(-log(1-a)).^(3/4); %A4

fa_10 = 3*(1-a).*(-log(1-a)).^(2/3); %A3

fa_11 = 2*(1-a).*(-log(1-a)).^(1/2); %A2

fa_12 = 1.5*(1-a).*(-log(1-a)).^(1/3); %A15


% geometrical contraction

fa_13 = 3*(1-a).^(2/3); % Contracting sphere

fa_14 = 2*(1-a).^(1/2); % contracting cylinder

%random scission

fa_15 = 2*(a.^(1/2)-a) ; % random scission


% order of reaction
fa_16 = 1-a; %1st order

fa_17 = (1-a).^2; %second order

fa_18 = (1-a).^3;  %third order

fa_19 = (1-a).^4;  %fourth order

fa_7(fa_7~=real(fa_7)) = 0;
fa_8(fa_8~=real(fa_8)) = 0;
fa_9(fa_9~=real(fa_9)) = 0;
fa_13(fa_13~=real(fa_13)) = 0;

fa = real([fa_1 fa_2 fa_3 fa_4 fa_5 fa_6 fa_7 fa_8 fa_9 fa_10 fa_11 fa_12 fa_13 fa_14 fa_15 fa_16 fa_17 fa_18 fa_19]);


R = 8.314E-3;

for i = 1:19
    
izq = da + log(b) - log(fa(:,i));
%izq = da  - log(fa(:,i));
M = [ones(size(da)),-1./(R*T)];
MM = M'*M;
izq = M'*izq;
EAA(:,i) = MM\izq;
end



   x = EAA(2,:);
   y = EAA(1,:);
   
   
   
    [p,S]=polyfit(x,y,1);
    % Calculate slope standard error and t-value
alpha = 0.05;
s = sqrt(S.normr / (S.df) * (1 / diff(x([1 end]))^2));
t = tinv(1 - alpha / 2, S.df);

logA = polyval(p,Ea);

se_logA = s * sqrt(1 / length(x) + (Ea - mean(x)).^2 / sum((x - mean(x)).^2));
CI_logA = [logA - t * se_logA; logA + t * se_logA];

% Calculate confidence interval for A
A = exp(logA);
CI_A = exp(CI_logA);    
    
 xx = [x  Ea];
    lx = [min(xx) max(xx)];
    ly = polyval(p,lx);
   hfig = figure;
   hold on
   plot(x,y,'o','LineWidth',1.5)
   plot(Ea,logA,'*','LineWidth',1.5)
   plot(lx,ly,'-','LineWidth',1.5)
   hold off
xlabel('E_{\alpha}')
ylabel('ln(A_{\alpha})')
legend('Modelos','Isoconversional','Línea de mínimos cuadrados','Location','southeast')

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional

set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%print(hfig,'grafico11','-dpdf','-r0')




   
   