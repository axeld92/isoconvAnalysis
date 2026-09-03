function lambdamasterplot(Ea,data)
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
b = (T(end)-T(1))/(t(end)-t(1));
R = 8.314/1000;



da_05 = interp1(a,da,0.5);
T_05 = interp1(a,T,0.5);


lambda_t = da.*exp(-Ea/R/T_05);
lambda_b = da_05.*exp(-Ea/R./T);
lambda = lambda_t./lambda_b;

% Calculate model from data
%fa = b*da./(A*exp(-Ea./(R*T)));


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
figure
hold on
for i = 1:19
    fa_05 = interp1(a,fa(:,i),0.5);
    lambda_f = fa(:,i)/fa_05;
    plot(a,lambda_f)
    MSD(i) = sum((lambda_f-lambda).^2);
end

plot(a,lambda,'--')
[minimum,mindex] = min(MSD);
names = {'PL4','PL3','PL2','PL2/3','1D','2D','3D','GB','A4','A3','A2','A3/2','R3','R2','RS','F1','F2','F3','F4'};
closest = names{mindex}