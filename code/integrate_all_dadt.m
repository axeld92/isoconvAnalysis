fs_1 = fs_fit(St01);
fs_2 = fs_fit(St06);
fs_3 = fs_fit(St11);
fs_4 = fs_fit(St16);

data_1 = St01;
data_2 = St06;
data_3 = St11;
data_4 = St16;


%Get temperatures
Temp_1 = data_1(:,1)+273.15;
Temp_2 = data_2(:,1)+273.15;
Temp_3 = data_3(:,1)+273.15;
Temp_4 = data_4(:,1)+273.15;

%Get times
time_1 = data_1(:,2)-data_1(1,2);
time_2 = data_2(:,2)-data_2(1,2);
time_3 = data_3(:,2)-data_3(1,2);
time_4 = data_4(:,2)-data_4(1,2);

%get massloss

massloss_1 = data_1(:,4);
massloss_2 = data_2(:,4);
massloss_3 = data_3(:,4);
massloss_4 = data_4(:,4);

%calculate conversion
alpha_1 = (massloss_1(1)-massloss_1)/(massloss_1(1)-massloss_1(end));
alpha_2 = (massloss_2(1)-massloss_2)/(massloss_2(1)-massloss_2(end));
alpha_3 = (massloss_3(1)-massloss_3)/(massloss_3(1)-massloss_3(end));
alpha_4 = (massloss_4(1)-massloss_4)/(massloss_4(1)-massloss_4(end));

ind = getind(alpha_1);
ind = logical(ind);
alpha_1 = alpha_1(ind);
Temp_1 = Temp_1(ind);
time_1 = time_1(ind);

ind = getind(alpha_2);
ind = logical(ind);
alpha_2 = alpha_2(ind);
Temp_2 = Temp_2(ind);
time_2 = time_2(ind);

ind = getind(alpha_3);
ind = logical(ind);
alpha_3 = alpha_3(ind);
Temp_3 = Temp_3(ind);
time_3 = time_3(ind);

ind = getind(alpha_4);
ind = logical(ind);
alpha_4 = alpha_4(ind);
Temp_4 = Temp_4(ind);
time_4 = time_4(ind);



aa_1 = zeros(length(Temp_1),4);

for ii = 2:length(Temp_1)
   
    aa_1(ii,:) = trapz(Temp_1(1:ii),fs_1(1:ii,:));
    
end

% figure;hold on;
% plot(Temp_1,aa_1(:,1))
% plot(Temp_1,aa_1(:,2))
% plot(Temp_1,aa_1(:,3))
% plot(Temp_1,aa_1(:,4))
% legend('mix','1','2','3')
% title('5 K/min')

aa_2 = zeros(length(Temp_2),4);

for ii = 2:length(Temp_2)
   
    aa_2(ii,:) = trapz(Temp_2(1:ii),fs_2(1:ii,:));
    
end

% figure;hold on;
% plot(Temp_2,aa_2(:,1))
% plot(Temp_2,aa_2(:,2))
% plot(Temp_2,aa_2(:,3))
% plot(Temp_2,aa_2(:,4))
% legend('mix','1','2','3')
% title('10 K/min')


aa_3 = zeros(length(Temp_3),4);

for ii = 2:length(Temp_3)
   
    aa_3(ii,:) = trapz(Temp_3(1:ii),fs_3(1:ii,:));
    
end

% figure;hold on;
% plot(Temp_3,aa_3(:,1))
% plot(Temp_3,aa_3(:,2))
% plot(Temp_3,aa_3(:,3))
% plot(Temp_3,aa_3(:,4))
% legend('mix','1','2','3')
% title('15 K/min')




aa_4 = zeros(length(Temp_4),4);

for ii = 2:length(Temp_4)
   
    aa_4(ii,:) = trapz(Temp_4(1:ii),fs_4(1:ii,:));
    
end

% figure;hold on;
% plot(Temp_4,aa_4(:,1))
% plot(Temp_4,aa_4(:,2))
% plot(Temp_4,aa_4(:,3))
% plot(Temp_4,aa_4(:,4))
% legend('mix','1','2','3')
% title('20 K/min')

% figure
% hold on
% plot(Temp_1,aa_1(:,3));
% plot(Temp_2,aa_2(:,3));
% plot(Temp_3,aa_3(:,3));
% plot(Temp_4,aa_4(:,3));
% legend('1','2','3','4')
% grid


alphas = 0.1:0.1:0.9;
datos_5 = [Temp_1, time_1, zeros(size(Temp_1)), aa_1(:,2)];
datos_10 = [Temp_2, time_2, zeros(size(Temp_2)), aa_2(:,2)];
datos_15 = [Temp_3, time_3, zeros(size(Temp_3)), aa_3(:,2)];
datos_20 = [Temp_4, time_4, zeros(size(Temp_4)), aa_4(:,2)];
vya1 = Vyazovkina(alphas,datos_5,datos_10,datos_15,datos_20);
% 
figure
plot(alphas,vya1(:,1));