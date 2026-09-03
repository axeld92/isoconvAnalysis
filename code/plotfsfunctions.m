data_1 = cleandata(St01);

T = data_1(:,1);
a = data_1(:,4);
 dadT = diff_alpha(T,a);


% a_s = lowpass(a,0.0001,'Steepness',0.99);

% dadT = diff_alpha(T,a);

dadT_s = lowpass(dadT,0.00001,'Steepness',0.99);
% ftdadt = fft(dadT_s);
% PSD = ftdadt.*conj(ftdadt);

% figure;plot(PSD)
figure
plot(T,dadT,T,dadT_s)
% figure;plot(T,a,T,a_s)





















% data_2 = cleandata(St06);
% data_3 = cleandata(St11);
% data_4 = cleandata(St16);
% 
% da_1 = data_1(:,3);
% da_2 = data_2(:,3);
% da_3 = data_3(:,3);
% da_4 = data_4(:,3);
% 
% T_1 = data_1(:,1);
% T_2 = data_2(:,1);
% T_3 = data_3(:,1);
% T_4 = data_4(:,1);
% 
% Temp = 390:0.1:1000;
% da_1 = interp1(T_1,da_1,Temp)';
% da_2 = interp1(T_2,da_2,Temp)';
% da_3 = interp1(T_3,da_3,Temp)';
% da_4 = interp1(T_4,da_4,Temp)';
% 
% X = [da_1 da_2 da_3 da_4];
% 
% 
% [U,S,V] = svd(X);

