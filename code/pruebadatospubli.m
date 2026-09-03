%try published data

  fs_param_1 = [0.0201 -0.0298 570 79.60 0.0280 -0.2886 605 31.84 0.0016 -8E-5 659 331.70];
% fs_param_2 = [0.0198 -0.0288 583 79.83 0.0237 -0.2322 618 31.94 0.0018 -8E-5 663 322.91];
% fs_param_3 = [0.01788 -0.0288 595 79.95 0.0204 -0.2417 631 32.69 0.0023 -9.3E-5 669 317.86];

%fs_param_1 = [0.0363    0.0682  560.0000   55.1830    0.0597   -0.2543  632.3738   35.0086    0.0118    0.9939  609.8258  160.0000];

Temp = (400:1:1000)';

fs_1 = fs_mixture2(fs_param_1,Temp);
% fs_2 = fs_mixture2(fs_param_2,Temp);
% fs_3 = fs_mixture2(fs_param_3,Temp);

% fs_1_1 = fs_function(Temp,fs_param_1(1:4));
% fs_1_2 = fs_function(Temp,fs_param_1(5:8));
% fs_1_3 = fs_function(Temp,fs_param_1(9:12));
% 
% figure
% plot(Temp,fs_1)
% hold on
% plot(Temp,fs_1_1)
% plot(Temp,fs_1_2)
% plot(Temp,fs_1_3)
% hold off


time_1 = (Temp-400)/5;
% time_2 = (Temp-400)/10;
% time_3 = (Temp-400)/20;

data_1 = [Temp, time_1 , fs_1, zeros(size(Temp))];
% data_2 = [Temp, time_2 , fs_2, zeros(size(Temp))];
% data_3 = [Temp, time_3 , fs_3, zeros(size(Temp))];
% 

for i=2:length(Temp)
    data_1(i,4) = trapz(data_1(1:i,1),data_1(1:i,3));
end
   
% for i=2:length(Temp)
%     data_2(i,4) = trapz(Temp(1:i),data_2(1:i,3));
% end
% 
% for i=2:length(Temp)
%     data_3(i,4) = trapz(Temp(1:i),data_3(1:i,3));
% end
% 
% data_1(:,4) = (data_1(1,4) - data_1(:,4))/(data_1(1,4) - data_1(end,4));
% 



[out_1_05,out_2_05,out_3_05] = fs_out2(data_1);
% [out_1_10,out_2_10,out_3_10] = fs_out2(data_2);
% [out_1_20,out_2_20,out_3_20] = fs_out2(data_3);

% %aa = 0.2:0.1:0.9;
% 
% figure
% hold on
% plot(out_1_05(:,1),out_1_05(:,4))
% plot(out_1_10(:,1),out_1_10(:,4))
% plot(out_1_15(:,1),out_1_15(:,4))
% plot(out_1_20(:,1),out_1_20(:,4))
% hold off
% grid
% legend('5','10','15','20')
% title('hemicellulose')
% % 
% 
% figure
% hold on
% plot(out_2_05(:,1),out_2_05(:,4))
% plot(out_2_10(:,1),out_2_10(:,4))
% plot(out_2_15(:,1),out_2_15(:,4))
% plot(out_2_20(:,1),out_2_20(:,4))
% hold off
% grid
% legend('5','10','15','20')
% title('cellulose')
% 
% 
% figure
% hold on
% plot(out_3_05(:,1),out_3_05(:,4))
% plot(out_3_10(:,1),out_3_10(:,4))
% plot(out_3_15(:,1),out_3_15(:,4))
% plot(out_3_20(:,1),out_3_20(:,4))
% hold off
% grid
% legend('5','10','15','20')
% title('lignin')
% 
% 
% figure
% hold on
% plot(out_1_05(:,1),out_1_05(:,3))
% plot(out_1_10(:,1),out_1_10(:,3))
% plot(out_1_15(:,1),out_1_15(:,3))
% plot(out_1_20(:,1),out_1_20(:,3))
% hold off
% grid
% legend('5','10','15','20')
% title('hemicellulose')
% % 
% 
% figure
% hold on
% plot(out_2_05(:,1),out_2_05(:,3))
% plot(out_2_10(:,1),out_2_10(:,3))
% plot(out_2_15(:,1),out_2_15(:,3))
% plot(out_2_20(:,1),out_2_20(:,3))
% hold off
% grid
% legend('5','10','15','20')
% title('cellulose')
% 
% 
% figure
% hold on
% plot(out_3_05(:,1),out_3_05(:,3))
% plot(out_3_10(:,1),out_3_10(:,3))
% plot(out_3_15(:,1),out_3_15(:,3))
% plot(out_3_20(:,1),out_3_20(:,3))
% hold off
% grid
% legend('5','10','15','20')
% title('lignin')
% 
% 
