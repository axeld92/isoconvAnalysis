
data1 = cleandatam(St01);
% data2 = cleandatam(St06);
% data3 = cleandatam(St11);
% data4 = cleandatam(St16);



[out_1_05,out_2_05,out_3_05] = fs_out(data1);
% [out_1_10,out_2_10,out_3_10] = fs_out(data2);
% [out_1_15,out_2_15,out_3_15] = fs_out(data3);
% [out_1_20,out_2_20,out_3_20] = fs_out(data4);

%aa = 0.2:0.1:0.9;

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
% % data_1 = out_1_05;
% % data_2 = out_1_10;
% % data_3 = out_1_15;
% % data_4 = out_1_20;
% 
% 
% 
% % Temp_1 = data_1(:,1);
% % Temp_2 = data_2(:,1);
% % Temp_3 = data_3(:,1);
% % Temp_4 = data_4(:,1);
% 
% %Get times
% % time_1 = data_1(:,2)-data_1(1,2);
% % time_2 = data_2(:,2)-data_2(1,2);
% % time_3 = data_3(:,2)-data_3(1,2);
% % time_4 = data_4(:,2)-data_4(1,2);
% 
% % alpha_1 = data_1(:,4);
% % alpha_2 = data_2(:,4);
% % alpha_3 = data_3(:,4);
% % alpha_4 = data_4(:,4);
% 
% 
% 
% % T_1 = Tfit(alpha_1,Temp_1);
% % T_2 = Tfit(alpha_2,Temp_2);
% % T_3 = Tfit(alpha_3,Temp_3);
% % T_4 = Tfit(alpha_4,Temp_4);
% 
% % t_1 = timefit(alpha_1,time_1);
% % t_2 = timefit(alpha_2,time_2);
% % t_3 = timefit(alpha_3,time_3);
% % t_4 = timefit(alpha_4,time_4);
% 
