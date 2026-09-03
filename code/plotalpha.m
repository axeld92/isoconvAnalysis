function plotalpha(data_1,data_2,data_3,data_4)

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

massloss_1 = data_1(:,3);
massloss_2 = data_2(:,3);
massloss_3 = data_3(:,3);
massloss_4 = data_4(:,3);

%calculate conversion
alpha_1 = (massloss_1(1)-massloss_1)/(massloss_1(1)-massloss_1(end-1));
alpha_2 = (massloss_2(1)-massloss_2)/(massloss_2(1)-massloss_2(end-1));
alpha_3 = (massloss_3(1)-massloss_3)/(massloss_3(1)-massloss_3(end-1));
alpha_4 = (massloss_4(1)-massloss_4)/(massloss_4(1)-massloss_4(end-1));


a_1 = smoothdata(alpha_1,'sgolay',50);
a_2 = smoothdata(alpha_2,'sgolay',50);
a_3 = smoothdata(alpha_3,'sgolay',50);
a_4 = smoothdata(alpha_4,'sgolay',50);

figure
hold on
plot(time_1)
plot(time_2)
plot(time_3)
plot(time_4)
hold off
legend('5','10','15','20')
grid