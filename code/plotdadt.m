function plotdadt(data_1,data_2,data_3,data_4)

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

dadT_1 = diff_alpha(Temp_1,alpha_1);
dadT_2 = diff_alpha(Temp_2,alpha_2);
dadT_3 = diff_alpha(Temp_3,alpha_3);
dadT_4 = diff_alpha(Temp_4,alpha_4);

sdadT_1 = smoothdata(dadT_1,'sgolay',500);
sdadT_2 = smoothdata(dadT_2,'sgolay',500);
sdadT_3 = smoothdata(dadT_3,'sgolay',500);
sdadT_4 = smoothdata(dadT_4,'sgolay',500);

figure
hold on
plot(Temp_1,sdadT_1)
plot(Temp_2,sdadT_2)
plot(Temp_3,sdadT_3)
plot(Temp_4,sdadT_4)
hold off
legend('5','10','15','20')
grid
