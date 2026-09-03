[Temp_1,fs_1] = fs_fit(St05);
[Temp_2,fs_2] = fs_fit(St10);
[Temp_3,fs_3] = fs_fit(St15);
[Temp_4,fs_4] = fs_fit(St20);


figure
hold on
plot(Temp_1,fs_1(:,2))
plot(Temp_2,fs_2(:,2))
plot(Temp_3,fs_3(:,2))
plot(Temp_4,fs_4(:,2))
hold off
legend('5','10','15','20')