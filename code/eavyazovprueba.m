aa = 0.2:0.01:0.78;

EA_v = Vyazovkinfor3(aa,St01,St06,St16);
figure
errorbar(aa,EA_v(:,1),EA_v(:,2))