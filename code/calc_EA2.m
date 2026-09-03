aa = 0.1:0.01:0.9;

EA_v1 = Vyazovkin(aa,St01,St06,St11,St16);
EA_v2 = Vyazovkin(aa,St02,St07,St12,St17);
EA_v3 = Vyazovkin(aa,St03,St08,St13,St18);
EA_v4 = Vyazovkin(aa,St04,St09,St14,St19);
EA_v5 = Vyazovkin(aa,St05,St10,St15,St20);

figure
hold on
errorbar(aa,EA_v1(:,1),EA_v1(:,2),'s-')
errorbar(aa,EA_v2(:,1),EA_v2(:,2),'o-')
errorbar(aa,EA_v3(:,1),EA_v3(:,2),'^-')
errorbar(aa,EA_v4(:,1),EA_v4(:,2),'v-')
errorbar(aa,EA_v5(:,1),EA_v5(:,2),'x-')
hold off
legend('0 % waste tires','25 % waste tires','50 % waste tires','75 % waste tires','100 % waste tires')
xlabel('Conversion')
ylabel('Activation energy [kJ/mol]')
title('Activation energy - Vyazovkin method')
grid



%pre exponential factor for 0% WT
p5_0 = comeffect3(St01,aa);
p10_0 = comeffect2(St06);
p15_0 = comeffect2(St11);
p20_0 = comeffect2(St16);


for i = 1:length(aa)
 y5_0(i) = polyval(p5_0(i,:),EA_v1(i,1));
end
y10_0 = polyval(p10_0,EA_v1(:,1));
y15_0 = polyval(p10_0,EA_v1(:,1));
y20_0 = polyval(p20_0,EA_v1(:,1));

figure
hold on
plot(aa,y5_0)
plot(aa,y10_0)
plot(aa,y15_0)
plot(aa,y20_0)
hold off
grid
legend('5 K/min','10 K/min','15 K/min','20 K/min')
xlabel('\alpha')
ylabel('ln(A_{\alpha})')
title('A_{\alpha} vs \alpha for 0% WT')
