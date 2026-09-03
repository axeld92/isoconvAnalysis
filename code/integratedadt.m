aa = zeros(size(Temp));
aa_1 = zeros(size(Temp));
aa_2 = zeros(size(Temp));
aa_3 = zeros(size(Temp));





for ii = 2:length(Temp)
    aa_1(ii) = trapz(Temp(1:ii),out_1(1:ii));
end

for ii = 2:length(Temp)
    aa_2(ii) = trapz(Temp(1:ii),out_2(1:ii));
end

for ii = 2:length(Temp)
    aa_3(ii) = trapz(Temp(1:ii),out_3(1:ii));
end
figure
plot(Temp,aa,Temp,aa_1,Temp,aa_2,Temp,aa_3)

legend('mix','1','2','3')