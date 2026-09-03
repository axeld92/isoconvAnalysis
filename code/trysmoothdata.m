dadT = diff_alpha(T,a);
figure;%plot(T,dadT);
sdadT = smoothdata(dadT,'sgolay',1000);
hold on
plot(T,sdadT)
b = smoothdadT2(T,a);
plot(T,b(T))
hold off