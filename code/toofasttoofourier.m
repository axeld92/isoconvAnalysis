data = St01;
ind = (data(:,1)>110).*(data(:,1)<780);
ind = logical(ind);
 
data = data(ind,:);

%eliminar ultima fila de nan
if isnan(data(end,:))
    data(end,:) = [];
end

T = data(:,1)+273.15;
t = data(:,2);
m = data(:,3);


%calcular conversión
a = (m(1)-m)/(m(1)-m(end));

%eliminar errores con remove error
%  a = smoothdata(a,'sgolay',100);




%calcular derivada da/dT
% figure
 dadt = diff_alpha(T,a);
 dadts = smoothdata(dadt,'sgolay',500);

 
 %figure
 %plot(T,dadt)
% hold on

%plot(T,dadt)
%hold off
t = t - t(1);
t = t*60;
n = length(t);
%fa = fft(dadt,n);
fas = fft(a,n);
%fa = fftshift(fa);
%PSD = fa.*conj(fa)/n;
PSDs = fas.*conj(fas)/n;




dt = t(1000)-t(999);
freq = 1/(dt*n)*(0:n-1);
figure
%subplot(2,1,1)
%plot(freq,PSD);
hold on
plot(PSDs)
hold off
%subplot(2,1,2)
%plot(freq,PSD-PSDs)

PSDs(1000:end) = 0;
fay = ifft(PSDs);
figure
plot(fay)




% indices = (freq>0.6).*(freq < 0.7);
% indices = logical(indices);
% fa = fa(indices);
% %fa = fftshift(fa);
% 
% %PSDc = PSD.*indices;
% %fac = indices.*fa;
% %fac(fac==0)= [];
% fafilt = ifft(fa);




% 
% figure
% hold on
% %plot(t,dadt);
% plot(fafilt)
% hold off



%dadtf = smoothdata(ffilt,'sgolay',300);

% figure
% hold on
% plot(t,ffilt)
% %plot(t,dadtf)
% hold off


