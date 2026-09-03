%Script para la limpieza y adecuación de datos
function data_out = cleandata2(data)

%remover datos anteriores a T = 110 °C
ind = (data(:,1)>100).*(data(:,1)<700);
ind = logical(ind);
 
data = data(ind,:);

%eliminar ultima fila de nan
if isnan(data(end,:))
    data(end,:) = [];
end

T = data(:,1)+273.15;
t = data(:,2)-data(1,2);
m = data(:,3);

p = polyfit(t,T,1);
beta = p(1);


%calcular conversión
a = (m(1)-m)/(m(1)-m(end));

%eliminar errores con remove error
%   [t,T,a] = removeerror(t,T,a);

%  a = smoothdata(a,'sgolay');




%calcular derivada da/dT
% figure
%dadt = [0;diff(a)./diff(T)];
 dadt = diff_alpha(t,a);
 %ddadt = diff_alpha(t,dadt);
  %dadt = lowpass(dadt,0.00001,'Steepness',0.99);
  %ddadt = lowpass(ddadt,0.00001,'Steepness',0.99);
 % dadt = smoothdata(dadt,'sgolay',500);
 
dadT = dadt/beta;
%ddadT = ddadt/beta/beta;
%  
% figure
% plot(T,a)
% yyaxis right
% plot(T,dadT)


 
 
 % plot(T,dadt)
% hold on


% dadt = smoothdata(dadt,'sgolay',500);
% plot(T,dadt)
% hold off

data_out = [T t dadt a ];%ddadT];




% figure
% yyaxis left
% plot(data_out(:,1),data_out(:,4))
% ylim([0,1])
% yyaxis right
% plot(data_out(:,1),data_out(:,3))


