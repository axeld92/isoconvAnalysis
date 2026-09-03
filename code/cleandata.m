%Script para la limpieza y adecuación de datos
function data_out = cleandata(data)


% figure
% plot(data(:,1),data(:,3),'LineWidth',1.5)

%remover datos anteriores a T = 110 °C
ind = (data(:,1)>150).*(data(:,1)<700);
ind = logical(ind);
 
data = data(ind,:);

%eliminar ultima fila de nan
if isnan(data(end,:))
    data(end,:) = [];
end

T = data(:,1)+273.15;
t = data(:,2)-data(1,2);
m = data(:,3);

% m = lowpass(m,0.000000000000001,'Steepness',0.99999)
%    figure
%    plot(T,m,'LineWidth',1.5)

 p = polyfit(t,T,1);
beta = p(1);


%calcular conversión
a = (m(1)-m)/(m(1)-m(end));
%  a = lowpass(a,0.000000000000001,'Steepness',0.99999);
%   figure
%   plot(T,a,'LineWidth',1.5)
%eliminar errores con remove error
%   [t,T,a] = removeerror(t,T,a);

 a = smoothdata(a,'sgolay',50);




%calcular derivada da/dT
% figure
%dadt = [0;diff(a)./diff(T)];
 dadt = diff_alpha(t,a);
dadt = lowpass(dadt,0.000000001,'Steepness',0.9999999); 
 % dadt = lowpass(dadt,0.000000001,'Steepness',0.9999999);
  % dadt = smoothdata(dadt,'sgolay',2000);
 
 dadT = dadt/beta;
 
%  ac = cumtrapz(T,dadT);
 
%  plot(T,a,T,ac)

 
%  figure
%   plot(T,dadT,'LineWidth',1.5)
% hold on
dadT = lowpass(dadT,0.000000001,'Steepness',0.9999999);
% figure
% plot(T,dadT,'LineWidth',1.5)




% dadt = smoothdata(dadt,'sgolay',500);
% plot(T,dadt)
% hold off

data_out = [T t dadT a];


ddadT = diff_alpha(T,dadT);
%   ddadT = lowpass(ddadT,0.000000001,'Steepness',0.9999999);

  
  
  
% 
%   figure
  %hold on
%   plot(T,ddadT,'LineWidth',1.5)
  
%   ddadT = lowpass(ddadT,0.000000001,'Steepness',0.9999999);
%   figure
%   plot(T,ddadT,'LineWidth',1.5)
  
  
  
%   plot([T(1) T(end)],[0 0])
%   yyaxis right
%   plot(T,dadT)
%   
  
  
  
%     plot([T(1) T(end)],[0 0])
%   
  
  
% %   
%   figure
% yyaxis left
% plot(data_out(:,1),data_out(:,4))
% ylim([0,1])
% yyaxis right
% plot(data_out(:,1),data_out(:,3))
% 
% figure
% yyaxis left
% plot(data_out(:,1),data_out(:,3))
% ylim([min(data_out(:,3)),max(data_out(:,3))])
% %yyaxis right
% %plot(data_out(:,1),ddadT)
% grid
% 
% a_int = cumtrapz(t,dadt);
% figure
% plot(t,a,t,a_int,'--','Linewidth',1.5)
% legend('Original','Filtrado')
% figure
% plot(t,a-a_int)

