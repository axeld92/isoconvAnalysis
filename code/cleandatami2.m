function data_out = cleandatami2(data)
%data = St01;

%remover datos anteriores a T = 110 °C
 ind = (data(:,1)>110).*(data(:,1)<795);
 ind = logical(ind);
% 
 data = data(ind,:);

%eliminar ultima fila de nan
if isnan(data(end,:))==1
    data(end,:) = [];
end

%plot(data(:,2),data(:,4))    
% % %calcular conversión
  %drift = -0.0010*data(:,2);%+data(1,4);
  %m = data(:,4)+drift;

 T = data(:,1) + 273.15;
 t = data(:,2) - data(1,2);
  
m = data(:,4)-data(1,4);
%a = (m(1)-m)/(m(1)-m(end));

%eliminar errores con remove error
% m = smoothdata(m,'sgolay',100);
% data(:,4) = m;

dmdt = diff_alpha(t,m);
dmdt = lowpass(dmdt,0.000000001,'Steepness',0.9999999);
data_out = [T t dmdt m]; 



% plot(t,T)
%  figure
% yyaxis left
%  plot(data_out(:,2),data_out(:,4))
%  yyaxis right
%  plot(data_out(:,2),data_out(:,3))


