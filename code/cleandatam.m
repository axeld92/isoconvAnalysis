function data_out = cleandatam(data)
%data = St01;


%remover datos anteriores a T = 110 °C
 ind = (data(:,1)>50).*(data(:,1)<799);
 ind = logical(ind);
% 
 data = data(ind,:);

%eliminar ultima fila de nan
if isnan(data(end,:))==1
    data(end,:) = [];
end

  
m = data(:,3);
m = m-m(1);
T = data(:,1)+273.15;
%a = (m(1)-m)/(m(1)-m(end));


%
m = smoothdata(m,'sgolay',100);
data(:,4) = m;

 dmdt = diff_alpha(T,m);
 dmdt = lowpass(dmdt,0.000000001,'Steepness',0.9999999);
 
%  ipt = findchangepts(m)
% 
% 
% T_onset = T(ipt)
% %  
% 
% [DTGmax,idmax] = max(-dmdt)
% T_dtgmax = T(idmax)
% 
% 
% idx10 = find(-m>10,'first')

data(:,3) = dmdt;
 %[data_out(:,2),data_out(:,1),data_out(:,4)] = removeerrorm(data(:,2),data(:,1),data(:,4));
data_out = data;
%calcular derivada da/dt

 %dmdt = diff_alpha(data_out(:,1),data_out(:,4));

 
 
%   
%  figure
% yyaxis left
%   plot(data_out(:,1),data_out(:,4))
% ylabel('Pérdida de masa [%]') 
%   yyaxis right
%  plot(data_out(:,1),data_out(:,3),'--')
%  ylabel('Velocidad de pérdida de masa [%/min]')
%  legend('TGA','DTG')
%  xlabel('Temperatura [°C]')
%  grid
%  


