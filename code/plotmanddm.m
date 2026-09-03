function plotmanddm(data)
%data = St01;

Temp = data(:,1);
ind = Temp>110;

time = data(:,2);
time = time(ind);
Temp = data(:,1);
Temp = Temp(ind);
mass = data(:,4);
mass = mass(ind);
a = (mass(1)-mass)/(mass(1)-mass(end-1));
%yyaxis left
plot(Temp, a)
