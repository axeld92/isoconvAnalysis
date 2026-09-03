function g = integralmodel(data,E,A)

%Get temperatures
Temp = data(:,1)+273.15;

%Get times
time = data(:,2)-data(1,2);

%get massloss

massloss = data(:,4);

%calculate conversion
alpha = (massloss(1)-massloss)/(massloss(1)-massloss(end));

%clean data
[time,Temp,alpha] = removeerror(time,Temp,alpha);

%get rid of data outside of alpha = 0.05 to alpha = 0.095

ind = getind(alpha);
ind = logical(ind);
alpha = alpha(ind);
Temp = Temp(ind);
time = time(ind);

g = A.*innerint2*