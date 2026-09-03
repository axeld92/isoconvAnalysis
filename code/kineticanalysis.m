function holas = kineticanalysis(data)

ind = (data(:,4)>0.05).*(data(:,4)<0.9);
ind = logical(ind);
data = data(ind,:);


T = data(:,1);
t = data(:,2);
da = data(:,3);
a = data(:,4);
lda = log(da);
b = (T(end)-T(1))/(t(end)-t(1));
%b = 1;

R = 8.314E-3;
    
izq = lda + log(b) ;
M = [ones(size(da)),-1./(R*T), ones(size(da)) , log(a) , log(1-a)];
MM = M'*M;
izq = M'*izq;
holas = MM\izq;

logA = holas(1)
E = holas(2)
C = exp(holas(3))
n = holas(4)
m = holas(5)

logda = logA - log(b) - E./(R*T) + log(C) + n*log(a) + m*log(1-a);
dac = exp(logda);
plot(T,da,T,dac)


   
   
   