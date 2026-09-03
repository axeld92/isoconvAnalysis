data1 = cleandata(St05);
data2 = cleandata(St10);
data3 = cleandata(St15);
data4 = cleandata(St20);

t1 = data1(:,2);
T1 = data1(:,1);
a1 = data1(:,4);
da1 = data1(:,3);

t2 = data2(:,2);
T2 = data2(:,1);
a2 = data2(:,4);

t3 = data3(:,2);
T3 = data3(:,1);
a3 = data3(:,4);

t4 = data4(:,2);
T4 = data4(:,1);
a4 = data4(:,4);

Tmin = min([T1;T2;T3;T4]);
Tmax = max([T1;T2;T3;T4]);

tmax = max([t1;t2;t3;t4]);

T1 = (T1-Tmin)/(Tmax-Tmin);
T2 = (T2-Tmin)/(Tmax-Tmin);
T3 = (T3-Tmin)/(Tmax-Tmin);
T4 = (T4-Tmin)/(Tmax-Tmin);

t1 = t1/tmax;
t2 = t2/tmax;
t3 = t3/tmax;
t4 = t4/tmax;

% figure
% hold on
% plot(T1,a1)
% plot(T2,a2)
% plot(T3,a3)
% plot(T4,a4)
% hold off

% figure
% 
% plot3(t1,T1,a1,t2,T2,a2,t3,T3,a3,t4,T4,a4)
% 
% xlabel('t')
% ylabel('T')
% zlabel('a')

% X = [a1 interp1(T2,a2,T1,'linear','extrap') interp1(T3,a3,T1,'linear','extrap') interp1(T4,a4,T1,'linear','extrap')];
%X = [a1 interp1(t2,a2,t1,'linear','extrap') interp1(t3,a3,t1,'linear','extrap') interp1(t4,a4,t1,'linear','extrap')];

%[U,S,V] = svd(X);

lT1 = log(T1);
plot(lT1,da1)
figure
plot(T1,da1)