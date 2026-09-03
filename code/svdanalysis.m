 data1 = St01;
 data2 = St06;
 data3 = St11;
 data4 = St16;

 data1 = cleandata2(data1);
 data2 = cleandata2(data2);
 data3 = cleandata2(data3);
 data4 = cleandata2(data4);
 
 a1 = data1(:,3);
 a2 = data2(:,4);
 a3 = data3(:,4);
 a4 = data4(:,4);
 
 
 T1 = data1(:,1);
 T2 = data2(:,1);
 T3 = data3(:,1);
 T4 = data4(:,1);
 
 a1 = interp1(T1,a1,T4);
 a2 = interp1(T2,a2,T4);
 a3 = interp1(T3,a3,T4);
% 
   [a1,a2] = prepareCurveData(a1,a2);
 [a3,a4] = prepareCurveData(a3,a4);

 M = length(T4);
 
 h = 1;
 
%  a1(isnan(a1))=[];
 x = [a1 a2 a3 a4];
%  x = a1;
 
xavg = mean(x,2);

xx = x-xavg;
 
 c_xx = xx'*xx;
 c_xx = c_xx*h/M;
 
 
 [U,S,V] = svd(x,0);
%%
 
 %l = diag(S);
 
 %at=sqrt(M)*V*sqrt(diag(S));
 %%
 
 fx=(x*at)\(at'*at);
 
 % coeff = pca(x');
 