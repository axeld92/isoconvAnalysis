function [EA,cin] = friedmannea(aa,a1,a2,a3,a4,T1,T2,T3,T4)



% [a1,T1] = removeerror2(a1,T1);
% [a2,T2] = removeerror2(a2,T2);
% [a3,T3] = removeerror2(a3,T3);
% [a4,T4] = removeerror2(a4,T4);

 Temp1 = Tfit(T1,a1)';
 Temp2 = Tfit(T2',a2)';
 Temp3 = Tfit(T3',a3)';
 Temp4 = Tfit(T4',a4)';

dT1 = smoothdadT2(T1,diff_alpha(T1,a1))';
dT2 = smoothdadT2(T2,diff_alpha(T2,a2))';
dT3 = smoothdadT2(T3,diff_alpha(T3,a3))';
dT4 = smoothdadT2(T4,diff_alpha(T4,a4))';

dadT1 = smoothdadT(a1,dT1)';
dadT2 = smoothdadT(a2,dT2)';
dadT3 = smoothdadT(a3,dT3)';
dadT4 = smoothdadT(a4,dT4)';




figure
hold on

%aa = 0.1:0.1:0.9;
for i=1:length(aa)
    TT1 = Temp1(aa(i));
    TT2 = Temp2(aa(i));
    TT3 = Temp3(aa(i));
    TT4 = Temp4(aa(i));
    dadTT1 = dadT1(aa(i));
    dadTT2 = dadT2(aa(i));
    dadTT3 = dadT3(aa(i));
    dadTT4 = dadT4(aa(i));
    
    x = [1/TT1 ; 1/TT2 ; 1/TT3 ;1/TT4];
    y = [log(5*dadTT1); log(10*dadTT2) ; log(15*dadTT3) ; log(20*dadTT4)];
    [p,S]=polyfit(x,y,1);
    ci = polyparci(p,S);
    y1 = polyval(p,x);
    %hold on
    plot(x,y,'o')
    plot(x,y1)
   
    %legend()
    %hold off
    EA(i) = -p(1)*8.314/1000;
    cin(:,i) = -ci(:,1)*8.314/1000;
end    
hold off
grid