function Eea = aicvy2(aa,t1,t2,t3,t4,T1,T2,T3,T4,a1,a2,a3,a4)

% [t1,a1,T1] = removeerror(t1,a1,T1);
% [t2,a2,T2] = removeerror(t2,a2,T2);
% [t3,a3,T3] = removeerror(t3,a3,T3);
% [t4,a4,T4] = removeerror(t4,a4,T4);



% %calcular las tasas de calentamiento
% b1 = betaprom(t1,T1);
% b2 = betaprom(t2,T2);
% b3 = betaprom(t3,T3);
% b4 = betaprom(t4,T4);
% bet = [b1 b2 b3 b4];
%crear objetos cfit para T=T(alpha)

TT1 = Tfit(T1,a1);
TT2 = Tfit(T2,a2);
TT3 = Tfit(T3,a3);
TT4 = Tfit(T4,a4);

tt1 = timefit(t1,a1);
tt2 = timefit(t2,a2);
tt3 = timefit(t3,a3);
tt4 = timefit(t4,a4);

delta_a = 0.1;
z=10;
tt = [tt1(aa-delta_a:delta_a/z:aa+delta_a),tt2(aa-delta_a:delta_a/z:aa+delta_a),tt3(aa-delta_a:delta_a/z:aa+delta_a),tt4(aa-delta_a:delta_a/z:aa+delta_a)];
TT = [TT1(aa-delta_a:delta_a/z:aa+delta_a),TT2(aa-delta_a:delta_a/z:aa+delta_a),TT3(aa-delta_a:delta_a/z:aa+delta_a),TT4(aa-delta_a:delta_a/z:aa+delta_a)];
%asumir valor de E

E = 250; %kJ/kmol
R = 8.31446261815324;

%calcular phi

phi =@(e) phiofe2(tt,TT,e);
  options = optimset('MaxFunEvals',100000,'TolFun',1e-10,'TolX',1e-8);
[Ea,fval] = fminsearch(phi,E,options);

errorp = abs(fval-12)*100/12;
Eea = [Ea, errorp];












