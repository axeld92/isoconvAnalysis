
%parameters to create data
pa = [0.0363    0.0682  560.0000-273   55.1830];
%independent variable vector
Temp = (220:0.1:350)';
%create data
fs_data = frasuz(Temp,pa);


Tempr = Temp;
Temp = (Tempr-Tempr(1))/(Tempr(end)-Tempr(1));


%[h s p w]
%[h p 1/w s]
pa2 = [pa(1) pa(3) 1/pa(4) pa(2)];

%  a0 = [0.03 555 0.017 55]';
rng default
a0 = pa2' + rand(size(pa2'));


hola = fs(Temp,pa2);
% figure
% plot(Temp,hola)
% hold on




% respuesta = minsolving(Temp,pa2,fs_data)

% aaaa = fs(Temp,a0);
data = fs_data;
a = a0
d = 1;
k=0;
% while d > 0.18315
%     
%     A = minsolving(Temp,a,data);
%     an = a+A;
%     d = abs(sum((an-a)./a));
%     a = an;
%     k=k+1;
% end




figure 
hold on
for i = 1:6
    
    plot(Temp,fs(Temp,a))
    A = minsolving(Temp,a,data);
    an = a+A;
    d(i) = abs(sum((an-a)./a))
    a = an
end

plot(Temp,fs(Temp,a))
hold off


  function A = minsolving(Temp,a,data)


P = zeros(4,4);

for l = 1:4
    for m = 1:4
        P(l,m) = sum(dfa(Temp,a,l).*dfa(Temp,a,m));
    end
end

P
detP = det(P)

Z = zeros(4,1);

for m = 1:4
    Z(m,1) = sum(dfa(Temp,a,m).*(data-fs(Temp,a)));
end
Z
 A = P\Z;
  end


function f = fs(x,a)

f = a(1)*exp(-log(2)*(2*a(3)*(x-a(2))).^2).*exp(log(2)*(2*a(3)*(x-a(2))).^3*a(4)).*exp(-log(2)*(11/12)*(2*a(3)*(x-a(2))).^4*a(4)^2).*exp(-(1/9)*(2*a(3)*(x-a(2))).^6*a(4)^4);

f(isnan(f)==1)=0;
end

function f = dfa(x,a,p)
if p==1
f = fs(x,a)./a(1);

elseif p==2
f = log(2)*(128*a(3)^6*a(4)^4*(x-a(2)).^5/3 - 160*a(3)^5*a(4)^3*(x-a(2)).^4/3 + ...
    176*a(3)^4*a(4)^2*(x-a(2)).^3/3 - 24*a(3)^3*a(4)*(x-a(2)).^2 + 8*a(3)^2*(x-a(2))).*fs(x,a);

elseif p==3
f = log(2)*(-128*a(3)^5*a(4)^4*(x-a(2)).^6/3 + 160*a(3)^4*a(4)^3*(x-a(2)).^5/3 - ...
    176*a(3)^3*a(4)^2*(x-a(2)).^4/3 + 24*a(3)^2*a(4)*(x-a(2)).^3 - 8*a(3)*(x-a(2)).^2).*fs(x,a);

elseif p==4
f = log(2)*(-256*a(3)^6*a(4)^3*(x-a(2)).^6/9 + 32*a(3)^5*a(4)^2*(x-a(2)).^5 ...
- 176*a(3)^4*a(4)*(x-a(2)).^4/9 + 8*a(3)^2*(x-a(2)).^2).*fs(x,a);
end

f(isnan(f)==1)=0;

end

function y = frasuz(T,p)

 in = 2 * p(2) * ((T-p(3)) / p(4));
 out = - log(2) / p(2)^2;
 
 y = p(1) * exp(out .* (log(1+in)).^2);
 y(y~=real(y))=0;
end


% % 
% % 
% % function f = fs(x,a)
% % 
% % f = a(1)*exp(-log(2)*(2*a(3)*(x-a(2))).^2).*exp(log(2)*(2*a(3)*(x-a(2))).^3*a(4)).*exp(-log(2)*(11/12)*(2*a(3)*(x-a(2))).^4*a(4)^2).*exp(-(1/9)*(2*a(3)*(x-a(2))).^6*a(4)^4);
% % 
% % end
% % 
% % function f = fa1(x,a)
% % 
% % f = fs(x,a)./a(1);
% % 
% % end
% % 
% % function f = fa2(x,a)
% % 
% % f = log(2)*(128*a(3)^6*a(4)^4*(x-a(2)).^5/3 - 160*a(3)^5*a(4)^3*(x-a(2)).^4/3 + ...
% %     176*a(3)^4*a(2)^2*(x-a(2)).^3/3 - 24*a(3)^3*a(4)*(x-a(2)).^2 + 8*a(3)^2*(x-a(2))).*fs(x,a);
% % end
% % 
% % function f = fa3()
% % 
% % f = log(2)*(-128*a(3)^5*a(4)^4*(x-a(2)).^6/3 + 160*a(3)^4*a(4)^3*(x-a(2)).^5/3 - ...
% %     176*a(3)^3*a(2)^2*(x-a(2)).^4/3 + 24*a(3)^2*a(4)*(x-a(2)).^3 - 8*a(3)*(x-a(2))).*fs(x,a);
% % end
% % 
% % function f = fa4(x,a)
% % 
% % f = log(2)*(-256*a(3)^6*a(4)^3*(x-a(2)).^6/9 + 32*a(3)^5*a(4)^2*(x-a(2)).^5 ...
% % - 176*a(3)^4*a(4)*(x-a(2)).^4/9 + 8*a(3)^2*(x-a(2)).^2).*fs(x,a);
% % end