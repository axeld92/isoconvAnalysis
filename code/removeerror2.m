function [aa,TT] = removeerror2(a,T)

%Newt(1) = Oldt(1);
TT(1) = T(1);
aa(1) = a(1);
i=1;
j=1;

i=i+1;

while i<=length(a)
    
    if aa(j)<a(i)
        j=j+1;
        %Newt(j) = Oldt(i);
        TT(j) = T(i);
        aa(j) = a(i);
        i=i+1;
    else i=i+1;
    end
end
        