function [Newt,NewT,Newa] = removeerror(Oldt,OldT,Olda)

Newt(1) = Oldt(1);
NewT(1) = OldT(1);
Newa(1) = Olda(1);
i=1;
j=1;

i=i+1;

while i<=length(Oldt)
    
    if Newa(j)<Olda(i)
        j=j+1;
        Newt(j) = Oldt(i);
        NewT(j) = OldT(i);
        Newa(j) = Olda(i);
        i=i+1;
    else i=i+1;
    end
end
Newt = Newt';
NewT = NewT';
Newa = Newa';