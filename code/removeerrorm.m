function [Newt,NewT,Newm] = removeerrorm(Oldt,OldT,Oldm)

Newt(1) = Oldt(1);
NewT(1) = OldT(1);
Newm(1) = Oldm(1);
i=1;
j=1;

i=i+1;

while i<=length(Oldt)
    
    if Newm(j)>Oldm(i)
        j=j+1;
        Newt(j) = Oldt(i);
        NewT(j) = OldT(i);
        Newm(j) = Oldm(i);
        i=i+1;
    else i=i+1;
    end
end
Newt = Newt';
NewT = NewT';
Newm = Newm';