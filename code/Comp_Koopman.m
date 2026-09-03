function [Comp residu]=Comp_Koopman(K,n);
    %return the companion matrix

    Km=K(:,1:n);

    %Proejction coefficient
    c=(Km'*Km)\Km'*K(:,n+1);


    %companion
    Comp=zeros(n,n);
    Comp(:,n)=c;
    for i=2:n,
     Comp(i,i-1)=1;
    end;

    cbis=Km*Comp;
    cbis2=cbis(:,end);
    residu=[K(:,n+1) cbis2];