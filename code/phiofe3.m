function phi = phiofe3(t,T,E)

n = 3;
for i = 1:n
TempInt(:,i) = innerint2(t(:,i),T(:,i),E);
end

for ii = 1:n
    for jj = 1:n
        if ii~=jj
            
            phi_num = TempInt(:,ii);
            phi_den = TempInt(:,jj);
            a = phi_num;
            b = phi_den;
            myphi(ii,jj) = a/b;
        end
    end
end
phi = sum(myphi(:));