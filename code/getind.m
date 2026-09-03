function ind = getind(a)

n = length(a);
ind = ones(size(a));
for i = 1:n
    if a(i)<0.01 || a(i)>0.99
        ind(i) = 0;
    end
end

