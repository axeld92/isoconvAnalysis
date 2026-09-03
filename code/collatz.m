function collatz(n)
i=1;
while n>1
    x(i) = n;
    if mod(n,2) == 0
        n = n/2;
    else
        n = 3*n+1;
    end
    i=i+1;
end
  plot(x)   
  N = length(x) + 1