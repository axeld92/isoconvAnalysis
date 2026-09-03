function da = diff_alpha(T,a)

n = length(a);
da = zeros(n,1);
for i = 1:n
    if i == 1 || i == 2
        da(i) = (-a(i+2)+4*a(i+1)-3*a(i))/(T(i+2)-T(i));
    elseif i == n || i == n-1
        da(i) = (3*a(i)-4*a(i-1)+a(i-2))/(T(i)-T(i-2));
    else
        da(i) = (a(i+1)-a(i-1))/(T(i+1)-T(i-1));
%         da(i) = (  -a(i+2) + 8*a(i+1) - 8*a(i-1) + a(i-2)  ) / 6*(T(i+1) - T(i-1));
    
    end
end
% da(1) = 0;
% da(end) = 0;