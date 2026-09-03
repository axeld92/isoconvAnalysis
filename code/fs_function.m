function answer = fs_function(temp,parameters)
%Calcula funcion de Fraser Suzuki
height = parameters(1);
skew = parameters(2);
position = parameters(3);
width = parameters(4);



interior = 2 * skew * ((temp-position) / width);
exterior = - log(2) / skew^2;

l = length(interior);
for i = 1:l
    if skew ~= 0
        if interior(i)>-1
            answer(i) = height*exp(-log(2)/skew^2*(log(1+2*skew*(temp(i)-position)/width))^2);
        else answer(i) = 0;
        end
    else answer(i) = height*exp(-4*log(2)*((temp(i)-position)/width)^2);
    end
end

%answer = height * exp(exterior .* (log(1+interior)).^2);
%answer(answer~=real(answer))=0;