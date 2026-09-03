function sda = smoothdadT2(Temp,alpha)

%Halla la derivada de conversion en funcion a la temperatura y suaviza la
%curva para reducir el ruido.
%la salida es una estructura cfit

n = length(alpha);
dalpha = zeros(n,1);
for i = 1:n
    if i == 1
        dalpha(i) = (-alpha(i+2)+4*alpha(i+1)-3*alpha(i))/(Temp(i+2)-Temp(i));
    elseif i == n
        dalpha(i) = (3*alpha(i)-4*alpha(i-1)+alpha(i-2))/(Temp(i)-Temp(i-2));
    else
        dalpha(i) = (alpha(i+1)-alpha(i-1))/(Temp(i+1)-Temp(i-1));
    end
end

ind = getind(alpha);
ind = logical(ind);
alpha = alpha(ind);
dalpha = dalpha(ind);
Temp = Temp(ind);





%fit in temperature
[T, dalpha] = prepareCurveData( Temp, dalpha );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.001;

% Fit model to data.
[fitresult, gof] = fit( T, dalpha, ft, opts );
sda =  fitresult;