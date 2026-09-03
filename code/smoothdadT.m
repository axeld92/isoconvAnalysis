function sda = smoothdadT(alpha,dalpha)
%fit in temperature
[a, da] = prepareCurveData( alpha, dalpha );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 1;

% Fit model to data.
[sda, gof] = fit( a, da, ft, opts );

