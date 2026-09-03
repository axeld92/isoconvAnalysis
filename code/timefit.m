function  fitresult = timefit(alpha,time)
%Creates functions from data


%fit in temperature
[a,t] = prepareCurveData(alpha,time);
% a = alpha;
% t = time;

% Set up fittype and options.
%ft = fittype( 'smoothingspline' );
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline');
opts.SmoothingParam = 1;

% Fit model to data.
 [fitresult, gof] = fit(a,t,ft,opts);
%[fitresult, gof] = fit(a,t,ft);
 %fita = fitresult(Temp);
