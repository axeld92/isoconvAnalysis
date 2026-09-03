function  fitresult = Tfit(alpha,Temp)
%Creates functions from data


%fit in temperature
[a,T] = prepareCurveData(alpha,Temp);

% if sum(isnan(a))>0
%     disp('hule a')
% elseif sum(isnan(T))>0
%     disp('hule T')
% else disp('No sé por qué pero hule')
% end
    

% Set up fittype and options.
 ft = fittype( 'smoothingspline' );
 opts = fitoptions( 'Method', 'SmoothingSpline' );
 opts.SmoothingParam = 0.999999;

%  Fit model to data.
 fitresult = fit(a,T,ft,opts);

