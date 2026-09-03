%%
R = 8.314E-3;

c = mean([C]);
% 

datareal = cleandatami(b75t550);

t0 = datareal(1,2);
tf = datareal(end,2);
T0 = datareal(1,1);
Tf = datareal(end,1);
b = 20;

tp = (Tf-T0)/b;
%         
dadt = @(t,aT)    [ A_1*exp(-E{1}/R/aT(5))*real(aT(1)^nmp_1(1)*(1-aT(1))^nmp_1(2)) ;          %1
                    A_2*exp(-E{2}/R/aT(5))*real(aT(2)^nmp_2(1)*(1-aT(2))^nmp_2(2)) ;          %2
                    A_3*exp(-E{3}/R/aT(5))*real(aT(3)^nmp_3(1)*(1-aT(3))^nmp_3(2)) ;          %3
                    A_4*exp(-E{4}/R/aT(5))*real(aT(4)^nmp_4(1)*(1-aT(4))^nmp_4(2)) ;          %4
                    b - b*heaviside(t-tp)];                
                    

  options = odeset('RelTol',1e-8,'AbsTol',1e-8) ;
%   [t,aT] = ode15s(dadt,[0 tf],[1E-3 1E-3 1E-3 T0],options);
 [t,aT] = ode15s(dadt,[0 tf],[1E-3 1E-3 1E-3 1E-3 T0],options);
%  aT(aT~=real(aT)) = 0;

% figure
% plot(t,aT(:,5),datareal(:,2),datareal(:,1))
dreal = interp1(datareal(:,2),datareal(:,4),t);
 %a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3);
a_calc = c(1)*aT(:,1) + c(2)*aT(:,2) + c(3)*aT(:,3) + c(4)*(aT(:,4));
%hold on
color = '#EDB120';
%plot(t,a_calc,'LineWidth',1.5,'Color',color)
%plot(datareal(:,2),datareal(:,4),'--','LineWidth',1.5,'Color',color)                
                
                
for i = 1:10:length(t)
    t_i = t(1:i);
    x_i = dreal(1:i);
    y_i = a_calc(1:i);
    
    plot(t_i,x_i,'-',t_i,y_i,'--','LineWidth',1.5);
    ylim([0 1]);
    xlim([0,max(t)])
    legend('Experimental','Modelo','Location','east')
       %hold on
    %plot(t,xx,t,yy);
    
    %grid on
    
    drawnow
    %pause(0.2)
    movieVector(i) = getframe;
end


%% save the movie

moviewr = VideoWriter('movietga','MPEG-4');
moviewr.FrameRate = 20;

open(moviewr);
writeVideo(moviewr, movieVector);
close(moviewr);

