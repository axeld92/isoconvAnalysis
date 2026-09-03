pa = [0.006487 0.0682 547.8 55.1830 0.01162 -0.2543 611.6 35.0086 0.0012 0.1 649 160];
Temp = (450:850)';
fs_data = fsm(Temp,pa);
fs_int = cumtrapz(Temp,fs_data);


hfig = figure;
hold on

plot(Temp,fs_data,'-','Linewidth',1.5)
plot(Temp,fs(Temp,pa(1:4)),'--','Linewidth',1.5)
plot(Temp,fs(Temp,pa(5:8)),'--','Linewidth',1.5)
plot(Temp,fs(Temp,pa(9:12)),'--','Linewidth',1.5)
hold off
ylim([-0.5e-3 13e-3])
ylabel('d\alpha/dt')
xlabel('T')
legend('Señal DTG','Pico 1','Pico 2','Pico 3')
picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
%%
% figure
% plot(Temp,fs_int)
% ylim([-0.5e-1 1.1])
% figure
% plot(Temp,fs_data)
% ylim([-0.5e-3 13e-3])
% figure
% plot(0,0)
% ylim([-0.5e-1 1.1])
% xlim([450 850])

par = [(0.1:0.1:1)' (-0.5:0.1:0.5)' (-0.5:0.1:0.5)'  (0.1:0.1:1)'];











%Define Fraser Suzuki Function
function y = fs(T,p)

 in = 2 * p(2) * ((T-p(3)) / p(4));
 out = - log(2) / p(2)^2;
 
 y = p(1) * exp(out .* (log(1+in)).^2);
 y(y~=real(y))=0;
end

%Define Fraser suzuki mixture

function y = fsm(T,pm)

y = fs(T,pm(1:4)) + fs(T,pm(5:8)) + fs(T,pm(9:12));

end

function y = intfsm(T,pm)

y = cumtrapz(T,fsm(T,pm));

end