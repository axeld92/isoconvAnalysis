aa = 0.1:0.01:0.9;

EA_v = Vyazovkin(aa,St01,St06,St11,St16);
% EA_f = friedman(aa,St01,St06,St11,St16);
% EA_FWO = FWO(aa,St01,St06,St11,St16);
% EA_KAS = KAS(aa,St01,St06,St11,St16);
%EA_Star = Starink(aa,St01,St06,St11,St16);

figure
hold on
errorbar(aa,EA_v(:,1),EA_v(:,2),'-')
grid
title('Vyazovkin')
xlabel('Conversión')
ylabel('E_A')

% figure
% errorbar(aa,EA_f(:,1),EA_f(:,2),'x')
% grid
% title('Friedman')
% 
% figure
% errorbar(aa,EA_FWO(:,1),EA_FWO(:,2),'^');
% grid
% title('FWO')
% figure
% errorbar(aa,EA_KAS(:,1),EA_KAS(:,2),'o');
% grid
% title('KAS')
% figure
% errorbar(aa,EA_Star(:,1),EA_Star(:,2),'*');
% grid
% title('Starink')
% %hold off

% legend('Vyazovkin','FWO','KAS','Starink')
% grid
% title('E_{\alpha} vs \alpha')
% xlabel('\alpha')
% ylabel('E_{\alpha}')



% A = preexpcomp(St01,EA_v);
% figure
% plot(aa,A)

%  EA_v = Vyazovkin(aa,St02,St07,St12,St17);
% EA_f = friedman(aa,St02,St07,St12,St17);
% figure
% hold on
% errorbar(aa,EA_v(:,1),EA_v(:,2))
% errorbar(aa,EA_f(:,1),EA_f(:,2))
% hold off
% grid
% 
% EA_v = Vyazovkin(aa,St03,St08,St13,St18);
% EA_f = friedman(aa,St03,St08,St13,St18);
% figure
% hold on
% errorbar(aa,EA_v(:,1),EA_v(:,2))
% errorbar(aa,EA_f(:,1),EA_f(:,2))
% hold off
% grid
% 
%  EA_v = Vyazovkin(aa,St04,St09,St14,St19);
% EA_f = friedman(aa,St04,St09,St14,St19);
%  figure
% hold on
%  errorbar(aa,EA_v(:,1),EA_v(:,2))
% errorbar(aa,EA_f(:,1),EA_f(:,2))
% hold off
%  grid
% 
% EA_v = Vyazovkin(aa,St05,St10,St15,St20);
% figure
% hold on
% errorbar(aa,EA_v(:,1),EA_v(:,2),'-')
% grid
% title('Vyazovkin')
% xlabel('Conversión')
% ylabel('E_A')
% 
% % %EA_f = friedman(aa,St05,St10,St15,St20);
% % EA_FWO = FWO(aa,St05,St10,St15,St20);
% % EA_KAS = KAS(aa,St05,St10,St15,St20);
% % EA_Star = Starink(aa,St05,St10,St15,St20);
% % figure
% % hold on
% % 
% % errorbar(aa,EA_v(:,1),EA_v(:,2),'s')
% % %errorbar(aa,EA_f(:,1),EA_f(:,2),'x')
% % errorbar(aa,EA_FWO(:,1),EA_FWO(:,2),'^')
% % errorbar(aa,EA_KAS(:,1),EA_KAS(:,2),'o')
% % errorbar(aa,EA_Star(:,1),EA_Star(:,2),'o')
% % hold off
% %  grid
% %  legend('Vyazovkin','FWO','KAS','Starink')