

[t5_0,T5_0,a5_0] = removeerror(t5_0,T5_0,a5_0);
[t10_0,T10_0,a10_0] = removeerror(t10_0,T10_0,a10_0);
[t15_0,T15_0,a15_0] = removeerror(t15_0,T15_0,a15_0);
[t20_0,T20_0,a20_0] = removeerror(t20_0,T20_0,a20_0);


alpha = 0.2:0.01:0.8;

%Activation energy for 0% WT

eaf_0 = friedmannea(alpha,a5_0,a10_0,a15_0,a20_0,T5_0,T10_0,T15_0,T20_0);

% for i=1:length(alpha)
%     eea_0(i,:) = aicvy2(alpha(i),t5_0,t10_0,t15_0,t20_0,T5_0,T10_0,T15_0,T20_0,a5_0,a10_0,a15_0,a20_0);
% end
%  ea_0= eea_0(:,1);
% error_0= eea_0(:,2);
% figure
% hold on
% plot(alpha,ea_0,'s');
% plot(alpha,eaf_0,'^');
% hold off
% grid
% legend('Vyazovkin','Friedman')
% title('E_{\alpha} vs \alpha for  0% WT')
% ylabel('E_{\alpha}')
% xlabel('\alpha')
% figure
% plot(alpha,error_0)
% title('error % 0')
% grid
