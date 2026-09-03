%%
% % Calculate activation energy of individual steps
% 
aa = 0.1:0.01:0.9;
EA_v = {};%%
% % Calculate activation energy of individual steps
% 
aa = 0.1:0.01:0.9;
EA_v = {};
ea_max = {};
for j = 1:3
EA_v{j} = Starinkv(aa,out{1}{j},out{2}{j},out{3}{j},out{4}{j});
ea_max{j} = max(EA_v{j}(:,1));
ea_min{j} = min(EA_v{j}(:,1));
E{j} = mean(EA_v{j}(:,1));
variation{j} = 100*(ea_max{j}-ea_min{j})/E{j};
end% 
%  
%% 
 data1c = cleandata(data1);
 data2c = cleandata(data2);
 data3c = cleandata(data3);
 data4c = cleandata(data4);

EA_t = Starinkv(aa,data1c,data2c,data3c,data4c);
% figure
% hold on
% yyaxis left
% errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
% yyaxis right
% plot(data1c(:,4),data1c(:,3))
% hold off
%%
%figure
%hold on
figure
hold on
errorbar(aa,EA_v{1}(:,1),EA_v{1}(:,2),'o-')
errorbar(aa,EA_v{2}(:,1),EA_v{2}(:,2),'s-')
errorbar(aa,EA_v{3}(:,1),EA_v{3}(:,2),'^-')
errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
hold off
grid
%title('3')
xlabel('\alpha')
ylabel('E_{\alpha}')
legend('1','2','3','Total')
ea_max = {};
for j = 1:3
EA_v{j} = Starinkv(aa,out{1}{j},out{2}{j},out{3}{j},out{4}{j});
ea_max{j} = max(EA_v{j}(:,1));
ea_min{j} = min(EA_v{j}(:,1));
E{j} = mean(EA_v{j}(:,1));
variation{j} = 100*(ea_max{j}-ea_min{j})/E{j};
end% 

 data1c = cleandata(data1);
 data2c = cleandata(data2);
 data3c = cleandata(data3);
 data4c = cleandata(data4);

EA_t = Starinkv(aa,data1c,data2c,data3c,data4c);
% figure
% hold on
% yyaxis left
% errorbar(aa,EA_t(:,1),EA_t(:,2),'*-')
% yyaxis right
% plot(data1c(:,4),data1c(:,3))
% hold off
%%
%figure
%hold on
figure
hold on
plot(aa,EA_v{1}(:,1),'o-')
plot(aa,EA_v{2}(:,1),'s-')
plot(aa,EA_v{3}(:,1),'^-')
plot(aa,EA_t(:,1),'*-')
hold off
grid
%title('3')
xlabel('\alpha')
ylabel('E_{\alpha}')
legend('1','2','3','Total')