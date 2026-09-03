function EA = Starinkv(alphas,varargin)

lower = 0.01;
upper = 0.99;

n = length(varargin);
data = varargin;
 for i = 1:n
     idx = (data{i}(:,4)>lower).*(data{i}(:,4)<upper);
     idx = logical(idx);
     data{i} = data{i}(idx,:);
 end
Temp = {};
time = {};
alpha = {};
 
 
for i = 1:n
    Temp{i} = data{i}(:,1);
    time{i} = data{i}(:,2);
    alpha{i} = data{i}(:,4);
end

T = {};
t = {};
a = {};

for i = 1:n
    T{i} = Tfit(alpha{i},Temp{i});
    t{i} = timefit(alpha{i},time{i});
    beta{i} = (Temp{i}(end)-Temp{i}(1))/(time{i}(end)-time{i}(1));

end

aa = alphas;
hfig = figure;
hold on
for i = 1:length(alphas)
    aa = alphas(i);

TT = [T{1}(aa) T{2}(aa) T{3}(aa) T{4}(aa)];
%betai = [dadT1(TT(1)) dadT2(TT(2)) dadT3(TT(3)) dadT4(TT(4))];

x = [1/TT(1) ; 1/TT(2) ; 1/TT(3) ;1/TT(4)];
    y = [log(beta{1}/TT(1)^(1.92)); log(beta{2}/TT(2)^(1.92)) ; log(beta{3}/TT(3)^(1.92)) ; log(beta{4}/TT(4)^(1.92))];
    [p,S]=polyfit(x,y,1);
    ci = polyparci(p,S);
    y1 = polyval(p,x);
    ys(:,i) = y;
    y1s(:,i)= y1;
    xs(:,i) = x;
     plot(x,y,'x','LineWidth',1.5)
     plot(x,y1,'LineWidth',1.5)
    Ee(i) = -p(1)*8.314/1000/1.0008;
    CI(i) = abs(-ci(2,1)*8.314/1000/1.0008-Ee(i));
end
EA = [Ee' CI'];



hold off
%grid
xlabel('1/T_{\alpha,i}')
ylabel('ln(\beta_i/T_{\alpha,i}^{1.92})')

picturewidth = 20; % set this parameter and keep it forever
set(findall(hfig,'-property','FontSize'),'FontSize',13)
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])