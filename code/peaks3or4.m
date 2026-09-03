function f = peaks3or4(data)

data = cleandata(data);

Temp = data(:,1);
time = data(:,2);
da = data(:,3);
a = data(:,4);

figure
plot(Temp,da)
grid

numpeaks = input('how many peaks do you want to fit?: ')
p0 = input('input your best initial guesses based on the graph: ')


if numpeaks == 3
    p0 = [0.5   600     100     0.5     700     100     0.5     800     100];
    lb = [0     500     10      0       600     10      0       700     50];
    ub = [1     700     200     1       800     200     1       900     200];
elseif numpeaks == 4
    p0 = [0.5   600     100     0.5     700     100     0.5     800     100     0.5     800     100];
    lb = [0     500     10      0       600     10      0       700     50      0       700     50];
    ub = [1     700     200     1       800     200     1       900     200     1       900     200];   
end

objfun = @(p) sum((da - multiplepeaks(Temp,p,numpeaks)).^2);

bfit = patternsearch(objfun,p0,[],[],[],[],lb,ub);


f = bfit;
figure
plot(Temp,da,Temp,multiplepeaks(Temp,bfit,numpeaks))

function y = multiplepeaks(Temp,p,numpeaks)
    
    if numpeaks == 3
        y = p(1)*normpdf(Temp,p(2),p(3)) + ...
            p(4)*normpdf(Temp,p(5),p(6)) + ...
            p(7)*normpdf(Temp,p(8),p(9));
    elseif numpeaks == 4
        y = p(1)*normpdf(Temp,p(2),p(3)) + ...
            p(4)*normpdf(Temp,p(5),p(6)) + ...
            p(7)*normpdf(Temp,p(8),p(9)) + ...
            p(10)*normpdf(Temp,p(11),p(12));
    end
end
end
