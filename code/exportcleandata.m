function exportcleandata(varargin)


n = nargin;

%Extract data

for i = 1:n
    data{i} = varargin{i};
    out{i} = cleandata(data{i});
    out{i}(:,1) = out{i}(:,1);%-273.15;
end
k = 0;
names = {}
for i = 1:n
names{i} = char('Data') + string(i);
end

for i=1:n
writematrix(out{i},names{i},'Delimiter',';')
end
