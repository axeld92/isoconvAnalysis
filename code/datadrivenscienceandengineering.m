%%
clear all;clc
x = 3;
a = [-2:.1:2]'; 
b = a*x + (1*randn(size(a)));

rng default
plot(a,x*a,'k','LineWidth',2);
hold on
plot(a,b,'rx','LineWidth',2);
hold off

[U,S,V] = svd(a,'econ');
xtilde = V*inv(S)*U'*b;
hold on
plot(a,xtilde*a,'b--','LineWidth',2)
hold off

%%
clear all;close all; clc
load hald;
A = ingredients;
b = heat ; 

[U,S,V] = svd(A,'econ');
x = V*inv(S)*U'*b;


plot(b,'k','LineWidth',2); hold on
plot(A*x,'r-o','LineWidth',1.,'Markersize',2);

%%
clear all;close all; clc
load housing
b = housing(:,14); % housing values in $1000s
A = housing(:,1:13); % other factors,
A = [A ones(size(A,1),1)]; % Pad with ones y-intercept
x = regress(b,A);
plot(b,'k-o');
hold on, plot(A*x,'r-o');
[b sortind] = sort(housing(:,14)); % sorted values
plot(b,'k-o')
hold on, plot(A(sortind,:)*x,'r-o')

%%

%PCA

xC = [2;1];
sig = [2;.5];

theta = pi/3;

R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
nP = 10000;

X = R*diag(sig)*randn(2,nP) + diag*xC*ones(2,nP);

X_avg = mean(X,2);
B = X-X_avg*ones(1,nP);
[U,S,V] = svd(B/sqrt(nP),'econ');


%%

clear all, close all, clc
n = 2^10;
w = exp(-i*2*pi/n);
% Slow
tic
DFT = zeros(n);
for i=1:n
for j=1:n
DFT(i,j) = w^((i-1)*(j-1));
end
end
toc
% Fast
tic
[I,J] = meshgrid(1:n,1:n);
DFT = w.^((I-1).*(J-1));
toc
imagesc(real(DFT))
