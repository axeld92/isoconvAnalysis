function y = fslogn(T,p)
y = lognpdf(T,p(1),p(2)) + lognpdf(T,p(3),p(4)) + normpdf(T,p(5),p(6));
plot(T,y)
end
