data = csvread('dump/dump4.csv', 1, 1);
plot(autocorrelation(data(1, :)))

data = csvread('dump/correlation4.csv');
lx = data(1,1);
ly = data(1, 2);
figure()
hold on
plot(data(2, 1:ly))
plot(data(2, (ly + 1):2*ly))

figure()
plot(integrate_corr(data(2:end), lx, ly))