%%for T < Tc, we want the magnetization to follow the given curve

gamma = 1/2.269:0.001:1;
m_th = 1 - (sinh(2*gamma)).^-4;
m_th = m_th .^(.125);
figure()
hold on
plot(gamma, m_th)


%%now plot simulation results alongside
%mag is 5th row of each dump file, don't forget to measure the absolute magnetization
dumppath = 'C:/Users/Matthew/Dropbox/Code/class_mc/results/2D_Ising/dump';
g = [.45, .5, .55, .6, .65];
d = 100*exp(-2*g);
J = 100*g;
mags = zeros(1, 5);
for i = 1:5
	data = csvread(strcat(dumppath, num2str(i), '.csv'));
	mag = mean(abs(data(5, 2:end)));
	plot(g(i), mag, 'rx');
end