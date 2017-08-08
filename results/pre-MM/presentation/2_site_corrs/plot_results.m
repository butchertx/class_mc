%n = 10;
%for i = 1:n
%	data = csvread(strcat('corr', num2str(i), '.csv'));
%	figure()
%	plot(data')
%end


corr1 = csvread('corr3.csv');
corr2 = csvread('corr4.csv');
corr3 = csvread('corr9.csv');
corr4 = csvread('corr10.csv');

figure()
hold on
plot(corr1', '-*', 'linewidth', 3)
plot(corr3', '-*', 'linewidth', 3)
hold off
figure()
hold on
plot( corr2', '-*', 'linewidth', 3)
plot( corr4', '-*', 'linewidth', 3)
