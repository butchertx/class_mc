%n = 21;
%figure()
%hold on
%for i = 1:n
%	data = csvread(strcat('corr', num2str(i), '.csv'));
%	plot(data(:, 1:.5*length(data))', 'linewidth', 4)
%end


corr1 = csvread('corr1.csv');
corr2 = csvread('corr2.csv');
corr3 = csvread('corr15.csv');
corr4 = csvread('corr13.csv');
corr5 = csvread('corr14.csv');
corr6 = csvread('corr21.csv');

figure()
hold on
plot(corr1', '-*', 'linewidth', 3)
plot(corr2', '-*', 'linewidth', 3)
plot(corr3', '-*', 'linewidth', 3)
hold off
figure()
hold on
plot( corr4', '-*', 'linewidth', 3)
plot( corr5', '-*', 'linewidth', 3)
plot( corr6', '-*', 'linewidth', 3)
