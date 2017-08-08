n = 84;
N = [10000, 1000, 100, 2000, 200, 4000, 400];
clear sx
for num = 1:length(N)
	for i = 1:12
		data = csvread(strcat('corr', num2str((num-1)*12 + i), '.csv'));
		corr(:, i) = data(1:.5*length(data));
		sx(12*(num - 1) + i) = .5*(1 - corr(2, i));
	end
	%figure()
	tau = linspace(1/(N(num)), .5, length(corr) - 1);
	%plot(tau, corr(2:end, :))
	clear corr
end


corr1 = csvread('corr3.csv');
corr2 = csvread('corr7.csv');
corr3 = csvread('corr8.csv');
corr4 = csvread('corr10.csv');
corrs = [corr1;corr2;corr3;corr4];
corrs = corrs(:, 2:.5*10000);
plot(corrs', '-x', 'linewidth', 2)


%alpha1 = [1.05, 1.1, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.30, 1.325, 1.35, 1.4];
%interval1 = 12;
%for i = 1:7
%	sx_mat(i, :) = sx(((i - 1)*interval1 + 1):(i*interval1));
%end
%figure()
%hold on
%plot(alpha_mat1', sx_mat')
%plot(alpha_mat1', mag_mat')