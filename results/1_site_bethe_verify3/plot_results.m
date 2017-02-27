%plot the results listed in the result files here
%126 files: 13 different alpha values, 3 different h values (listed by constant h first) and 2 temperatures
for i = 1:78
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(6, 2);
	mag_err(i) = data(6, 3);
  
end

alpha = 0:0.05:0.6;
%split mags
mag_mat = [mags(1:13); mags(14:26); mags(27:39); mags(40:52); mags(53:65); mags(66:78)];
mag_err_mat = [mag_err(1:13); mag_err(14:26); mag_err(27:39); mag_err(40:52); mag_err(53:65); mag_err(66:78)];
alpha_mat = [alpha; alpha; alpha; alpha; alpha; alpha];


%make 1 plot with 6 different lines on it
errorbar(2*alpha_mat', -mag_mat', mag_err_mat')
legend('h = .001, beta = 100', 'h = .001, beta = 10', 'h = .01, beta = 100', 'h = .01, beta = 10', 'h = .05, beta = 100', 'h = .05, beta = 10')

%plot a correlation function or two
data = csvread('results53.csv');
corr = data(11, 2:end);
plot(corr);