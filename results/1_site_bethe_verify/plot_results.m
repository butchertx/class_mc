%plot the results listed in the result files here
%126 files: 21 different alpha values, 6 different h values (listed by constant h first)
for i = 1:126
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(6, 2);
	mag_err(i) = data(6, 3);
end

alpha = 0:0.05:0.5;
%split mags
mag_mat = [mags(1:11); mags(22:32); mags(43:53); mags(64:74); mags(85:95); mags(106:116)];
mag_err_mat = [mag_err(1:11); mag_err(22:32); mag_err(43:53); mag_err(64:74); mag_err(85:95); mag_err(106:116)];
alpha_mat = [alpha; alpha; alpha; alpha; alpha; alpha];


%make 1 plot with 6 different lines on it
errorbar(alpha_mat', -mag_mat', mag_err_mat')
legend('h = 1e-10', 'h = 1e-8', 'h = 1e-6', 'h = 1e-4', 'h = 1e-3', 'h = 5e-3')