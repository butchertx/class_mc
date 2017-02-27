%plot the results listed in the result files here
%52 files: 13 different alpha values, 4 different delta values
for i = 1:52
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(6, 2);
	mag_err(i) = data(6, 3);
  
end

alpha = 2*[.3, .4, .45, .47, .48, .49, .50, .51, .52, .53, .55, .60, .70];
%split mags
mag_mat = [mags(1:13); mags(14:26); mags(27:39); mags(40:52)];
mag_err_mat = [mag_err(1:13); mag_err(14:26); mag_err(27:39); mag_err(40:52)];
alpha_mat = [alpha; alpha; alpha; alpha];


%make 1 plot with 6 different lines on it
errorbar(alpha_mat', mag_mat', mag_err_mat')
legend('delta = .01h', 'delta = .1h', 'delta = h', 'delta = 10h')