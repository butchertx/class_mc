%plot the results listed in the result files here
%68 files: 17 different alpha values, 4 different delta values
for i = 1:68
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(6, 2);%1 - data(8, 2) / (3*data(7, 2)*data(7,2));
	mag_err(i) = data(6, 3);
  
end

alpha = 2*[.05, .1, .14, .18, .2, .22, .23, .24, .25, .26, .27, .28, .30, .32, .36, .40, .45];
%split mags
mag_mat = [mags(1:17); mags(18:34); mags(35:51); mags(52:68)];
mag_err_mat = [mag_err(1:17); mag_err(18:34); mag_err(35:51); mag_err(52:68)];
alpha_mat = [alpha; alpha; alpha; alpha];


%make 1 plot with 6 different lines on it
errorbar(alpha_mat', mag_mat', mag_err_mat')
legend('delta = .001', 'delta = .01', 'delta = .1', 'delta = 1')