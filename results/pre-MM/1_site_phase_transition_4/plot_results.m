%plot the results listed in the result files here
%81 files: 9 different alpha values
n = 81;
for i = 1:n
	try
		data = csvread(strcat('results', num2str(i), '.csv'));
		mags(i) = 1 - data(8, 2) / (3*data(7, 2)*data(7,2));
		mag_err(i) = data(6, 3);
	catch
		mags(i) = 0;
		mag_err(i) = 0;
	end_try_catch
end

alpha = 2*[0.30, 0.40, 0.45, 0.48, 0.50, 0.52, 0.55, 0.60, 0.70];
%split mags
interval = 9;
ind = 1:interval:n
for i = 1:length(ind)
	mag_mat(i, :) = mags(ind(i):ind(i) + interval - 1);
	mag_err_mat(i, :) = mag_err(ind(i):ind(i) + interval - 1);
	alpha_mat(i, :) = alpha;
end


%mag_mat = [mags(1:17); mags(18:34); mags(35:51); mags(52:68)];
%mag_err_mat = [mag_err(1:17); mag_err(18:34); mag_err(35:51); mag_err(52:68)];
%alpha_mat = [alpha; alpha; alpha; alpha];
%
%
%%make 1 plot with 6 different lines on it
errorbar(alpha_mat', mag_mat', mag_err_mat')
xlabel('alpha')
ylabel('Binder Cumulant')
legend('', 'delta = .1, Nt = 1000', 'delta = .1, Nt = 100', 'delta = .1, Nt = 1250', 'delta = .1, Nt = 500', 'delta = 1, Nt = 1000', 'delta = 1, Nt = 100', 'delta = 1, Nt = 1250', 'delta = 1, Nt = 500', 'location', 'southeast')