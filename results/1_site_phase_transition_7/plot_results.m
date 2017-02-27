%plot the results listed in the result files here
n = 84;
for i = 1:n
	try
		data = csvread(strcat('results', num2str(i), '.csv'));
		mags(i) = 1 - data(8, 2) / (3*data(7, 2)*data(7, 2));
		mag_err(i) = (1 - mags(i))*(data(8, 3)/data(8, 2) + 2*data(7, 3)/data(7, 2));
	catch
		mags(i) = 0;
		mag_err(i) = 0;
	end_try_catch
end

alpha2 = 2*[.45, .475, .5, .525, .55, .575, .6];
alpha3 = 2*[.55, .575, .6, .625, .65, .675, .7];
alpha1 = 2*[.35, .375, .4, .425, .45, .475, .5];
%split mags
interval1 = 7;
interval2 = 7;
interval3 = 7;
ind1 = 1:interval1:28;
ind2 = 29:interval2:56;
ind3 = 57:interval3:84;
for i = 1:length(ind1)
	mag_mat1(i, :) = mags(ind1(i):ind1(i) + interval1 - 1);
	mag_err_mat1(i, :) = mag_err(ind1(i):ind1(i) + interval1 - 1);
	alpha_mat1(i, :) = alpha1;
end

for i = 1:length(ind2)
	mag_mat2(i, :) = mags(ind2(i):ind2(i) + interval2 - 1);
	mag_err_mat2(i, :) = mag_err(ind2(i):ind2(i) + interval2 - 1);
	alpha_mat2(i, :) = alpha2;
end

for i = 1:length(ind3)
	mag_mat3(i, :) = mags(ind3(i):ind3(i) + interval3 - 1);
	mag_err_mat3(i, :) = mag_err(ind3(i):ind3(i) + interval3 - 1);
	alpha_mat3(i, :) = alpha3;
end


%mag_mat = [mags(1:17); mags(18:34); mags(35:51); mags(52:68)];
%mag_err_mat = [mag_err(1:17); mag_err(18:34); mag_err(35:51); mag_err(52:68)];
%alpha_mat = [alpha; alpha; alpha; alpha];
%
%
%%make 3 plots
figure()
errorbar(alpha_mat1', mag_mat1', mag_err_mat1')
xlabel('alpha')
ylabel('Binder Cumulant')
legend('delta = .1, Nt = 1000', 'delta = .1, Nt = 100', 'delta = .1, Nt = 1250', 'delta = .1, Nt = 500', 'location', 'southeast')

figure()
errorbar(alpha_mat2', mag_mat2', mag_err_mat2')
xlabel('alpha')
ylabel('Binder Cumulant')
legend('delta = .5, Nt = 1000', 'delta = .5, Nt = 100', 'delta = .5, Nt = 1250', 'delta = .5, Nt = 500', 'location', 'southeast')

figure()
errorbar(alpha_mat3', mag_mat3', mag_err_mat3')
xlabel('alpha')
ylabel('Binder Cumulant')
legend('delta = 1, Nt = 1000', 'delta = 1, Nt = 100', 'delta = 1, Nt = 1250', 'delta = 1, Nt = 500', 'location', 'southeast')
