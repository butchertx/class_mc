%plot the results listed in the result files here
n = 84;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag(i) = data(2, 2);
	mag2(i) = data(3, 2);
	mag2_err(i) = data(3, 3);
	binders(i) = 1 - data(4, 2) / (3*data(3, 2)*data(3, 2));
	binder_err(i) = (1 - binders(i))*(data(4, 3)/data(4, 2) + 2*data(3, 3)/data(3, 2));
end

alpha1 = [1.05, 1.1, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.30, 1.325, 1.35, 1.4];

%split mags
interval1 = 12;

for i = 1:7
	mag_mat(i, :) = mag(((i - 1)*interval1 + 1):(i*interval1));
	mag2_mat(i, :) = alpha1.*mag2(((i - 1)*interval1 + 1):(i*interval1));
	mag2_err_mat(i, :) = alpha1.*mag2_err(((i - 1)*interval1 + 1):(i*interval1));
	binder_mat1(i, :) = binders(((i - 1)*interval1 + 1):(i*interval1));
	binder_err_mat1(i, :) = binder_err(((i - 1)*interval1 + 1):(i*interval1));
	alpha_mat1(i, :) = alpha1;
end
figure()
plot(alpha_mat1', mag_mat', '-*', 'linewidth', 4, 'markersize', 10)
figure()
plot(alpha_mat1', binder_mat1', '-*', 'linewidth', 4, 'markersize', 10)
%%%make 2 plots
figure()
plot(alpha_mat1', mag2_mat', '-*', 'linewidth', 4, 'markersize', 10)
%xlabel('alpha')
%ylabel('Scaled Order Parameter')
%legend('Nt = 100','Nt = 200','Nt = 400','Nt = 1000','Nt = 2000','Nt = 4000','Nt = 10000','location', 'southeast')
%
%figure()
%plot(alpha_mat1', binder_mat1')
%xlabel('alpha')
%ylabel('Binder Cumulant')
%legend('Nt = 10000',  'Nt = 1000', 'Nt = 100', 'Nt = 2000', 'Nt = 200', 'Nt = 4000',  'Nt = 400', 'location', 'southeast')


%%find the critical alpha
%%for each row, find the alpha where mag2_mat crosses 0.5
alpha_c = zeros(7, 1);
for i = 1:7
	row = mag2_mat(i, :) - 0.5;
	last_neg_index = 0;
	for j = 1:length(row)
		if (row(j) < 0)
			last_neg_index = j;
		end
	end
	slope = (row(last_neg_index + 1) - row(last_neg_index))/(alpha1(last_neg_index + 1) - alpha1(last_neg_index));
	dy = - row(last_neg_index);
	alpha_c(i) = alpha1(last_neg_index) + dy / slope;
end
N = [10000, 1000, 100, 2000, 200, 4000, 400];

%figure()
%plot(N.^-1, alpha_c, '*')

%%Try to do least squares fit to find alpha_c
%%K = 1.3, i = 9
psi_k_inf = zeros(1, length(alpha1));
for i = 1:length(alpha1)
	L = fliplr([10000, 4000, 2000, 1000, 400, 200, 100]);
	psi_L = fliplr([mag2_mat(1, i), mag2_mat(6, i), mag2_mat(4, i), mag2_mat(2, i), mag2_mat(7, i), mag2_mat(5, i), mag2_mat(3, i)]);
	f = @(x, p) p(1) * (1 + p(2)*x.^(2 - 2*p(1)/p(4)) + p(3)*x.^(4 - 4*p(1)/p(4)));
	[vals, new_p] = leasqr(L, psi_L, [1.0, 1.0, 1.0, .5], f, .0001, 200);
	new_p(4)
	psi_k_inf(i) = new_p(1);
end
figure()
hold on
plot(alpha_mat1', mag2_mat', '-*', 'linewidth', 4, 'markersize', 10)
plot(alpha1(7:end), psi_k_inf(7:end), 'linewidth', 4)
%legend('Nt = 100','Nt = 200','Nt = 400','Nt = 1000','Nt = 2000','Nt = 4000','Nt = 10000', 'Fitted', 'location', 'southeast')

