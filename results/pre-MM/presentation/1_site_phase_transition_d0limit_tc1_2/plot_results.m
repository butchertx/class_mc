%plot the results listed in the result files here
n = 75;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag2(i) = data(3, 2);
	binders(i) = 1 - data(4, 2) / (3*data(3, 2)*data(3, 2));
	binder_err(i) = (1 - binders(i))*(data(4, 3)/data(4, 2) + 2*data(3, 3)/data(3, 2));
end

alpha_mat = [.5, .6, .7, .8, .9;.5, .6, .7, .8, .9;.950, .975, 1.0, 1.025, 1.05;1.075, 1.1, 1.125, 1.15, 1.175;1.125, 1.15, 1.175, 1.2, 1.225];
delta = [0.01, 0.1, 0.5, 0.76, 1.0];

%split mags
interval1 = 5;

for i = 1:15
	mag2_mat(i, :) = alpha_mat(ceil(i/3), :).*mag2(((i - 1)*interval1 + 1):(i*interval1));
	binder_mat1(i, :) = binders(((i - 1)*interval1 + 1):(i*interval1));
	binder_err_mat1(i, :) = binder_err(((i - 1)*interval1 + 1):(i*interval1));
	alpha_mat1(i, :) = alpha_mat(ceil(i/3), :);
end

%%%make 2 plots
figure()
plot(alpha_mat1', mag2_mat')
xlabel('alpha')
ylabel('Scaled Order Parameter')

figure()
plot(alpha_mat1', binder_mat1')
xlabel('alpha')
ylabel('Binder Cumulant')


%%find the critical alpha
%%for each row, find the alpha where mag2_mat crosses 0.5
alpha_c = zeros(15, 1);
for i = 1:15
	alpha1 = alpha_mat1(i, :);
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
N = [10000, 1000, 3000];

figure()
plot(1./log(N), alpha_c(13:15), '*', 1./log(N), alpha_c(10:12), '*', 1./log(N), alpha_c(7:9), '*', 1./log(N), alpha_c(4:6), '*', 1./log(N), alpha_c(1:3), '*')
legend('delta 1.0', 'delta .76', 'delta .5')
