%plot the results listed in the result files here
n = 77;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag2(i) = data(7, 2);
	binders(i) = 1 - data(8, 2) / (3*data(7, 2)*data(7, 2));
	binder_err(i) = (1 - binders(i))*(data(8, 3)/data(8, 2) + 2*data(7, 3)/data(7, 2));
end

alpha1 = [.7, .75, .8, .85, .9, .95, 1.0, 1.05, 1.1, 1.15, 1.2];

%split mags
interval1 = 11;

for i = 1:7
	mag2_mat(i, :) = alpha1.*mag2(((i - 1)*interval1 + 1):(i*interval1));
	binder_mat1(i, :) = binders(((i - 1)*interval1 + 1):(i*interval1));
	binder_err_mat1(i, :) = binder_err(((i - 1)*interval1 + 1):(i*interval1));
	alpha_mat1(i, :) = alpha1;
end

%%make 2 plots
figure()
plot(alpha_mat1', mag2_mat')
xlabel('alpha')
ylabel('Scaled Order Parameter')
legend('Nt = 10000',  'Nt = 1000', 'Nt = 100', 'Nt = 2000', 'Nt = 200', 'Nt = 4000',  'Nt = 400')

figure()
plot(alpha_mat1', binder_mat1')
xlabel('alpha')
ylabel('Binder Cumulant')
legend('Nt = 10000',  'Nt = 1000', 'Nt = 100', 'Nt = 2000', 'Nt = 200', 'Nt = 4000',  'Nt = 400')