%plot the results listed in the result files here
n = 220;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag2(i) = data(3, 2);
	%binders(i) = 1 - data(4, 2) / (3*data(3, 2)*data(3, 2));
	%binder_err(i) = (1 - binders(i))*(data(4, 3)/data(4, 2) + 2*data(3, 3)/data(3, 2));
end

alpha = [2.0, 2.5, 3.0, 5.0];
a = [10, 25, 50, 75, 100];
J = 0:.01:.1;

for j = 1:(length(alpha)*length(a))
	%split mags
	interval1 = 20;


	mag2_mat = mag2(j : interval1 : end);
%	binder_mat1(i, :) = binders((24*(j - 1) +(i - 1)*interval1 + 1):(24*(j - 1) +i*interval1));
%	binder_err_mat1(i, :) = binder_err((24*(j - 1) +(i - 1)*interval1 + 1):24*(j - 1) +(i*interval1));
	alpha_mat1 = J;

	%%%make 2 plots
	figure()
	plot(alpha_mat1', mag2_mat')
	xlabel('J')
	ylabel('<m^2>')
%
%	figure()
%	plot(alpha_mat1', binder_mat1')
%	xlabel('alpha')
%	ylabel('Binder Cumulant')
%	title(strcat('J', num2str(J(j))))
end