%plot the results listed in the result files here
n = 330;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag2(i) = data(3, 2);
	binders(i) = 1 - data(4, 2) / (3*data(3, 2)*data(3, 2));
	binder_err(i) = (1 - binders(i))*(data(4, 3)/data(4, 2) + 2*data(3, 3)/data(3, 2));
end

alpha1 = [.2, .4, .6, .8, .85, .875, .9, .925, .95,  1.0, 1.2];
alpha2 = [.2, .4, .6, .675, .7, .725, .75, .775, .8, 1.0, 1.2];
alpha3 = [.2, .4, .6, .65, .675, .7, .725, .75, .8,  1.0, 1.2];
alpha_j = [alpha1;alpha2;alpha3;alpha3;alpha3];
J = [0.2, 0.4, 0.6, 0.8, 1.0];

for j = 1:length(J)
	%split mags
	interval1 = 11;

	for i = 1:6
		mag2_mat(i, :) = alpha_j(j, :).*mag2((66*(j - 1) + (i - 1)*interval1 + 1):(66*(j - 1) +i*interval1));
		binder_mat1(i, :) = binders((66*(j - 1) +(i - 1)*interval1 + 1):(66*(j - 1) +i*interval1));
		binder_err_mat1(i, :) = binder_err((66*(j - 1) +(i - 1)*interval1 + 1):66*(j - 1) +(i*interval1));
		alpha_mat1(i, :) = alpha_j(j, :);
	end

%	%%%make 2 plots
	figure()
	plot(alpha_mat1', mag2_mat')
	xlabel('alpha')
	ylabel('Scaled Order Parameter')
	title(strcat('J', num2str(J(j))))

	figure()
	plot(alpha_mat1', binder_mat1')
	xlabel('alpha')
	ylabel('Binder Cumulant')
	title(strcat('J', num2str(J(j))))

%	%%find the critical alpha
%	%%for each row, find the alpha where mag2_mat crosses 0.5
%	alpha_c = zeros(6, 1);
%	for i = 1:6
%		row = mag2_mat(i, :) - 0.5;
%		last_neg_index = 0;
%		for j = 1:length(row)
%			if (row(j) < 0)
%				last_neg_index = j;
%			end
%		end
%		slope = (row(last_neg_index + 1) - row(last_neg_index))/(alpha1(last_neg_index + 1) - alpha1(last_neg_index));
%		dy = - row(last_neg_index);
%		alpha_c(i) = alpha1(last_neg_index) + dy / slope;
%	end
%	N = [1000, 100, 2000, 200, 4000, 400];
%
%	figure()
%	plot(1./log(N), alpha_c, '*')
end




%%%plot some correlation functions
%data = csvread('results4.csv');
%corr = 0.5*(data(11,3:end) + fliplr(data(11, 3:end)));
%figure()
%plot(corr)
%figure()
%plot(data(11, 2:end))