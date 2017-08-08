%plot the results listed in the result files here
n = 270;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag2(i) = data(3, 2);
	binders(i) = 1 - data(4, 2) / (3*data(3, 2)*data(3, 2));
	binder_err(i) = (1 - binders(i))*(data(4, 3)/data(4, 2) + 2*data(3, 3)/data(3, 2));
end

alpha = [0, .2, .4, .6, .8, 1.0];
J = [0, .5, 10, 1, 1.5, 2, 2.5, 3, 5];
delta = [.01, .03, .1, .3, 1];

for j = 1:length(J)
	%split mags
	interval1 = length(alpha);

	for i = 1:length(delta)
		mag2_mat(i, :) = mag2((j - 1)*30 + (i - 1)*interval1 + 1 : (j - 1)*30 + i*interval1);
%		binder_mat1(i, :) = binders((24*(j - 1) +(i - 1)*interval1 + 1):(24*(j - 1) +i*interval1));
%		binder_err_mat1(i, :) = binder_err((24*(j - 1) +(i - 1)*interval1 + 1):24*(j - 1) +(i*interval1));
		alpha_mat1(i, :) = alpha;
	end

	%%%make 2 plots
	figure()
	plot(alpha_mat1', mag2_mat')
	xlabel('alpha')
	ylabel('<m^2>')
	title(strcat('J', num2str(J(j))))
%
%	figure()
%	plot(alpha_mat1', binder_mat1')
%	xlabel('alpha')
%	ylabel('Binder Cumulant')
%	title(strcat('J', num2str(J(j))))
end

%
%%%find the critical alpha
%%%for each row, find the alpha where mag2_mat crosses 0.5
%alpha_c = zeros(7, 1);
%for i = 1:7
%	row = mag2_mat(i, :) - 0.5;
%	last_neg_index = 0;
%	for j = 1:length(row)
%		if (row(j) < 0)
%			last_neg_index = j;
%		end
%	end
%	slope = (row(last_neg_index + 1) - row(last_neg_index))/(alpha1(last_neg_index + 1) - alpha1(last_neg_index));
%	dy = - row(last_neg_index);
%	alpha_c(i) = alpha1(last_neg_index) + dy / slope;
%end
%N = [10000, 1000, 100, 2000, 200, 4000, 400];
%
%figure()
%plot(1./log(N), alpha_c, '*')

%%%plot some correlation functions
%data = csvread('results4.csv');
%corr = 0.5*(data(11,3:end) + fliplr(data(11, 3:end)));
%figure()
%plot(corr)
%figure()
%plot(data(11, 2:end))