%plot the results listed in the result files here
n = 100;
for i = 1:n
	data = csvread(strcat('results', num2str(i), '.csv'));
	mag2(i) = data(3, 2);
	binders(i) = 1 - data(4, 2) / (3*data(3, 2)*data(3, 2));
	binder_err(i) = (1 - binders(i))*(data(4, 3)/data(4, 2) + 2*data(3, 3)/data(3, 2));
end


alpha_J0 = [1.0, 1.33, 1.66, 2.0, 2.5];
alpha_J1 = [.5, .75, 1.0, 1.25, 2.0];

%first extract alpha = 0 data
%binders_a0 = [mag2(33:5:48),mag2(69:76), binders(77:5:92), mag2(97:100)];
%plot(alpha_J0, binders_a0(1:4:end), alpha_J0, binders_a0(2:4:end), alpha_J0, binders_a0(3:4:end))

%extract data for each distinct alpha and plot as a function of J
for i = 1:4
	temp_mag = [mag2(i:4:32), mag2((33 + i):5:52), mag2((52 + i):4:68), mag2((77 + i):5:96)];
	temp_binder = [binders(i:4:32), binders((33 + i):5:52), binders((52 + i):4:68), binders((77 + i):5:96)];
	figure()
	plot(alpha_J1, temp_mag(1:4:end), alpha_J1, temp_mag(2:4:end), alpha_J1, temp_mag(3:4:end), alpha_J1, temp_mag(4:4:end))
	figure()
	plot(alpha_J1, temp_binder(1:4:end), alpha_J1, temp_binder(2:4:end), alpha_J1, temp_binder(3:4:end), alpha_J1, temp_binder(4:4:end))
end

%%split mags
%interval1 = 12;
%
%for i = 1:7
%	mag2_mat(i, :) = alpha1.*mag2(((i - 1)*interval1 + 1):(i*interval1));
%	binder_mat1(i, :) = binders(((i - 1)*interval1 + 1):(i*interval1));
%	binder_err_mat1(i, :) = binder_err(((i - 1)*interval1 + 1):(i*interval1));
%	alpha_mat1(i, :) = alpha1;
%end
%
%figure()
%plot(alpha_mat1', binder_mat1')
%xlabel('J')
%ylabel('Binder Cumulant')
%
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