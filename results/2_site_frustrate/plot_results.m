%plot the results listed in the result files here
n = 126;
for i = 1:n
	if exist(strcat('results', num2str(i), '.csv'))
		data = csvread(strcat('results', num2str(i), '.csv'));
		mag2(i) = data(3, 2);
		cluster(i) = data(5, 2);
		a_mag(i) = data(6, 2);
		a_mag2(i) = data(7, 2);
		sx(i) = data(8, 2);
		loc(i) = data(9, 2);
	else
		mag2(i) = 0;
		cluster(i) = 0;
		a_mag(i) = 0;
		a_mag2(i) = 0;
		sx(i) = 0;
		loc(i) = 0;
	end
end

alpha = 0:0.2:4;
Nt = [10000, 1000, 100, 3000, 300, 50000];
beta = 10;

%split mags
interval = length(alpha);

for i = 1:length(Nt)
	mag2_mat(i, :) = mag2(((i - 1)*interval + 1):(i*interval));
	cluster_mat(i, :) = 0.5*cluster(((i - 1)*interval + 1):(i*interval))/Nt(i);
	a_mag_mat(i, :) = a_mag(((i - 1)*interval + 1):(i*interval));
	a_mag2_mat(i, :) = a_mag2(((i - 1)*interval + 1):(i*interval));
	a_mag_sus_mat(i, :) = a_mag2_mat(i, :) - a_mag_mat(i, :).*a_mag_mat(i, :);
	sx_mat(i, :) = sx(((i - 1)*interval + 1):(i*interval));
	loc_mat(i, :) = loc(((i - 1)*interval + 1):(i*interval));
	alpha_mat(i, :) = alpha;
end

%%%make 2 plots
figure()
plot(alpha_mat', loc_mat', '-*', 'markersize', 10, 'linewidth', 4)
xlabel('alpha')
ylabel('localization')
legend('Nt = 10000', '1000', '100', '3000', '300', '50000', 'location', 'northeast')

figure()
plot(alpha_mat', sx_mat', '-*', 'markersize', 10, 'linewidth', 4)
xlabel('alpha')
ylabel('sx')
legend('Nt = 10000', '1000', '100', '3000', '300', '50000', 'location', 'southeast')

figure()
plot(alpha_mat', mag2_mat', '-*', 'markersize', 10, 'linewidth', 4)
xlabel('alpha')
ylabel('mag2')
legend('Nt = 10000', '1000', '100', '3000', '300', '50000', 'location', 'southeast')

figure()
plot(alpha_mat', a_mag2_mat', '-*', 'markersize', 10, 'linewidth', 4)
xlabel('alpha')
ylabel('a\_mag2')
legend('Nt = 10000', '1000', '100', '3000', '300', '50000', 'location', 'southeast')

%figure()
%plot(alpha_mat', cluster_mat')
%xlabel('alpha')
%ylabel('cluster fraction')
%legend('Nt = 10000', '1000', '100', '3000', '300')



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
