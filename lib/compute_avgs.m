%%read a list of dump and param files and display measurements
imax = 36;
params = textread('C:/Users/Matthew/Dropbox/Code/class_mc/results/2_site_sb_run2/params.txt', '%s')
for i = 1:imax
	data = csvread(strcat('C:/Users/Matthew/Dropbox/Code/class_mc/results/2_site_sb_run2/dump',num2str(i),'.csv'));
	params(i, :)
	time = data(2, end)
	sx_avg = mean(data(3, 2:end))
	sz_avg = mean(data(4, 2:end))
	mag_avg = mean(data(5, 2:end))
	binder =  mean(data(6, 2:end)) * mean(data(6, 2:end)) / mean(data(7, 2:end))
	ising_order = mean(data(8, 2:end))
	display('')
	%param_temp = csvread(strcat('C:/Users/Matthew/Dropbox/Code/class_mc/results/2_site_sb/params',num2str(i),'.csv'));
	%params(i, :, :) = param_temp;
end

%for i = 1:imax
%%	A0 = params(i, 6, 2)
%%	tcJ = params(i, 4, 2)
%%	gamma = params(i, 4, 3)
%%	a_over_beta_v = params(i, 3, 2) / params(i, 5, 2) / params(i, 8, 2)
%
%end
