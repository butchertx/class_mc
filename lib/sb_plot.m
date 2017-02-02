imax = 16;
params = textread('C:/Users/Matthew/Dropbox/Code/class_mc/results/1_site_sb/params.txt', '%s');
%read in results and compute averages
for i = 1:imax
	data = csvread(strcat('C:/Users/Matthew/Dropbox/Code/class_mc/results/2_site_sb_run3/dump',num2str(i),'.csv'));
	sx_avg(i) = mean(data(3, 2:end));
	sz_avg(i) = mean(data(4, 2:end));
	mag_avg(i) = mean(data(5, 2:end));
	binder(i) =  mean(data(6, 2:end)) * mean(data(6, 2:end)) / mean(data(7, 2:end));
	ising_order(i) = mean(data(8, 2:end));
	corr_t{i} = data(9, 2:100);
	[J(i), D(i), a(i), A0(i)] = parse_sb_params(params{i});
end

%run through param lists and plot points
%4 different quantum ising couplings, 3 seperations, then 4 A0 couplings
%->4 different figures for each observable with 4 lines on each

%sx
for i = 1:4:12
		plot(A0(i : i + 3), sx_avg(i: i + 3),...
			A0(i + 12: i + 15), sx_avg(i + 12: i + 15),...
			A0(i + 24: i + 27), sx_avg(i + 24: i + 27),...
			A0(i + 36: i + 39), sx_avg(i + 36: i + 39));
		legend(strcat('tcJ = ', num2str(J(i) / 100)), strcat('tcJ = ', num2str(J(i + 12) / 100)), strcat('tcJ = ', num2str(J(i + 24) / 100)), strcat('tcJ = ', num2str(J(i + 36) / 100)));
		title(strcat('a/(beta*v) = ', num2str(a(i))));
		xlabel('A0');
		ylabel('<sx> / Nx');
		saveas(1, strcat('../plots/2site_sb_run3/Sx_a', num2str(a(i)), '.png'));
end

%sz
for i = 1:4:12
		plot(A0(i : i + 3), sz_avg(i: i + 3),...
			A0(i + 12: i + 15), sz_avg(i + 12: i + 15),...
			A0(i + 24: i + 27), sz_avg(i + 24: i + 27),...
			A0(i + 36: i + 39), sz_avg(i + 36: i + 39));
		legend(strcat('tcJ = ', num2str(J(i) / 100)), strcat('tcJ = ', num2str(J(i + 12) / 100)), strcat('tcJ = ', num2str(J(i + 24) / 100)), strcat('tcJ = ', num2str(J(i + 36) / 100)));
		title(strcat('a/(beta*v) = ', num2str(a(i))));
		xlabel('A0');
		ylabel('<sz> / Nx');
		saveas(1, strcat('../plots/2site_sb_run3/Sz_a', num2str(a(i)), '.png'));
end


%ising order
for i = 1:4:12
		plot(A0(i : i + 3), ising_order(i: i + 3),...
			A0(i + 12: i + 15), ising_order(i + 12: i + 15),...
			A0(i + 24: i + 27), ising_order(i + 24: i + 27),...
			A0(i + 36: i + 39), ising_order(i + 36: i + 39));
		legend(strcat('tcJ = ', num2str(J(i) / 100)), strcat('tcJ = ', num2str(J(i + 12) / 100)), strcat('tcJ = ', num2str(J(i + 24) / 100)), strcat('tcJ = ', num2str(J(i + 36) / 100)));
		title(strcat('a/(beta*v) = ', num2str(a(i))));
		xlabel('A0');
		ylabel('<C(a, 0)>_t');
		saveas(1, strcat('../plots/2site_sb_run3/corrx_a', num2str(a(i)), '.png'));
end


%binder cumulant
for i = 1:4:12
		plot(A0(i : i + 3), binder(i: i + 3),...
			A0(i + 12: i + 15), binder(i + 12: i + 15),...
			A0(i + 24: i + 27), binder(i + 24: i + 27),...
			A0(i + 36: i + 39), binder(i + 36: i + 39));
		legend(strcat('tcJ = ', num2str(J(i) / 100)), strcat('tcJ = ', num2str(J(i + 12) / 100)), strcat('tcJ = ', num2str(J(i + 24) / 100)), strcat('tcJ = ', num2str(J(i + 36) / 100)));
		title(strcat('a/(beta*v) = ', num2str(a(i))));
		xlabel('A0');
		ylabel('<m^2>^2 / <m^4>');
		saveas(1, strcat('../plots/2site_sb_run3/binder_a', num2str(a(i)), '.png'));
end


%correlation of t
%
%for i = 1:4
%	figure()
%	hold on
%	for j = 1:12
%		plot(1:length(corr_t{(i - 1)*12 + j}), corr_t{(i - 1)*12 + j})
%	end
%end