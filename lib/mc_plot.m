%analyze results of mc runs
i = 1:5;
j = 1:6;
beta = [1, 3, 10, 15];
times = zeros(length(beta), length(i), length(j));
sxs = zeros(length(beta), length(i), length(j));
szs = zeros(length(beta), length(i), length(j));
for J = i
	for D = j
		for B = 1:length(beta)
			dumppath = strjoin({'C:\Users\Matthew\Dropbox\Code\class_mc\results\Beta', num2str(beta(B)),'\8x128_J',num2str(J),'_D',num2str(D),'\dump\dump9.csv'}, '');
			data = csvread(dumppath);
			%get time - end of row 2
			times(B, J, D) = data(2, end);
			%get sz - average over row 3
			szs(B, J, D) = mean(data(3, 2:end));
			%get sx - average over row 4
			sxs(B, J, D) = mean(data(4, 2:end));
		end
	end
end

%plot - a new line for each BJ, plot 
figure()
plot(j./5, szs(3, 5, :), 'linewidth', 2, j./4, szs(4, 4, :), 'linewidth', 2, j./3, szs(2, 3, :), 'linewidth', 2)
h = legend('Beta*J = 50', 'Beta*J = 60', 'Beta*J = 9')
set(h, "fontsize", 15)
h = xlabel('Delta / J')
set(h, "fontsize", 15)
h = ylabel('<Sz> (per site)')
set(h, "fontsize", 15)
h = gca()
set(h, "fontsize", 15)

%check some of the lower delta values
dumppath1 = 'C:\Users\Matthew\Dropbox\Code\class_mc\results\Beta1\8x128_J1_D05\dump\dump9.csv';
	data = csvread(dumppath1);
	sz(1) = mean(data(3, 2:end));
	sx(1) = mean(data(4, 2:end));
dumppath2 = 'C:\Users\Matthew\Dropbox\Code\class_mc\results\Beta1\8x128_J1_D01\dump\dump9.csv';
	data = csvread(dumppath2);
	sz(2) = mean(data(3, 2:end));
	sx(2) = mean(data(4, 2:end));
dumppath3 = 'C:\Users\Matthew\Dropbox\Code\class_mc\results\Beta1\8x128_J1_D005\dump\dump8.csv';
	data = csvread(dumppath3);
	sz(3) = mean(data(3, 2:end));
	sx(3) = mean(data(4, 2:end));
dumppath4 = 'C:\Users\Matthew\Dropbox\Code\class_mc\results\Beta1\8x128_J1_D001\dump\dump2.csv';
	data = csvread(dumppath4);
	sz(4) = mean(data(3, 2:end));
	sx(4) = mean(data(4, 2:end));
dumppath5 = 'C:\Users\Matthew\Dropbox\Code\class_mc\results\Beta1\8x128_J1_D0001\dump\dump1.csv';
	data = csvread(dumppath5);
	sz(5) = mean(data(3, 2:end));
	sx(5) = mean(data(4, 2:end));
figure()
delta_over_J = [0.5, 0.1, 0.05, 0.01, 0.001];
loglog(delta_over_J, sz, 'linewidth', 2, delta_over_J, sx, 'linewidth', 2)
h = legend('Sz', 'Sx')
set(h, "fontsize", 15)
h = xlabel('Delta / J')
set(h, "fontsize", 15)
h = ylabel('<O> (per site)')
set(h, "fontsize", 15)
h = gca()
set(h, "fontsize", 15)
