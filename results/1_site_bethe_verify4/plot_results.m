%plot the results listed in the result files here
%126 files: 13 different alpha values, 3 different h values (listed by constant h first) and 2 temperatures
for i = 1:11
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(6, 2);
	mag_err(i) = data(6, 3);
end

alpha = 0:0.05:0.5;


%make 1 plot
errorbar(2*alpha, mags, mag_err)