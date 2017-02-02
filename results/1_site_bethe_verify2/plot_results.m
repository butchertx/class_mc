%plot the results listed in the result files here
for i = 1:11
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(5, 2);
	mag_err(i) = data(5, 3);
end

alpha = 0:0.1:1.0;


%make 1 plot with 6 different lines on it
errorbar(alpha, mags, mag_err)
