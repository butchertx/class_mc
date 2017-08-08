%plot the results listed in the result files here
%15 different alpha values
for i = 1:15
	data = csvread(strcat('results', num2str(i), '.csv'));
	mags(i) = data(6, 2);
	mag_err(i) = data(6, 3);
end

alpha = 0:0.1:1.4;


%make 1 plot
errorbar(alpha, mags, mag_err)