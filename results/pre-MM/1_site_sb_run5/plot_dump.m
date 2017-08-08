for i = 1:120
	file = strcat('results', num2str(i), '.csv');
	data = csvread(file);
	mag(i) = data(6, 2);
	mag_err(i) = data(6, 3);
	mag2(i) = data(7, 2);
	mag4(i) = data(8, 2);
end
binder = 1 - (1/3)*mag4 ./ mag2./mag2;

h = [0.01, 0.1, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4];
A0 = [0.00, 0.01, 0.1, 1.0];
beta = [1.0, 10.0, 100.0];
for i = 1:3
%	figure()
%	plot(h, binder(i:12:end), h, binder(i+3:12:end), h, binder(i+6:12:end), h, binder(i + 9:12:end))
%	title(strcat('binder, beta = ', num2str(beta(i))))
%	legend('A0 = 0.0', 'A0 = 0.01', 'A0 = 0.1', 'A0 = 1.0');
	figure()
	plot(h, mag(i:12:end), h, mag(i+3:12:end), h, mag(i+6:12:end), h, mag(i + 9:12:end))
	title(strcat('magnetization, beta = ', num2str(beta(i))))
	legend('A0 = 0.0', 'A0 = 0.01', 'A0 = 0.1', 'A0 = 1.0');
	ylabel('<m>')
	xlabel('h/delta')
end