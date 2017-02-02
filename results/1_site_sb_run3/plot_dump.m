for i = 1:40
	file = strcat('dump', num2str(i), '.csv');
	data = csvread(file);
	mag2(i) = mean(data(6, 2:end));
	mag4(i) = mean(data(7, 2:end));
end
binder = mag2.*mag2 ./ mag4

h = [0.01, 0.1, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4];
beta = [0.1, 1.0, 10.0, 100.0];
plot(h, binder(1:4:end), h, binder(2:4:end), h, binder(3:4:end), h, binder(4:4:end))