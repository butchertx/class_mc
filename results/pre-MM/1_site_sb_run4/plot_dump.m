for i = 1:90
	file = strcat('dump', num2str(i), '.csv');
	data = csvread(file);
	if (length(data(6, :)) ~= 10001)
		display(i)
	end
	mag2(i) = mean(data(6, 2:end));
	mag4(i) = mean(data(7, 2:end));
end
binder = mag2.*mag2 ./ mag4;

h = [0.01, 0.1, 0.4, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4];
A0 = [0.01, 0.1, 1.0];
for i = 1:3
	figure()
	plot(h, binder(i:9:end), h, binder(i+3:9:end), h, binder(i+6:9:end))
end