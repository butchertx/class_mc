for i = 1:16
	file = strcat('dump', num2str(i), '.csv');
	data = csvread(file);
	mag(i) = mean(data(5, 2:end))
end