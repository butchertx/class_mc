for i = 1:35
	file = strcat('dump', num2str(i), '.csv');
	data = csvread(file);
	mag(i) = mean(data(5, 2:end));
end
mag'