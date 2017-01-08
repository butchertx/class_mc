data = csvread('autocorrelation_results.csv');
plot_data = data(4, 2:end)/ max(abs(data(4, 2:end)));
for i = 8:4:400
	plot_data = [plot_data; data(i, 2:end)];
	plot_data(end, :) = plot_data(end, :) / max(abs(plot_data(end, :)));
end
figure()
hold on
for i = 90
	plot(plot_data(i, :))
end