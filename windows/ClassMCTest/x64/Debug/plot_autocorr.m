data = csvread('autocorrelation_results.csv');
plot_data = data(4, 2:end);
plot(plot_data)
