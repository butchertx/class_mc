%%Read full_input.txt to obtain parameters
N = 12;%number of different param sets
params = SB_Params('full_input.txt', N)

%%Read all output files to obtain and organize data
%read all dump files
dump0 = csvread('dump0.csv', 0, 1);
full_meas = zeros(N, size(dump0, 1), size(dump0, 2));
size(full_meas)
for i = 1:N
    full_meas(i, :, :) = csvread(strcat('dump',num2str(i - 1),'.csv'), 0, 1);
end

%read results file
fid = fopen('results.dat');
data_columns = textscan(fgetl(fid), '%s');
column_names = data_columns{1}';
results = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f');

%%Plot results
figure()
hold on
plot(params.alpha, results{3})
xlabel('alpha')
ylabel('loc2')

%%Plot results with resampling
%try just one interval with two different data sets, and interpolate
[S1,S2] = split_action(full_meas(4, 7, :), params.alpha(4), params.gamma, full_meas(4, 8, :), params.ly, params.lx);
resample_loc2_alpha02 = swendsen_resample_point(params.alpha(4), S2, 0.2:0.005:0.25, full_meas(4, 2, :));

[S1,S2] = split_action(full_meas(5, 7, :), params.alpha(5), params.gamma, full_meas(5, 8, :), params.ly, params.lx);
resample_loc2_alpha025 = swendsen_resample_point(params.alpha(5), S2, 0.2:0.005:0.25, full_meas(5, 2, :));

plot(0.2:0.005:0.25, resample_loc2_alpha02, '-x')
plot(0.2:0.005:0.25, resample_loc2_alpha025, '-x')
plot(0.2:0.005:0.25, 0.5*(resample_loc2_alpha02 + resample_loc2_alpha025), '-x')

%Now try resampling over the whole interval
[S1,S2] = split_action(full_meas(4, 7, :), params.alpha(4), params.gamma, full_meas(4, 8, :), params.ly, params.lx);
resample_loc2_alpha02 = swendsen_resample_point(params.alpha(4), S2, 0.05:0.01:0.6, full_meas(4, 2, :));

[S1,S2] = split_action(full_meas(5, 7, :), params.alpha(5), params.gamma, full_meas(5, 8, :), params.ly, params.lx);
resample_loc2_alpha025 = swendsen_resample_point(params.alpha(5), S2, 0.05:0.01:0.6, full_meas(5, 2, :));

plot(0.05:0.01:0.6, resample_loc2_alpha02, '-x')
plot(0.05:0.01:0.6, resample_loc2_alpha025, '-x')
plot(0.05:0.01:0.6, 0.5*(resample_loc2_alpha02 + resample_loc2_alpha025), '-x')



%%Fit correlations 