function success = check_rel_freq(datapath)
display(datapath)
beta = 1.0;
N = 64;
mags = 0:2:N;
mrange = .5*(N - mags);
erange = -.5*(mags.*mags - N) / N;
energies = csvread(datapath, [2 1 2 1000]);
zvals = zeros(size(mags));
for i = 1:length(erange)
	zvals(i) = bincoeff(N, mrange(i)) * exp(-beta*erange(i));
end
min(zvals)
max(zvals)
zvals;
figure()
hold on
hist(energies, 15*N, 1);
plot(erange, .055*zvals / max(zvals));

hold off
%check evolution of energies as well
%for i = 1:10
%	energychunk = energies((i - 1) * 100 + 1 : i*100);
%	figure();
%	hist(energychunk, 20, 1);
%end

success = 1;
end