data = csvread('dump/dump9.csv');

figure()
loglog(abs(data(9, 2:6)))
title('loglog')

figure()
semilogy(abs(data(9, 2:6)))
title('semilogy')