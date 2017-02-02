%%plot the interaction strengths as a function of x and imaginary time in the quantum
%%to classical mapping of the spin boson model

%first just plot 1d: 2 sites
%time goes from 0 to pi
%x = pi/v/beta * r can take any value - will want to choose different values for a / v / beta
avbeta = pi*0.001;
time = 0 : 0.001 : pi/2;
coeff = pi * pi / length(time) / length(time) .* (sin(time).*sin(time).*cosh(avbeta).*cosh(avbeta) - sinh(avbeta).*sinh(avbeta).*cos(time).*cos(time))...
	./(sin(time).*sin(time).*cosh(avbeta).*cosh(avbeta) + sinh(avbeta).*sinh(avbeta).*cos(time).*cos(time))./(sin(time).*sin(time).*cosh(avbeta).*cosh(avbeta) + sinh(avbeta).*sinh(avbeta).*cos(time).*cos(time));
%figure()
plot(time, coeff);

%many sites
site = 1:10;
[t, x] = meshgrid (time, avbeta * site);
z = pi * pi / length(time) / length(time) .* (sin(t).*sin(t).*cosh(x).*cosh(x) - sinh(x).*sinh(x).*cos(t).*cos(t))...
	./(sin(t).*sin(t).*cosh(x).*cosh(x) + sinh(x).*sinh(x).*cos(t).*cos(t))./(sin(t).*sin(t).*cosh(x).*cosh(x) + sinh(x).*sinh(x).*cos(t).*cos(t));
%figure()
%mesh(time, avbeta*site, z);