%%integrating spin-spin interaction

%%test
%quadgk(@(x) 2*x, 10, 13)

%%calculate the "known" quantity: case for ohmic dispersion at same site, different tau
beta = 100.0;
%s = [.5, 1.0, 2.0];
%tau = 0.01:0.001:.5*beta;
%chi = ones(length(s), length(tau));
%for i = 1:length(tau)
%  for j = 1:length(s)
%	  chi(i, j) = quadcc(@(x) cos(s(j)*x).*x.*cosh((.5*beta - tau(i))*x)./sinh(.5*beta*x), 0, Inf);
%   end
%end
%figure()
%plot(tau, chi(:, 1), tau, chi(:, 2), tau, chi(:, 3))
%legend('s = .5')

%plot the integrand
%f = @(x) -besselj(0, 1*x).*x.*x.*cosh(.45*x.*sqrt(1 + x.*x))./sinh(.5*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x);
%figure()
%plot(0:.01:10, f(0:.01:10))
%

%%%start with chi(tau) for different s
%beta = 100.0;
%s = beta*[.01, .1, .5, 1.0, 2.0];
%tau = 0.01*beta:(0.01*beta):.5*beta;
%chi = ones(length(s), length(tau));
%for i = 1:length(s)
%	for j = 1:length(tau)
%		chi(i, j) = quadcc(@(x) -sin(s(i) * x).*x.*cosh((.5*beta - tau(j))*x.*sqrt(1 + x.*x))./sinh(.5*beta*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x), 0, 100);
%	end
%end
%chi_max = 1;
%figure()
%%plot(tau, chi(1, :), tau, chi(2, :) , tau, chi(3,:), tau, chi(4,:), tau, chi(5,:))
%plot(tau, chi(1, :) / max(abs(chi(1, :))), tau, chi_max / max(abs(chi(2, :)))* chi(2, :) , tau, chi_max /max(abs(chi(3, :)))*chi(3,:), tau, chi_max / max(abs(chi(4, :)))*chi(4,:), tau, chi_max / max(abs(chi(5, :)))*chi(5,:))
%title('1D')
%legend('s = .1', 's = 1', 's = 5', 's = 10', 's = 20')
%xlabel('tau-squiggle')
%ylabel('chi')
%
%
%s = beta*[.01, .1, .5, 1.0, 2.0];
%tau = 0.01*beta:(0.01*beta):.5*beta;
%chi = ones(length(s), length(tau));
%for i = 1:length(s)
%	for j = 1:length(tau)
%		chi(i, j) = quadcc(@(x) -besselj(0, s(i)*x).*x.*x.*cosh((.5*beta - tau(j))*x.*sqrt(1 + x.*x))./sinh(.5*beta*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x), 0, 100);
%	end
%end
%figure()
%%plot(tau, chi(1, :), tau, chi(2, :) , tau, chi(3,:), tau, chi(4,:), tau, chi(5,:))
%plot(tau, chi(1, :) / max(abs(chi(1, :))), tau, chi_max / max(abs(chi(2, :)))* chi(2, :) , tau, chi_max /max(abs(chi(3, :)))*chi(3,:), tau, chi_max / max(abs(chi(4, :)))*chi(4,:), tau, chi_max / max(abs(chi(5, :)))*chi(5,:))
%title('2D')
%legend('s = .1', 's = 1', 's = 5', 's = 10', 's = 20')
%xlabel('tau-squiggle')
%ylabel('chi')
%
%
%s = beta*[.01, .1, .5, 1.0, 2.0];
%chi = ones(length(s), length(tau));
%for i = 1:length(s)
%	for j = 1:length(tau)
%		chi(i, j) = 1/s(i) * quadcc(@(x) -sin(s(i)*x).*x.*x.*cosh((.5*beta - tau(j))*x.*sqrt(1 + x.*x))./sinh(.5*beta*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x), 0, 100);
%	end
%end
%figure()
%%plot(tau, chi(1, :), tau, chi(2, :) , tau, chi(3,:), tau, chi(4,:), tau, chi(5,:))
%plot(tau, chi(1, :) / max(abs(chi(1, :))), tau, chi_max / max(abs(chi(2, :)))* chi(2, :) , tau, chi_max /max(abs(chi(3, :)))*chi(3,:), tau, chi_max / max(abs(chi(4, :)))*chi(4,:), tau, chi_max / max(abs(chi(5, :)))*chi(5,:))
%title('3D')
%legend('s = .1', 's = 1', 's = 5', 's = 10', 's = 20')
%xlabel('tau-squiggle')
%ylabel('chi')
%
%%
%
%
%%%%now do chi(s) for different tau
tau = beta*[0, .05, .1, .25, .5];
s = 0.01:0.5:beta;
chi = ones(length(s), length(tau));
for i = 1:length(s)
	for j = 1:length(tau)
		chi(i, j) = quadcc(@(x) -sin(s(i)*x).*x.*cosh((.5*beta - tau(j))*x.*sqrt(1 + x.*x))./sinh(.5*beta*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x), 0, .1);
	end
end
figure()
plot(s, chi(:, 1), s, chi(:, 2), s, chi(:, 3), s, chi(:, 4), s, chi(:, 5))
%plot(s, chi(:, 1)/max(abs(chi(:, 1))), s, chi(:, 2)/max(abs(chi(:, 2))), s, chi(:, 3)/max(abs(chi(:, 3))), s, chi(:, 4)/max(abs(chi(:, 4))), s, chi(:, 5)/max(abs(chi(:, 5))))
title('1D')
legend('tau = .01', 'tau = .05', 'tau = .1', 'tau = .25', 'tau = .5')
xlabel('s')
ylabel('chi')


s = 0.01:0.5:beta;
chi = ones(length(s), length(tau));
for i = 1:length(s)
	for j = 1:length(tau)
		chi(i, j) = quadcc(@(x) -besselj(0, s(i)*x).*x.*x.*cosh((.5*beta - tau(j))*x.*sqrt(1 + x.*x))./sinh(.5*beta*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x), 0, .1);
	end
end
figure()
plot(s, chi(:, 1), s, chi(:, 2), s, chi(:, 3), s, chi(:, 4), s, chi(:, 5))
%plot(s, chi(:, 1)/max(abs(chi(:, 1))), s, chi(:, 2)/max(abs(chi(:, 2))), s, chi(:, 3)/max(abs(chi(:, 3))), s, chi(:, 4)/max(abs(chi(:, 4))), s, chi(:, 5)/max(abs(chi(:, 5))))
title('2D')
legend('tau = .01', 'tau = .05', 'tau = .1', 'tau = .25', 'tau = .5')
xlabel('s')
ylabel('chi')



s = 0.01:0.5:beta;
chi = ones(length(s), length(tau));
for i = 1:length(s)
	for j = 1:length(tau)
		chi(i, j) = 1/s(i) * quadcc(@(x) -sin(s(i)*x).*x.*x.*cosh((.5*beta - tau(j))*x.*sqrt(1 + x.*x))./sinh(.5*beta*x.*sqrt(1 + x.*x))./sqrt(1 + x.*x), 0, .1);
	end
end
figure()
plot(s, chi(:, 1), s, chi(:, 2), s, chi(:, 3), s, chi(:, 4), s, chi(:, 5))
%plot(s, chi(:, 1)/max(abs(chi(:, 1))), s, chi(:, 2)/max(abs(chi(:, 2))), s, chi(:, 3)/max(abs(chi(:, 3))), s, chi(:, 4)/max(abs(chi(:, 4))), s, chi(:, 5)/max(abs(chi(:, 5))))
title('3D')
legend('tau = .01', 'tau = .05', 'tau = .1', 'tau = .25', 'tau = .5')
xlabel('s')
ylabel('chi')

%%%%now do chi(s) for the linearly dispersing case
tau = beta*[0, .05, .1, .25, .5];
s = 0.01:0.5:beta;
chi = ones(length(s), length(tau));
for i = 1:length(s)
	for j = 1:length(tau)
		chi(i, j) = quadcc(@(x) -cos(s(i)*x).*x.*cosh((.5*beta - tau(j))*x)./sinh(.5*beta*x), 0, 1);
	end
end
figure()
plot(s, chi(:, 1), s, chi(:, 2), s, chi(:, 3), s, chi(:, 4), s, chi(:, 5))
%plot(s, chi(:, 1)/max(abs(chi(:, 1))), s, chi(:, 2)/max(abs(chi(:, 2))), s, chi(:, 3)/max(abs(chi(:, 3))), s, chi(:, 4)/max(abs(chi(:, 4))), s, chi(:, 5)/max(abs(chi(:, 5))))
title('1D')
legend('tau = .01', 'tau = .05', 'tau = .1', 'tau = .25', 'tau = .5')
xlabel('s')
ylabel('chi')

%%plot chi(s) from what we calculated analytically in Jed's notes
for j = 1:length(tau)
	chi(:, j) = 1/beta*(sin(pi*tau(j)/beta)*sin(pi*tau(j)/beta)*cosh(pi*s).*cosh(pi*s) - cos(pi*tau(j)/beta)*cos(pi*tau(j)/beta)*sinh(pi*s).*sinh(pi*s))...
				./(sin(pi*tau(j)/beta)*sin(pi*tau(j)/beta)*cosh(pi*s).*cosh(pi*s) + cos(pi*tau(j)/beta)*cos(pi*tau(j)/beta)*sinh(pi*s).*sinh(pi*s))...
				./(sin(pi*tau(j)/beta)*sin(pi*tau(j)/beta)*cosh(pi*s).*cosh(pi*s) + cos(pi*tau(j)/beta)*cos(pi*tau(j)/beta)*sinh(pi*s).*sinh(pi*s));
end
figure()
plot(s, chi(:, 1), s, chi(:, 2), s, chi(:, 3), s, chi(:, 4), s, chi(:, 5))


%s = 0:0.01:2;
%mrho = .1;
%chi_1D = -.5*mrho*exp(-s);
%chi_2D = -2/pi*mrho*mrho * besselk(0, s);
%chi_3D = -1*mrho*mrho/(10*pi)./s.*exp(-s);
%plot(s, chi_1D, s, chi_2D, s, chi_3D)
%xlabel('s')
%ylabel('chi')
%title('Long-Time Interactions')
%legend('1D', '2D', '3D', 'location', 'southeast')