t = 0:0.01:pi;
coshx = [1.01, 1.1, 1.5];
figure()
hold on
for i = coshx
	plot(t, sin(2*t)./(cos(2*t) - i))
end

q1 = 12*[0.01, 0.1, 0.5];
q2 = 1/3*q1;
q3 = 1/12*q1;
t1 = acos(sqrt(.5./(1 + q1.*q1) .* (1 + q1.*q1*(1 - 1.01) + sqrt(1 + q1.*q1*(1 - 1.01*1.01)))));
t2 = acos(sqrt(.5./(1 + q2.*q2) .* (1 + q2.*q2*(1 - 1.1) + sqrt(1 + q2.*q2*(1 - 1.1*1.1)))));
t3 = acos(sqrt(.5./(1 + q3.*q3) .* (1 + q3.*q3*(1 - 1.5) + sqrt(1 + q3.*q3*(1 - 1.5*1.5)))));

plot(t1 + pi/2, q1, 'r*', t2 + pi/2, q2, 'r*', t3 + pi/2, q3, 'r*')