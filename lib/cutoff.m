%%find the crossover distance for different values of omega_c (outside which the interaction
%%between same-time neighbors becomes antiferromagnetic)
y = [.001, .01, .1];
x = 0:.0001:.5;
zeta2 = pi*pi/6;
zeta3 = 1.202056903159594285399738161511449990764986292340498881792;
zeta4 = pi*pi*pi*pi/90;
a = -12*zeta4*ones(size(y));
b = 2*zeta2 - 8*zeta3*y - 12*zeta4*y.^2;
c = -2 + 2*zeta2*y.^2 - 16*zeta3*y.^3 + 12*zeta4*y.^4;
d = 2*y.^2 + 2*zeta2*y.^4 - 8*zeta3*y.^5 + 12*zeta4*y.^6;

mat = kron(x.^6, a') + kron(x.^4, b') + kron(x.^2, c') + kron(ones(size(x)), d');
x_mat = [x;x;x;x];

plot(x', mat')