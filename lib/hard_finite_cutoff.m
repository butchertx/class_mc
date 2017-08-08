%%Calculate interactions in the spin-boson model for a finite, hard cutoff
%%Want to calculate as a function of r for some given values of tau
%%Use Simpson's 3/8 rule for integration
M = 10;
Nt = 100;
Nu = 100;
alpha = 1.0;
beta = 1.0;
omega_c = 100;
rbv = linspace(0, beta*omega_c, M);
taub = linspace(0, beta, Nt);
kernel = zeros(M, Nt);

for i = 1:Nt
    for j = 1:M
        integrand = @(u) 2*alpha/Nt/Nt * u * cos(rbv(j)*u)*(exp(u*(1 - taub(i))) + exp(u*taub(i)))/(exp(u) - 1);
        kernel(j, i) = simpson3_8(integrand, 0, beta*omega_c, 4*alpha/Nt/Nt, integrand(beta*omega_c), Nu);
    end
end
kernel