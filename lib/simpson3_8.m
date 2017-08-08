function result = simpson3_8(f_integrand, a, b, f_a, f_b, N)
%%Perform simpson's 3/8 rule integration on the function f_integrand
%%a: starting point, b: ending point, f_a: f(a), f_b: f(b), N: number of
%%intervals (multiply this by 3)
h = (b - a)/(3*N);
result = 3*h*(f_a + f_b)/8;
for i = 1:(3*N - 1)
    if(mod(i, 3) == 0)
        result = result + 6*h*f_integrand(a + i*h)/8;
    else
        result = result + 9*h*f_integrand(a + i*h)/8;
    end
end