%find the roots of the trigamma function
real_x = linspace(0, 1, 100);
imag_y = linspace(0, 1, 100);
[xx, yy] = meshgrid(real_x, imag_y);
zz = gsl_sf_psi_l(