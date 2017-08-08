#pragma once

#include <gsl/gsl_sf_zeta.h>
#include <cmath>

double real_trigamma_cmplxarg(double x, double y, int m, int l){
    //compute the real part of the trigamma function at z = x + iy using
    //my accelerated series formula
    //argument should have positive real part (x > 0)
    //m gives the number of zeta function terms to use
    //l gives the number of terms to use in the residual series
    double result = 0;

    //if x is small, the result will be large and thus, inaccurate.  Use the
    //polygamma shift formula to give a larger value of x for the computation
    if(x < 1000){
        return ((x*x - y*y)/(x*x + y*y)/(x*x + y*y)) + real_trigamma_cmplxarg(x + 1, y, m, l);
    }
    else{
        double phase, ypow = 1.0, xpow, denom;
        //compute finite sum of Hurwitz zeta functions
        for (int i = 1; i < m; ++i){
            phase = 2*(i/2) == i ? -1.0 : 1.0;
            result += phase*(2*i - 1)*ypow*gsl_sf_hzeta(2*i, x);
            ypow = ypow*y*y;
        }
        //compute the infinite sum of residuals
        phase = 2*(m/2) == m ? 1.0 : -1.0;
        for (int n = 0; n < l; ++n){
            xpow = pow(x + n, 2*m);
            denom = ((x + n)*(x + n) + y*y)*((x + n)*(x + n) + y*y);
            result += phase*ypow*((2*m - 1)*y*y + (2*m + 1)*(x + n)*(x + n)) / xpow / denom;
        }
        return result;
    }
}