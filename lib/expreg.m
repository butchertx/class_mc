function [ A, xi ] = expreg( xvals, yvals )
%EXPREG perform an exponential regression on vals
%   return A and xi from the form y = Ae^-x/xi
%   Assume length(xvals) = length(yvals)
modelfun = @(b, x) b(1)*exp(-x/b(2));
[A, xi] = nlinfit(xvals, yvals, modelfun, [yvals(1), xvals(2) - xvals(1)]); 

end

