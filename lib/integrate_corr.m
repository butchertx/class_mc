function [ space_corr ] = integrate_corr( corr, lx, ly)
%INTEGRATE_CORR integrate the im. time dimension of the correlation
%function

space_corr = zeros(lx, 1);
for i = 1:lx
    space_corr(i) = sum(corr(((i - 1)*ly + 1):i*ly));
end


end

