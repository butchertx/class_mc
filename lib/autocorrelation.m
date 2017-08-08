function [ autocorr ] = autocorrelation( vals )
%AUTOCORRELLATION calculate the autocorrelation function of a list of
%values

autocorr = fft(vals);
autocorr = abs(autocorr);
autocorr = ifft(autocorr) - mean(vals)*ones(1,length(vals));

end

