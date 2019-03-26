function [cdf,logcdfminushalf] = logNormCDF(x,terms)

s = x;
value = x;
for i = 1:terms
  value = (value*x*x/(2*i+1));
  s = s+value;
end
cdf  = 0.5+(s/sqrt(2*pi))*exp(-(x*x)/2);

logcdfminushalf = log(s) - 0.5 * log(2 * pi) + -1*(x*x)/2;
