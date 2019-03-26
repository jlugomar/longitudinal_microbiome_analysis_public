function data = NormCont(data)
% data = NormCont(data)
% 
%  normalize a column of continuous data:
%
%  subtract mean and divide by std dev

d = data - mean(data);
d = d ./ std(data);
data = d;

