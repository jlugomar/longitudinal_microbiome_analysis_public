function ddata = simpleDisc(data,k)
% ddata = simpleDisc(data,k)
%
% Extremely simple discretization scheme that splits each continuous data 
% column into k discrete bins
%
% INPUT
% DATA: data matrix
% K: number of bins to discretize continuous columns into.
%
% OUTPUT:
% DDATA: discretized version of DATA matrix.


disc = IsDiscrete(data);
ddata = data;
for i = 1:length(disc)
    if (~disc(i))
        low = min(data(:,i));
        high = max(data(:,i));
        step = (high-low)/k;
        for j = 1:k-1
            ddata(data(:,i) >= low + (j-1) * step & data(:,i) < low + j * step,i) = j;
        end
        ddata(data(:,i) >= low + (k-1) * step & data(:,i) <= high, i) = k;    
    end
end
    
    
    