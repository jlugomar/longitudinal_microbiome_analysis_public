function auc = AUCfromBNfinderFile(fname, trueclass)
%auc = AUCfromBNfinderFile(fname, trueclass)
%
% load a .cls file output from the BNfinder bnc command
% compute a AUC for the predictions contained therein
%
% 

% load file:
fid = fopen(fname);
% skip header line
% need to skip the %[^\n] construction, which is aparently broken in R2013a
textscan(fid, '%s', 1, 'bufsize', 1000000, 'Delimiter', '\n');

textscan(fid, '%*s', 1, 'delimiter','\t');
data = textscan(fid, '%f', 'delimiter','\t');
data = data{1};

auc = AUCWorker(.5, data > .5, data, trueclass, true, true);

fclose(fid);

