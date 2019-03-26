function OutputForBNFinder(ofile, data, names, disc, backwardsK2)
%OutputForBNFinder(ofile, data, names, disc, backwardsK2)
%
% Writes a dataset to files that are readable by the BNfinder2.0 program.
% Writes out node-parent search constraints for each of the variables
% listed in NAMES in order, assuming a K2-ordered list of nodes.  Works
% best if NAMES is a topological ordering of the nodes in the true network
% that generated the DATA.
%
% Writes two files out, one suitable for training with the BNF command, and
% one suitable for testing/prediction, with the BNC command.
%
% INPUT
% OFILE: output file name
% DATA: data matrix, variables stored in columns
% NAMES: column names for the variables in DATA.
% DISC: boolean array parallel to NAMES indicating if the data column is
%   discrete or continuous
% BACKWARDSK2: if true, will output the NAMES in node-parent constraints in
%   the reverse order. 
%
% (c) Michael McGeachie, 2013, all rights reserved.

if (nargin < 5)
    backwardsK2 = false;
end
  
if (backwardsK2)
    fid = fopen([ofile,'_revK2.txt'],'w');
    fidtest = fopen([ofile,'_revK2_test.txt'],'w');
else
    fid = fopen([ofile,'.txt'],'w');
    fidtest = fopen([ofile,'_test.txt'],'w');
end
fids = [fid, fidtest];



% nominate node 1 and 2 as regulators
%fprintf(fid, '#regulators %s %s\n', names{1}, names{2});

fprintf(fidtest, '#regulators');
for i = 2:length(disc)
    fprintf(fidtest, ' N%s', names{i});
end
fprintf(fidtest, '\n');

% make a list of allowable parents, using the K2 system:
% node number N can have parents numbered 1...N-1
if (~backwardsK2)
    for i = 1:length(disc)
        fprintf(fid, '#parents N%s', names{i});
        for j = i-1:-1:1
            fprintf(fid, ' N%s', names{j}); 
        end
        fprintf(fid,'\n');
    end
end

% backwards K2 listing
if (backwardsK2)
    for i = length(disc):-1:1
        fprintf(fid, '#parents N%s', names{i});
        for j = i+1:length(disc)
            fprintf(fid, ' N%s', names{j}); 
        end
        fprintf(fid,'\n');
    end
end

% write out which values (variables) are continuous
for f = 1:length(fids)
    fprintf(fids(f), '#continuous');
    for i= 1:length(disc)
        if (~disc(i))
            fprintf(fids(f),' N%s', names{i});
        end
    end

    % write out a list of all data points, which are named
    fprintf(fids(f), '\nclasses');
    for i = 1:size(data,1)
        fprintf(fids(f), ['\tEXP', num2str(i)]);
    end
    fprintf(fids(f),'\n');

    % write out the actual data, which is transposed
    % trick here: don't want to print phenotype out on test data file:
    for i = f:size(data,2)
        fprintf(fids(f), 'N%s', names{i});
        for j = 1:size(data,1)
            fprintf(fids(f), ['\t', num2str(data(j,i))]);
        end
        fprintf(fids(f),'\n');
    end

    fclose(fids(f));
end
