function outputBayesNetGraph(adjMatrix,nodeNames, fname)
%outputBayesNetGraph(adjMatrix,nodeNames, fname)
%
% This function outputs Bayesian network structure into .tgf format
% 
% INPUT: 
% ADJMATRIX: adjaceny matrix representing a network
% NODENAMES: names of variable, ordered as in ADJMATRIX
% FNAME: file name to output too.
%
% Copyright Hsun-Hsien Chang, 2010.  MIT license. See cgbayesnets_license.txt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%% print network structure
if (nargin < 3)
    fname = 'Result_BayesNetGraph.tgf';
else
    fname = [fname, '.tgf'];
end
fileID = fopen(fname,'wt'); 

% print class node
for g=1:length(nodeNames)
    fprintf(fileID,'%s\t%s\n',num2str(g),nodeNames{g});
end

% print #
fprintf(fileID, '#\n');

% print edges
for g=1:length(nodeNames)
    childNodes = find(adjMatrix(g,:)>0);
    if ~isempty(childNodes)
        for c = childNodes
            fprintf(fileID,'%s\t%s\n',num2str(g),num2str(c));
        end
    end
end

fclose(fileID);
