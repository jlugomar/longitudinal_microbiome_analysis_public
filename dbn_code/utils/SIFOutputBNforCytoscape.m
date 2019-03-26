function SIFOutputBNforCytoscape(adjmat, nodenames, fname)
%SIFOutputBNforCytoscape(adjmat, nodenames, fname)
%
% function to output a bayesian network to the SIF (simple interaction
% format) format for importing into Cytoscape.
%
% INPUT:
% ADJMAT: adjacency matrix representing network structure.
% NODENAMES: names of each node (corresponding to variable / column names)
% FNAME: filename for output.
%
% SIF format has one line per edge and is written:
% NODEPARENT <relation> NODECHILD
%
% where <relation> is an identifier indicating the type of edge between the
% two nodes. Nodes with no edges appear by themselves on lines.  Should
% have a ".sif" file extension.
%
% Cytoscape says of delimiters and node names:
%   "Whitespace (space or tab) is used to delimit the names in the simple 
%    interaction file format. However, in some cases spaces are desired in 
%    a node name or edge type. The standard is that, if the file contains 
%    any tab characters, then tabs are used to delimit the fields and 
%    spaces are considered part of the name. If the file contains no tabs, 
%    then any spaces are delimiters that separate names (and names cannot 
%    contain spaces)."
%
%
% Copyright Michael McGeachie, 2012.  MIT license. See cgbayesnets_license.txt.

if (nargin < 3)
    fname = 'Result_BayesNetGraph.sif';
else
    fname = [fname, '.sif'];
end
fid = fopen(fname,'wt'); 

edgesfound = false(1,length(nodenames));
% print edges
[rowinds,colinds] = find(adjmat);
for i = 1:length(rowinds)
    fprintf(fid,'%s\tlink\t%s\n',nodenames{rowinds(i)}, nodenames{colinds(i)});
    edgesfound(rowinds(i)) = true;
    edgesfound(colinds(i)) = true;
end

% now print out all the nodes that weren't found
for i = 1:length(nodenames)
    if (~edgesfound(i))
        fprintf(fid,'%s\n',nodenames{i});
    end
end

fclose(fid);







