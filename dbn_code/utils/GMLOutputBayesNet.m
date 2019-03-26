function GMLOutputBayesNet(BN, fname, add_nodeatts, add_nodeatt_names, add_edgeatts, add_edgeatt_names, weightlimit)
%GMLOutputBayesNet(BN, fname, add_nodeatts, add_nodeatt_names, add_edgeatts, add_edgeatt_names, weightlimit)
%
% Output a Bayes Net in Graph ML format.
%
% INPUT:
%   BN: A class BayesNet object encoding a network
%   FNAME: filename to output, will have '.graphml' added to the end.
%       optional.
%   ADD_NODEATTS: additional node attributes.  a cell array of matrices the
%       same size as # of nodes in the BN.  Can be empty.  Optional.
%   ADD_NODEATT_NAMES: names of node attributes, parallel to ADD_NODEATTS
%   ADD_EDGEATTS: additional edge attributes.  a cell array of matrices the
%       same size as # adjacency matrix of the BN.  Can be empty.  Optional.
%   ADD_EDGEATT_NAMES: names of edge attributes, parallel to ADD_EDGEATTS
%   WEIGHTLIMIT: maximum value for edge weights.  Can be used to limit the
%       range of possible weights assigned to edges for future display in
%       cytoscape.
%   
%
% Copyright Michael McGeachie, 2013.  MIT license. See cgbayesnets_license.txt.

if (nargin < 2)
    fname = 'Result_BayesNetGraph.graphml';
else
    fname = [fname, '.graphml'];
end
if (nargin < 3)
    add_nodeatts = {};
end
if (nargin < 4)
    add_nodeatt_names = {};
end
if (nargin < 5)
    add_edgeatts = {};
end
if (nargin < 6)
    add_edgeatt_names = {};
end
if (nargin < 7)
    weightlimit = inf;
end

fid = fopen(fname,'wt'); 

fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<graphml xmlns="http://graphml.graphdrawing.org/xmlns"\n');
fprintf(fid,'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n');
fprintf(fid,'xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns\n');
fprintf(fid,'http://graphml.graphdrawing.org/xmlns/1.1/graphml.xsd">\n');
fprintf(fid,'<key id="key_weight" for="edge" attr.name="weight" attr.type="double"/>\n');
for i = 1:length(add_nodeatt_names)
    fprintf(fid,'<key id="key_%s" for="node" attr.name="%s" attr.type="double"/>\n',...
        add_nodeatt_names{i}, add_nodeatt_names{i});
end
for i = 1:length(add_edgeatt_names)
    fprintf(fid,'<key id="key_%s" for="edge" attr.name="%s" attr.type="double"/>\n',...
        add_edgeatt_names{i}, add_edgeatt_names{i});
end
fprintf(fid,'<graph id="%s" edgedefault="directed">\n',BN.title);

% now print each node:
for i = 1:length(BN.cols)
    if (length(add_nodeatt_names) == 0)
        fprintf(fid,'<node id="%s"/>\n',BN.cols{i});
    else
        fprintf(fid,'<node id="%s">\n',BN.cols{i});
        for j = 1:length(add_nodeatt_names)
            fprintf(fid,'<data key="key_%s">%f</data>\n',add_nodeatt_names{j},...
                add_nodeatts{j}(i));
        end
        fprintf(fid,'</node>\n');
    end
end

% and now print each edge:
edgenum = 0;
[rowinds,colinds] = find(BN.adjmat);
for i = 1:length(rowinds)
    edgenum = edgenum + 1;
    fprintf(fid,'<edge id="e%d" source="%s" target="%s">\n', edgenum,...
        BN.cols{rowinds(i)},BN.cols{colinds(i)});
    wlm = min(full(BN.weightMatrix(rowinds(i),colinds(i))),weightlimit);
    fprintf(fid,'<data key="key_weight">%f</data>\n',wlm);
    for k = 1:length(add_edgeatt_names)
        fprintf(fid,'<data key="key_%s">%f</data>\n',...
            add_edgeatt_names{k}, add_edgeatts{k}(rowinds(i),colinds(i)));
    end
    fprintf(fid,'</edge>\n');
end

% now close the elements:
fprintf(fid,'</graph>\n</graphml>\n');

fclose(fid);







