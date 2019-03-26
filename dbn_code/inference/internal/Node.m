classdef Node
%
% An array of NODE() represents a hybrid discrete-continuous bayesian
% network.  The connectivity is maintained in the parents an children links
% of the NODE structures.
%
% the first Node(1) should always be the root of the Bayes Net.
%
% node(1).filename: bdn file name
% node(n).self: name of the node (just a string)
% node(n).parents: name(s) of parent nodes(s) (cell array of strings)
% node(n).children: name(s) of children nodes(s) (cell array of strings)
% node(n).values: list of the values the node variable can assume (array)
% node(n).xy: x,y location of node (array) in bayesware discoverer GUI
%   layout
% node(n).index: index into node array; redundant so that nodes(i) =
%   node.index
% node(n).cindex: indices of children into the node array, such that
%   nodes(node(i).cindex(j)) is a child of NODE(i).
% node(n).pindex: indices of parent nodes into the node array, as above
% node(n).discrete: true/false indicating if the node is a discrete
%   datatype
% node(n).colind: index into the datafile such that node.self is the column
%   name for node.colind
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

    properties
        filename
        self        % a string
        parents
        children
        values
        xy
        index
        cindex
        pindex
        discrete
        colind      % index into DATA array such that Node.self is the column name
                    % for that column of data.
    end
    
    methods
        % constructor
        function node = Node(filename, self, parents, children, values, ...
                xy, index, cindex, pindex, discrete, colind)
            if (nargin > 0)
                node.filename = filename;
                node.self = self;
                node.parents = parents;
                node.children = children;
                node.values = values;
                node.xy = xy;
                node.index = index;
                node.cindex = cindex;
                node.pindex = pindex;
                node.discrete = discrete;
                node.colind = colind;
            else
                node.filename = '';
                node.self = '';
                node.parents = {};
                node.children = {};
                node.values = {};
                node.xy = [];
                node.index = [];
                node.cindex = [];
                node.pindex = [];
                node.discrete = [];
                node.colind = [];
            end
        end
        
    end

end
    