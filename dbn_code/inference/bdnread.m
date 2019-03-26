function node = bdnread(filename)
%node = bdnread(filename)
%
% bdnread(filename) read a network from Bayesware Discoverer .bdn file
%
% node: a cell array of network node structures
% node(1) is the root (highest parent) node node(1).parents: {NIL}
% node(1).filename: bdn file name
% node(n).self: name of the node (cell containing a string)
% node(n).parents: name(s) of parent nodes(s) (cell array of strings)
% node(n).values: list of the values the node variable can assume (array)
% node(n).xy: x,y location of node (array)
%
% node(n).children: determined by crawling the network
%                   name(s) of child node(s) (cell array of strings)
% 
% Note: uses only .bdn file "defnode" entries, other data is ignored.

  node = struct('filename',{},'self',{},'parents',{},...
      'children',{},'values',{},'xy',{},'discrete',{});

node(1).filename = filename;
i=0;
fid=fopen(char(filename));

tline = fgetl(fid);

while ischar(tline);
    
    [cmd, remain] = strtok(tline);
    
    switch cmd
        
        % parse defnode commands to get network nodes and edges
        case 'defnode'
            
            i=i+1;
            [self, remain] = strtok(remain);
            node(i).self = textscan(self,'%s',1);
            node(i).self = node(i).self{1};
            % mmcgeach: all BWD nodes are discrete:
            node(i).discrete = true;
            %node(i).self =cparse(self);
            remain = strtrim(remain); % pull off left over space
            [nodevalues, remain] = strtok(remain,'()');
            node(i).values = str2num(nodevalues);
            [discard, remain] = strtok(remain);
            remain = strtrim(remain); % pull off left over space
            
            % asymmetry in .bdn, "(parents)" vs "NIL" for no parents
            if ~isempty(findstr(remain,'('))
                [parents, remain] = strtok(remain,'()');
                node(i).parents = textscan(parents, '%s');
                node(i).parents = node(i).parents{1};
                %node(i).parents = cparse(parents);
                [discard, remain] = strtok(remain);
            else
                [parents, remain] = strtok(remain);
                node(i).parents = textscan(parents, '%s');
                node(i).parents = node(i).parents{1};
                %node(i).parents = cparse(parents);
                remain = strtrim(remain);
            end
            
            node(i).xy = str2num(remain);
            %disp(sprintf('%s', node(i).name, ' <-- ', node(i).parents));
            
            %ignore the rest
        otherwise
            
    end
    
    tline = fgetl(fid);
end
fclose(fid);

node = bdnfindchildren(node);

function node = bdnfindchildren(node)

% traverse the network and add children to each node
% this is just a convenience to make the network easier to work with

nodes=size(node,2);
for i=1:nodes
    if strcmp(node(i).parents{1}, 'NIL')
        % it is the root, do nothing
    else
        % make self a child on each parent
        jj = length(node(i).parents);
        for j= 1: jj
            node = bdnaddchild(node, node(i).parents{j}, node(i).self{1});
            % child to add: node(i).self
            % where to add: node(i).parents{j}
        end
    end
end

function node = bdnaddchild(node, name, childname)

% add a child to node named name

nodes=size(node,2);

for i=1:nodes
    % find node named name
    if strcmp(node(i).self, name)
        if isempty(node(i).children)
            node(i).children{1} = childname;
            return
        else
            node(i).children{length(node(i).children)+1} = childname;
            return
        end
    end
end


