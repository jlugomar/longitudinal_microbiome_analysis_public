function [edgemat, cellmat, regressions] = ContEdgeStrengths(nodes, tree, cols)
% function to find the parameter strengths of all continuous edges in the
% network.


tables = GetProbTables(nodes, tree);

% tables is a cell array of either CondProbTables (discrete) or LPPotential
% (continuous)

edgemat = zeros(length(cols),length(cols));
cellmat = cell(size(edgemat));
regressions = {};
for i = 1:length(tables)
    lpp = tables{i};
    nodeclass = class(lpp);
    if (strcmpi(nodeclass, 'CondProbTable'))
        % discrete node:
        continue;
    end
    
    % could be a longer cell of these if it's conditioned on many disc
    % vars
    for k = 1:length(lpp)
        klpp = lpp{k};

        % find edge address:
        % indexes of klpp.head and klpp.tail index into nodes()
        sink = nodes(klpp.head).colind; 

        regress = [cols{sink}, ' = ', num2str(klpp.const)];
        for j = 1:length(klpp.tail)
            source = nodes(klpp.tail(j)).colind;
            if (~isempty(cellmat{source,sink}))
                cellmat{source, sink} = {cellmat{source, sink}{:}, klpp.params(j)};
            else
                cellmat{source, sink} = {klpp.params(j)};
            end
            edgemat(source, sink) = klpp.params(j);
            regress = [regress, ' + ', num2str(klpp.params(j)),'*', cols{source}];
        end
        regressions{end+1} = regress;
    end
    
end



        
        