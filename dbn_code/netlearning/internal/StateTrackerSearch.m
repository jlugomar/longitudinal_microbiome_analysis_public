classdef StateTrackerSearch
% class for implementing network search algorithms that know which
% states need to be re-evalutated and which don't.
% Keeps lists of all possible edges and whether or not they are: legal,
% illegal, possible, and need to be reevaluated.
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.
    
    properties
        contData        % required input
        discData        % required input
        priorPrecision  % required input
        phencol         % required input
        backtracking    % optional input
        bf_thresh       % optional input
        nophenotype     % optional input
        checkRepeats    % optional input
        numContNode     % computed from input
        numDiscNode     % computed from input
        n               % computed from input
        bn              % initialized to default
        weightMatrix    % initialized to default
        lldiff          % initialized to default
        llscore         % initialized to default
        llrecompute     % initialized to default
        numparents      % initialized to default
        edges           % initialized to default
        closed          % initialized to default
        never           % initialized to default
        done            % initialized to default
        remflag         % initialized to default
        prevBNmap       % initialized to default
        cycles          % initialized to default
        self            % allow self-loops? true for dynamic bayes nets 
    end
    
    methods
        % constructor
        function sts = StateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats)
            sts.contData = contData;
            sts.discData = discData;
            sts.priorPrecision = priorPrecision;
            if (nargin < 4)
                % default is last discrete var is the phenotype, which is
                % the bottom-est row of data.
                sts.phencol = size(contData,1) + size(discData,1);
            else
                sts.phencol = phencol;
            end
            if (nargin < 5)
                sts.backtracking = false;
            else
                sts.backtracking = backtracking;
            end
            if (nargin < 6)
                sts.bf_thresh = 0;
            else
                sts.bf_thresh = bf_thresh;
            end
            if (nargin < 7)
                sts.nophenotype = false;
            else
                sts.nophenotype = nophenotype;
            end
            if (nargin < 8)
                sts.checkRepeats = true;
            else
                sts.checkRepeats = checkRepeats;
            end
            sts.numContNode = size(contData,1);
            sts.numDiscNode = size(discData,1);
            sts.n = sts.numContNode+sts.numDiscNode;
            sts.bn = zeros(sts.n);
            sts.weightMatrix = zeros(sts.n);
            sts.lldiff = zeros(sts.n);
            sts.llscore = zeros(sts.n);
            sts.llrecompute = true(sts.n);
            sts.numparents = zeros(sts.n,1);
            sts.edges = 0;
            sts.closed = false(sts.n);
            sts.never = false(sts.n);
            sts.done = false;
            sts.remflag = false;
            sts.prevBNmap = containers.Map('KeyType','char','ValueType','logical');
            sts.cycles = false;
            sts.self = false;
        end
        
        function obj = Init(obj)
            % set up closed edges:
            % any edge from a cont node to a disc node is closed:
            for i = 1:obj.numContNode
                obj.closed(i,obj.numContNode+1:end) = ones(1,obj.numDiscNode);
            end
            % phenotype node can't have any edges going in:
            if (~obj.nophenotype)
                obj.closed(:,obj.phencol) = ones(obj.n,1);
            end
            for i = 1:obj.n
                % also close edges from a node to itself...
                obj.closed(i,i) = true;
            end

            % set up edges we can never add:
            obj.never = obj.closed;
            
            % also close edges where the max number of edges is already
            % present
            for i = 1:length(obj.bn)
                if (sum(obj.bn(:,i)) >= obj.priorPrecision.maxParents)
                    obj.closed(:,i) = ones(obj.n,1);
                end
                % also close all existing edges and the reverse arc of
                % those edges:
                for j = 1:length(obj.bn)
                    if (obj.bn(i,j))
                        obj.closed(i,j) = true;
                        obj.closed(j,i) = true;
                    end
                end
            end
        end

        function obj = SetVariables2Learn(obj, ccols2learn, dcols2learn)
            % close any edges going into nodes that are not learned (~ccols2learn)
            for i = 1:length(ccols2learn)
                if (~ccols2learn(i))
                    obj.closed(:,i) = ones(obj.n,1);
                end
            end
            for i = 1:length(dcols2learn)
                choseni = i + obj.numContNode;
                if (~dcols2learn(i))
                    obj.closed(:,choseni) = ones(obj.n,1);
                end
            end
        end

        function obj = ReOpenSearch(obj)
            % re-init this object to search again, starting from the
            % existing BN, weightmatrix, etc.
            obj.lldiff = zeros(obj.n);
            obj.llscore = zeros(obj.n);
            obj.llrecompute = true(obj.n);
            sts.done = false;
            sts.remflag = false;
            sts.prevBNmap = containers.Map('KeyType','char','ValueType','logical');
            
            % then re-init:
            obj = obj.Init();
        end
  
        function [contParentData, discParentData] = GetParentData(obj, parents)
            % parent data can be different by dynamic bayes net or static
            % bayes net
            allData = [obj.contData;obj.discData];
            contParentData = allData(parents(parents <= obj.numContNode),:);
            discParentData = allData(parents(parents > obj.numContNode),:);            
        end
        
        function childData = GetChildData(obj, child)
            allData = [obj.contData;obj.discData];
%            if (child > size(obj.contData,2))
                childData = allData(child,:);
%            else
%                childData = double(allData(child,:));
%            end
        end
        
        function [obj, evalq] = EvalOneEdge(obj, choseni, chosenj)
            % continuous data comes before discrete data; phenotype comes last in disc
            % data.
            evalq = false;

            if (~obj.self && choseni == chosenj)
                return;
            end
            % if we don't have to recompute this edge, don't bother:
            if (~obj.llrecompute(choseni,chosenj))
                return;
            end

            if (~obj.closed(choseni,chosenj) && ~obj.bn(choseni,chosenj))
                % we're adding this edge to the network
                adding = true;
            elseif (obj.backtracking && obj.bn(choseni,chosenj))
                % we're removing an existing edge
                adding = false;
            else 
                % no other conditions to worry about
                return;
            end
            %% consider making an edge from node i to j.
            % data of the child node
            childData = obj.GetChildData(chosenj);

            % initial model of node j 
            parents = find(obj.bn(:,chosenj)' ~= 0);
            [contParentData, discParentData] = obj.GetParentData(parents);
            if (chosenj <= obj.numContNode)
                % is a continuous node
                model_ini=learnLocalBN_MixToCont(contParentData,discParentData,childData,obj.priorPrecision);
            else
                model_ini=learnLocalBN_DiscToDisc(discParentData,childData,obj.priorPrecision);
            end

            % now add edge (choseni,chosenj)
            if (adding)
                parents2 = unique([parents,choseni],'legacy');
            else
                parents2 = setdiff(parents,choseni,'legacy');
            end
            [contParentData, discParentData] = obj.GetParentData(parents2);
            if (chosenj <= obj.numContNode)
                % is a continuous node
                model_addedge=learnLocalBN_MixToCont(contParentData,discParentData,childData,obj.priorPrecision);
            else
                model_addedge=learnLocalBN_DiscToDisc(discParentData,childData,obj.priorPrecision);
            end

            evalq = true;
            % update the search matrices:
            obj.lldiff(choseni,chosenj) = model_addedge.logLLH - model_ini.logLLH;
            obj.llscore(choseni,chosenj) = model_addedge.logLLH;
            obj.llrecompute(choseni,chosenj) = false;
        end
    
        function [obj, numevals] = EvalAllEdges(obj)
            numevals = 0;
            for j = 1:obj.n
                % Consider node j the potential child:
                for i = 1:obj.n
                    [obj, evalq] = obj.EvalOneEdge(i, j);
                    if (evalq)
                        numevals = numevals + 1;
                    end
                end
            end
        end

        function [obj, numevals] = EvalSelfEdges(obj)
            numevals = 0;
            if (~obj.self)
                return;
            end
            for i = 1:obj.n
                % Consider self-edge (i, i)
                [obj, evalq] = obj.EvalOneEdge(i, i);
                if (evalq)
                    numevals = numevals + 1;
                end
            end
        end

        function [obj, numParameters] = GetNumParameters(obj, choseni, chosenj, validOperation)
            tempAdjMat = obj.bn;
            if (validOperation)
                if (~obj.bn(choseni, chosenj))
                    if (~obj.bn(chosenj, choseni))
                        % Edge operation is adding (choseni, chosenj) 
                        tempAdjMat(choseni, chosenj) = 1;
                    else
                        % Edge operation is reversing (chosenj, choseni) to (choseni, chosenj)
                        tempAdjMat(chosenj, choseni) = 0;
                        tempAdjMat(choseni, chosenj) = 1;
                    end
                else
                    % Edge operation is removing (choseni, chosenj) 
                    tempAdjMat(choseni, chosenj) = 0;
                end
            end
            numParameters = 0;
            % Compute number of free parameters for each node
            for nodei = 1:obj.n
                contParents = 0;
                discParents = 0;
                discParentsConfig = 1;
                nodeParameters = 0;
                parents = find(tempAdjMat(:,nodei)' ~= 0);
                for i = 1:length(parents)
                    parenti = parents(i);
                    if (parenti <= obj.numContNode)
                        % Parent is continuous
                        contParents = contParents + 1;
                    else
                        % Parent is discrete
                        [~, discParentData] = obj.GetParentData(parenti); % data of the child node
                        parentConfig = length(unique(discParentData, 'legacy'));
                        discParentsConfig = discParentsConfig * parentConfig;
                        discParents = discParents + 1;
                    end
                end
                if (~isempty(parents))
                    if (nodei <= obj.numContNode)
                        % Node is a continuous 
                        nodeParameters = (discParentsConfig * contParents) + discParentsConfig;
                    else
                        % Node is discrete 
                        childData = obj.GetChildData(nodei); % Data of the child node
                        childConfig = length(unique(childData, 'legacy'));
                        nodeParameters = discParentsConfig * (childConfig - 1);
                    end
                end
                numParameters = numParameters + nodeParameters;
            end
        end
        
        function obj = UpdateCoefficients(obj)
            %JLM: Should be obj.numContNode but we might add functionality for discrete nodes as well
            for chosenj = 1:obj.n 
                if (chosenj > obj.numContNode)
                    continue;
                end
            
                % Data of the child node
                childData = obj.GetChildData(chosenj);

                parents = find(obj.bn(:,chosenj)' ~= 0);
                [contParentData, discParentData] = obj.GetParentData(parents);
                numDiscParent = size(discParentData,1);
                if (numDiscParent > 0)
                    continue;
                end
                % Node is a continuous without any discrete parent
                model_addedge = learnLocalBN_MixToCont(contParentData, discParentData, childData, obj.priorPrecision);
                regreCoeff = model_addedge.regreCoeff;
                for i = 2:length(regreCoeff)
                    choseni = parents(i - 1);
                    obj.weightMatrix(choseni, chosenj) = regreCoeff(i);
                    %obj.weightMatrix(choseni, chosenj) = abs(obj.weightMatrix(choseni, chosenj)) * sign(regreCoeff(i));
                end
            end
        end
        
        function [obj, success, adding] = DoEdgeAdd(obj, choseni, chosenj, lldiffchosen)        
            % debugging:
            if (obj.numparents(chosenj) > obj.priorPrecision.maxParents)
                outline = ['maxP reached;', num2str(choseni), ';', num2str(chosenj), ';', num2str(obj.lldiff(choseni, chosenj))];
                disp(outline);
                error('Error: exceeded number of parents for node %d\n', chosenj);
            end
            
            success = false;
            % we are adding/removing edge (choseni, chosenj):        
            adding = ~obj.bn(choseni,chosenj); 
            if (~adding)
                % disallow backtracking if that's the setting
                if (~obj.backtracking)
                    % may also want to update the object (closed, etc)
                    return;
                end
                if (obj.checkRepeats)
                    obj.remflag = true;
                end
            end
            % decide if this edge results in a cycle before deciding to add it:
            if (~obj.cycles && adding)
                tempadjmat = obj.bn;
                tempadjmat(choseni, chosenj) = 1;
                iscycle = dfsCycleCheck(tempadjmat, choseni);
            else
                % if we don't care if there's cycles, we can pretend we didn't
                % find one:
                iscycle = false;
            end

            if (~iscycle)
                if (adding)
                    obj.bn(choseni,chosenj) = 1;
                    obj.weightMatrix(choseni,chosenj) = lldiffchosen;
                    obj.numparents(chosenj) = obj.numparents(chosenj) + 1;
                    obj.edges = obj.edges + 1;
                    % only want to prevent a reverse edge if we're adding this edge
                    % A reverse edge is actually warranted if this edge causes a
                    % cycle.
                    if (~obj.cycles)
                        obj.closed(chosenj,choseni) = true;
                    end
                    % debugging:
                    if (obj.numparents(chosenj) > obj.priorPrecision.maxParents)
                        error('Error: exceeded number of parents for node %d\n', chosenj);
                    end
                else
                    obj.bn(choseni,chosenj) = 0;
                    obj.weightMatrix(choseni,chosenj) = 0;
                    obj.numparents(chosenj) = obj.numparents(chosenj) - 1;
                    obj.edges = obj.edges - 1;
                    % also open up the reverse edge as a possiblity again:
                    if (~obj.never(chosenj,choseni) && ~obj.cycles)
                        obj.closed(chosenj,choseni) = ...
                            obj.numparents(choseni) >= obj.priorPrecision.maxParents;
                    end
                end
                obj.lldiff(chosenj,choseni) = 0;
                obj.llscore(chosenj,choseni) = 0;
                obj.llrecompute(chosenj,choseni) = true;
                % fix up indicators
                obj.llrecompute(:,chosenj) = true(obj.n,1);
                obj.lldiff(:,chosenj) = zeros(obj.n,1);
                obj.llscore(:,chosenj) = zeros(obj.n,1);
                if (obj.numparents(chosenj) >= obj.priorPrecision.maxParents)
                    obj.closed(:,chosenj) = true(obj.n,1);
                elseif (~adding)
                    % re-open all the incoming edges, except those that
                    % already have edges:
                    obj.closed(:,chosenj) = obj.never(:,chosenj) | obj.bn(:,chosenj);
                end
                success = true;
                if (obj.remflag && obj.checkRepeats)
                    % use container.Map() for string representations of the whole BN matrix.
                    if (isKey(obj.prevBNmap,mat2str(obj.bn)))
                        obj.done = true;
                    end
                end
                if (obj.checkRepeats)
                    obj.prevBNmap(mat2str(obj.bn)) = true;
                end
            else
                outline = ['Cycle created;', num2str(choseni), ';', num2str(chosenj), ';', num2str(obj.lldiff(choseni, chosenj))];
                disp(outline);
            end
            % want to prevent this edge from being chosen again, whether or not
            % it was added or disallowed for cycle-making
            obj.lldiff(choseni,chosenj) = 0;
            obj.llscore(choseni,chosenj) = 0;
            obj.closed(choseni,chosenj) = true;
        end
        
    end
    
    
    
end

    