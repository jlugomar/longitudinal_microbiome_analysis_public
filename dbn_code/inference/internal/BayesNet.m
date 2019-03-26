classdef BayesNet
%
% A BayesNet structure represents a conditional Guassian Bayesian Network. 
%
%     nodes       % array of NODES representing the BN
%     title       % a string that names the BN
%     adjmat      % the adjacency matrix representing the BN
%     weightMatrix %the weight matrix representing the LLH improvement for each edge
%     mb          % indices into nodes that are the markov blanket of the network
%     isMB        % boolean: indicator of if the NODES represent a Markov Blanket
%     disc        % boolean array: is each column discrete?
%     data        % actual data that the BN was built on        
%     cols        % column names for the data.
%     pheno       % string matching one of the cols
%     priorPrecision % structure with BN bayesian hyperparameters:
%       priorPrecision.nu; % prior sample size for prior variance estimate
%       priorPrecision.sigma2; % prior variance estimate
%       priorPrecision.alpha; % prior sample size for discrete nodes
%       priorPrecision.maxParents; % maximum number of parents for any node
%               in the network
%     discvals    % cell array of cells of values for discrete variables.
%                   Useful for when discrete variables have values that do
%                   not occur in all datasets.
%     tree        % array of ClusterSetTree objects representing the 
%                 % parameters learned for the bayes net
%
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

    properties
        nodes       % array of NODES representing the BN
        title       % a string that names the BN
        adjmat      % the adjacency matrix representing the BN
        weightMatrix %the weight matrix representing the LLH improvement for each edge
        mb          % indices into nodes that are the markov blanket of the network
        isMB        % boolean: indicator of if the NODES represent a Markov Blanket
        disc        % boolean array: is each column discrete?
        data        % actual data that the BN was built on        
        cols        % column names for the data.
        pheno       % string matching one of the cols
        priorPrecision % structure for parameters of bayesian learning
        discvals    % optional inclusion of discrete value arrays
        tree        % array of ClusterSetTree objects representing the 
                    % parameters learned for the bayes net
        bootStrapMatrix %BootsScore edge matrix

    end
    
    methods
        % constructor
        function bn = BayesNet(nodes, title, adjmat, weightMatrix, mb, ...
                isMB, disc, data, cols, pheno, priorPrecision, discvals, tree,bootStrapMatrix)
            if (nargin < 14)
                bn.bootStrapMatrix = [];
            else
                bn.bootStrapMatrix = bootStrapMatrix;
            end
            if (nargin < 13)
                bn.tree = [];
            else
                bn.tree = tree;
            end
            if (nargin < 12)
                bn.discvals = {};
            else
                bn.discvals = discvals;
            end
            if (nargin < 11)
                bn.priorPrecision.sigma2 = 1;
                bn.priorPrecision.nu = 10;
                bn.priorPrecision.alpha = 10;
                bn.priorPrecision.maxParents = 2;
            else
                bn.priorPrecision = priorPrecision;
            end
            if (nargin < 10)
                bn.pheno = '';
            else
                bn.pheno = pheno;
            end
            if (nargin < 9)
                bn.cols = {};
            else
                bn.cols = cols;
            end
            if (nargin < 8)
                bn.data = [];
            else
                bn.data = data;
            end
            if (nargin < 7)
                bn.disc = [];
            else
                bn.disc = disc;
            end
            if (nargin < 6)
                bn.isMB = false;
            else
                bn.isMB = isMB;
            end
            if (nargin < 5)
                bn.mb = [];
            else
                bn.mb = mb;
            end
            if (nargin < 4)
                bn.weightMatrix = [];
            else
                bn.weightMatrix = weightMatrix;
            end
            if (nargin < 3)
                bn.adjmat = [];
            else
                bn.adjmat = adjmat;
            end
            if (nargin < 2)
                bn.title = '';
            else
                bn.title = title;
            end
            if (nargin < 1)
                bn.nodes = [];
            else
                bn.nodes = nodes;
            end
            if (isempty(bn.disc) && ~isempty(bn.data))
                bn.disc = IsDiscrete(bn.data,5);
            end
            if (isempty(bn.nodes) && ~isempty(bn.disc) && ~isempty(bn.cols) && ~isempty(bn.adjmat))
                bn.nodes = BNfromAdjMat(bn.adjmat, bn.disc, bn.cols);
            end
            if (isempty(bn.adjmat) && ~isempty(bn.nodes) && ~isempty(bn.cols))
                bn.nodes = AddColIndstoNodes(bn.nodes,bn.cols);
                bn.nodes = nodetree(bn.nodes, bn.data, bn.cols);
                bn.adjmat = AdjMatfromBN(bn.nodes, bn.cols);
            end
            if (isempty(bn.discvals) && ~isempty(bn.disc))
                bn.discvals = ComputeDiscVals(bn);
            end
        end
        
        % convert the BayesNet into a Markov Blanket for the BN:
        function obj = MakeIntoMB(obj, newpheno)
            if (nargin > 1)
                obj.pheno = newpheno;
            end
            targetNode = strcmp(obj.pheno, obj.cols);

            % get nodes under Markov blanket
            parentNodes = find(obj.adjmat(:,targetNode));
            childNodes = find(obj.adjmat(targetNode,:));
            parentOfChildNodes = [];
            for c = childNodes
                parentOfChildNodes = union(parentOfChildNodes, find(obj.adjmat(:,c)),'legacy');    
            end

            % get all nodes under the Markov blanket
            inds = find(targetNode);
            inds = union(inds,parentNodes,'legacy');
            inds = union(inds,childNodes,'legacy');
            inds = union(inds,parentOfChildNodes,'legacy');
            
            % now whack everything down to just the inds:
            obj.isMB = true;
            obj.mb = inds;
            if (~isempty(obj.data))
                obj.data = obj.data(:,inds);
            end
            if (~isempty(obj.adjmat))
                obj.adjmat = obj.adjmat(inds,inds);
            end
            if (~isempty(obj.weightMatrix))
                obj.weightMatrix = obj.weightMatrix(inds,inds);
            end
            if (~isempty(obj.discvals))
                obj.discvals = obj.discvals(inds);
            end
            % obj.priorPrecision doesn't change
            obj.disc = obj.disc(inds);
            obj.cols = obj.cols(inds);
            % and recompute the NODES structure:
            obj.nodes = BNfromAdjMat(obj.adjmat, obj.disc, obj.cols);
            % have to kill the tree here:
            obj.tree = [];
        end
        
        function obj = InflateFromMBtoFull(obj, fullcols, fulldata, fulldisc, fulldiscvals)
            % find the mapping from the BayesNets columns to the full
            % columns:
            inds = 1:length(fullcols);
            map = zeros(1,length(obj.cols));
            for i = 1:length(obj.cols)
                matches = strcmp(obj.cols(i),fullcols);
                map(i) = inds(matches);
            end
            obj.cols = fullcols;
            newadjmat = zeros(length(fullcols));
            newadjmat(map,map) = obj.adjmat;
            obj.adjmat = newadjmat;
            newweightmat = zeros(length(fullcols));
            newweightmat(map,map) = obj.weightMatrix;
            obj.weightMatrix = newweightmat;
            obj.isMB = false;
            obj.mb = map;
            if (nargin >= 3)
                obj.data = fulldata;
            else
                obj.data = [];
            end
            if (nargin >= 4)
                obj.disc = fulldisc;
            else
                obj.disc = [];
            end
            if (nargin >= 5)
                obj.discvals = fulldiscvals;
            else
                obj.discvals = {};
            end
            % rework nodes array:
            obj.nodes = AddColIndstoNodes(obj.nodes,fullcols);
            % test this: 
            obj.nodes = nodetree(obj.nodes, obj.data, fullcols);
            % kill tree here: prediction no longer makes sense
            obj.tree = [];
        end
        
        function obj = ReorderByNames(obj, colnames)
            % take an input list of column names and reorder the elements
            % of this BayesNet to match that ordering of columns.
            % colnames must be a subset of the existing columns of the
            % BayesNet object.
            map = zeros(size(colnames));
            inds = 1:length(colnames);
            for i = 1:length(colnames)
                matches = strcmp(colnames(i),obj.cols);
                if (sum(matches) > 1)
                    warning('column names not unique!\n');
                end
                if (sum(matches) == 0)
                    warning('unable to find column %s!\n',colnames(i));
                end
                map(i) = inds(matches);
            end
            % map is a mapping of indices of existing columns to the new
            % columns.  
            obj.cols = obj.cols(map);
            obj.adjmat = obj.adjmat(map,map);
            obj.weightMatrix = obj.weightMatrix(map,map);
            % obj.isMB doesn't change
            if (~isempty(obj.mb))
                obj.mb = obj.mb(map);
            end
            if (~isempty(obj.data))
                obj.data = obj.data(:,map);
            end
            if (~isempty(obj.disc))
                obj.disc = obj.disc(map);
            end
            if (~isempty(obj.discvals))
                obj.discvals = obj.discvals(map);
            end
            % rework nodes array:
            obj.nodes = AddColIndstoNodes(obj.nodes,obj.cols);
            % test this: 
            obj.nodes = nodetree(obj.nodes, obj.data, obj.cols);
            % and kill the tree:
            obj.tree = [];
        end
        
        function obj = ReplaceData(obj, newdata, newcols)
            % take a new dataset and reorder it to match the columns in
            % this BayesNet 
            % keep the current COLS intact.
            % use input NEWCOLS to identify columns in the NEWDATA
            
            obj.discvals = ComputeDiscVals(obj);
            
            map = zeros(size(obj.cols));
            inds = 1:length(newcols);
            for i = 1:length(obj.cols)
                matches = strcmp(obj.cols(i),newcols);
                if (sum(matches) > 0)
                    map(i) = inds(matches);
                end
            end
            obj.data = newdata(:,map);
            
            % check that new discrete columns don't exceed the values of
            % old discrete columns:
            for i = 1:length(obj.cols)
                if (obj.disc(i))
                    % data(:,i) is the new data, discvals{i} are the values
                    % of the discrete data for previous data
                    newvals = unique(obj.data(:,i),'legacy');
                    rarevals = setdiff(newvals,cell2mat(obj.discvals{i}),'legacy');
                    count = zeros(1,length(obj.discvals{i}));
                    for j = 1:length(obj.discvals{i})
                        count(j) = sum(obj.data(:,i) == obj.discvals{i}{j});
                    end
                    [~,rarestind] = min(count);
                    for j = 1:length(rarevals)
                        % re-assign these to the least likely occurrences in
                        % data
                        downval = obj.discvals{i}{rarestind};
                        obj.data(obj.data(:,i) == rarevals(j),i) = downval;
                    end
                end
            end
        end
        
        function obj = ConvertToDBN(obj,newdata)
            % assume all edges are from T0 to Tn;
            % all edges go from previous time point of a data variable to
            % the current time point of a data variable.
            
            % first take each variable and suffix it with '_t0' or '_tn'
            cols0 = obj.cols;
            colsn = obj.cols;
            for i = 1:length(obj.cols)
                cols0{i} = [obj.cols{i},'_t0'];
                colsn{i} = [obj.cols{i},'_tn'];
            end
            
            % then convert the adjacency matrix to one with twice as many
            % variables
            n = length(obj.cols);
            adjmatnew = zeros(n * 2);
            weightMatrixnew = zeros(n * 2);
            for i = 1:length(obj.cols)
                for j = 1:length(obj.cols)
                    if (obj.adjmat(i,j))
                        adjmatnew(i,n + j) = obj.adjmat(i,j);
                        weightMatrixnew(i,n+j) = obj.weightMatrix(i,j);
                    end                    
                end
            end
            
            % fix up disc : 
            discnew = [obj.disc, obj.disc];
            discvalsnew = {obj.discvals{:}, obj.discvals{:}};

            % set data to be split data by timeslice:
            obj.data = newdata;

            % then fix up everything else to have twice as many variables:
            obj.disc = discnew;
            obj.discvals = discvalsnew;
            obj.weightMatrix = weightMatrixnew;
            obj.adjmat = adjmatnew;
            obj.cols = {cols0{:},colsn{:}};
            % rework nodes array:
            obj.nodes = BNfromAdjMat(obj.adjmat, obj.disc, obj.cols);
            obj.nodes = AddColIndstoNodes(obj.nodes,obj.cols);
            obj.nodes = nodetree(obj.nodes, obj.data, obj.cols,obj.discvals,obj.disc);
            % and kill the tree:
            obj.tree = [];
        end
               
        function obj = AddDiscVals(obj, newdiscvals)
            if (length(newdiscvals) == length(obj.cols))
                obj.discvals = newdiscvals;
            else
                if (~isempty(obj.mb))
                    obj.discvals = newdiscvals(obj.mb);
                else
                    error('Discvals must match columns of BayesNet');
                end
            end
        end
        
        function discvals = ComputeDiscVals(obj)
            if (isempty(obj.disc))
                obj.disc = IsDiscrete(obj.data);
            end
            discvals = cell(1,length(obj.cols));
            for i = 1:length(obj.cols)
                if (obj.disc(i))
                    discvals{i} = num2cell(unique(obj.data(:,i),'legacy'));
                end
            end
            
        end
        
        function datacol = GetPhenoCol(obj, newpheno)
            if (nargin < 2)
                newpheno = obj.pheno;
            end
            match = strcmp(newpheno, obj.cols);
            datacol = obj.data(:,match);
        end
        
        function WriteToSIF(obj, fname)
            if (nargin < 2)
                fname = obj.title;
            end
            SIFOutputBNforCytoscape(obj.adjmat,obj.cols,fname);
        end
        
         function WriteToGML(obj, fname)
            if (nargin < 2)
                fname = obj.title;
            end
            if size(obj.bootStrapMatrix) >0  attrNames = {'bootScore'};, else attrNames = {};, end
            GMLOutputBayesNet(obj, fname,{},{},{obj.bootStrapMatrix},attrNames);
        end
        
        function WriteToTGF(obj, fname)
            if (nargin < 2)
                fname = obj.title;
            end
            outputBayesNetGraph(obj.adjmat,obj.cols, fname);
        end
        
        function obj = RebuildNodes(obj)
        	obj.nodes = BNfromAdjMat(obj.adjmat, obj.disc, obj.cols);
        end
        
        function obj = RebuildAdjmat(obj)
           obj.adjmat = AdjMatfromBN(obj.nodes, obj.cols); 
        end

        function obj = LearnParams(obj)
            obj = LearnParamsBN(obj);
        end
        
        function [acc, p, z] = Predict(obj, verbose)
            if (nargin < 2)
                verbose = false;
            end
            [acc, p, z] = PredictPheno(obj, verbose);
        end
    end

end
    