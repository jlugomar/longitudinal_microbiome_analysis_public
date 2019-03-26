classdef TSDBNStateTrackerSearch < StateTrackerSearch
% class for implementing network search algorithms that know which
% states need to be re-evalutated and which don't on Two-Stage Dynamic BNs
% Keeps lists of all possible edges and whether or not they are: legal,
% illegal, possible, and need to be reevaluated.
%
% Michael McGeachie (c) 2014. MIT license. See cgbayesnets_license.txt.
    
    properties
        timeseries      % dynamic bayes nets
        t0cont          % dynamic bayes nets  
        tncont          % dynamic bayes nets  
        t0disc          % dynamic bayes nets  
        tndisc          % dynamic bayes nets  
    end
    
    methods
        
        % constructor
        function dbnsts = TSDBNStateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats)
            if (nargin < 5)
                backtracking = false;
            end
            if (nargin < 6)
                bf_thresh = 0;
            end
            if (nargin < 7)
                % change default here for dynamic bayes net:
                nophenotype = true;
            end
            if (nargin < 8)
                checkRepeats = true;
            end
            dbnsts@StateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats);
            dbnsts.cycles = false;
            dbnsts.self = false;
            dbnsts.timeseries = true;
            dbnsts.t0cont = [];
            dbnsts.tncont = [];
            dbnsts.t0disc = [];
            dbnsts.tndisc = [];
        end

        function obj = SetTimeSeries(obj, t0cont, t0disc, tncont, tndisc)
            % call after dbnsts.Init()
        	% this sets up the search to be a dynamic bayesian network for
        	% time series data
            obj.timeseries = true;
            
            % set up the two datasets: t0 dataset and tn dataset
            obj.t0cont = t0cont;
            obj.tncont = tncont;
            obj.t0disc = t0disc;
            obj.tndisc = tndisc;
            
            % dynamic bayes nets don't care if there's cycles
            obj.cycles = true;
            
            % dynamic bayes nets don't care if there's self-edges
            obj.self = true;
            for i = 1:obj.n
                % re-open edges from a node to itself...
                obj.closed(i,i) = false;
            end
            obj.never = obj.closed;
            
            % also re-open reverse edges of existing edges 
            for i = 1:length(obj.bn)
                for j = 1:length(obj.bn)
                    if (obj.bn(i,j))
                        obj.closed(i,j) = true;
                        obj.closed(j,i) = false;
                    end
                end
            end
        end
                
        function [contParentData, discParentData] = GetParentData(obj, parents)
            % override base class function:
            % parent data can be different by dynamic bayes net or static
            % bayes net
            allData = [obj.t0cont;obj.t0disc];
            contParentData = allData(parents(parents <= obj.numContNode),:);
            discParentData = allData(parents(parents > obj.numContNode),:);
        end
        
        function childData = GetChildData(obj, child)
            % override base class function:
            % child data can be different by dynamic bayes net or static
            % bayes net
            allData = [obj.tncont;obj.tndisc];
            childData = allData(child,:);       
        end
        
    end
end
