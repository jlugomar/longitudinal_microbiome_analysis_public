classdef PCStateTrackerSearch < StateTrackerSearch
% class for implementing network search algorithms that know which
% states need to be re-evalutated and which don't on Pheno-Centric Networks
% Keeps lists of all possible edges and whether or not they are: legal,
% illegal, possible, and need to be reevaluated.
%
% BETA, UNFINISHED
%
% Michael McGeachie (c) 2015. MIT license. See cgbayesnets_license.txt.
    
    methods
        
        % constructor
        function pcsts = PCStateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats)
            if (nargin < 5)
                backtracking = false;
            end
            if (nargin < 6)
                bf_thresh = 0;
            end
            if (nargin < 7)
                nophenotype = false;
            end
            if (nargin < 8)
                checkRepeats = true;
            end
            pcsts@StateTrackerSearch(contData, discData, priorPrecision, ...
                phencol, backtracking, bf_thresh, nophenotype, checkRepeats);
        end

        function obj = Init(obj)
            % Should be possible to do PC search this way by starting with
            % everything closed except edges from the phenotype; then
            % opening edges into those nodes.
            
            obj = Init@StateTrackerSearch(obj);
            
            % set closed: everything except edges from phenotype
            colinds = true(1,obj.n);
            colinds(obj.phencol) = false;
            obj.closed(colinds,:) = true(obj.n-1,obj.n);
        end

        function [obj, success] = DoEdgeAdd(obj, choseni, chosenj, lldiffchosen)
            [obj, success, adding] = DoEdgeAdd@StateTrackerSearch(obj, choseni, chosenj, lldiffchosen);
            
            if (success)
                if (adding)
                    % nothing if the node has >= maxparents already
                    if (obj.numparents(chosenj) >= obj.priorPrecision.maxParents)
                        % close everything going into this edge:
                        obj.closed(:,chosenj) = true(obj.n,1);
                        return;
                    end
                    % if we added that edge, now open up edges into that one:
                    % 1) open all edges into chosenj except:
                    closedvec = true(obj.n,1);
                    % 2) edges that already exist
                    closedvec(obj.bn(:,chosenj)) = false;
                    % 3) edges who's reverse already exists
                    closedvec(obj.bn(chosenj,:)') = false;

                    obj.closed(closedvec,chosenj) = false;
                    % also for an incoming edge being added, we need to
                    % update what's computed:
                    obj.llrecompute(closedvec,chosenj) = true;
                    
                else
                    % removed the edge;
                    % check if chosenj still has edge from phenotype
                    % or an edge into a child from the phenotype
                    phenoconnection = obj.bn(obj.phencol, chosenj);
                    children = obj.bn(chosenj,:);
                    phenchild = sum(obj.bn(obj.phencol, children)) > 0;
                    phenoconnection = phenoconnection | phenchild;
                    
                    if (phenoconnection)
                        % don't remove this node from further consideration
                    else
                        % close everything going INTO this node
                    end
                        
                end
                
            end
            
        end
    end
end
