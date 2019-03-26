classdef LPPotential
% This class defines the structure for an LP-POTENTIAL 
%
% MEMBERS: 
% lppot.head : single node index into NODES
% lppot.tail : array of conditioning node indices (in NODES)
% lppot.params : parallel array of weights for linear combination of TAIL
%               variables.
% lppot.const : weight for the empty parameter (constant)
% lppot.sigma : parameter representing SIGMA SQUARED for the distribution
% lppot.conditionvars : array of indices to discrete vars that condition this LPPot
% lppot.conditionvalinds : arracy of indices to VALUES of those same VARs that
%                          condition this LPPot
%
% DEBUGGING status: throws error on zero-variance evidence.  this handling
% could be improved.
%
% Copyright Michael McGeachie, 2010.  MIT license. See cgbayesnets_license.txt.

    properties
        head % index into NODES that indicates what conditional ...
             % distribution this represents
        tail % array of indices into NODES representing conditioning variables
        params % weights assigned to elements of tail;
        const % weight assigned to a hypothetical constant (1) tail element
        sigma % the std deviation of the distribution
        conditionvars % array of indices to discrete vars that condition this LPPot
        conditionvalinds % arracy of indices to VALUES of those same VARs that
                         % condition this LPPot
    end
    
    methods
        % constructor
        function lppot = LPPotential(head, tail, params, const, sigma, vars, vals)
            if (nargin > 0)
                lppot.head = head;
                lppot.tail = tail;
                lppot.params = params;
                lppot.const = const;
                lppot.sigma = sigma;
                lppot.conditionvars = vars;
                lppot.conditionvalinds = vals;
                if (length(lppot.params) ~= length(lppot.tail))
                    error('missing parameters for lp potential dependencies!');
                end
            else
                lppot.head = 0;
                lppot.tail = [];
                lppot.params = [];
                lppot.const = 0;
                lppot.sigma = 0;
                lppot.conditionvars = [];
                lppot.conditionvalinds = [];
            end
        end
        
        % enter evidence into the LPPotential that continuous variables VAR have
        % obtained values VALS
        function obj = AddEvidence(obj, vars, vals)
            % find the vars to update
            varinds = zeros(length(vars));
            for i = 1: length(vars)
                for j = 1: length(obj.tail)
                    if (vars(i) == obj.tail(j))
                        varinds(i) = j;
                    break;
                    end
                end
            end
            if (sum(varinds) == 0)
                return;
                % Evidence on irrelevant vars just doesn't matter.  
                %error('Adding Irrelevant Evidence: LP Potential doesn''t depend on VAR');
            end
            obj.const = obj.const + obj.params(varinds) * vals';
            obj.params(varinds) = [];
            obj.tail(varinds) = [];
            if (length(obj.params) ~= length(obj.tail))
                error('missing parameters for lp potential dependencies!');
            end
        end
        
        % enter evidence on the HEAD var when the lppot has an empty tail
        function [w,cdf] = MakeWeight(obj, var, varval)
            % [w,cdf] = MakeWeight(obj, var, varval)
            %
            % OUTPUT:
            %   CDF: the (log) probability of observing a value more
            %   extreme for variable VAR than VARVAL
            %
            if (obj.head ~= var || ~isempty(obj.tail))
                error('LP Potential should be unconditional when making weight.');
            end
            if (obj.sigma == 0)
                % in this case the potential is a deterministic function of 
                % its variables, so the observed outcome is the only possible
                % outcome, and has likelihood 1:
                w = 0; % log (1) = 0.
                cdf = 0;
            else
                % previous non-log weight computation:
                % w = exp( -1 * (varval - obj.const)^2 / (2 * obj.sigma)) / sqrt(2 * pi * obj.sigma);
                % new log-weight computation:
                % the weight is NOT a probability, its the value of the PDF:
                w = -1 * (varval - obj.const)^2 / (2 * obj.sigma) - ...
                    0.5 * log(2 * pi * obj.sigma);

                % also use the CDF to get the probability of a more extreme
                % value :
                % flip varval so it's on the lower end of the distribution:
                if (varval > obj.const)
                    varval = obj.const - (varval - obj.const);
                end
                % do a two-sided test here (add log(2) to multiply by 2):
                cdf = log(2)  + log(normcdf(varval, obj.const, obj.sigma));
                if (isinf(cdf))
                    % use crude approximation to log CDF based on sigma:
                    z = abs((varval - obj.const)/obj.sigma);
                    cdf = -1 * z^2/2;
                end
            end
        end
            
        
    end
end
