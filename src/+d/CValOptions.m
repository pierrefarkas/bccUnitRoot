classdef CValOptions < handle
% to store a concrete value as well as the list of options
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
    
    properties(SetAccess = immutable)
        pRestrictions   % restrictions
    end
        
    properties
        pV              % values
        pOptions        % options
    end
    
    methods
        
        function cValOptions = CValOptions(v)
            if nargin > 0
                cValOptions.pRestrictions = v;
            end
        end
        
        function set.pV(cValOptions, v)
            cValOptions.pV = cValOptions.mCheckVal(v);
        end
            
        function set.pOptions(cValOptions, v)    
            cValOptions.pOptions = cValOptions.mCheckVal(v);
        end
        
    end
    
    methods(Access = protected)
        
        function v = mCheckVal(cValOptions, v)
            if isempty(cValOptions.pRestrictions) == 0
                if u.isAllMember(v, cValOptions.pRestrictions) == 0
                    if isnumeric(v)
                        error([num2str(v), ' is not in ', num2str(cValOptions.pRestrictions)]);
                    else
                        error([v, ' is not in ', cValOptions.pRestrictions]);
                    end
                end
            end
        end
        
    end
    
end