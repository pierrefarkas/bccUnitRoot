classdef(Abstract) CWithTest < handle
% For example to call all methods starting with mTest
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
    
    methods
        
        function res = mRunMethodsStarting(obj, strStart)
            M = methods(obj);
            % res = cell(0,4); % toDo: update with fullObjName
            res = cell(0, 2);
            rw  = 0;
            for i = 1:numel(M)
                if length(M) > (numel(strStart)-1)
                    if strfind(M{i}, strStart) == 1
                        disp(' ');
                        disp(['Running ', M{i}])
                        rw = rw + 1;
                        res{rw, 1} = M{i};
                        res{rw, 2} = obj.(M{i});
                    end
                end
            end        
        end
        
    end

end