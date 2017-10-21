function [str] = cell2Str(cellIn)
% [str] = cell2Str(cellIn)
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        	
    str = '';
    for i = 1:numel(cellIn)
        if isnumeric(cellIn{i})
            str = [str, ' ', num2str(cellIn{i})];
        else
            str = [str, ' ', cellIn{i}]; %#ok<AGROW>
        end
    end
end