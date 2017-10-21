function v = isAllMember(MSmaller, MBigger)
% v = isAllMember(MSmaller, MBigger)
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        		
    if sum(ismember(MSmaller, MBigger)) == numel(MSmaller)
        v = true;
    else
        v = false;
    end
end