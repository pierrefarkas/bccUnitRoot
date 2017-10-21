function [mStr] = matrixSizeStr(M)
% 	[mStr] = matrixSizeStr(M)
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        		
    mStr = [num2str(u.rw(M)), 'x', num2str(u.cl(M))];
end