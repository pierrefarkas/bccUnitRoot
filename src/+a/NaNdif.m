function [YD] = NaNdif(Y)
% [YD] = NaNdif(Y)
% calculates difference in case of missing data
% Y is the incoming data matrix
% some observation may be missing, NaN
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com

    [rw,cl] = size(Y);
    YD = zeros(rw,cl);
    [~, ~, y0ind] = a.NaNInd(Y);
    % disp('here')
    for i=1:cl
        d1 = Y( y0ind(i), i);
        for j = y0ind(i):(rw-1)
            if isnan( Y(j+1,i) ) == 1
                % missing data go to the next
            else
                YD(j+1,i) = Y(j+1,i) - d1;
                d1 = Y(j+1,i);
            end   
        end
    end

end