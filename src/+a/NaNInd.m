function [de,dm,y0ind] = NaNInd_f(Y)
% [de,dm,y0ind] = NaNInd_f(Y)
%  Y: data which has missing values
%  output:
        % de: 1 where data exists, 0 otherwise
        % dm: 0 where data exists, 1 otherwise
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
  
    % init
    [rw,cl] = size(Y);
    Yind = repmat([1:1:rw]',1,cl);
    de = zeros(rw,cl);
    dm = zeros(rw,cl);
    %
    for i=1:cl_f(Y)
      dm(:,i) = isnan(Y(:,i));  
    end
    de = ~dm;
    % first valid observation
    tmp1          = Yind .* de;
    tmp1(tmp1==0) = rw+1;
    y0ind         = min(tmp1);
    y0ind(y0ind>rw) = rw;
end