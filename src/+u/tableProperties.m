function [res] = tableProperties(T)
    res = properties(T)';
    res = setdiff(res, {'Properties'});
end