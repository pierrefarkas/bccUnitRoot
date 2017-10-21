function [ r ] = randBetween(min, max, rw, cl)
% [ r ] = randBetween(min, max, rw, cl)
% generate uniform random number between min and max, rw, cl 
% \copyright Peter Farkas
% peter.farkas.dr@gmail.com        
    r = repmat(min,rw,cl) + (max-min) * rand(rw,cl);
end
