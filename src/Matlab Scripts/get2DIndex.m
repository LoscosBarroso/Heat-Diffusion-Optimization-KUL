
function [i, j] = get2DIndex(nr, nc, ind1d)
    i = floor((ind1d)/nc)+1;
    j = ind1d+1-(i-1)*nc;
end
