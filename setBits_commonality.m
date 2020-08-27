function indexing = setBits_commonality(col,indexing)
% out = setBits_commonality(col,input)
%   
%   Matlab implementation of R-function setBits by Kim Nimon
%
%   2018 - Paolo Papale fecit

row = 1;
val = col;
while val ~= 0
    if mod(val,2)
        indexing(row,col) = 1;
    end
    val = floor(val/2);
    row = row+1;
end
