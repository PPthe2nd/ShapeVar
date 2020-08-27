function newlist = genList_commonality(list,value)
% newlist = genList_commonality(list,value)
%   
%   Matlab implementation of R-function genList by Kim Nimon
%
%   2018 - Paolo Papale fecit

for i = 1:length(list)
    newlist(i) = abs(list(i))+abs(value);
    if (((list(i) < 0) && (value >= 0)) || ((list(i) >= 0) && (value < 0)))
        newlist(i) = newlist(i)*-1;
    end
end

end