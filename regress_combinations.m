function [r2 indexing B_all] = regress_combinations(y,predictors)
% [r2_temp all_p indexing B_all] = regress_combinations(y,predictors)
%
%   returns r2 for each combination of predictors
%
%   Input
%       y >>> dependent variable [1 n]
%       predictors >>> matrix containing all the predictors [n 1+p]
%                      N.B.: include a column of ones!
%
%   Output:
%       r2 >>> r2 for each combination [2p-1 1]
%       indexing >>> binary matrix reporting if a predictor contributes
%                    to each partition [2p-1 p]
%       B_all >>> matrix reporting betas of each predictor in each
%                 combination [2p-1 p+1]
%
%   2018 - Paolo Papale fecit

all_p = 1:size(predictors,2)-1;
indexing = zeros([2^length(all_p)-1 length(all_p)]);
y = zscore(y);
B_all = zeros([2^length(all_p)-1 size(predictors,2)]);
r2 = zeros([2^length(all_p)-1 1]);
ind = 1;
for i = 1:length(all_p)
    clear index_temp
    index_temp = combnk(all_p,i);
      
    for j = 1:size(index_temp,1)
        clear pred_temp B predicted
        pred_temp = [ones([size(predictors,1) 1]) predictors(:,1+index_temp(j,:))];
        B = regress(y,pred_temp);
        predicted = B'*pred_temp';
        r2(ind) = (corr(y,predicted'))^2;
        indexing(ind,index_temp(j,:)) = 1;
        B_all(ind,[1 1+index_temp(j,:)]) = B;
        ind = ind+1;
    end
end

end