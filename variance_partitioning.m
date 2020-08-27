function [r2 indexing] = variance_partitioning(y,predictors,varargin)
% [r2 indexing] = variance_partitioning(y,predictors,B_all)
%
%   returns r2 for each unique partition - this method is best known as
%   commonality analysis outside the neuroimaging field
%
%   Input
%       y >>> dependent variable [1 n]
%       predictors >>> matrix containing all the predictors [n p]
%                      N.B.: do not include a column of ones!
%       B_all >>> [optional] output of regress_combinations: matrix [2p-1 p+1]
%                 reporting betas of each predictor in each combination 
%                 N.B.: it is needed for crossvalidation!
%   Output:
%       r2 >>> r2 for each partion [2p 1]
%              N.B.: the first element is the total r2
%       indexing >>> binary matrix reporting if a predictor contributes
%                    to each partition [2p-1 p]
%
%   2018 - Paolo Papale fecit
%
%       heavily based on R-script 'commonalityCoefficients' by Kim Nimon
%       Reference: Nimon et al., 2008, Behavior Research Methods.
%       FYI: if you have 2 or 3 models and want to plott a Venn diagram, the 
%            r2(2:end) perfectly fits the input format requested by venn.m 
%            by Darix on MathExchange
%

%%% preprocess data
all_p = 1:size(predictors,2);
y = zscore(y);
predictors = [ones([size(predictors,1) 1]) zscore(predictors)];

%%% compute the total r2
B = regress(y,predictors);
predicted = B'*predictors';
r2(1) = (corr(y,predicted'))^2;

%%% compute the r2 for each combination of the models
if isempty(varargin)
    [r2_temp indexing] = regress_combinations(y,predictors);
else
    B_all = varargin{1};
    indexing = zeros([2^length(all_p)-1 length(all_p)]);
    r2_temp = zeros([2^length(all_p)-1 1]);
    ind = 1;
    for i = 1:length(all_p)
        clear index_temp
        index_temp = combnk(all_p,i);
        
        for j = 1:size(index_temp,1)
            clear pred_temp B predicted
            indexing(ind,index_temp(j,:)) = 1;
            pred_temp = predictors(:,[1 1+index_temp(j,:)]);
            B = B_all(ind,[1 1+index_temp(j,:)])';
            predicted_t = B'*pred_temp';
            r2_temp(ind) = (corr(y,predicted_t'))^2;
            ind = ind+1;
        end
    end
end

%%% define wich combination contribute to each partition - as in Nimon R-script
%%% as in Nimon R-script

% format data
for i = 0:length(all_p)-1
    all_p_ind(i+1) = 2^i;
end

indexing_nimon = zeros(size(indexing'));
for i = 1:size(indexing,1)
    indexing_nimon = setBits_commonality(i,indexing_nimon);
end

for i = 1:size(indexing_nimon,2)
    row_temp = indexing_nimon(:,i);
    index_good(i) = find(ismember(indexing,row_temp','rows'));
end

r2_temp_nimon = r2_temp(index_good);

for i = 1:size(indexing_nimon,2)
    clear el_temp list_temp
    el_temp = indexing_nimon(1,i);
    if el_temp == 1
        list_temp = [0 -all_p_ind(1)];
    else
        list_temp = all_p_ind(1);
    end
    for j = 2:length(all_p_ind)
        clear el_temp a_list list_temp2
        el_temp = indexing_nimon(j,i);
        if el_temp == 1
            a_list = list_temp;
            list_temp2 = genList_commonality(list_temp,-all_p_ind(j));
            list_temp = [a_list list_temp2];
        else
            list_temp = genList_commonality(list_temp,all_p_ind(j));
        end
    end
    list_temp = list_temp*-1;
    list_comps{i} = list_temp;
end

for i = 1:size(indexing_nimon,2)
    clear r2_list_temp list_l_temp r2_sum_temp
    r2_list_temp = list_comps{i};
    list_l_temp = length(r2_list_temp);
    r2_sum_temp = 0;
    for j = 1:list_l_temp
        clear index_temp index_abs_temp el_temp
        index_temp = r2_list_temp(j);
        index_abs_temp = abs(index_temp);
        if index_abs_temp ~= 0
            el_temp = r2_temp_nimon(index_abs_temp);
            if (index_temp < 0)
                el_temp = el_temp*-1;
            end
            r2_sum_temp = r2_sum_temp+el_temp;
        end
    end
    r2_nimon(i) = r2_sum_temp;
end

r2_good(index_good) = r2_nimon;
r2 = [r2 r2_good];

end