function stats = react_rdm_stats(rdms, category, varargin)
%REACT_RDM_STATS   Calculate averages for subject RDMs.
%
%  stats = react_rdm_stats(rdms, category, ...)
%
%  INPUTS
%  rdms : [items x items x subjects] numeric array
%      Representational dissimilarity matrices for each subject,
%      comparing pre-exposure representations to study representations.
%
%  category : [items] numeric array
%      Category of each item.
%
%  OUTPUTS
%  stats : struct
%      Contains item, within, and between fields with average
%      fisher z scores for those parts of the RDM.
%
%  OPTIONS
%  subj_ind : [subjects x items] boolean array ([])
%      Indicates items to include for each subject
%      (e.g. subsequently correct trials)

def.subj_ind = [];
opt = propval(varargin, def);

n_item = size(rdms, 1);
n_subj = size(rdms, 3);

item_mat = eye(n_item);
cat_mat = NaN(n_item);
for i = 1:n_item
    for j = 1:n_item
        if category(i) == category(j)
            cat_mat(i,j) = 1;
        else
            cat_mat(i,j) = 2;
        end
    end
end

ucat = unique(category);
n_cat = length(ucat);
stats.item = NaN(n_subj, 1);
stats.within = NaN(n_subj, 1);
stats.item_cat = NaN(n_subj, n_cat);
stats.within_cat = NaN(n_subj, n_cat);
stats.between = NaN(n_subj, 1);
for i = 1:n_subj
    if ~isempty(opt.subj_ind)
        % allow for subject indices that apply to all 180 pairs,
        % although we will only ever have 120 items in the
        % RDM. This assumes that the subj_ind matrix is sorted by
        % group, which it should be anyway to match the RDMs
        ind = opt.subj_ind(i,1:n_item);
        rdm = rdms(ind, ind, i);
        item_subj = item_mat(ind, ind);
        cat_subj = cat_mat(ind, ind);
        labels_subj = category(ind);
    else
        rdm = rdms(:,:,i);
        item_subj = item_mat;
        cat_subj = cat_mat;
        labels_subj = category;
    end
    
    stats.item(i) = mean(rdm(item_subj==1));
    stats.within(i) = mean(rdm(item_subj==0 & cat_subj==1));
    stats.between(i) = mean(rdm(cat_subj==2));
    for j = 1:n_cat
        ind = labels_subj == ucat(j);
        item_submat = item_subj(ind, ind);
        rdm_submat = rdm(ind, ind);
        stats.item_cat(i,j) = mean(rdm_submat(item_submat==1));
        stats.within_cat(i,j) = mean(rdm_submat(item_submat==0));
    end
end
