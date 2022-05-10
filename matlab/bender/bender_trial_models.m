function model = bender_trial_models(model_names, fnames, subjects)
%BENDER_TRIAL_MODELS   Load trial-level semantic models for all subjects.
%
%  model = bender_trial_models(model_names, fnames, subjects)
%
%  INPUTS
%  model_names - cell array of strings
%      Trial-level models to load (e.g., {'a' 'b' 'c'}).
%
%  fnames - cell array of strings
%      Field names for output struct.
%
%  subjects - numeric array.
%      All subject numbers.
%
%  OUTPUTS
%  model - struct
%      Struct with a field for each model. Each will contain an
%      [items x items x subjects] matrix of dissimilarity values.

model = struct();
for i = 1:length(model_names)
    rdm = bender_load_rdms('models3', 'wiki_w2v', model_names{i}, ...
                           'model', subjects);
    model.(fnames{i}) = rdm(1:120,1:120,:);
end
