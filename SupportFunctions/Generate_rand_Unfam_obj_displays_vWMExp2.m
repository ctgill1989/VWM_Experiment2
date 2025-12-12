function [Unfamiliar_display_objs]=Generate_rand_Unfam_obj_displays_vWMExp2(TF_pairs_unfam)
%By CTGill
%Last updated: 5/17/22
%
% This Version of Generate_rand_obj_displays is used in vWM_Experiment2 and
% ensures that each target location is selected an equal number of times for
% each condition.
%
% This versiion of Generate_rand_obj_displays is the functional equivalent
% of Generate_rand_obj_displays_V2 used in vWM_Experiment_V1 (i.e.
% Experiment1)
%
% Function generates sets of real-world object stimuli for each trial for
% the Unfamiliar Object Conidtion. No images are displayed more than once.
% To find the maximally dissimilar Target-Foil pairs I extradted the
% features of all objects in my Familiar and Unfamiliar image sets
% from the top max-pooling layer of ImageNet-pretrained VGG16 model and
% then computed similarity metric between all pairs of objects within each
% set using Cosine Similarity (i.e. length-normalized dot product).
% Subsequently I chose the 80 most dissimilar pairs from each set to use as
% TF pairs (see CosineSimilarity_vWM_Exp2.mat). All Target-Foil pairs are
% categorically, semantically, and perceptually dissimilar.
%
%
% INPUTS
% TF_pairs_unfam     --> a struct containing both objects in an unfamiliar
%                        pair and the cosine similarity measure.
%
% Num_stim_per_trial --> The number of stimuli in each display of the
%                        encoding phase of the vWM task. (e.g. 6)
%
%
% Num_obj_trials     --> The total number of real-world-object trials per familiarity condition.
%
% OUTPUTS
%
% display_objs       --> a struct containing obj display info.
%                        display_objs.Target            --> contains the Target object filename.
%                        display_objs.Target_Cat        --> contains the Object category of the Target object.
%                        display_objs.Foil              --> contains the Foil object filename.
%                        display_objs.Foil_Cat          --> contains the Object category of the Foil object.
%                        display_objs.Study_stimuli     --> contains filename of all stim in study display.
%                        display_objs.Target_study_pos  --> contains the Target object position in study display (Pos 1 is always at the 9'oclock position and subseqent positions go clockwise).
%                        display_objs.TF_Sim            --> contains Cosine Similarity measure of the Target-Foil pair   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_stim_per_trial = 6;
Num_obj_trials = 80;

%% Load Stim Files
%Set path to stim folders
Unfamiliar_object_StimFolder = '/Users/christophergill/Documents/MATLAB/vWM_Experiment2/Stimuli/Objects/Unfamiliar';

%Find Object Image Categories
temp = dir(Unfamiliar_object_StimFolder);
categories = {temp(~ismember({temp.name},{'.','..','.DS_Store'})).name}';
if strcmp(Unfamiliar_object_StimFolder(1),'C')
    categories = categories(2:end);
end

Unfamiliar_imgFiles = [];
for cat = 1:length(categories)
    filePattern = fullfile([Unfamiliar_object_StimFolder filesep categories{cat}], '*.jpg');
    Unfamiliar_imgFiles.(categories{cat}) = dir(filePattern);
    Unfamiliar_imgFiles.(categories{cat}) = Unfamiliar_imgFiles.(categories{cat})(~startsWith({Unfamiliar_imgFiles.(categories{cat}).name},'.'));
end

repeat_gate=1;
while repeat_gate==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Unfamiliar_Already_Used_Images = {};
    Unfamiliar_img_count = 1;
    Unfamiliar_display_objs = struct('Target',{},'Target_Cat',{},'Foil',{},'Foil_Cat',{},'Study_stimuli',{},'Study_stimuli_Cat',{},'Target_study_pos',[],'TF_Sim',[]);
    
    %create shuffled target Location Index such that in each block the target
    %will be in each potential location as close as possible to an equal number
    %of times. Because there are 40 trials per condition, the target will be shown
    %in 4 of the study locations 7 times and in 2 of the study locations the
    %target will only be shown 6 times.
    temp = repmat([1:Num_stim_per_trial],1,ceil(Num_obj_trials/Num_stim_per_trial/2));
    %randomly select which of the 2 study locations will only be used 6 times
    locs2remove = randi([1,6],1,2);
    for i = 1:2
        idx2remove = find(temp==locs2remove(i),1,'last');
        temp(idx2remove) = [];
    end
    temp = temp(1:(Num_obj_trials/2));
    %shuffle index separately for 1st block and 2nd block then add together
    temp1 = temp(randperm(Num_obj_trials/2));
    temp2 = temp(randperm(Num_obj_trials/2));
    shuffled_target_loc_idx_Unfamiliar = [temp1 temp2];
    
    %% Deal out Unfam Target Foil Pairs, randomizing the trial order and also which obj in each pair is target and which is foil
    unfam_shuffled_order = randperm(Num_obj_trials);
    
    for trl = 1:Num_obj_trials
        tmp = randperm(2);
        if tmp(1)==1
            Unfam_Target{trl} = TF_pairs_unfam.Obj1{unfam_shuffled_order(trl)};
            Unfam_Target_Cat{trl} = TF_pairs_unfam.Obj1Cat{unfam_shuffled_order(trl)};
            Unfam_Foil{trl} = TF_pairs_unfam.Obj2{unfam_shuffled_order(trl)};
            Unfam_Foil_Cat{trl} = TF_pairs_unfam.Obj2Cat{unfam_shuffled_order(trl)};
        else
            Unfam_Target{trl} = TF_pairs_unfam.Obj2{unfam_shuffled_order(trl)};
            Unfam_Target_Cat{trl} = TF_pairs_unfam.Obj2Cat{unfam_shuffled_order(trl)};
            Unfam_Foil{trl} = TF_pairs_unfam.Obj1{unfam_shuffled_order(trl)};
            Unfam_Foil_Cat{trl} = TF_pairs_unfam.Obj1Cat{unfam_shuffled_order(trl)};
        end
        Unfam_TF_Sim(trl) = TF_pairs_unfam.SimRanking(unfam_shuffled_order(trl));
    end
    
    %find number of Targets per category
    for cat = 1:length(categories)
        unfam_target_cat_counts(cat) = sum(strcmp(categories{cat},Unfam_Target_Cat));
    end
    
    %Order Target from Most to least common Target_Cat placing the "Other"
    %category last regardless of count
    [Unfam_sorted,unfam_idx] = sort(unfam_target_cat_counts);
    unfam_idx = fliplr(unfam_idx);
    %remove "other" category
    idx2remove = find(unfam_idx == 12);
    unfam_idx(idx2remove) = [];
    %add "other" at end
    unfam_cat_order = [unfam_idx, 12];
    
    
    %Now deal out the Target_Foil pairs using the order just found
    unfam_trl_count = 0;
    for cat = 1:length(categories)
        idx = find(strcmp(categories{unfam_cat_order(cat)},Unfam_Target_Cat)==1);
        for i = 1:numel(idx)
            unfam_trl_count = unfam_trl_count + 1;
            temp_Unfamiliar_display_objs(unfam_trl_count).Target = Unfam_Target(idx(i));
            temp_Unfamiliar_display_objs(unfam_trl_count).Target_Cat = Unfam_Target_Cat(idx(i));
            temp_Unfamiliar_display_objs(unfam_trl_count).Foil = Unfam_Foil(idx(i));
            temp_Unfamiliar_display_objs(unfam_trl_count).Foil_Cat = Unfam_Foil_Cat(idx(i));
            temp_Unfamiliar_display_objs(unfam_trl_count).TF_Sim = Unfam_TF_Sim(idx(i));
            %Add selected Foil image to Used Image array for tracking
            Unfamiliar_Already_Used_Images(Unfamiliar_img_count) = Unfam_Target(idx(i));
            Unfamiliar_img_count = Unfamiliar_img_count + 1;
            Unfamiliar_Already_Used_Images(Unfamiliar_img_count) = Unfam_Foil(idx(i));
            Unfamiliar_img_count = Unfamiliar_img_count + 1;
        end
    end
    
    %Randomly select five study display images for each trial that do not
    %overlap with Target Category (unless they are in the "Other
    %category). First list the categories in order of how many objects are in each.
    %Then, working from the largest to smallest category, distribute the
    %objects into 80 groups of 5. Deal out all objects in each category such
    %that "Other" category repeats only occur once all of the 80 groups have at least
    %1 object from that category. Check each object to make sure it hasnt been
    %used as a target or foil before including.
    unfam_cat_sizes = zeros(1,length(fieldnames(Unfamiliar_imgFiles)));
    for c = 1:length(unfam_cat_sizes)
        unfam_cat_sizes(c) = numel(Unfamiliar_imgFiles.(categories{c}));
    end
    [~,unfam_cat_ranks] = sort(unfam_cat_sizes);
    unfam_cat_ranks = fliplr(unfam_cat_ranks);
    
    unfam_obj_count = 0;
    for i = 1:length(unfam_cat_ranks)
        shuffled_idx = randperm(length(Unfamiliar_imgFiles.(categories{unfam_cat_ranks(i)})));
        for j = 1:length(Unfamiliar_imgFiles.(categories{unfam_cat_ranks(i)}))
            if sum(strcmp(Unfamiliar_imgFiles.(categories{unfam_cat_ranks(i)})(shuffled_idx(j)).name,Unfamiliar_Already_Used_Images)) == 0
                unfam_obj_count = unfam_obj_count + 1;
                Unfam_Objects_name{unfam_obj_count} = Unfamiliar_imgFiles.(categories{unfam_cat_ranks(i)})(shuffled_idx(j)).name;
                Unfam_Objects_cat{unfam_obj_count} = categories{unfam_cat_ranks(i)};
            end
        end
    end
    tmp = reshape(Unfam_Objects_name,[80,5]);
    tmp2 = reshape(Unfam_Objects_cat,[80,5]);
    
    try
        %select set of 5 filler study display objects that dont overlap with the Target
        %category
        shuffled_idx = randperm(Num_obj_trials);
        tmp = tmp(shuffled_idx,:);
        tmp2 = tmp2(shuffled_idx,:);
        for i = 1:Num_obj_trials
            stim(1) = temp_Unfamiliar_display_objs(i).Target;
            stim_Cat(1) = temp_Unfamiliar_display_objs(i).Target_Cat;
            
            for j = 1:length(tmp)
                if isempty(intersect(temp_Unfamiliar_display_objs(i).Target_Cat,tmp2(j,:))) | strcmp(temp_Unfamiliar_display_objs(i).Target_Cat,'other')==1
                    stim(2:6) = tmp(j,:);
                    stim_Cat(2:6) = tmp2(j,:);
                    
                    temp_Unfamiliar_display_objs(i).filler_stim = stim;
                    temp_Unfamiliar_display_objs(i).filler_stim_Cat = stim_Cat;
                    
                    %Add selected Foil image to Used Image array for tracking
                    Unfamiliar_Already_Used_Images(Unfamiliar_img_count:Unfamiliar_img_count+4) = stim(2:6);
                    Unfamiliar_img_count = Unfamiliar_img_count + 5;
                    
                    tmp(j,:) = [];
                    tmp2(j,:) = [];
                    break;
                end
                
                if j+1==length(tmp)
                    error(['Unable to find suitable Study Display filler images on trial num = ' num2str(i)])
                end
            end
        end
        repeat_gate = 0;
    catch
        repeat_gate=1;
    end
end

%Now Randomize trial order
shuffled_trial_idx = randperm(Num_obj_trials);
for trl = 1:Num_obj_trials
    Unfamiliar_display_objs(trl).Target = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).Target;
    Unfamiliar_display_objs(trl).Target_Cat = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).Target_Cat;
    Unfamiliar_display_objs(trl).Foil = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).Foil;
    Unfamiliar_display_objs(trl).Foil_Cat = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).Foil_Cat;
    Unfamiliar_display_objs(trl).TF_Sim = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).TF_Sim;
    
    %randomize Target Study Position such that each target
    %location occurs an equal number of times per condition
    stim_locs = zeros(1,Num_stim_per_trial);
    stim_locs(shuffled_target_loc_idx_Unfamiliar(trl)) = 1;
    shuffled_other_stim_locs = randperm(Num_stim_per_trial);
    shuffled_other_stim_locs = shuffled_other_stim_locs(shuffled_other_stim_locs>1);
    count = 0;
    for j = 1:length(stim)
        if stim_locs(j) == 0
            count = count + 1;
            stim_locs(j) = shuffled_other_stim_locs(count);
        end
    end
    stim = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).filler_stim;
    stim_cat = temp_Unfamiliar_display_objs(shuffled_trial_idx(trl)).filler_stim_Cat;
    %Save Encoding phase stimuli and Target Position
    Unfamiliar_display_objs(trl).Study_stimuli = stim(stim_locs);
    Unfamiliar_display_objs(trl).Study_stimuli_Cat = stim_cat(stim_locs);
    Unfamiliar_display_objs(trl).Target_study_pos = find(stim_locs == 1);
end

