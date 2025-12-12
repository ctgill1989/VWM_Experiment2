function [Unfamiliar_display_faces1,Unfamiliar_display_faces2]=Generate_rand_Unfam_face_displays_vWMExp2(TF_pairs_unfam_faces_RS_binned,nRepsPerCondition,Unfamiliar_faces_StimFolder )
%By CTGill
%Last updated: 11/3/22
%
% This Version of Generate_rand_obj_displays is used in vWM_Experiment2 and
% ensures that each target location is selected an equal number of times for
% each condition.
%
% This versiion of Generate_rand_obj_displays is the functional equivalent
% of Generate_rand_obj_displays_V2 used in vWM_Experiment_V1 (i.e.
% Experiment1)
%
% Function generates sets of face stimuli for each trial for
% the Familiar Face Conidtion. No Faces are displayed more than once.
% To find the maximally dissimilar Target-Foil pairs I extradted the
% features of all objects in my Famous and Unfamiliar face sets
% from the top max-pooling layer of VGGFace2-pretrained ResNet-50 model and
% then computed similarity metric between all pairs of objects within each
% set using Cosine Similarity (i.e. length-normalized dot product).
%
%
% INPUTS
% TF_pairs_fam_faces --> a struct containing both objects in a familiar
%                        pair and the cosine similarity measure.
%
% Num_stim_per_trial --> The number of stimuli in each display of the
%                        encoding phase of the vWM task. (e.g. 6)
%
%
% Num_obj_trials_per_block     --> The total number of real-world-object trials per familiarity condition.
%
% OUTPUTS
%
% Unfamiliar_display_faces       --> a struct containing face display info.
%                        display_faces.Target            --> contains the Target object filename.
%                        display_faces.Foil              --> contains the Foil object filename.
%                        display_faces.Study_stimuli     --> contains filename of all stim in study display.
%                        display_faces.Target_study_pos  --> contains the Target object position in study display (Pos 1 is always at the 9'oclock position and subseqent positions go clockwise).
%                        display_faces.Cat               --> contains the Race-Sex category of the pair
%                        display_faces.TF_Sim            --> contains Cosine Similarity measure of the Target-Foil pair   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_stim_per_trial = 6;
Num_obj_trials_per_block = nRepsPerCondition/2; %40 trials in each of 2 blocks, each block should have 20 female face trials and 20 male face trials

%% Load Stim Files
%Set path to stim folders
% Unfamiliar_faces_StimFolder = '/Users/christophergill/Documents/MATLAB/vWM_Experiment2/Stimuli/Faces/Unfamiliar';

%6 Face Image Categories (Race and Sex)
categories = {'AF','AM','BF','BM','WF','WM'}';

Unfamiliar_imgFiles = [];
filePattern = fullfile(Unfamiliar_faces_StimFolder , '*.jpg');
Unfamiliar_imgFiles = dir(filePattern);
Unfamiliar_imgFiles = struct2cell(Unfamiliar_imgFiles)';
Unfamiliar_imgFiles(:,2:end) = [];
%find Cat associated with each Filename and merge these together
for i = 1:numel(Unfamiliar_imgFiles)
    cat_idx = strfind(Unfamiliar_imgFiles{i},'-');
    if strcmp(Unfamiliar_imgFiles{i}(cat_idx(1)+1),'L')
        imgFiles_cats{i} = ['W' Unfamiliar_imgFiles{i}(cat_idx(2)-1)];
    else
        imgFiles_cats{i} = Unfamiliar_imgFiles{i}(cat_idx(1)+1:cat_idx(2)-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Unfamiliar_Already_Used_Images_female = {};
Unfamiliar_Already_Used_Images_male = {};
Unfamiliar_img_count_female = 1;
Unfamiliar_img_count_male = 1;
Unfamiliar_display_faces1 = struct('Target',{},'Foil',{},'Cat',{},'Study_stimuli',{},'Target_study_pos',[],'TF_Sim',[]);
Unfamiliar_display_faces2 = struct('Target',{},'Foil',{},'Cat',{},'Study_stimuli',{},'Target_study_pos',[],'TF_Sim',[]);

%create shuffled target Location Index such that in each block the target
%will be in each potential location as close as possible to an equal number
%of times. Because there are 40 trials per block, the target will be shown
%in 4 of the study locations 7 times and in 2 of the study locations the
%target will only be shown 6 times.
temp1 = repmat([1:Num_stim_per_trial],1,ceil(Num_obj_trials_per_block/Num_stim_per_trial));
temp1(randi([1,6],1,2)) = [];
temp2 = repmat([1:Num_stim_per_trial],1,ceil(Num_obj_trials_per_block/Num_stim_per_trial));
temp2(randi([1,6],1,2)) = [];
%shuffle indices for Block 1 and Block 2 separately to approximately
%balance target location across the 2 blocks
shuffled_target_loc_idx_Unfamiliar1 = temp1(randperm(Num_obj_trials_per_block));
shuffled_target_loc_idx_Unfamiliar2 = temp2(randperm(Num_obj_trials_per_block));

%% Deal out Unfam Target Foil Pairs, randomizing the trial order and also which Face in each pair is target and which is foil
unfam_female_trl_count = 0;
unfam_male_trl_count = 0;
for cat = 1:numel(categories)
    if sum(strcmp(categories{cat},{'AF','BF','WF'}))
        for i = 1:numel(TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1)
            unfam_female_trl_count = unfam_female_trl_count + 1;
            tmp = randperm(2);
            if tmp(1)==1
                temp_Unfamiliar_display_faces_female(unfam_female_trl_count).Target = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1{i};
                temp_Unfamiliar_display_faces_female(unfam_female_trl_count).Foil = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face2{i};
            else
                temp_Unfamiliar_display_faces_female(unfam_female_trl_count).Target = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face2{i};
                temp_Unfamiliar_display_faces_female(unfam_female_trl_count).Foil = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1{i};
            end
            temp_Unfamiliar_display_faces_female(unfam_female_trl_count).TF_Sim = TF_pairs_unfam_faces_RS_binned.(categories{cat}).SimRanking(i);
            temp_Unfamiliar_display_faces_female(unfam_female_trl_count).Cat = categories{cat};
            %Add selected Foil image to Used Image array for tracking
            Unfamiliar_Already_Used_Images_female{Unfamiliar_img_count_female} = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1{i};
            Unfamiliar_img_count_female = Unfamiliar_img_count_female + 1;
            Unfamiliar_Already_Used_Images_female{Unfamiliar_img_count_female} = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face2{i};
            Unfamiliar_img_count_female = Unfamiliar_img_count_female + 1;
        end
    else
        for i = 1:numel(TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1)
            unfam_male_trl_count = unfam_male_trl_count + 1;
            tmp = randperm(2);
            if tmp(1)==1
                temp_Unfamiliar_display_faces_male(unfam_male_trl_count).Target = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1{i};
                temp_Unfamiliar_display_faces_male(unfam_male_trl_count).Foil = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face2{i};
            else
                temp_Unfamiliar_display_faces_male(unfam_male_trl_count).Target = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face2{i};
                temp_Unfamiliar_display_faces_male(unfam_male_trl_count).Foil = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1{i};
            end
            temp_Unfamiliar_display_faces_male(unfam_male_trl_count).TF_Sim = TF_pairs_unfam_faces_RS_binned.(categories{cat}).SimRanking(i);
            temp_Unfamiliar_display_faces_male(unfam_male_trl_count).Cat = categories{cat};
            %Add selected Foil image to Used Image array for tracking
            Unfamiliar_Already_Used_Images_male{Unfamiliar_img_count_male} = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face1{i};
            Unfamiliar_img_count_male = Unfamiliar_img_count_male + 1;
            Unfamiliar_Already_Used_Images_male{Unfamiliar_img_count_male} = TF_pairs_unfam_faces_RS_binned.(categories{cat}).face2{i};
            Unfamiliar_img_count_male = Unfamiliar_img_count_male + 1;
        end
    end
end


%Select five study display filler faces from the same Race-Sex category for
%each trial. Check each face to make sure it hasnt been already been
%used as a target or foil before including.
unfam_female_count = 0;
unfam_male_count = 0;
for i = 1:length(temp_Unfamiliar_display_faces_female)
    for cat = 1:length(categories)
        idx = find(strcmp(categories{cat},imgFiles_cats)==1);
        shuffled_idx = randperm(length(idx));
        
        Filenames.(categories{cat}) = Unfamiliar_imgFiles(idx(shuffled_idx));
        %find already used Filenames
        if sum(strcmp(categories{cat},{'AF','BF','WF'}))
            [~,idx2remove,~] = intersect(Filenames.(categories{cat}),Unfamiliar_Already_Used_Images_female);
        else
            [~,idx2remove,~] = intersect(Filenames.(categories{cat}),Unfamiliar_Already_Used_Images_male);
        end
        Filenames.(categories{cat})(idx2remove) = [];
        shuff_idx = randperm(length(Filenames.(categories{cat})));
        names = Filenames.(categories{cat})(shuff_idx);
        tmp = reshape(names,[(size(Filenames.(categories{cat}),1)/5),5]);
        
        
        for trl = 1:size(tmp,1)
            if sum(strcmp(categories{cat},{'AF','BF','WF'}))
                unfam_female_count = unfam_female_count + 1;
                stim{1} = temp_Unfamiliar_display_faces_female(unfam_female_count).Target;
                stim(2:6) = tmp(trl,:);
                temp_Unfamiliar_display_faces_female(unfam_female_count).filler_stim = stim;
                 %Add selected Foil image to Used Image array for tracking
                Unfamiliar_Already_Used_Images_female(Unfamiliar_img_count_female:Unfamiliar_img_count_female+4) = stim(2:6);
                Unfamiliar_img_count_female = Unfamiliar_img_count_female + 5;
            else
                unfam_male_count = unfam_male_count + 1;
                stim{1} = temp_Unfamiliar_display_faces_male(unfam_male_count).Target;
                stim(2:6) = tmp(trl,:);
                temp_Unfamiliar_display_faces_male(unfam_male_count).filler_stim = stim;
                 %Add selected Foil image to Used Image array for tracking
                Unfamiliar_Already_Used_Images_male(Unfamiliar_img_count_male:Unfamiliar_img_count_male+4) = stim(2:6);
                Unfamiliar_img_count_male = Unfamiliar_img_count_male + 5;
            end
        end
    end
end


%Now take half the males from each race and half the females from each race
%creating two groups, one for Group1 and one for Group2.
%Group1 and one for Group2 will always have the same number of male and female faces
%and will always have the same number of faces in each race group. However,
%the number of trials using BM, BF, WM, and WF will differ slightly.
%Group1 will always have 3 AM, 2 AF, 5 BM, 4 BF, 12 WM, and 14 WF trials,
%whereas Group2 will have 3 AM, 2 AF, 4 BM, 5 BF, 13 WM, and 13 WF trials.

Group1_count = 0;
Group2_count = 0;
for i = 1:length(temp_Unfamiliar_display_faces_male)
    if mod(i,2)
        Group1_count = Group1_count + 1;
        Group1(Group1_count).Target = temp_Unfamiliar_display_faces_male(i).Target;
        Group1(Group1_count).Foil = temp_Unfamiliar_display_faces_male(i).Foil;
        Group1(Group1_count).TF_Sim = temp_Unfamiliar_display_faces_male(i).TF_Sim;
        Group1(Group1_count).Cat = temp_Unfamiliar_display_faces_male(i).Cat;
        Group1(Group1_count).filler_stim = temp_Unfamiliar_display_faces_male(i).filler_stim;
    else
        Group2_count = Group2_count + 1;
        Group2(Group2_count).Target = temp_Unfamiliar_display_faces_male(i).Target;
        Group2(Group2_count).Foil = temp_Unfamiliar_display_faces_male(i).Foil;
        Group2(Group2_count).TF_Sim = temp_Unfamiliar_display_faces_male(i).TF_Sim;
        Group2(Group2_count).Cat = temp_Unfamiliar_display_faces_male(i).Cat;
        Group2(Group2_count).filler_stim = temp_Unfamiliar_display_faces_male(i).filler_stim;
    end
end

for i = 1:length(temp_Unfamiliar_display_faces_female)
    if mod(i,2)
        Group2_count = Group2_count + 1;
        Group2(Group2_count).Target = temp_Unfamiliar_display_faces_female(i).Target;
        Group2(Group2_count).Foil = temp_Unfamiliar_display_faces_female(i).Foil;
        Group2(Group2_count).TF_Sim = temp_Unfamiliar_display_faces_female(i).TF_Sim;
        Group2(Group2_count).Cat = temp_Unfamiliar_display_faces_female(i).Cat;
        Group2(Group2_count).filler_stim = temp_Unfamiliar_display_faces_female(i).filler_stim;
    else
        Group1_count = Group1_count + 1;
        Group1(Group1_count).Target = temp_Unfamiliar_display_faces_female(i).Target;
        Group1(Group1_count).Foil = temp_Unfamiliar_display_faces_female(i).Foil;
        Group1(Group1_count).TF_Sim = temp_Unfamiliar_display_faces_female(i).TF_Sim;
        Group1(Group1_count).Cat = temp_Unfamiliar_display_faces_female(i).Cat;
        Group1(Group1_count).filler_stim = temp_Unfamiliar_display_faces_female(i).filler_stim;
    end
end


%Now randomly select which, Group1 or Group2 will be used in the Long
%condition and which will be used in the Short condition and Randomize
%trial order

if randi([0,1])
    
    shuffled_trial_idx = randperm(Num_obj_trials_per_block);
    for trl = 1:Num_obj_trials_per_block
        Unfamiliar_display_faces1(trl).Target = Group1(shuffled_trial_idx(trl)).Target;
        Unfamiliar_display_faces1(trl).Foil = Group1(shuffled_trial_idx(trl)).Foil;
        Unfamiliar_display_faces1(trl).Cat = Group1(shuffled_trial_idx(trl)).Cat;
        Unfamiliar_display_faces1(trl).TF_Sim = Group1(shuffled_trial_idx(trl)).TF_Sim;
        
        %randomize Target Study Position such that each target
        %location occurs an equal number of times per condition
        stim_locs = zeros(1,Num_stim_per_trial);
        stim_locs(shuffled_target_loc_idx_Unfamiliar1(trl)) = 1;
        shuffled_other_stim_locs = randperm(Num_stim_per_trial);
        shuffled_other_stim_locs = shuffled_other_stim_locs(shuffled_other_stim_locs>1);
        count = 0;
        for j = 1:length(stim)
            if stim_locs(j) == 0
                count = count + 1;
                stim_locs(j) = shuffled_other_stim_locs(count);
            end
        end
        stim = Group1(shuffled_trial_idx(trl)).filler_stim;
        %Save Encoding phase stimuli and Target Position
        Unfamiliar_display_faces1(trl).Study_stimuli = stim(stim_locs);
        Unfamiliar_display_faces1(trl).Target_study_pos = find(stim_locs == 1);
    end
    
    shuffled_trial_idx = randperm(Num_obj_trials_per_block);
    for trl = 1:Num_obj_trials_per_block
        Unfamiliar_display_faces2(trl).Target = Group2(shuffled_trial_idx(trl)).Target;
        Unfamiliar_display_faces2(trl).Foil = Group2(shuffled_trial_idx(trl)).Foil;
        Unfamiliar_display_faces2(trl).Cat = Group2(shuffled_trial_idx(trl)).Cat;
        Unfamiliar_display_faces2(trl).TF_Sim = Group2(shuffled_trial_idx(trl)).TF_Sim;
        
        %randomize Target Study Position such that each target
        %location occurs an equal number of times per condition
        stim_locs = zeros(1,Num_stim_per_trial);
        stim_locs(shuffled_target_loc_idx_Unfamiliar2(trl)) = 1;
        shuffled_other_stim_locs = randperm(Num_stim_per_trial);
        shuffled_other_stim_locs = shuffled_other_stim_locs(shuffled_other_stim_locs>1);
        count = 0;
        for j = 1:length(stim)
            if stim_locs(j) == 0
                count = count + 1;
                stim_locs(j) = shuffled_other_stim_locs(count);
            end
        end
        stim = Group2(shuffled_trial_idx(trl)).filler_stim;
        %Save Encoding phase stimuli and Target Position
        Unfamiliar_display_faces2(trl).Study_stimuli = stim(stim_locs);
        Unfamiliar_display_faces2(trl).Target_study_pos = find(stim_locs == 1);
    end
else
    shuffled_trial_idx = randperm(Num_obj_trials_per_block);
    for trl = 1:Num_obj_trials_per_block
        Unfamiliar_display_faces1(trl).Target = Group2(shuffled_trial_idx(trl)).Target;
        Unfamiliar_display_faces1(trl).Foil = Group2(shuffled_trial_idx(trl)).Foil;
        Unfamiliar_display_faces1(trl).Cat = Group2(shuffled_trial_idx(trl)).Cat;
        Unfamiliar_display_faces1(trl).TF_Sim = Group2(shuffled_trial_idx(trl)).TF_Sim;
        
        %randomize Target Study Position such that each target
        %location occurs an equal number of times per condition
        stim_locs = zeros(1,Num_stim_per_trial);
        stim_locs(shuffled_target_loc_idx_Unfamiliar1(trl)) = 1;
        shuffled_other_stim_locs = randperm(Num_stim_per_trial);
        shuffled_other_stim_locs = shuffled_other_stim_locs(shuffled_other_stim_locs>1);
        count = 0;
        for j = 1:length(stim)
            if stim_locs(j) == 0
                count = count + 1;
                stim_locs(j) = shuffled_other_stim_locs(count);
            end
        end
        stim = Group2(shuffled_trial_idx(trl)).filler_stim;
        %Save Encoding phase stimuli and Target Position
        Unfamiliar_display_faces1(trl).Study_stimuli = stim(stim_locs);
        Unfamiliar_display_faces1(trl).Target_study_pos = find(stim_locs == 1);
    end
    
    shuffled_trial_idx = randperm(Num_obj_trials_per_block);
    for trl = 1:Num_obj_trials_per_block
        Unfamiliar_display_faces2(trl).Target = Group1(shuffled_trial_idx(trl)).Target;
        Unfamiliar_display_faces2(trl).Foil = Group1(shuffled_trial_idx(trl)).Foil;
        Unfamiliar_display_faces2(trl).Cat = Group1(shuffled_trial_idx(trl)).Cat;
        Unfamiliar_display_faces2(trl).TF_Sim = Group1(shuffled_trial_idx(trl)).TF_Sim;
        
        %randomize Target Study Position such that each target
        %location occurs an equal number of times per condition
        stim_locs = zeros(1,Num_stim_per_trial);
        stim_locs(shuffled_target_loc_idx_Unfamiliar2(trl)) = 1;
        shuffled_other_stim_locs = randperm(Num_stim_per_trial);
        shuffled_other_stim_locs = shuffled_other_stim_locs(shuffled_other_stim_locs>1);
        count = 0;
        for j = 1:length(stim)
            if stim_locs(j) == 0
                count = count + 1;
                stim_locs(j) = shuffled_other_stim_locs(count);
            end
        end
        stim = Group1(shuffled_trial_idx(trl)).filler_stim;
        %Save Encoding phase stimuli and Target Position
        Unfamiliar_display_faces2(trl).Study_stimuli = stim(stim_locs);
        Unfamiliar_display_faces2(trl).Target_study_pos = find(stim_locs == 1);
    end
end

