%% vWM_Exp2
% by CTGill
% Last updated: 11/3/22
%
% Gaze-contingent visual working memory task with EyeLink integration
%
% - 6-item VWM display (objects, faces, or colors) arranged on a circle
% - Gaze-contingent study phase: stimuli appear only after central fixation
%   is held within a window for a minimum duration
% - Simultaneous EEG + EyeLink recording with event triggers
%License: GNU-3.0-or-later (see LICENSE)
%
% NOTE: To Quit/Exit/End experiment press 'q' key
%
% **Note: Before running Experiment make sure to set the StimFolder path
% to the correct folder, and adjust monitor resolution and dimensions!
%
% User data is saved as a .mat file to Subject Data folder as "SubjectNum" (e.g. Sbj001.mat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
rng('default');
rng('shuffle');

params.port = 1;                                    % 0 = off, 1 = on
params.eyeTrack = 1;                                % 0 = off, 1 = on
params.eyeMode = 0;                                 % 0 = chin rest, 1 = remote
params.debug = 0;                                   % 1 = debug mode
params.skip_practice = 0;                           % 1 = Skipping practice trials
params.windowed = 0;                                % 1 = reduced window size for debugging

%save current directory
root = pwd;
%add folder with called functions to path
SupportFunctionsDir = [root, filesep 'SupportFunctions' filesep];
addpath(SupportFunctionsDir);
%add stim folders to path
Famous_faces_StimFolder = [root, filesep 'Stimuli' filesep 'Faces' filesep 'Famous' filesep];
Unfamiliar_faces_StimFolder = [root, filesep 'Stimuli' filesep 'Faces' filesep 'Unfamiliar' filesep];
Familiar_objects_StimFolder = [root, filesep 'Stimuli' filesep 'Objects' filesep 'Familiar' filesep];
Unfamiliar_objects_StimFolder = [root, filesep 'Stimuli' filesep 'Objects' filesep 'Unfamiliar' filesep];
addpath(Famous_faces_StimFolder,Unfamiliar_faces_StimFolder,Familiar_objects_StimFolder,Unfamiliar_objects_StimFolder);
%add folders to practice stimuli
Practice_Long_FamFace_StimFolder = [root, filesep 'Stimuli' filesep 'Practice' filesep 'Long_Fam_Face_Practice' filesep];
Practice_Long_UnfamFace_StimFolder = [root, filesep 'Stimuli' filesep 'Practice' filesep 'Long_Unfam_Face_Practice' filesep];
Practice_Long_FamObj_StimFolder = [root, filesep 'Stimuli' filesep 'Practice' filesep 'Long_Fam_Obj_Practice' filesep];
Practice_Long_UnfamObj_StimFolder = [root, filesep 'Stimuli' filesep 'Practice' filesep 'Long_Unfam_Obj_Practice' filesep];
addpath(Practice_Long_FamFace_StimFolder,Practice_Long_UnfamFace_StimFolder,Practice_Long_FamObj_StimFolder,Practice_Long_UnfamObj_StimFolder);

if ~params.debug
    %User input prompt
    prompt = {'Subject Number:'}; %description of fields
    title = 'Enter Subject Info';
    dims = [1 50];
    defaults = {'000'};
    answer = inputdlg(prompt,title,dims,defaults);
    SubjectNum = answer{1,:}; %Gets Subject ID Number
    output.ExperimentalSession_start_time = clock;
    
    % Build an output directory if it doesn't already exist
    if ~exist([root filesep 'Subject Data' filesep],'dir')
        mkdir([root filesep 'Subject Data' filesep]);
    end
    filename = [root filesep 'Subject Data' filesep 'vWM2_' num2str(SubjectNum) '.mat'];
    
    while exist(filename) > 1
        %User input prompt
        prompt = {'Subject Number:'}; %description of fields
        title = 'Entered Subject Number Already In Use Select New Subject Number';
        dims = [1 100];
        defaults = {'000'};
        answer = inputdlg(prompt,title,dims,defaults);
        SubjectNum = answer{1,:}; %Gets Subject ID Number
        
        filename = [root filesep 'Subject Data' filesep 'vWM_' num2str(SubjectNum) '.mat'];% num2str(SubjectNum);%strcat(SubjectNum,'_',num2str(c(2)),num2str(c(3)),num2str(c(1)),num2str(c(4)),num2str(c(5)));
    end
    
    %%Counterbalance Conditions
    Sbj_count = numel(dir([root filesep 'Subject Data' filesep '*.mat']));
    [Cond_Order] = rand_Cond_Order_vWMExp2();
    output.SubjectNum = SubjectNum;
else
    prompt = {'Debug Filename'}; %description of fields
    title = 'Running In Debug Mode';
    dims = [1 50];
    defaults = {'Testing'};
    answer = inputdlg(prompt,title,dims,defaults);
    SubjectNum = answer{1,:}; %Gets Subject ID Number
    output.ExperimentalSession_start_time = clock;
    
    % Build an output directory if it doesn't already exist
    if ~exist([root filesep 'Debug Mode Data' filesep],'dir')
        mkdir([root filesep 'Debug Mode Data' filesep]);
    end
    filename = [root filesep 'Debug Mode Data' filesep SubjectNum '.mat'];
    
    [Cond_Order] = rand_Cond_Order_vWMExp2();
    output.SubjectNum = SubjectNum;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Params
params.Presentation_Dur_conditions = 1;             % 1=Long(2000ms)
params.Stim_type_conditions = [1 2 3 4 5];          % 1=Fam Obj; 2=Unfam Obj; 3=Fam Face; 4=Unfam Face; 5=Color
params.Num_Blocks = 10;                             % Number of Blocks
params.Num_conditions = 5;                          % total number of conditions
params.nRepsPerCondition = 80;                      % Number of trials per condition
params.Num_Trials = 40;                             % Number of trials per block
params.Long_StimDuration = 2;                       % Stim display duration in seconds for Long Condition (2000ms)
% params.Long_StimDuration = .25;                    % Stim display duration in seconds for Short Condition (250ms)
params.DelayPeriod = .8;                            % WM delay period (800ms)
params.preTrialDelay = .9;                          % fixed Delay Interval between offset of Articulatory suppression task stim and onset of vWM task stimuli (100ms additional fixation is then required to initiate vWM study stim onset)
params.Blink_time = 1;                              % This is the time period during the ITI with the red fixation cross drawn to the screen (Red fixation cross is a sign for the subject that it is ok to blink)
params.ITI = 2;                                     % Intertial Interval in seconds (time between trials) (seconds)
params.rehearsal_digits_duration = 1;               % rehearsal digits display duration, 1000ms
params.RectFrame_color = [.75 .75 .75 .5];          % Color of rectangle outlines (light gray)
params.Target_RectFrame_color = [0 0 0];            % color of target position rectangle (black)
params.Background_color = [1 1 1];                  % background color (white)
params.calib_target_colors = [0 0 0];               % calibration target colors  (black)
params.trialSetSize = 6;                            % Number of stimuli displayed per trial
params.stimSize = 2.5;%2.1;                              % size of individual stimuli in visual degrees
params.rectSize = 2.6;%2.2;                              % size of each rectangle in visual degrees
params.familiarity_ranking_stimSize = 4.2;          % size of the familiarity ranking stimulus in visual degrees
params.fixationSize = 0.15;                         % size of fixation cross in visual degrees
params.fixWinSize = 1.5;                            % size of fixation window in visual degrees (used in fixation contingent vWM study stim onset)
params.fixateTime = .1;                             % fixation for 100ms is needed to initiate onset of VWM study phase
params.arrow_size = 1;                              % length of arrow image in rating phase of experiment
params.fixationColor = [0 0 0];                     % color of fixation cross when it is not ok for sbj to blink (black)
params.Blink_fixationColor = [1 0 0];               % color of fixation cross when it is ok for sbj to blink (red)
params.circRad_studyStim = 4.2;                     % distance from fixation to center of each rect in visual degrees for stimuli in study phase
params.circRad_testStim = (params.rectSize/2)+0.2;  % distance from fixation to center of each rect in visual degrees for Test and and Foil stimuli
params.monitor_distance = 940;                      % Distance from screen in mm
params.rehearsal_digits_font_size = 40;             % Rehearsal Digits Font size
params.text_font_size = 26;                         % Text Instructions Font size
params.NumPracticeTrials = 5;                      % Number of practice trials - 1 practice trial for each of the 10 conditions
params.color_radius = 49;                           % In the CIE L*a*b space based on Brady & Stormer (2020)
params.color_L = 54;
params.color_a = 21.5;
params.color_b = 11.5;
params.color_min_distance = 45;                     % All displayed colors in each trial must be a minimum of 15deg apart on CIE color circle (6 study colors and 1 Foil).
params.ET_sampleing_rate = 1000;                    % eye-tracker sampling rate (samples per sec (i.e. Hz))
%% Event Descriptors
params.event.Exp_st = 1;                              %Experiment start
params.event.Exp_end = 2;                             %Experiment end
params.event.Long_Fam_Obj_Blk_st = 3;                 %Long Familiar Object Cond. block start
params.event.Long_Unfam_Obj_Blk_st = 4;               %Long Unfamiliar Object Cond. block start
params.event.Long_Fam_Face_Blk_st = 5;                %Long Famous Face Cond. block start
params.event.Long_Unfam_Face_Blk_st = 6;              %Long Unfamiliar Face Cond. block start
params.event.Long_Color_Blk_st = 7;                   %Long Color Cond. block start
% params.event.Short_Fam_Obj_Blk_st = 8;                %Short Familiar Object Cond. block start
% params.event.Short_Unfam_Obj_Blk_st = 9;              %Short Unfamiliar Object Cond. block start
% params.event.Short_Fam_Face_Blk_st = 10;              %Short Famous Face Cond. block start
% params.event.Short_Unfam_Face_Blk_st = 11;            %Short Unfamiliar Face Cond. block start
% params.event.Short_Color_Blk_st = 12;                 %Short Color Cond. block start
params.event.Blk_End = 13;                            %Block End
params.event.Trl_st = 15;                             %Start of each trial
params.event.Trl_end = 16;                            %End of each trial
params.event.Stim_Onset = 20;                         %Onset of Study stim display
params.event.Delay_Period_Onset = 40;                 %Onset of wm delay period between final stimulus offset and test probe
params.event.Test_Stim_Onset = 50;                    %Onset of Target Probe and Foil
params.event.Resp = 55;                               %Response given for vWM task
params.event.ArticSupressTask_Study_Stim_Onset = 70;  %Onset of stimulus (digits) for study phase of silent articulatory supression task (AKA rehearsal task)
params.event.ArticSupressTask_Study_Stim_Offset = 71; %Offset of stimulus (digits) for study phase of silent articulatory supression task (AKA rehearsal task)
params.event.ArticSupressTask_Test_Stim_Onset = 72;   %Onset of Target and Foil (digits) for test phase of silent articulatory supression task (AKA rehearsal task)
params.event.ArticSupressTask_Resp = 75;              %Response given for articulatory supression task
params.event.EyeTracking_calib_st = 80;               %Onset of Eyetacking calibration before each block
params.event.EyeTracking_calib_end = 81;              %Offset of Eyetacking calibration before each block
params.event.FamRank_Phase_st = 90;                   %Start of FamRank phase of experiment
params.event.FamRank_Phase_end = 91;                  %End of FamRank phase of experiment
params.event.FamRank_Face_Blk_st = 92;
params.event.FamRank_Object_Blk_st = 93;
params.event.FamRank_Blk_end = 94;
params.event.FamRank_Trl_st = 95;
params.event.FamRank_Stim_Onset = 96;
params.event.FamRank_Resp = 97;
params.event.FamRank_Trl_end = 98;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make a Trial Info Struct
% column 1 = Duration_condition (1=Long)
% column 2 = Stimulus_type_condition (1=Fam Obj; 2=Unfam Obj; 3=Fam Face; 4=Unfam Face; 5=Color)
% column 3 = block number
% column 4 = trial number
% column 5 = stim file names, or colors, for each of the 6 stimulus positions
% column 6 = Target stimulus (filename for obj and rgb for color)
% column 7 = Target study position (positions 1-6: position 1 is always at 9 o'clock position and then goes clockwise (as in Brady 2020))
% column 8 = Target test position ('left' or 'right')
% column 9 = Foil stimulus (filename for obj and rgb for color)
% column 10 = subject response ('left' (i.e. button #7) or 'right' (i.e. button #8))
% column 11 = RT in ms
% column 12 = whether the response is correct(1), incorrect(0)
% column 13 = rehearsal digits
% column 14 = foil digits
% column 15 = digits target location ('left' or 'right' of fixation)
% column 16 = Rehearsal task RT
% column 17 = whether the response to the rehearsal digits is correct(1), incorrect(0)
% column 18 = theta angles for all color stimuli
% column 19 = Target-Foil pair Cosine Similarity
Data = struct('Duration_cond',[],'Stimulus_type_cond',[],'Blk_num',[],'Trl_num',[],...
    'Stim',[],'Target',{},'Target_study_pos',[],'Target_test_pos',[],'Foil',{},'Resp',[],'RT',[],...
    'Is_Correct',[],'rehearsal_digits',[],'foil_digits',[],'digit_target_loc',[],'rehearsal_Resp',[],...
    'rehearsal_task_RT',[],'rehearsal_resp_Is_Correct',[],'Color_theta_Info',[],'TF_Sim',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Timing vectors
time.ITI.vblstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.ITI.onset = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.ITI.flipstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.ITI.missed = nan(1,params.nRepsPerCondition*params.Num_conditions);

time.digits_study.vblstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_study.onset = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_study.flipstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_study.missed = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_test.vblstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_test.onset = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_test.flipstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.digits_test.missed = nan(1,params.nRepsPerCondition*params.Num_conditions);

time.stim_study.Long.vblstamp = nan(params.nRepsPerCondition*(params.Num_conditions/2),params.trialSetSize);
time.stim_study.Long.onset = nan(params.nRepsPerCondition*(params.Num_conditions/2),params.trialSetSize);
time.stim_study.Long.flipstamp = nan(params.nRepsPerCondition*(params.Num_conditions/2),params.trialSetSize);
time.stim_study.Long.missed = nan(params.nRepsPerCondition*(params.Num_conditions/2),params.trialSetSize);
time.stim_test.Long.vblstamp = nan(1,params.nRepsPerCondition*(params.Num_conditions/2));
time.stim_test.Long.onset = nan(1,params.nRepsPerCondition*(params.Num_conditions/2));
time.stim_test.Long.flipstamp = nan(1,params.nRepsPerCondition*(params.Num_conditions/2));
time.stim_test.Long.missed = nan(1,params.nRepsPerCondition*(params.Num_conditions/2));

time.WM_delay.vblstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.WM_delay.onset = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.WM_delay.flipstamp = nan(1,params.nRepsPerCondition*params.Num_conditions);
time.WM_delay.missed = nan(1,params.nRepsPerCondition*params.Num_conditions);

%% Load Stim Files and select all stim for each trial
load('TF_pairs_fam.mat')
load('TF_pairs_unfam.mat')
load('TF_pairs_fam_faces_RS_binned.mat')
load('TF_pairs_unfam_faces_RS_binned.mat')
% UpArrow = imread('up_arrow.jpg');
LeftArrow = imread('left_arrow.jpg');
RightArrow = imread('right_arrow.jpg');

Instruction_Message_1 = imread('Instructions_Message_1.jpg');
Instruction_Message_2 = imread('Instructions_Message_2.jpg');
Instruction_Message_3 = imread('Instructions_Message_3.jpg');
Instruction_Message_4 = imread('Instructions_Message_4.jpg');

[color_stim_displays] = Generate_rand_color_displays_minSep45(params.trialSetSize,params.nRepsPerCondition);
[Familiar_display_objs] = Generate_rand_Fam_obj_displays_vWMExp2(TF_pairs_fam);
[Unfamiliar_display_objs] = Generate_rand_Unfam_obj_displays_vWMExp2(TF_pairs_unfam);
[Famous_display_faces1,Famous_display_faces2] = Generate_rand_Fam_face_displays_vWMExp2(TF_pairs_fam_faces_RS_binned,params.nRepsPerCondition,Famous_faces_StimFolder);
[Unfamiliar_display_faces1,Unfamiliar_display_faces2] = Generate_rand_Unfam_face_displays_vWMExp2(TF_pairs_unfam_faces_RS_binned,params.nRepsPerCondition,Unfamiliar_faces_StimFolder );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup portcodes
if params.debug == 0
    if params.port
        config_io;
        params.portcode = hex2dec('0x00003FF8'); 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    %% Setting Up Experiment
    rand('state', sum(100*clock));
    Screen('Preference', 'SkipSyncTests', 0); %Should be set to zero when running experiment
    PsychDefaultSetup(2);
    KbName('UnifyKeyNames'); %used for cross-platform compatibility of keynaming
    ListenChar(2);
    
    % make a dummy call to GetSecs to load the .dll before we need it
    dummy = GetSecs; clear dummy;
    
    %define allowed keyboard responses
    params.buttons.right_bumper = 8;
    params.buttons.left_bumper = 7;
    params.buttons.Familiar = 1;
    params.buttons.Unfamiliar = 3;
    quitKey = KbName('q');
    
    Devices.Logitech = 2;
    Devices.Keyboard = [];
    
    %Get the screen numbers
    screens = Screen('Screens'); %Get the screen numbers
    screenNumber = max(screens);
    %     screenNumber = 0;
    
    %Colors
    black = BlackIndex(screenNumber);
    white = WhiteIndex(screenNumber);
    gray = white / 2;
    
    % Open an on screen window and color it white
    if params.windowed %|| params.debug
        %if we're debugging, open up a small window
        [window,windowRect] = PsychImaging('OpenWindow', screenNumber, params.Background_color,[120 120 1024,768]);
    else
        [window, windowRect] = PsychImaging('OpenWindow', screenNumber, params.Background_color);
        Priority(MaxPriority(window));
        HideCursor();
    end
    
    %Query vertical refreshrate
    params.refreshRate = Screen('GetFlipInterval', window);
    params.refreshCycle = 1/params.refreshRate;
    
    %set to max priority level
    topPriorityLevel = MaxPriority(window);
    
    
    %Find monitor visual angle and pixelangle for stim sizing and placement
    [params.monitor_width,  params.monitor_height] = Screen('DisplaySize', screenNumber);  %monitor width and height in mm
    params.width_pixels = windowRect(3);
    params.x_center = params.width_pixels/2;
    params.height_pixels = windowRect(4);
    params.y_center = params.height_pixels/2;
    params.MonitorVisualAngle = rad2deg(2*atan(params.monitor_width/(2*params.monitor_distance)));
    params.PixelAngle = params.MonitorVisualAngle/params.width_pixels; %Number of visual Degrees per pixel
    params.PixelAngle_arcmin = 60*(params.MonitorVisualAngle/params.width_pixels); %Number of arcmin per pixel
    params.PixelsPerDeg = 1/params.PixelAngle; % number of pixels per visual deg
    
    
    %Convert circRad, stimSize, rectSize, and fixationSize from visual degrees to pixels
    params.circRad_studyStim = round(params.circRad_studyStim * params.PixelsPerDeg);
    params.circRad_testStim = round(params.circRad_testStim * params.PixelsPerDeg);
    params.stimSize = round(params.stimSize * params.PixelsPerDeg);
    params.familiarity_ranking_stimSize = round(params.familiarity_ranking_stimSize *  params.PixelsPerDeg);
    params.rectSize = round(params.rectSize * params.PixelsPerDeg);
    params.fixationSize = round(params.fixationSize * params.PixelsPerDeg);
    params.fixWinSize = round(params.fixWinSize * params.PixelsPerDeg);
    params.arrow_size = round(params.arrow_size * params.PixelsPerDeg);
    
    %set up fixation cross
    fix_xCoords = [-params.fixationSize params.fixationSize 0 0];
    fix_yCoords = [0 0 -params.fixationSize params.fixationSize];
    fixation_Coords = [fix_xCoords; fix_yCoords];
    lineWidthPix = 3;
    
    %set locations for arrow images in familiarity ranking phase
    arrows_theta = [180, 0];%[180, 270, 0];
    baseRect_arrow = [0 0 params.arrow_size params.arrow_size];
    baseRect_text = [0 0 params.rectSize params.rectSize];
    arrow_positions = nan(4,length(arrows_theta));
    Ranking_phase_text_positions = nan(4,length(arrows_theta));
    for i=1:length(arrows_theta)
        X_arrows(i) = round(params.x_center + (params.familiarity_ranking_stimSize/2 + params.arrow_size) .* cosd(arrows_theta(i)));
        Y_arrows(i) = round((params.y_center + params.y_center/4) + (params.familiarity_ranking_stimSize/2 + params.rectSize/1.5) .* sind(arrows_theta(i)));
        arrow_positions(:, i) = CenterRectOnPointd(baseRect_arrow, X_arrows(i), Y_arrows(i));
        
        X_Text(i) = round(params.x_center + (params.familiarity_ranking_stimSize/2 + params.arrow_size + params.rectSize) .* cosd(arrows_theta(i)));
        Y_Text(i) = round((params.y_center + params.y_center/4) + (params.familiarity_ranking_stimSize/2 + params.arrow_size) .* sind(arrows_theta(i)));
        
        Ranking_phase_text_positions(:, i) = CenterRectOnPointd(baseRect_text, X_Text(i), Y_Text(i));
    end
    
    
    %set angles for rect placement
    change_theta = 360/params.trialSetSize;
    rect_theta = 0:change_theta:300;
    rect_theta = circshift(rect_theta,[0 3]); %rearrage so that stim in Sequential condition always start at 9 o'clock position and go clockwise
    test_rect_theta = [180 ,0]; %test position 1 is on the left of fixation, test position 2 is on the right
    
    %set locations for stimuli and surrounding rectangles (study phase)
    baseRect = [0 0 params.rectSize params.rectSize];
    baseRect_stim = [0 0 params.stimSize params.stimSize];
    study_stim_positions = nan(4,params.trialSetSize);
    study_rect_positions = nan(4,params.trialSetSize);
    for i=1:params.trialSetSize
        X_study(i) = round(params.x_center + params.circRad_studyStim .* cosd(rect_theta(i)));
        Y_study(i) = round(params.y_center + params.circRad_studyStim .* sind(rect_theta(i)));
        study_stim_positions(:, i) = CenterRectOnPointd(baseRect_stim, X_study(i), Y_study(i));
        study_rect_positions(:, i) = CenterRectOnPointd(baseRect, X_study(i), Y_study(i));
    end
    
    %set locations for Target and Foil stimuli and surrounding rectangles (test phase)
    test_stim_positions = nan(4,2);
    test_rect_positions = nan(4,2);
    for i=1:2
        X_test(i) = round(params.x_center + params.circRad_testStim .* cosd(test_rect_theta(i)));
        Y_test(i) = round(params.y_center + params.circRad_testStim .* sind(test_rect_theta(i)));
        test_stim_positions(:, i) = CenterRectOnPointd(baseRect_stim, X_test(i), Y_test(i));
        test_rect_positions(:, i) = CenterRectOnPointd(baseRect, X_test(i), Y_test(i));
    end
    
    %initialize counters
    Trl_count = 0;
    FamObj_Trl_count = 0;
    UnfamObj_Trl_count = 0;
    FamFace_Trl_count_long1 = 0;
    FamFace_Trl_count_long2 = 0;
    UnfamFace_Trl_count_long1 = 0;
    UnfamFace_Trl_count_long2 = 0;
    Color_Trl_count = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initiate the Eyetracker.
    edfFile = [];
    
    if params.eyeTrack
        % initialization of the connection with the Eyelink Gazetracker
        if EyelinkInit()~= 1;return;end
        
        % load in default settings
        EyeLinkDefaults=EyelinkInitDefaults(window);
        
        % update defaults
        EyeLinkDefaults.backgroundcolour = params.Background_color;        % match calibration background color trial background color
        EyeLinkDefaults.calibrationtargetcolour = params.calib_target_colors;
        EyeLinkDefaults.cal_target_beep=[600 0 0.05];
        EyeLinkDefaults.drift_correction_target_beep=[600 0 0.05];
        EyeLinkDefaults.calibration_failed_beep=[400 0 0.25];
        EyeLinkDefaults.calibration_success_beep=[800 0 0.25];
        EyeLinkDefaults.drift_correction_failed_beep=[400 0 0.25];
        EyeLinkDefaults.drift_correction_success_beep=[800 0 0.25];
        EyelinkUpdateDefaults(EyeLinkDefaults);                            % update the eyetrack settings
        
        % setup eye tracker
        if params.eyeMode
            % remote mode
            Eyelink('Command', 'elcl_select_configuration = RTABLER'); % remote mode
            Eyelink('Command', 'calibration_type = HV5'); % 5-pt calibration
            Eyelink('Command', 'sample_rate = 500'); % sample at 500 Hz for remote mode
        else
            % chin rest
            Eyelink('Command', 'elcl_select_configuration = BTABLER'); % chin rest
            Eyelink('Command', 'calibration_type = HV9'); % 9-pt calibration
            Eyelink('Command', 'sample_rate = 1000'); % sample at 1000 Hz for chin rest
            Eyelink('Command', 'binocular_enabled = YES');
            Eyelink('Command', 'enable_automatic_calibration = NO'); %Force manual calibration acceptance
            Eyelink('Command', 'heuristic_filter == OFF'); %Turn off Heuristic Filter
            Eyelink('command', 'use_ellipse_fitter = no'); %use Centroid tracking algorithm
        end
        
        % stamp header in EDF file
        Eyelink('Command', 'add_file_preamble_text','vWM_Exp2'); %
        % Setting the proper recording resolution, proper calibration type,
        % as well as the data file content;
        [width, height]=Screen('WindowSize',window);
        Eyelink('Command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, width-1, height-1);
        Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, width-1, height-1);
        % make sure that we get gaze data from the Eyelink
        % set EDF file contents using the file_sample_data and
        % file-event_filter commands
        Eyelink('Command', 'file_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('Command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,HTARGET,GAZERES,STATUS,INPUT');
        % set link data through link_sample_data and link_event_filter
        Eyelink('Command', 'link_event_filter = LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT');
        Eyelink('Command', 'link_sample_data  = LEFT,RIGHT,GAZE,GAZERES,AREA,HTARGET,STATUS,INPUT');
        % proportion commands to adjust size of calibrated area
        Eyelink('Command', 'calibration_area_proportion 0.5 0.5');                                           %Do these two lines need to be changed?
        Eyelink('Command', 'validation_area_proportion 0.5 0.5');
        
        % get host tracker version
        [v,vs]=Eyelink('GetTrackerVersion');
        
        fprintf('Running experiment on a ''%s'' tracker.\n', vs );
        fprintf('Running experiment on version ''%d''.\n', v );
        
        % open file to record data to (filename must be fewer than 8 characters)
        edfFile = ['vWM2_', num2str(SubjectNum),'.edf'];
        Eyelink('Openfile', edfFile);
        
        %make sure we're still connected
        if Eyelink('IsConnected')~=1
            fprint('not connected');
            Eyelink('Shutdown');
            Screen('CloseAll');
            return;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Preexperiment Instructions for Subject
    %Instruction message 1
    Instructions_1 = Screen('MakeTexture', window, Instruction_Message_1);
    
    % Get the size of the image
    [s1, s2, s3] = size(Instruction_Message_1);
    
    % Get the aspect ratio of the image.
    aspectRatio = s2 / s1;
    
    % set the height of each image to the height of the screen
    imageHeight = params.height_pixels * .9;
    imageWidth = imageHeight .* aspectRatio;
    
    % Make the destination rectangles for instructions message image.
    baseRect = [0 0 imageWidth imageHeight];
    InstruRect = CenterRectOnPointd(baseRect, params.x_center,params.y_center);
    
    Screen('DrawTextures', window, Instructions_1,[],InstruRect)
    Screen('DrawingFinished', window);
    Screen('Flip', window);
    
    [Resp, secs, output, Data] = getGamepadResp([1,2,3,4,5,6,7,8,9,10], Devices,quitKey,params,output,Data,1,filename);
    Screen('Close', Instructions_1);
    
    %Instruction message 2
    Instructions_2 = Screen('MakeTexture', window, Instruction_Message_2);
    [s1, s2, s3] = size(Instruction_Message_2);
    aspectRatio = s2 / s1;
    imageHeight = params.height_pixels * .85;
    imageWidth = imageHeight .* aspectRatio;
    baseRect = [0 0 imageWidth imageHeight];
    InstruRect = CenterRectOnPointd(baseRect, params.x_center,params.y_center);
    
    Screen('DrawTextures', window, Instructions_2,[],InstruRect)
    Screen('DrawingFinished', window);
    Screen('Flip', window);
    
    [Resp, secs, output, Data] = getGamepadResp([1,2,3,4,5,6,7,8,9,10], Devices,quitKey,params,output,Data,1,filename);
    Screen('Close', Instructions_2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Practice Trials
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    if params.skip_practice == 0
        [practice_color_stim_displays] = Generate_rand_color_displays_minSep45(params.trialSetSize,2);
        
        for prac_trl = 1:params.NumPracticeTrials
            if prac_trl == 1 || prac_trl == 2
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.ITI);
                
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
            else
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.ITI);
                
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
            end
            
            %%setting up digit rehearsal task
            [~,rehearsal_digits,foil_digits,digit_target_rect,digit_foil_rect,~,~] = setup_Digit_Rehearsal(params,Data,Trl_count,test_stim_positions,window,1,time);
            WaitSecs(params.rehearsal_digits_duration);
            % clear rehearsal digits from screen
            if prac_trl == 1 || prac_trl == 2
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
            else
                % Draw the rects and fixation cross to the screen
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
            end
            % wait 1 second between articulatory suppression/rehearsal task and onset of vWM stimuli
            WaitSecs(params.preTrialDelay);
            
            if prac_trl == 1 %Long color practice trial
                Target_study_pos = practice_color_stim_displays(1).Target_study_pos;
                Screen('FillOval', window, practice_color_stim_displays(1).Study_stimuli, study_stim_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.Long_StimDuration);
                
                % clear stim from display
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.DelayPeriod);
                
                %%2AFC - test phase
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                Test_Colors = zeros(3,2);
                if target_loc == 0
                    Test_Colors(:,1) = practice_color_stim_displays(1).Target';
                    Test_Colors(:,2) = practice_color_stim_displays(1).Foil';
                else
                    Test_Colors(:,1) = practice_color_stim_displays(1).Foil';
                    Test_Colors(:,2) = practice_color_stim_displays(1).Target';
                end
                
                Screen('FillOval', window, Test_Colors, test_stim_positions);
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameOval', window, params.Target_RectFrame_color, study_rect_positions(:,Target_study_pos));
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to vWM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                Screen('Close')
            end
            
            if prac_trl == 2  %Long Familiar Object practice trial
                filePattern = fullfile(Practice_Long_FamObj_StimFolder, '*.jpg');
                practice_imgFiles = dir(filePattern);
                practice_imgFiles(strncmp({practice_imgFiles.name}, '.', 1)) = [];
                Target_stim = imread([practice_imgFiles(6).folder filesep practice_imgFiles(6).name]);
                Foil_stim = imread([practice_imgFiles(3).folder filesep practice_imgFiles(3).name]);
                for j = 1:length(practice_imgFiles)
                    stim{j} = imread([Practice_Long_FamObj_StimFolder filesep practice_imgFiles(j).name]);
                end
                
                %Draw all stim to display
                imageTexture1 = Screen('MakeTexture', window, stim{1});
                imageTexture2 = Screen('MakeTexture', window, stim{6});
                imageTexture3 = Screen('MakeTexture', window, stim{4});
                imageTexture4 = Screen('MakeTexture', window, stim{5});
                imageTexture5 = Screen('MakeTexture', window, stim{2});
                imageTexture6 = Screen('MakeTexture', window, stim{7});
                imageTexture = [imageTexture1; imageTexture2; imageTexture3; imageTexture4; imageTexture5; imageTexture6];
                
                %show display
                Screen('DrawTextures', window, imageTexture, [], study_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                WaitSecs(params.Long_StimDuration);
                
                % clear stim from display
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.DelayPeriod);
                
                %%2AFC - test phase
                target = 2;
                imageTexture_Target = Screen('MakeTexture', window, Target_stim);
                imageTexture_Foil = Screen('MakeTexture', window, Foil_stim);
                
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                if target_loc == 0
                    imageTexture_Test = [imageTexture_Target imageTexture_Foil];
                else
                    imageTexture_Test = [imageTexture_Foil imageTexture_Target];
                end
                Screen('DrawTextures', window, imageTexture_Test, [], test_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameRect', window, params.Target_RectFrame_color, study_rect_positions(:,target));
                Screen('FrameRect', window, params.RectFrame_color, test_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to vWM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                Screen('Close')
            end
            
            if prac_trl ==  3 %Long Unfamiliar Object practice trial
                filePattern = fullfile(Practice_Long_UnfamObj_StimFolder, '*.jpg');
                practice_imgFiles = dir(filePattern);
                practice_imgFiles(strncmp({practice_imgFiles.name}, '.', 1)) = [];
                Target_stim = imread([practice_imgFiles(2).folder filesep practice_imgFiles(2).name]);
                Foil_stim = imread([practice_imgFiles(4).folder filesep practice_imgFiles(4).name]);
                for j = 1:length(practice_imgFiles)
                    stim{j} = imread([Practice_Long_UnfamObj_StimFolder filesep practice_imgFiles(j).name]);
                end
                
                %Draw all stim to display
                imageTexture1 = Screen('MakeTexture', window, stim{2});
                imageTexture2 = Screen('MakeTexture', window, stim{1});
                imageTexture3 = Screen('MakeTexture', window, stim{3});
                imageTexture4 = Screen('MakeTexture', window, stim{5});
                imageTexture5 = Screen('MakeTexture', window, stim{6});
                imageTexture6 = Screen('MakeTexture', window, stim{7});
                imageTexture = [imageTexture1; imageTexture2; imageTexture3; imageTexture4; imageTexture5; imageTexture6];
                
                %show display
                Screen('DrawTextures', window, imageTexture, [], study_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                WaitSecs(params.Long_StimDuration);
                
                % clear stim from display
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.DelayPeriod);
                
                %%2AFC - test phase
                target = 1;
                imageTexture_Target = Screen('MakeTexture', window, Target_stim);
                imageTexture_Foil = Screen('MakeTexture', window, Foil_stim);
                
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                if target_loc == 0
                    imageTexture_Test = [imageTexture_Target imageTexture_Foil];
                else
                    imageTexture_Test = [imageTexture_Foil imageTexture_Target];
                end
                Screen('DrawTextures', window, imageTexture_Test, [], test_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameRect', window, params.Target_RectFrame_color, study_rect_positions(:,target));
                Screen('FrameRect', window, params.RectFrame_color, test_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to vWM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                Screen('Close')
            end
            
            if prac_trl == 4  %Long Familiar Faces practice trial
                filePattern = fullfile(Practice_Long_FamFace_StimFolder, '*.jpg');
                practice_imgFiles = dir(filePattern);
                practice_imgFiles(strncmp({practice_imgFiles.name}, '.', 1)) = [];
                Target_stim = imread([practice_imgFiles(3).folder filesep practice_imgFiles(3).name]);
                Foil_stim = imread([practice_imgFiles(7).folder filesep practice_imgFiles(7).name]);
                for j = 1:length(practice_imgFiles)
                    stim{j} = imread([Practice_Long_FamFace_StimFolder filesep practice_imgFiles(j).name]);
                end
                
                %Draw all stim to display
                imageTexture1 = Screen('MakeTexture', window, stim{1});
                imageTexture2 = Screen('MakeTexture', window, stim{2});
                imageTexture3 = Screen('MakeTexture', window, stim{4});
                imageTexture4 = Screen('MakeTexture', window, stim{5});
                imageTexture5 = Screen('MakeTexture', window, stim{6});
                imageTexture6 = Screen('MakeTexture', window, stim{3});
                imageTexture = [imageTexture1; imageTexture2; imageTexture3; imageTexture4; imageTexture5; imageTexture6];
                
                %show display
                Screen('DrawTextures', window, imageTexture, [], study_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                WaitSecs(params.Long_StimDuration);
                
                % clear stim from display
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.DelayPeriod);
                
                %%2AFC - test phase
                target = 6;
                imageTexture_Target = Screen('MakeTexture', window, Target_stim);
                imageTexture_Foil = Screen('MakeTexture', window, Foil_stim);
                
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                if target_loc == 0
                    imageTexture_Test = [imageTexture_Target imageTexture_Foil];
                else
                    imageTexture_Test = [imageTexture_Foil imageTexture_Target];
                end
                Screen('DrawTextures', window, imageTexture_Test, [], test_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameRect', window, params.Target_RectFrame_color, study_rect_positions(:,target));
                Screen('FrameRect', window, params.RectFrame_color, test_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to vWM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                Screen('Close')
            end
            
            if prac_trl == 5 %Long Unfamiliar Face practice trial
                filePattern = fullfile(Practice_Long_UnfamFace_StimFolder, '*.jpg');
                practice_imgFiles = dir(filePattern);
                practice_imgFiles(strncmp({practice_imgFiles.name}, '.', 1)) = [];
                Target_stim = imread([practice_imgFiles(6).folder filesep practice_imgFiles(6).name]);
                Foil_stim = imread([practice_imgFiles(5).folder filesep practice_imgFiles(5).name]);
                for j = 1:length(practice_imgFiles)
                    stim{j} = imread([Practice_Long_UnfamFace_StimFolder filesep practice_imgFiles(j).name]);
                end
                
                %Draw all stim to display
                imageTexture1 = Screen('MakeTexture', window, stim{2});
                imageTexture2 = Screen('MakeTexture', window, stim{1});
                imageTexture3 = Screen('MakeTexture', window, stim{6});
                imageTexture4 = Screen('MakeTexture', window, stim{3});
                imageTexture5 = Screen('MakeTexture', window, stim{6});
                imageTexture6 = Screen('MakeTexture', window, stim{7});
                imageTexture = [imageTexture1; imageTexture2; imageTexture3; imageTexture4; imageTexture5; imageTexture6];
                
                %show display
                Screen('DrawTextures', window, imageTexture, [], study_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                WaitSecs(params.Long_StimDuration);
                
                % clear stim from display
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                WaitSecs(params.DelayPeriod);
                
                %%2AFC - test phase
                target = 3;
                imageTexture_Target = Screen('MakeTexture', window, Target_stim);
                imageTexture_Foil = Screen('MakeTexture', window, Foil_stim);
                
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                if target_loc == 0
                    imageTexture_Test = [imageTexture_Target imageTexture_Foil];
                else
                    imageTexture_Test = [imageTexture_Foil imageTexture_Target];
                end
                Screen('DrawTextures', window, imageTexture_Test, [], test_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameRect', window, params.Target_RectFrame_color, study_rect_positions(:,target));
                Screen('FrameRect', window, params.RectFrame_color, test_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to vWM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                Screen('Flip', window);
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,1,filename);
                Screen('Close')
            end
        end
    end
    %Instruction message 3
    Instructions_3 = Screen('MakeTexture', window, Instruction_Message_3);
    [s1, s2, s3] = size(Instruction_Message_2);
    aspectRatio = s2 / s1;
    imageHeight = params.height_pixels * .85;
    imageWidth = imageHeight .* aspectRatio;
    baseRect = [0 0 imageWidth imageHeight];
    InstruRect = CenterRectOnPointd(baseRect, params.x_center,params.y_center);
    
    Screen('DrawTextures', window, Instructions_3,[],InstruRect)
    Screen('DrawingFinished', window);
    Screen('Flip', window);
    [Resp, secs, output, Data] = getGamepadResp([1,2,3,4,5,6,7,8,9,10], Devices,quitKey,params,output,Data,1,filename);
    Screen('Close', Instructions_3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Block Loops
    
    % send start-of-experiment trigger to EEG
    if params.port
        outp(params.portcode,params.event.Exp_st);
    end
    
    for Blk_count = 1:params.Num_Blocks
        Durations = 'LONG';
        StimDur = params.Long_StimDuration;
        
        
        if Blk_count > 1
            %Inter-block Break Screen
            Screen('TextSize',window, params.text_font_size + 30);
            break_message = 'Break';
            DrawFormattedText(window, break_message, 'center', params.y_center/3, [0 0 0]);
            Screen('TextSize',window, params.text_font_size);
            message = ['Next Block is ' num2str(Blk_count) ' of ' num2str(params.Num_Blocks)...
                '\n\n\n\nPress any button to continue'];
            DrawFormattedText(window, message, 'center', 'center', [0 0 0]);
            Screen('DrawingFinished', window);
            Screen('Flip', window);
            
            % Wait for a button press to continue
            [Resp, secs, output, Data] = getGamepadResp([1,2,3,4,5,6,7,8,9,10], Devices,quitKey,params,output,Data,1,filename);
            Screen('Close');
        else
            %Inform participant block is ablout to begin
            Screen('TextSize',window, params.text_font_size);
            message = 'Block 1 is about to begin.';
            DrawFormattedText(window, message, 'center', 'center', [0 0 0]);
            Screen('DrawingFinished', window);
            Screen('Flip', window);
            WaitSecs(2);
            Screen('Close');
        end
        
        % Setup eyetracker/Calibration
        if params.port
            outp(params.portcode,params.event.EyeTracking_calib_st);% send start-of-calibraation trigger to EEG
        end
        if params.eyeTrack
            EyelinkDoTrackerSetup(EyeLinkDefaults);
        end
        if params.port
            outp(params.portcode,params.event.EyeTracking_calib_end);% send end-of-calibration trigger to EEG
        end
        
        %Send Block Condition type triggers to EEG
        if params.port
            if Cond_Order(Blk_count) == 1
                outp(params.portcode,params.event.Long_Fam_Obj_Blk_st);% send LongFamObj Block Start trigger to EEG
            elseif Cond_Order(Blk_count) == 2
                outp(params.portcode,params.event.Long_Unfam_Obj_Blk_st);% send LongUnfamObj Block Start trigger to EEG
            elseif Cond_Order(Blk_count) == 3
                outp(params.portcode,params.event.Long_Fam_Face_Blk_st);% send LongFamFace Block Start trigger to EEG
            elseif Cond_Order(Blk_count) == 4
                outp(params.portcode,params.event.Long_Unfam_Face_Blk_st);% send LongUnfamFace Block Start trigger to EEG
            else 
                outp(params.portcode,params.event.Long_Color_Blk_st);% send LongColor Block Start trigger to EEG
            end
        end
        
        
        %%Trial Loops
        if params.debug
            NumTrials = 1;
        else
            NumTrials = params.Num_Trials;
        end
        for trl = 1:NumTrials
            Trl_count = Trl_count + 1;
            
            % tell eye tracker that the trial is beginning.
            if params.eyeTrack
                %Send a 'TRIALID' message to mark the start of a trial in
                %DataViewer. The viewer will not parse any messages,
                %events, or samples in the datafile prior to this message.
                Eyelink('Message', 'TRIALID %d', Trl_count);
                
                % set idle mode before start of recording
                Eyelink('Command', 'set_idle_mode');
                WaitSecs(0.05);
                % Put block, trial number at the bottom of operater display
                Eyelink('Command', 'record_status_message "BLOCK %d TRIAL %d"', Blk_count, Trl_count)
                % start recording
                Eyelink('StartRecording');
                Eyelink('Message', 'BLOCK %d ', Blk_count);
                % mark zero-plot time in EDF file
                Eyelink('Message', 'TrialStart');  % may be a lag in recording (not a problem because the baseline screen is much longer than needed.
                eye_used = Eyelink('EyeAvailable'); % get eye that's tracked
                % returns 0 (LEFT_EYE), 1 (RIGHT_EYE) or 2 (BINOCULAR) depending on what data is
                if eye_used == 2
                    eye_used = 1; % use the right_eye data
                end
            end
            % send start-of-trial trigger to EEG
            if params.port
                outp(params.portcode,params.event.Trl_st);
            end
            
            %Set up display depending on whether curent block is color or obj/Face block
            if Cond_Order(Blk_count,2) ~= 5 %Obj/Face block
                if trl == 1
                    % Draw the rects and fixation to the screen
                    Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    Screen('DrawingFinished',window);
                    [time.ITI.vblstamp(Trl_count), time.ITI.onset(Trl_count), time.ITI.flipstamp(Trl_count), time.ITI.missed(Trl_count)]= Screen('Flip', window);
                    Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    %%setting up rehearsal task
                    [Data,rehearsal_digits,foil_digits,digit_target_rect,digit_foil_rect,digit_target_pos,time] = setup_Digit_Rehearsal(params,Data,Trl_count,test_stim_positions,window,0,time,(time.ITI.onset(Trl_count)+params.ITI-0.5*params.refreshRate));
                else
                    Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                    Screen('DrawingFinished',window);
                    [time.ITI.vblstamp(Trl_count), time.ITI.onset(Trl_count), time.ITI.flipstamp(Trl_count), time.ITI.missed(Trl_count)]= Screen('Flip', window);
                    Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    Screen('DrawingFinished',window);
                    Screen('Flip', window,(time.ITI.onset(Trl_count)+params.Blink_time-0.5*params.refreshRate));
                    Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    %%setting up rehearsal task
                    [Data,rehearsal_digits,foil_digits,digit_target_rect,digit_foil_rect,digit_target_pos,time] = setup_Digit_Rehearsal(params,Data,Trl_count,test_stim_positions,window,0,time,(time.ITI.onset(Trl_count)+(params.ITI+params.Blink_time)-0.5*params.refreshRate));
                end
                
            else %Color Block
                if trl == 1
                    Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    Screen('DrawingFinished',window);
                    [time.ITI.vblstamp(Trl_count), time.ITI.onset(Trl_count), time.ITI.flipstamp(Trl_count), time.ITI.missed(Trl_count)]= Screen('Flip', window);
                    Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    %%setting up rehearsal task
                    [Data,rehearsal_digits,foil_digits,digit_target_rect,digit_foil_rect,digit_target_pos,time] = setup_Digit_Rehearsal(params,Data,Trl_count,test_stim_positions,window,0,time,(time.ITI.onset(Trl_count)+params.ITI-0.5*params.refreshRate));
                else
                    Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                    Screen('DrawingFinished',window);
                    [time.ITI.vblstamp(Trl_count), time.ITI.onset(Trl_count), time.ITI.flipstamp(Trl_count), time.ITI.missed(Trl_count)]= Screen('Flip', window);
                    Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    Screen('DrawingFinished',window);
                    Screen('Flip', window,(time.ITI.onset(Trl_count)+params.Blink_time-0.5*params.refreshRate));
                    Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                    Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                    %%setting up rehearsal task
                    [Data,rehearsal_digits,foil_digits,digit_target_rect,digit_foil_rect,digit_target_pos,time] = setup_Digit_Rehearsal(params,Data,Trl_count,test_stim_positions,window,0,time,(time.ITI.onset(Trl_count)+(params.ITI+params.Blink_time)-0.5*params.refreshRate));
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clear rehearsal digits from screen
            if Cond_Order(Blk_count,2) ~= 5 %Obj/Face block
                % Draw the rects and fixation cross to the screen
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [~, fix_onset] = Screen('Flip', window,(time.digits_study.onset(Trl_count)+(params.rehearsal_digits_duration)-0.5*params.refreshRate));
            else %Color block
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [~, fix_onset] = Screen('Flip', window,(time.digits_study.onset(Trl_count)+(params.rehearsal_digits_duration)-0.5*params.refreshRate));
            end
            if params.port
                outp(params.portcode,params.event.ArticSupressTask_Study_Stim_Offset);% send SupressTask Stim Offset trigger to EEG
            end
            if params.eyeTrack
                Eyelink('message','ArticSupressTask_Study_Stim_Offset');
            end
            
            % wait ~1 second between articulatory suppression/rehearsal task and onset of vWM stimuli
            WaitSecs('UntilTime',fix_onset+(params.preTrialDelay)-0.5*params.refreshRate);
            gazeWinStart = GetSecs; %  Start timer for fixation: 50ms fixation needed to start vWM study phase
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % wait for gaze to be within central fixation window for 50ms before starting vWM portion of trial
            if params.eyeTrack
                while 1 % loop till error or space bar is pressed
                    % Check if a new sample is available online via the link.
                    if Eyelink('NewFloatSampleAvailable') > 0
                        evt = Eyelink('NewestFloatSample');
                        x = evt.gx(eye_used+1);
                        y = evt.gy(eye_used+1);
                        if x >= (params.x_center - params.fixWinSize/2) && x <= (params.x_center + params.fixWinSize/2)...
                                && y >= (params.y_center - params.fixWinSize/2) && y <= (params.y_center + params.fixWinSize/2)
                            if (GetSecs - gazeWinStart)*1000 >= params.fixateTime % If gaze duration >= minimum fixation window time (fxateTime)
                                break; % break while loop to show stimulus
                            end
                        elseif ~(x >= (params.x_center - params.fixWinSize/2) && x <= (params.x_center + params.fixWinSize/2)...
                                && y >= (params.y_center - params.fixWinSize/2) && y <= (params.y_center + params.fixWinSize/2))
                            gazeWinStart = GetSecs; % Reset fixation window timer
                        end
                    end
                    % Give option to Start trial with quit key just in case code gets stuck due to lost eye-tracking
                    [~, ~, keyCode] = KbCheck(Devices.Keyboard);
                    if keyCode(quitKey)
                        % Write message to EDF file to mark the space bar press time
                        Eyelink('Message', 'FIXATION_Aborted');
                        ShowCursor;
                        sca;
                        Priority(0);
                        Eyelink('StopRecording');
                        Eyelink('CloseFile');
                        Eyelink('Shutdown');
                        
                        output.Data = Data;
                        output.params = params;
                        output.ExperimentalSession_end_time = clock;
                        output.ExperimentalSession_dur = etime(clock,output.ExperimentalSession_start_time);
                        save(filename,'output')
                        
                        % Retrieve edf file from tracker
                        try
                            fprintf('Receiving data file ''%s''\n', edfFile );
                            status=Eyelink('ReceiveFile');
                            if status > 0
                                fprintf('ReceiveFile status %d\n', status);
                            end
                            if 2==exist(edfFile, 'file')
                                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
                            end
                        catch
                            fprintf('Problem receiving data file ''%s''\n', edfFile );
                        end
                        
                        % send end-of-experiment trigger to EEG
                        if params.port
                            outp(params.portcode,params.event.Exp_end);
                        end
                        error('User quit');
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Start vWM portion of Trial            
            if Cond_Order(Blk_count,2) == 5 %Color Cond
                Color_Trl_count = Color_Trl_count + 1;
                Target_study_pos = color_stim_displays(Color_Trl_count).Target_study_pos;
                Data(Trl_count).Target = color_stim_displays(Color_Trl_count).Target;
                Data(Trl_count).Foil = color_stim_displays(Color_Trl_count).Foil;
                Data(Trl_count).Target_study_pos = Target_study_pos;
                Data(Trl_count).Stim = num2cell(color_stim_displays(Color_Trl_count).Study_stimuli',2)';
                Data(Trl_count).Color_theta_Info = color_stim_displays(Color_Trl_count).stim_theta_angles;
                Screen('FillOval', window, color_stim_displays(Color_Trl_count).Study_stimuli, study_stim_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [time.stim_study.Long.vblstamp(Trl_count), time.stim_study.Long.onset(Trl_count), time.stim_study.Long.flipstamp(Trl_count), time.stim_study.Long.missed(Trl_count)]=Screen('Flip', window);
                if params.eyeTrack
                    Eyelink('Message','Stim_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.Stim_Onset);
                end
                
                % clear stim from display
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [time.WM_delay.vblstamp(Trl_count), time.WM_delay.onset(Trl_count), time.WM_delay.flipstamp(Trl_count), time.WM_delay.missed(Trl_count)]=Screen('Flip', window,(time.stim_study.Long.onset(Trl_count)+StimDur)-0.5*params.refreshRate);
                if params.eyeTrack
                    Eyelink('Message','Delay_Period_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.Delay_Period_Onset);
                end
                
                %%2AFC - test phase
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                Test_Colors = zeros(3,2);
                if target_loc == 0
                    Test_Colors(:,1) = color_stim_displays(Color_Trl_count).Target';
                    Test_Colors(:,2) = color_stim_displays(Color_Trl_count).Foil';
                    Data(Trl_count).Target_test_pos = 'Left';
                else
                    Test_Colors(:,1) = color_stim_displays(Color_Trl_count).Foil';
                    Test_Colors(:,2) = color_stim_displays(Color_Trl_count).Target';
                    Data(Trl_count).Target_test_pos = 'Right';
                end
                
                Screen('FillOval', window, Test_Colors, test_stim_positions);
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameOval', window, params.Target_RectFrame_color, study_rect_positions(:,Target_study_pos));
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [time.stim_test.Long.vblstamp(Trl_count), time.stim_test.Long.onset(Trl_count), time.stim_test.Long.flipstamp(Trl_count), time.stim_test.Long.missed(Trl_count)]=Screen('Flip', window,(time.WM_delay.onset(Trl_count)+params.DelayPeriod)-0.5*params.refreshRate);
                if params.eyeTrack
                    Eyelink('Message','Test_Stim_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.Test_Stim_Onset)
                end
                
                %Get Response to WM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,0,filename);
                Data(Trl_count).RT = round((secs - time.stim_test.Long.onset(Trl_count)) * 1000,3);
                [Data] = CheckResp(Data,Resp, Trl_count,target_loc,1);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameOval', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                [time.digits_test.vblstamp(Trl_count), time.digits_test.onset(Trl_count), time.digits_test.flipstamp(Trl_count), time.digits_test.missed(Trl_count)]=Screen('Flip', window);
                if params.eyeTrack
                    Eyelink('Message','ArticSupressTask_Test_Stim_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.ArticSupressTask_Test_Stim_Onset);
                end
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,0,filename);
                
                % send end-of-trial trigger to EEG
                if params.port
                    outp(params.portcode,params.event.Trl_end);
                end
                
                Data(Trl_count).rehearsal_task_RT = round((secs - time.digits_test.onset(Trl_count)) * 1000,3);
                [Data] = CheckResp(Data,Resp, Trl_count,digit_target_pos,0);
                Screen('Close');
            else
                %Set variables for the 4 non-color conditions
                if Cond_Order(Blk_count,2) == 1 %FamObj Cond
                    FamObj_Trl_count = FamObj_Trl_count + 1;
                    count = FamObj_Trl_count;
                    StimFolder = Familiar_objects_StimFolder;
                    display_stim = Familiar_display_objs;
                    
                elseif Cond_Order(Blk_count,2) == 2 %UnfamObj Cond
                    UnfamObj_Trl_count = UnfamObj_Trl_count + 1;
                    count = UnfamObj_Trl_count;
                    StimFolder = Unfamiliar_objects_StimFolder;
                    display_stim = Unfamiliar_display_objs;
                    
                elseif Cond_Order(Blk_count,2) == 3 %FamFace Cond
                    if Cond_Order(Blk_count,1) == 1
                        FamFace_Trl_count_long1 = FamFace_Trl_count_long1 + 1;
                        count = FamFace_Trl_count_long1;
                        display_stim = Famous_display_faces1;
                    else
                        FamFace_Trl_count_long2 = FamFace_Trl_count_long2 + 1;
                        count = FamFace_Trl_count_long2;
                        display_stim = Famous_display_faces2;
                    end
                    StimFolder = Famous_faces_StimFolder;
                    
                else %UnfamFace Cond
                    if Cond_Order(Blk_count,1) == 1
                        UnfamFace_Trl_count_long1 = UnfamFace_Trl_count_long1 + 1;
                        count = UnfamFace_Trl_count_long1;
                        display_stim = Unfamiliar_display_faces1;
                    else
                        UnfamFace_Trl_count_long2 = UnfamFace_Trl_count_long2 + 1;
                        count = UnfamFace_Trl_count_long2;
                        display_stim = Unfamiliar_display_faces2;
                    end
                    StimFolder = Unfamiliar_faces_StimFolder;
                end
                
                %All non-Color conditions
                for j = 1:params.trialSetSize
                    if Cond_Order(Blk_count,2)==1 | Cond_Order(Blk_count,2)==2
                        stim{j} = imread([StimFolder display_stim(count).Study_stimuli_Cat{j} filesep display_stim(count).Study_stimuli{j}]);
                        Data(Trl_count).Stim{j} = display_stim(count).Study_stimuli{j};
                        Foil = imread([StimFolder display_stim(count).Foil_Cat{:} filesep display_stim(count).Foil{:}]);
                    else
                        stim{j} = imread([StimFolder display_stim(count).Study_stimuli{j}]);
                        Data(Trl_count).Stim{j} = display_stim(count).Study_stimuli{j};
                        Foil = imread([StimFolder filesep display_stim(count).Foil]);
                    end
                end
                
                %Draw all stim to display
                imageTexture1 = Screen('MakeTexture', window, stim{1});
                imageTexture2 = Screen('MakeTexture', window, stim{2});
                imageTexture3 = Screen('MakeTexture', window, stim{3});
                imageTexture4 = Screen('MakeTexture', window, stim{4});
                imageTexture5 = Screen('MakeTexture', window, stim{5});
                imageTexture6 = Screen('MakeTexture', window, stim{6});
                imageTexture = [imageTexture1; imageTexture2; imageTexture3; imageTexture4; imageTexture5; imageTexture6];
                
                %show display
                Screen('DrawTextures', window, imageTexture, [], study_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [time.stim_study.Long.vblstamp(Trl_count), time.stim_study.Long.onset(Trl_count), time.stim_study.Long.flipstamp(Trl_count), time.stim_study.Long.missed(Trl_count)]=Screen('Flip', window);
                if params.eyeTrack
                    Eyelink('Message','Stim_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.Stim_Onset);
                end
                
                % clear stim from display
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [time.WM_delay.vblstamp(Trl_count), time.WM_delay.onset(Trl_count), time.WM_delay.flipstamp(Trl_count), time.WM_delay.missed(Trl_count)]=Screen('Flip', window,(time.stim_study.Long.onset(Trl_count)+StimDur)-0.5*params.refreshRate);
                if params.eyeTrack
                    Eyelink('Message','Delay_Period_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.Delay_Period_Onset);
                end
                
                %%2AFC - test phase
                target = display_stim(count).Target_study_pos;
                imageTexture_Target = Screen('MakeTexture', window, stim{target});
                Data(Trl_count).Target = display_stim(count).Target;
                Data(Trl_count).Target_study_pos = target;
                Data(Trl_count).TF_Sim = display_stim(count).TF_Sim;
                Data(Trl_count).Foil = display_stim(count).Foil;
                imageTexture_Foil = Screen('MakeTexture', window, Foil);
                
                %Randomize position of Target and Foil
                target_loc = randi([0 1]); % 0=position Left of fixation, 1=position Right of fixation
                if target_loc == 0
                    imageTexture_Test = [imageTexture_Target imageTexture_Foil];
                    Data(Trl_count).Target_test_pos = 'Left';
                else
                    imageTexture_Test = [imageTexture_Foil imageTexture_Target];
                    Data(Trl_count).Target_test_pos = 'Right';
                end
                Screen('DrawTextures', window, imageTexture_Test, [], test_stim_positions);
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('FrameRect', window, params.Target_RectFrame_color, study_rect_positions(:,target));
                Screen('FrameRect', window, params.RectFrame_color, test_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.fixationColor, [params.x_center params.y_center], 2);
                Screen('DrawingFinished',window);
                [time.stim_test.Long.vblstamp(Trl_count), time.stim_test.Long.onset(Trl_count), time.stim_test.Long.flipstamp(Trl_count), time.stim_test.Long.missed(Trl_count)]=Screen('Flip', window,(time.WM_delay.onset(Trl_count)+params.DelayPeriod)-0.5*params.refreshRate);
                if params.eyeTrack
                    Eyelink('Message','Test_Stim_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.Test_Stim_Onset)
                end
                
                %%Get response to vWM task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,0,filename);
                Data(Trl_count).RT = round((secs - time.stim_test.Long.onset(Trl_count)) * 1000,3);
                [Data] = CheckResp(Data,Resp, Trl_count,target_loc,1);
                
                % Draw digits for test phase of digit rehearsal task to the display
                Screen('TextSize', window, params.rehearsal_digits_font_size);
                Screen('TextFont', window, 'Times');
                Screen('FrameRect', window, params.RectFrame_color, study_rect_positions);
                Screen('DrawLines', window, fixation_Coords, lineWidthPix, params.Blink_fixationColor, [params.x_center params.y_center], 2);
                DrawFormattedText(window, rehearsal_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_target_rect);
                DrawFormattedText(window, foil_digits, 'center', 'center', [0 0 0],[],[],[],[],[],digit_foil_rect);
                Screen('DrawingFinished',window);
                [time.digits_test.vblstamp(Trl_count), time.digits_test.onset(Trl_count), time.digits_test.flipstamp(Trl_count), time.digits_test.missed(Trl_count)]=Screen('Flip', window);
                if params.eyeTrack
                    Eyelink('Message','ArticSupressTask_Test_Stim_Onset');
                end
                if params.port
                    outp(params.portcode,params.event.ArticSupressTask_Test_Stim_Onset);
                end
                
                %%Get response to digit rehearsal task
                [Resp, secs, output, Data] = getGamepadResp([params.buttons.left_bumper, params.buttons.right_bumper], Devices,quitKey,params,output,Data,0,filename);
                
                % send end-of-trial trigger to EEG
                if params.port
                    outp(params.portcode,params.event.Trl_end);
                end
                
                Data(Trl_count).rehearsal_task_RT = round((secs - time.digits_test.onset(Trl_count)) * 1000,3);
                [Data] = CheckResp(Data,Resp, Trl_count,digit_target_pos,0);
                Screen('Close');
            end
            
            % stop recording eye data for current trial
            if params.eyeTrack
                Eyelink('Message', 'TrialEnd');
                Eyelink('StopRecording');
            end
            
            % send end-of-trial trigger to EEG
            if params.port
                outp(params.portcode,params.event.Trl_end);
            end
            
            Data(Trl_count).Duration_cond = Cond_Order(Blk_count,1);
            Data(Trl_count).Stimulus_type_cond = Cond_Order(Blk_count,2);
            Data(Trl_count).Blk_num = Blk_count;
            Data(Trl_count).Trl_num = Trl_count;
        end
        
        % send end-of-block trigger to EEG and message to EyeLink
        if params.port
            outp(params.portcode,params.event.Blk_End);
        end
        if params.eyeTrack
            Eyelink('Message', 'BlockEnd');
        end
    end
    WaitSecs(.01);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Familiarity Ranking Phase
    if params.port
        outp(params.portcode,params.event.FamRank_Phase_st);
    end
    if params.eyeTrack
        Eyelink('Message', 'FamRank_Phase_Start');
    end
    
    %% Familiarity Ranking Phase Params
    FamRank_params.Stim_type = [1 2];                           % 1=Faces, 2=Objects
    FamRank_params.Stim_familiarity = [1 2];                    % 1=Familiar/Famous, 2=Unfamiliar
    FamRank_params.Num_Blocks = 2;                              % Faces and Objects
    FamRank_params.nRepsPerBlock = 200;                         % Number of trials per block (100 Familiar stim and 100 Unfamiliar randomly mixed in each block)
    FamRank_params.stimSize = 5;                                % size of the stimulus in visual degrees
    %% Make a Trial Info Struct
    FamRank_Data = struct('Blk_num',[],'Trl_num',[],'Stim',{},'Stim_type',[],'Stim_familiarity',[],'Stim_Cat',{},'Resp',[],'RT',[]);
    
    %% Timing Info
    FamRank_time.stim.vblstamp = nan(1,FamRank_params.nRepsPerBlock*FamRank_params.Num_Blocks);
    FamRank_time.stim.onset = nan(1,FamRank_params.nRepsPerBlock*FamRank_params.Num_Blocks);
    FamRank_time.stim.flipstamp = nan(1,FamRank_params.nRepsPerBlock*FamRank_params.Num_Blocks);
    FamRank_time.stim.missed = nan(1,FamRank_params.nRepsPerBlock*FamRank_params.Num_Blocks);
    
    %% Randomly select 100 stim from each stimulus type and familiarity and then shuffle together the Faces trials and Object trials
    %giving 2 blocks of 200 trials each
    
    Fam_Objects_Filenames = dir([Familiar_objects_StimFolder filesep '*' filesep '*.jpg']);
    Fam_Objects_Filenames = Fam_Objects_Filenames(~startsWith({Fam_Objects_Filenames.name},'.'));
    Unfam_Objects_Filenames = dir([Unfamiliar_objects_StimFolder filesep '*' filesep '*.jpg']);
    Unfam_Objects_Filenames = Unfam_Objects_Filenames(~startsWith({Unfam_Objects_Filenames.name},'.'));
    Fam_Faces_Filenames = dir([Famous_faces_StimFolder '*.jpg']);
    Unfam_Faces_Filenames = dir([Unfamiliar_faces_StimFolder '*.jpg']);
    
    Fam_Obj_idx = datasample([1:560],100,'Replace',false);
    Unfam_Obj_idx = datasample([1:560],100,'Replace',false);
    Fam_Face_idx = datasample([1:560],100,'Replace',false);
    Unfam_Face_idx = datasample([1:560],100,'Replace',false);
    
    if ismac
        slash = '/';
    else
        slash = '\';
    end
    
    %Randomly select 100 stim per cond
    Fam_Obj2use = table();
    Unfam_Obj2use = table();
    Fam_Face2use = table();
    Unfam_Face2use = table();
    
    
    for i = 1:100
        Fam_Obj2use.stim_path{i} = fullfile(Fam_Objects_Filenames(Fam_Obj_idx(i)).folder, Fam_Objects_Filenames(Fam_Obj_idx(i)).name);
        Fam_Obj2use.cat{i} = Fam_Objects_Filenames(Fam_Obj_idx(i)).folder(max(strfind(Fam_Objects_Filenames(Fam_Obj_idx(i)).folder,slash))+1:end);
        
        Unfam_Obj2use.stim_path{i} = fullfile(Unfam_Objects_Filenames(Unfam_Obj_idx(i)).folder, Unfam_Objects_Filenames(Unfam_Obj_idx(i)).name);
        Unfam_Obj2use.cat{i} = Unfam_Objects_Filenames(Unfam_Obj_idx(i)).folder(max(strfind(Unfam_Objects_Filenames(Unfam_Obj_idx(i)).folder,slash))+1:end);
        
        Fam_Face2use.stim_path{i} = fullfile(Fam_Faces_Filenames(Fam_Face_idx(i)).folder, Fam_Faces_Filenames(Fam_Face_idx(i)).name);
        if strcmp(Fam_Faces_Filenames(Fam_Face_idx(i)).name(end-5),'L')
            Fam_Face2use.cat{i} = ['W' Fam_Faces_Filenames(Fam_Face_idx(i)).name(end-4)];
        else
            Fam_Face2use.cat{i} = Fam_Faces_Filenames(Fam_Face_idx(i)).name(end-5:end-4);
        end
        
        Unfam_Face2use.stim_path{i} = fullfile(Unfam_Faces_Filenames(Unfam_Face_idx(i)).folder, Unfam_Faces_Filenames(Unfam_Face_idx(i)).name);
        cat_idx = strfind(Unfam_Faces_Filenames(Unfam_Face_idx(i)).name,'-');
        if strcmp(Unfam_Faces_Filenames(Unfam_Face_idx(i)).name(cat_idx(1)+1),'L')
            Unfam_Face2use.cat{i} = ['W' Unfam_Faces_Filenames(Unfam_Face_idx(i)).name(cat_idx(2)-1)];
        else
            Unfam_Face2use.cat{i} = Unfam_Faces_Filenames(Unfam_Face_idx(i)).name(cat_idx(1)+1:cat_idx(2)-1);
        end
    end
    
    Face_Stim = [Fam_Face2use;Unfam_Face2use];
    Obj_Stim = [Fam_Obj2use;Unfam_Obj2use];
    
    rand_face_idx = randperm(200);
    rand_obj_idx = randperm(200);
    total_count = 0;
    face_count = 0;
    obj_count = 0;
    FamRank_Cond_Order = randi([0,1]);
    if FamRank_Cond_Order == 1 %Block1 is objects Block2 is faces
        for Blk = 1:FamRank_params.Num_Blocks
            for trl = 1:FamRank_params.nRepsPerBlock
                total_count = total_count + 1;
                FamRank_Data(total_count).Blk_num = Blk;
                FamRank_Data(total_count).Trl_num = trl;
                
                if Blk == 1
                    obj_count = obj_count + 1;
                    FamRank_Data(total_count).Stim = Obj_Stim.stim_path{rand_obj_idx(obj_count)};
                    FamRank_Data(total_count).Stim_type = 2;
                    if rand_obj_idx(obj_count) <=100
                        FamRank_Data(total_count).Stim_familiarity = 1;
                    else
                        FamRank_Data(total_count).Stim_familiarity = 2;
                    end
                    FamRank_Data(total_count).Stim_Cat = Obj_Stim.cat{rand_obj_idx(obj_count)};
                else
                    face_count = face_count + 1;
                    FamRank_Data(total_count).Stim = Face_Stim.stim_path{rand_face_idx(face_count)};
                    FamRank_Data(total_count).Stim_type = 1;
                    if rand_face_idx(face_count) <=100
                        FamRank_Data(total_count).Stim_familiarity = 1;
                    else
                        FamRank_Data(total_count).Stim_familiarity = 2;
                    end
                    FamRank_Data(total_count).Stim_Cat = Face_Stim.cat{rand_face_idx(face_count)};
                end
            end
        end
        
    else%Block1 is faces Block2 is objects
        for Blk = 1:FamRank_params.Num_Blocks
            for trl = 1:FamRank_params.nRepsPerBlock
                total_count = total_count + 1;
                FamRank_Data(total_count).Blk_num = Blk;
                FamRank_Data(total_count).Trl_num = trl;
                
                if Blk == 1
                    face_count = face_count + 1;
                    FamRank_Data(total_count).Stim = Face_Stim.stim_path{rand_face_idx(face_count)};
                    FamRank_Data(total_count).Stim_type = 1;
                    if rand_face_idx(face_count) <=100
                        FamRank_Data(total_count).Stim_familiarity = 1;
                    else
                        FamRank_Data(total_count).Stim_familiarity = 2;
                    end
                    FamRank_Data(total_count).Stim_Cat = Face_Stim.cat{rand_face_idx(face_count)};
                else
                    obj_count = obj_count + 1;
                    FamRank_Data(total_count).Stim = Obj_Stim.stim_path{rand_obj_idx(obj_count)};
                    FamRank_Data(total_count).Stim_type = 2;
                    if rand_obj_idx(obj_count) <=100
                        FamRank_Data(total_count).Stim_familiarity = 1;
                    else
                        FamRank_Data(total_count).Stim_familiarity = 2;
                    end
                    FamRank_Data(total_count).Stim_Cat = Obj_Stim.cat{rand_obj_idx(obj_count)};
                end
            end
        end
    end
    
    
    %% Setting Up FamRanking Phase of Experiment
    
    %Convert stimSize from visual degrees to pixels
    FamRank_params.stimSize = round(FamRank_params.stimSize * params.PixelsPerDeg);
    baseRect_stim_familiarity_ranking = [0 0 FamRank_params.stimSize FamRank_params.stimSize];
    ranking_stim_position = CenterRectOnPointd(baseRect_stim_familiarity_ranking, params.x_center, (params.y_center+(params.y_center/4)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Start FamRank Trials
    
    %Instruction message 4
    Instructions_4 = Screen('MakeTexture', window, Instruction_Message_4);
    [s1, s2, s3] = size(Instruction_Message_2);
    aspectRatio = s2 / s1;
    imageHeight = params.height_pixels * .85;
    imageWidth = imageHeight .* aspectRatio;
    baseRect = [0 0 imageWidth imageHeight];
   
    InstruRect = CenterRectOnPointd(baseRect, params.x_center,params.y_center);
    
    Screen('DrawTextures', window, Instructions_4,[],InstruRect)
    Screen('DrawingFinished', window);
    Screen('Flip', window);
    
    [Resp, secs, output, Data] = getGamepadResp([1,2,3,4,5,6,7,8,9,10], Devices,quitKey,params,output,Data,1,filename);
    Screen('Close',Instructions_4);
    
    %     i=0;
    %     while i < length(FamRank_Data)
    %         i = i + 1;
    FamRank_trl_count = 0;
    for Blk = 1:2
        if params.port
            if FamRank_Cond_Order == 1
                if Blk == 1
                    outp(params.portcode,params.event.FamRank_Object_Blk_st);% send Object Block Trigger
                else
                    outp(params.portcode,params.event.FamRank_Face_Blk_st);% send Face Block Trigger
                end
            else
                if Blk == 1
                    outp(params.portcode,params.event.FamRank_Face_Blk_st);% send Face Block Trigger
                else
                    outp(params.portcode,params.event.FamRank_Object_Blk_st);% send Object Block Trigger
                end
            end
        end
        
        for trl = 1:200
            FamRank_trl_count = FamRank_trl_count + 1;
            
            % send start-of-trial trigger to EEG
            if params.port
                outp(params.portcode,params.event.FamRank_Trl_st);
            end
            
            stim = imread(FamRank_Data(FamRank_trl_count).Stim);
            
            %arrows
            imageTexture1 = Screen('MakeTexture', window, LeftArrow);
            imageTexture2 = Screen('MakeTexture', window, RightArrow);
            imageTexture_arrows = [imageTexture1; imageTexture2];
            
            %Draw object image to screen
            Screen('TextSize',window, params.text_font_size + 30);
            message1 = 'Are you Familiar or Unfamiliar with this Person/Object?';
            DrawFormattedText(window, message1, 'center', params.y_center/3, [0 0 0]);
            Screen('TextSize',window, params.text_font_size);
            message2 = 'Familiar';
            message3 = 'Unfamiliar';
            DrawFormattedText(window, message2, 'center', 'center', [0 0 0],[],[],[],[],[],Ranking_phase_text_positions(:,1)');
            DrawFormattedText(window, message3, 'center', 'center', [0 0 0],[],[],[],[],[],Ranking_phase_text_positions(:,2)');
            imageTexture = Screen('MakeTexture', window, stim);
            Screen('DrawTextures', window, imageTexture,[],ranking_stim_position);
            Screen('DrawTextures', window, imageTexture_arrows,[], arrow_positions);
            Screen('DrawingFinished',window);
            [FamRank_time.stim.vblstamp(FamRank_trl_count), FamRank_time.stim.onset(FamRank_trl_count), FamRank_time.stim.flipstamp(FamRank_trl_count), FamRank_time.stim.missed(FamRank_trl_count)] = Screen('Flip', window);
                
            if params.port
                outp(params.portcode,params.event.FamRank_Stim_Onset);
            end
            
            %%Get response to Familiary rating
            [Resp, secs, output, Data] = getGamepadResp([params.buttons.Familiar, params.buttons.Unfamiliar], Devices,quitKey,params,output,Data,0,filename);
            FamRank_Data(FamRank_trl_count).RT = round((secs - FamRank_time.stim.onset(FamRank_trl_count)) * 1000,3);
            Screen('Close');
            
            if Resp == 1 %Familiar
                FamRank_Data(FamRank_trl_count).Resp = 1;
            elseif Resp == 3 %Unfamiliar
                FamRank_Data(FamRank_trl_count).Resp = 2;
            end

            % send end-of-trial trigger to EEG
            if params.port
                outp(params.portcode,params.event.FamRank_Trl_end);
            end
            WaitSecs(.01);
        end
    end
    
    if params.port
        outp(params.portcode,params.event.FamRank_Phase_end);
    end
    if params.eyeTrack
        Eyelink('Message', 'FamRank_Phase_Start');
    end
    
    WaitSecs(.01);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % send end-of-experiment trigger to EEG
    if params.port
        outp(params.portcode,params.event.Exp_end);
    end
    
    ShowCursor;
    fclose('all');
    Priority(0);
    Screen('CloseAll');
    ListenChar(0);
    output.ExperimentalSession_end_time = clock;
    
    %Save Data
    output.Data = Data;
    ExperimentalSession_dur = seconds(etime(output.ExperimentalSession_end_time,output.ExperimentalSession_start_time));
    output.ExperimentalSession_dur = ExperimentalSession_dur;
    output.params = params;
    output.timing_info = time;
    output.FamRank_Data = FamRank_Data;
    output.FamRank_params = FamRank_params;
    output.FamRank_timing_info = FamRank_time;
    save(filename,'output')
    
    if params.eyeTrack
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(.1);
        Eyelink('CloseFile');
        
        % download data file
        try
            fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        catch
            fprintf('Problem receiving data file ''%s''\n', edfFile );
        end
        Eyelink('Shutdown');
    end
catch
    Screen('CloseAll');
    ShowCursor;
    fclose('all');
    Priority(0);
    ListenChar(0);
    
    if params.eyeTrack
        Eyelink('Command', 'set_idle_mode');
        WaitSecs(.1);
        Eyelink('CloseFile');
        
        % download data file
        try
            fprintf('Receiving data file ''%s''\n', edfFile );
            status=Eyelink('ReceiveFile');
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(edfFile, 'file')
                fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
            end
        catch
            fprintf('Problem receiving data file ''%s''\n', edfFile );
        end
        Eyelink('Shutdown');
    end
    
    psychrethrow(psychlasterror);
    return;
end