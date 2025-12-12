function [Data,rehearsal_digits,foil_digits,digit_target_rect,digit_foil_rect,digit_target_pos,time] = setup_Digit_Rehearsal(params,Data,Trl_count,test_stim_positions,window,is_Practice,time,Wait_time)
%CTGill
%last updated 7/21/21
%
% is_Practice = 1 indicates this is a practice trial
% is_Practice = 0 indicates that this is a real trial

% Select random digits for rehearsal task
rand_digits = datasample([1:9],2,'Replace',false);
study_digit_1 = rand_digits(1);
study_digit_2 = rand_digits(2);
rehearsal_digits = num2str([study_digit_1 study_digit_2]);

if study_digit_1 == 1 & study_digit_2 == 2
    test_digit_1 = study_digit_1;
    test_digit_2 = study_digit_2 + 1;
elseif study_digit_1 == 2 & study_digit_2 == 1
    test_digit_1 = study_digit_1 + 1;
    test_digit_2 = study_digit_2;
elseif study_digit_1 == 8 & study_digit_2 == 9
    test_digit_1 = study_digit_1 - 1;
    test_digit_2 = study_digit_2;
elseif study_digit_1 == 9 & study_digit_2 == 8
    test_digit_1 = study_digit_1;
    test_digit_2 = study_digit_2 - 1;
else
    %decide which digit is getting changed 1st or 2nd
    temp = randi([0 1]);
    
    if temp == 0 %change 1st digit
        if minus(study_digit_1,study_digit_2) == -1 | study_digit_1 == 9
            test_digit_1 = study_digit_1 - 1;
        elseif minus(study_digit_1,study_digit_2) == 1 | study_digit_1 == 1
            test_digit_1 = study_digit_1 + 1;
        else
            %randomly choose if 1st digit should increase or decrease by 1
            temp2 = randi([0 1]);
            if temp2 == 0
                test_digit_1 = study_digit_1 + 1;
            else
                test_digit_1 = study_digit_1 - 1;
            end
        end
        test_digit_2 = study_digit_2;
        
    else %change 2nd digit
        if minus(study_digit_1,study_digit_2) == -1 | study_digit_2 == 1
            test_digit_2 = study_digit_2 + 1;
        elseif minus(study_digit_1,study_digit_2) == 1 | study_digit_2 == 9
            test_digit_2 = study_digit_2 - 1;
        else
             %randomly choose if 1st digit should increase or decrease by 1
            temp2 = randi([0 1]);
            if temp2 == 0
                test_digit_2 = study_digit_2 + 1;
            else
                test_digit_2 = study_digit_2 - 1;
            end
        end
        test_digit_1 = study_digit_1;
    end
end

foil_digits = num2str([test_digit_1 test_digit_2]);
if ~is_Practice 
    Data(Trl_count).rehearsal_digits = rehearsal_digits;
    Data(Trl_count).foil_digits = foil_digits;
end

%randomly select position of target and foil digits for test
digit_target_pos = randi([0 1]);
if digit_target_pos == 0 %target digits on left of fixation
    digit_target_rect = test_stim_positions(:,1)';
    digit_foil_rect = test_stim_positions(:,2)';
    if ~is_Practice 
        Data(Trl_count).digit_target_loc = 'Left';
    end
else %target digits on right of fixation
    digit_target_rect = test_stim_positions(:,2)';
    digit_foil_rect = test_stim_positions(:,1)';
    if ~is_Practice 
        Data(Trl_count).digit_target_loc = 'Right';
    end
end
% Draw digits for study phase of rehearsal task to the display
params.rehearsal_digits_pretrial_yloc = params.y_center*1.1;
Screen('TextSize', window, params.rehearsal_digits_font_size);
Screen('TextFont', window, 'Times');
DrawFormattedText(window, rehearsal_digits, 'center', params.rehearsal_digits_pretrial_yloc, [0 0 0]);
Screen('DrawingFinished',window);
if ~is_Practice
    [time.digits_study.vblstamp(Trl_count), time.digits_study.onset(Trl_count), time.digits_study.flipstamp(Trl_count), time.digits_study.missed(Trl_count)]=Screen('Flip', window,Wait_time);
else
    Screen('Flip', window);
end
if ~is_Practice
    if params.eyeTrack
        Eyelink('message','ArticSupressTask_Study_Stim_Onset');
    end
    if params.port
        outp(params.portcode,params.event.ArticSupressTask_Study_Stim_Onset);% send SupressTask Stim Onet trigger to EEG
    end
end
