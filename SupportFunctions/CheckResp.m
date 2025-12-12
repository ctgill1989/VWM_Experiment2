function [Data] = CheckResp(Data,Resp, Trl_count,target_loc,Task_type)
%CTGill
%Last updated 6/18/2021
%
% set Task_type=0 for Articulatory suppression task (i.e. digit rehearal) responses
%     Task_type=1 for 2AFC vWM task responses
%     Task_type=2 for Change Detection vWM task responses

if Task_type == 1 %2AFC WM task responses
    if Resp == 7 %Left bumper
        Data(Trl_count).Resp = 'Left';
        if target_loc == 0
            Data(Trl_count).Is_Correct = 1;
        else
            Data(Trl_count).Is_Correct = 0;
        end
    elseif Resp == 8 %Right bumper
        Data(Trl_count).Resp = 'Right';
        if target_loc == 1
            Data(Trl_count).Is_Correct = 1;
        else
            Data(Trl_count).Is_Correct = 0;
        end
    end
elseif Task_type == 0 %Articulatory suppression task responses
    if Resp == 7 %Left bumper
        Data(Trl_count).rehearsal_Resp = 'Left';
        if target_loc == 0
            Data(Trl_count).rehearsal_resp_Is_Correct = 1;
        else
            Data(Trl_count).rehearsal_resp_Is_Correct = 0;
        end
    elseif Resp == 8 %Right bumper
        Data(Trl_count).rehearsal_Resp = 'Right';
        if target_loc == 1
            Data(Trl_count).rehearsal_resp_Is_Correct = 1;
        else
            Data(Trl_count).rehearsal_resp_Is_Correct = 0;
        end
    end
    
    elseif Task_type == 2 %Change Detection vWM task responses
    if Resp == 7 %Left bumper
        Data(Trl_count).Resp = 'Same';
        if target_loc == 0 %No Change
            Data(Trl_count).Is_Correct = 1;
        else
            Data(Trl_count).Is_Correct = 0;
        end
    elseif Resp == 8 %Right bumper
        Data(Trl_count).Resp = 'Different';
        if target_loc == 1
            Data(Trl_count).Is_Correct = 1;
        else
            Data(Trl_count).Is_Correct = 0;
        end
    end
end