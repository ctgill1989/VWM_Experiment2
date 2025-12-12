function [Resp, secs, output, Data] = getGamepadResp(Resp_buttons,Devices,quitKey,params,output,Data,practice,filename)
%CTGill
%last updated 8/23/21
%
% set practice=1 for practice trial  and practice=0 for experimental trial
%
%

respToBeMade = true;
while respToBeMade == true
    [ButtonIsDown,secs, keyCode] = KbCheck(Devices.Logitech);
    if ButtonIsDown
        if length(find(keyCode>0)) == 1
            if sum(find(keyCode>0) == Resp_buttons)
                if params.eyeTrack
                    Eyelink('message','Resp'); %send Resp trigger to EyeLink
                end
                if params.port
                    outp(params.portCode,params.event.Resp); %send Resp trigger to EEG
                end
                
                Resp = Resp_buttons(find(keyCode(Resp_buttons))); %Determine Response
                
                while 1
                    [ButtonIsDown,~, ~] = KbCheck(Devices.Logitech);
                    
                    if ~ButtonIsDown
                        respToBeMade = false;
                        break;
                    end
                end
            end
        end
    end
    
    [KeyIsDown,~, keyCode] = KbCheck(Devices.Keyboard);
    if KeyIsDown
        if keyCode(quitKey)
            ShowCursor;
            Screen('CloseAll');
            Priority(0);
            if params.eyeTrack
                Eyelink('StopRecording');
                Eyelink('CloseFile');
                Eyelink('Shutdown');
            end
            
            if practice == 0
                output.Data = Data;
                output.params = params;
                output.ExperimentalSession_end_time = clock;
                output.ExperimentalSession_dur = etime(clock,output.ExperimentalSession_start_time);
                output.timing_info = time;
                save(filename,'output')
                
                % send end-of-experiment trigger to EEG
                if params.port
                    outp(params.portcode,params.event.Exp_end);
                end
            end
            error('User quit');
        end
    end
end

