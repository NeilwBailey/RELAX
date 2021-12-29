%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2022) Neil Bailey

%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see https://www.gnu.org/licenses/.

%% RELAX_epoching:

function [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg)

    %% This function splits the data into 1s epochs with 0.5s overlap, enabling the next functions to check for extreme outlying periods and to create the MWF cleaning template

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'OnlyIncludeTaskRelatedEpochs')==0
            RELAX_cfg.OnlyIncludeTaskRelatedEpochs=0; % if this is set, epochs are only created if they are within 5s of a task trigger are 
        end
        if isfield(RELAX_cfg, 'ms_per_sample')==0
            RELAX_cfg.ms_per_sample=(1000/EEG.srate);
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.OnlyIncludeTaskRelatedEpochs=0;  
        RELAX_cfg.ms_per_sample=(1000/EEG.srate);
    end

    %% Take a copy of the continuous EEG data for use in later steps:
    continuousEEG=EEG;
    
    %% Epoch 1 second of data every 500ms, and optionally exclude epochs that aren't within 5 seconds of a task trigger:    
    EEG.RELAXProcessing.Details.NaNsForNonEvents=NaN(1,EEG.pnts);
    EEG=eeg_regepochs(EEG,'recurrence',0.500,'eventtype','X','extractepochs','off');
    
    for e=1:size(EEG.event,2)
        EEG.event(e).originallatency=EEG.event(e).latency;
    end

    % if only including task related epochs, change X triggers (which will
    % be tested for artifacts) to Y triggers (which will be ignored), when
    % the trigger is not within 5s of a task related trigger:
    if RELAX_cfg.OnlyIncludeTaskRelatedEpochs==1
        event=EEG.event;
        for e=11:size(EEG.event,2)-11
            if (strcmp(EEG.event(e-10).type, 'X'))&&(strcmp(EEG.event(e-9).type, 'X'))&&(strcmp(EEG.event(e-8).type, 'X'))&&(strcmp(EEG.event(e-7).type, 'X'))...
                    &&(strcmp(EEG.event(e-6).type, 'X'))&&(strcmp(EEG.event(e-5).type, 'X'))&&(strcmp(EEG.event(e-4).type, 'X'))&&(strcmp(EEG.event(e-3).type, 'X'))...
                    &&(strcmp(EEG.event(e-2).type, 'X'))&&(strcmp(EEG.event(e-1).type, 'X'))&&(strcmp(EEG.event(e+1).type, 'X'))&&(strcmp(EEG.event(e+2).type, 'X'))...
                    &&(strcmp(EEG.event(e+3).type, 'X'))&&(strcmp(EEG.event(e+4).type, 'X'))&&(strcmp(EEG.event(e+5).type, 'X'))&&(strcmp(EEG.event(e+6).type, 'X'))...
                    &&(strcmp(EEG.event(e+7).type, 'X'))&&(strcmp(EEG.event(e+8).type, 'X'))&&(strcmp(EEG.event(e+9).type, 'X'))&&(strcmp(EEG.event(e+10).type, 'X'))
                event(e).type = 'Y';
            end
        end
        EEG.event=event;
    end
    EEG = pop_selectevent( EEG, 'type', 'X', 'deleteevents','on');
    % Change X's to Y's if they're the first or last 11 events, because you can't check if they're in a task related
    % section of the data, and they're often noisy and will contaminate potential MWF masks:
    for e=1:11
        if (strcmp(EEG.event(e).type, 'X'))
            EEG.event(e).type = 'Y';
        end
    end
    for e=size(EEG.event,2)-11:size(EEG.event,2)
        if (strcmp(EEG.event(e).type, 'X'))
            EEG.event(e).type = 'Y';
        end
    end
    %Epoch data
    EEG = pop_epoch( EEG, {'X'}, [0 1.0], 'epochinfo', 'yes');
    EEG = eeg_checkset( EEG );
    EEG = pop_selectevent( EEG, 'type', 'X', 'omitlatency', '1<=1999', 'deleteevents','on');
    
    % Create the template for the mask, with NaNs outside of the task
    % related periods, and 0's in the task related periods (the 0's are
    % interpreted as clean data in the mask for MWF, and are replaced with
    % 1's when artifacts are detected in the artifact detection and marking
    % functions:
    OneSecondOf0s=zeros(1,round(1000/RELAX_cfg.ms_per_sample));        
    for e=1:size(EEG.event,2)
        if (strcmp(EEG.event(e).type, 'X'))
            EEG.RELAXProcessing.Details.NaNsForNonEvents(EEG.event(e).originallatency:EEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf0s;
        end
    end
    EEG.RELAXProcessing.Details.NoiseMaskFullLength=EEG.RELAXProcessing.Details.NaNsForNonEvents;
    continuousEEG.RELAXProcessing.Details=EEG.RELAXProcessing.Details;
    epochedEEG=EEG;
    if isfield(epochedEEG, 'RELAX')==1
        continuousEEG.RELAX=epochedEEG.RELAX;
    end    
end