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

%% RELAX_blinks_IQR_method:

function [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg)
% This function uses an inter-quartile range threshold to create a mask specifying the location of
% eye blinks in the continuous EEG trace

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'BlinkMaskFocus')==0
            RELAX_cfg.BlinkMaskFocus=150; % this value decides how much data before and after the right and left base of the eye blink to mark as part of the blink artifact window.
        end
        if isfield(RELAX_cfg, 'BlinkElectrodes')==0
            RELAX_cfg.BlinkElectrodes={'FP1';'FPZ';'FP2';'AF3';'AF4';'F3';'F1';'FZ';'F2';'F4'}; % sets the electrodes to average for blink detection using the IQR method
        end
        if isfield(RELAX_cfg, 'LowPassFilterAt_6Hz_BeforeDetectingBlinks')==0
            RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks='no'; % low pass filters the data @ 6Hz prior to blink detection (helps if high power alpha is disrupting blink detection, not necessary in the vast majority of cases)
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.BlinkMaskFocus=150; % this value decides how much data before and after the right and left base of the eye blink to mark as part of the blink artifact window. 
        RELAX_cfg.BlinkElectrodes={'FP1';'FPZ';'FP2';'AF3';'AF4';'F3';'F1';'FZ';'F2';'F4'}; % sets the electrodes to average for blink detection using the IQR method
        RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks='no'; % low pass filters the data @ 6Hz prior to blink detection (helps if high power alpha is disrupting blink detection, not necessary in the vast majority of cases)
    end
           
    % Robust average re-reference EEG data for more reliable blink detection 
    EEG=continuousEEG;
    [EEG] = RELAX_average_rereference(EEG); % (only used for blink detection in this step, data output from this function is left with the reference montage previously set)
    
    % Create mask full of zeros:
    continuousEEG.RELAXProcessing.Details.eyeblinkmask=zeros(1,EEG.pnts);
    
    %% BLINK DETECTION VIA THRESHOLDING WITH 75TH% + 3*INTERQUARTILE INTERVAL METHOD.
    continuousEEG.RELAX.IQRmethodDetectedBlinks=0;
    % Obtain list of non-blink affected channels, then remove them before
    % averaging blink affected channels together to detect blinks:
    NonBlinkChannelList=struct2cell(EEG.chanlocs);
    NonBlinkChannelList=squeeze(NonBlinkChannelList(1,1,:));
    for s=1:size(RELAX_cfg.BlinkElectrodes,1)
        NonBlinkChannelList(strcmpi(NonBlinkChannelList(:),RELAX_cfg.BlinkElectrodes(s)))=[]; % RELAX v1.1.3 NWB adjusted to allow flexibility in the case of strings referring to electrode names
    end
    EEGEyeOnly=pop_select(EEG,'nochannel',NonBlinkChannelList);
    Message = ['electrodes removed here only to average blink affected electrodes to enable blink detection, data still contains ', num2str(continuousEEG.nbchan), ' electrodes'];
    disp(Message);
    % Use TESA to apply butterworth filter: 
    if strcmp(RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks,'no')
        EEGEyeOnly = RELAX_filtbutter( EEGEyeOnly, 1, 25, 4, 'bandpass' );
    elseif strcmp(RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks,'yes')
        EEGEyeOnly = RELAX_filtbutter( EEGEyeOnly, 1, 6, 4, 'bandpass' );
    end

    Message = ['Filtering here only performed to better detect blinks, output data will still be bandpass filtered from ', num2str(RELAX_cfg.HighPassFilter), ' to ', num2str(RELAX_cfg.LowPassFilter)];
    disp(Message);
    % Make values in extreme outlying periods = 0 to help blink detection perform
    % better (not done in final output data, just in data fed to blinker or the IQR blink detection method):
    if isfield(continuousEEG, 'RELAX')==1
        if isfield(continuousEEG.RELAX, 'ExtremelyBadPeriodsForDeletion')==1
            for x=1:size(continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion,1)
                EEGEyeOnly.data(:,continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion(x,1):continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion(x,2))=0;
            end
        end
    end
    
    EEGEyeOnly.data=mean(EEGEyeOnly.data,1); % This could be changed to median, which would increase robustness against bad channels or outliers (but wouldn't work so well if including electrodes barely affected by blinks)
    InterQuartileRangeAllTimepointsAndEpochs=iqr(EEGEyeOnly.data(:,:),2);
    Upper25 = prctile(EEGEyeOnly.data,75,2);
    UpperBound=squeeze(Upper25+(3*InterQuartileRangeAllTimepointsAndEpochs)); % sets the threshold, above which a period is assumed to have a blink

    for x=1:size(EEGEyeOnly.data,2)
        if EEGEyeOnly.data(1,x)<UpperBound
            BlinkIndexMetric(1,x)=0; 
        elseif EEGEyeOnly.data(1,x)>UpperBound
            BlinkIndexMetric(1,x)=1; 
        end
    end

    % Check that blinks exceed the IQR threshold for more
    % than 50ms, and to detect the blink peak within the period that
    % exceeds the threshold:
    ix_blinkstart=find(diff(BlinkIndexMetric)==1)+1;  % indices where BlinkIndexMetric goes from 0 to 1
    ix_blinkend=find(diff(BlinkIndexMetric)==-1);  % indices where BlinkIndexMetric goes from 1 to 0

    if ~isempty(ix_blinkstart)
        if ix_blinkend(1,1)<ix_blinkstart(1,1); ix_blinkend(:,1)=[]; end % if the first downshift occurs before the upshift, remove the first value in end
        if ix_blinkend(1,size(ix_blinkend,2))<ix_blinkstart(1,size(ix_blinkstart,2)); ix_blinkstart(:,size(ix_blinkstart,2))=[];end % if the last upshift occurs after the last downshift, remove the last value in start
        BlinkThresholdExceededLength=ix_blinkend-ix_blinkstart; % length of consecutive samples where blink threshold was exceeded
        BlinkRunIndex = find(BlinkThresholdExceededLength>round(50/RELAX_cfg.ms_per_sample)); % find locations where blink threshold was exceeded by more than 50ms
        % find latency of the max voltage within each period where the blink
        % threshold was exceeded:
        if size(BlinkRunIndex,2)>0
            continuousEEG.RELAX.IQRmethodDetectedBlinks=1;
            epochedEEG.RELAX.IQRmethodDetectedBlinks=1;
            for x=1:size(BlinkRunIndex,2)
                o=ix_blinkstart(BlinkRunIndex(x));
                c=ix_blinkend(BlinkRunIndex(x));
                [~,I]=max(EEGEyeOnly.data(1,o:c),[],2);
                BlinkMaxLatency(1,x)=o+I;
            end
        end
    end

    %%

    % Mark the blink peak into the event list:
    if continuousEEG.RELAX.IQRmethodDetectedBlinks==1
        for x=1:size(BlinkRunIndex,2) % v1.1.3 update, removed an erroneous "- 1" after size(BlinkRunIndex,2), thankyou Mana Biabani for the update
            if (BlinkMaxLatency(1,x)-400>0) && (BlinkMaxLatency(1,x)+400)<size(EEG.data,2)
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkLeftBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkMax';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkRightBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample);
                continuousEEG.RELAXProcessing.Details.eyeblinkmask(round(BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample)):round(BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample)))=1;
            else
                continuousEEG.RELAX.IQRmethodDetectedBlinks=0; % v1.1.3 update, so later functions do not search for blinks that aren't marked, if the blink was only present right on the edge of the file. Thankyou Mana Biabani for the update
            end
        end
    end    

    continuousEEG=eeg_checkset(continuousEEG,'eventconsistency');        
    continuousEEG.RELAX.eyeblinkmask=continuousEEG.RELAXProcessing.Details.eyeblinkmask;
    epochedEEG.RELAX.eyeblinkmask=continuousEEG.RELAX.eyeblinkmask;
    epochedEEG.RELAX.IQRmethodDetectedBlinks=continuousEEG.RELAX.IQRmethodDetectedBlinks;
end