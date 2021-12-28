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
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.BlinkMaskFocus=150; % this value decides how much data before and after the right and left base of the eye blink to mark as part of the blink artifact window. 
        RELAX_cfg.BlinkElectrodes={'FP1';'FPZ';'FP2';'AF3';'AF4';'F3';'F1';'FZ';'F2';'F4'}; % sets the electrodes to average for blink detection using the IQR method
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
        NonBlinkChannelList(strcmp(NonBlinkChannelList(:),RELAX_cfg.BlinkElectrodes(s)))=[];
    end
    EEGEyeOnly=pop_select(EEG,'nochannel',NonBlinkChannelList);
    Message = ['electrodes removed here only to average blink affected electrodes to enable blink detection, data still contains ', num2str(continuousEEG.nbchan), ' electrodes'];
    disp(Message);
    % Use TESA to apply butterworth filter: 
    EEGEyeOnly = pop_tesa_filtbutter( EEGEyeOnly, 1, 25, 4, 'bandpass' );
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

    % Use RunLength to check that blinks exceed the IQR threshold for more
    % than 50ms, and to detect the blink peak within the period that
    % exceeds the threshold
    [b, n, bi] = RunLength(BlinkIndexMetric);
    BlinkRunIndex = find(b>0 & n>round(50/RELAX_cfg.ms_per_sample));
    if size(BlinkRunIndex,2)>0
        continuousEEG.RELAX.IQRmethodDetectedBlinks=1;
        epochedEEG.RELAX.IQRmethodDetectedBlinks=1;
        for x=1:size(BlinkRunIndex,2)-1
            o=bi(BlinkRunIndex(x));
            c=bi(BlinkRunIndex(x)+1);
            [~,I]=max(EEGEyeOnly.data(1,o:c),[],2);
            BlinkMaxLatency(1,x)=o+I;
        end
    end
    %%

    % Mark the blink peak into the event list:
    if continuousEEG.RELAX.IQRmethodDetectedBlinks==1
        for x=1:size(BlinkRunIndex,2)-1
            if (BlinkMaxLatency(1,x)-400>0) && (BlinkMaxLatency(1,x)+400)<size(EEG.data,2)
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkLeftBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkMax';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x);
                continuousEEG.event(size(continuousEEG.event,2)+1).type='EyeBlinkRightBase';
                continuousEEG.event(size(continuousEEG.event,2)).latency=BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample);
                continuousEEG.RELAXProcessing.Details.eyeblinkmask(round(BlinkMaxLatency(1,x)-(400/RELAX_cfg.ms_per_sample)):round(BlinkMaxLatency(1,x)+(400/RELAX_cfg.ms_per_sample)))=1;
            end
        end
    end    
    
    continuousEEG=eeg_checkset(continuousEEG,'eventconsistency');        
    continuousEEG.RELAX.eyeblinkmask=continuousEEG.RELAXProcessing.Details.eyeblinkmask;
    epochedEEG.RELAX.eyeblinkmask=continuousEEG.RELAX.eyeblinkmask;
    epochedEEG.RELAX.IQRmethodDetectedBlinks=continuousEEG.RELAX.IQRmethodDetectedBlinks;
end