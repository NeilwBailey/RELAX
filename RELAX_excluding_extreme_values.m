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

%% RELAX_excluding_extreme_values:

function [continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg)

    %% This function marks epochs/periods that exceed the extreme outlier detection thresholds 

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'ExtremeVoltageShiftThreshold')==0
            RELAX_cfg.ExtremeVoltageShiftThreshold=25; % Threshold MAD from the median all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe voltage shifts within the epoch
        end
        if isfield(RELAX_cfg, 'ExtremeImprobableVoltageDistributionThreshold')==0
            RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; % Threshold SD from the mean of all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe improbable data
        end
        if isfield(RELAX_cfg, 'ms_per_sample')==0
            RELAX_cfg.ms_per_sample=(1000/epochedEEG.srate); 
        end
        if isfield(RELAX_cfg, 'ExtremeSingleChannelKurtosisThreshold')==0
            RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=8; % Threshold kurtosis of each electrode against the same electrode in different epochs. This could be set lower and would catch less severe kurtosis 
        end
        if isfield(RELAX_cfg, 'ExtremeAllChannelKurtosisThreshold')==0
            RELAX_cfg.ExtremeAllChannelKurtosisThreshold=8; % Threshold kurtosis across all electrodes. This could be set lower and would catch less severe kurtosis  
        end
        if isfield(RELAX_cfg, 'ExtremeAbsoluteVoltageThreshold')==0
            RELAX_cfg.ExtremeAbsoluteVoltageThreshold=500; % microvolts max or min above which will be excluded from cleaning and deleted from data
        end
        if isfield(RELAX_cfg, 'ExtremeBlinkShiftThreshold')==0
            RELAX_cfg.ExtremeBlinkShiftThreshold=8; % How many MAD from the median across blink affected epochs to exclude as extreme data 
            % (applies the higher value out of this value and the
            % RELAX_cfg.ExtremeVoltageShiftThreshold above as the
            % threshold, which caters for the fact that blinks don't affect
            % the median, so without this, if data is clean and blinks are
            % large, blinks can get excluded as extreme outliers)
        end
        if isfield(RELAX_cfg, 'ExtremeDriftSlopeThreshold')==0
            RELAX_cfg.ExtremeDriftSlopeThreshold=-4; % slope of log frequency log power below which to reject as drift without neural activity
        end
        if isfield(RELAX_cfg, 'LineNoiseOrNotchFilterAffectedData')==0
            RELAX_cfg.LineNoiseOrNotchFilterAffectedData=1; % Removes influence of line noise (or line noise filtering) from the muscle slope computations.
        end
    elseif exist('RELAX_cfg', 'var')==0     
        RELAX_cfg.ExtremeVoltageShiftThreshold=25; %MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
        RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; %SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
        RELAX_cfg.ms_per_sample=(1000/epochedEEG.srate);
        RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=8; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis 
        RELAX_cfg.ExtremeAllChannelKurtosisThreshold=8; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis  
        RELAX_cfg.ExtremeAbsoluteVoltageThreshold=500; % microvolts max or min above which will be excluded from cleaning and deleted from data
        RELAX_cfg.ExtremeBlinkShiftThreshold=3; % How many MAD from the median of blink affected epochs to exclude as extreme data
        RELAX_cfg.ExtremeDriftSlopeThreshold=-4; % slope of log frequency log power below which to reject as drift without neural activity
        RELAX_cfg.LineNoiseOrNotchFilterAffectedData=1;
    end
    %% Detect voltage shift in Epoch to identify outlying channels:
    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        % Mark blinks so the function is informed if an epoch contains
        % a blink (which will increase it's amplitude shift, but not
        % reflect an extreme artifact that should be deleted:
        [continuousEEG_initialblinkmarking, ~] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg);
        % Eliminate non-blink triggers (just for this step) and order the event struct by latency
        BlinkShiftUpperBound=zeros(size(continuousEEG_initialblinkmarking.data,1),1);
        if isfield(continuousEEG_initialblinkmarking, 'RELAX')
            if isfield(continuousEEG_initialblinkmarking.RELAX, 'IQRmethodDetectedBlinks')
                if continuousEEG_initialblinkmarking.RELAX.IQRmethodDetectedBlinks==1
                    EEGEyeOnly = pop_selectevent( continuousEEG_initialblinkmarking, 'type',{'EyeBlinkMax'},'deleteevents','on');
                    EEGEyeOnly = eeg_checkset(EEGEyeOnly, 'eventconsistency');
                    % Epoch data around blinks only:
                    EEGEyeOnly = pop_epoch( EEGEyeOnly, {'EyeBlinkMax'}, [-0.5 0.5], 'epochinfo', 'yes');
                    EEGEyeOnly = eeg_checkset( EEGEyeOnly );
                    % Compute blink affected epoch amplitude shifts and
                    % outlier threshold for blink affected epochs:
                    BlinkAmplitudeShiftWithinEachEpoch=range(EEGEyeOnly.data(:,:,:),2);
                    BlinkMedianAmplitudeShift=median(BlinkAmplitudeShiftWithinEachEpoch(:,:,:),3);
                    BlinkMADAmplitudeShift=mad(BlinkAmplitudeShiftWithinEachEpoch(:,:,:),1,3);
                    BlinkShiftUpperBound=BlinkMedianAmplitudeShift+(RELAX_cfg.ExtremeBlinkShiftThreshold*BlinkMADAmplitudeShift);
                end
            end
        end
    end

    % The following extracts both the (max - min) in the epoch for each
    % electrode, compares that absolute (max - min) voltage shift to all
    % other epochs for that electrode, and marks as an extreme outlier
    % voltage epochs with shifts more than the threshold:
    epochedEEG.RELAXProcessing.Details.EpochsExceedingShiftThreshold=zeros;
    epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold=zeros;
    epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections=zeros;
    AmplitudeShiftWithinEachEpoch=range(epochedEEG.data(:,:,:),2);
    % The following is an alternative method to create an upper bound
    % that accounts for potential blinks, if the IQR blink detection
    % method has failed. It calculates the median and MAD of the 20%
    % largest amplitude shifts within an epoch (the outcome of which
    % generally matches well to using the IQR method to detect blink
    % epochs)
    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        if continuousEEG_initialblinkmarking.RELAX.IQRmethodDetectedBlinks~=1
            ExtractingProbableBlinkAffectedEpochs=AmplitudeShiftWithinEachEpoch;
            ExtractingProbableBlinkAffectedEpochs(prctile(ExtractingProbableBlinkAffectedEpochs,80,3)>ExtractingProbableBlinkAffectedEpochs)=NaN;
            ProbableBlinkMedians=median(ExtractingProbableBlinkAffectedEpochs,3,'omitnan');
            ProbableBlinkMAD=mad(ExtractingProbableBlinkAffectedEpochs,1,3);
            BlinkShiftUpperBound=ProbableBlinkMedians+(RELAX_cfg.ExtremeBlinkShiftThreshold*ProbableBlinkMAD);
        end
    end
    % Obtain the median and MAD of all epoch amplitude shifts:
    MedianAmplitudeShift=median(AmplitudeShiftWithinEachEpoch(:,:,:),3);
    MADAmplitudeShift=mad(AmplitudeShiftWithinEachEpoch(:,:,:),1,3);
    % Use these values and user settings to establish the threshold:
    ShiftUpperBoundAllEpochs=MedianAmplitudeShift+(RELAX_cfg.ExtremeVoltageShiftThreshold*MADAmplitudeShift);
    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        % Use the maximum value out of the all epoch upper bound
        % threshold and the blink upper bound threshold:
        ShiftUpperBound=max(BlinkShiftUpperBound,ShiftUpperBoundAllEpochs);
    elseif RELAX_cfg.ProbabilityDataHasNoBlinks==2
        ShiftUpperBound=ShiftUpperBoundAllEpochs;
    end
    % If voltage shift within an epoch exceeds the threshold, mark for rejection:
    for epoch = 1:size(epochedEEG.data,3)
        if sum(AmplitudeShiftWithinEachEpoch(:,:,epoch)>ShiftUpperBound(:,1),1)>0
            epochedEEG.RELAXProcessing.Details.EpochsExceedingShiftThreshold(1,epoch)=1;
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1;
        end
    end
    % If a channel doesn't show a voltage shift of more than 2 microvolts,
    % mark as extreme outlier (assumed to indicate no neural activity):
    ShiftLowerBound=2*ones(size(epochedEEG.data,1),1);
    for epoch = 1:size(epochedEEG.data,3)
        if sum(AmplitudeShiftWithinEachEpoch(:,:,epoch)<ShiftLowerBound(:,1),1)>0
            epochedEEG.RELAXProcessing.Details.EpochsExceedingShiftThreshold(1,epoch)=1;% values in first column (channels) exceed upper bound
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1;% values in first column (channels) exceed upper bound
        end
    end
        
    %% Absolute threshold to identify absolute amplitude extreme values:
    MaxInEpoch=squeeze(max(epochedEEG.data,[],2));
    MinInEpoch=squeeze(min(epochedEEG.data,[],2));
    ExtremeAbsoluteVoltageThreshold = RELAX_cfg.ExtremeAbsoluteVoltageThreshold*ones(size(MaxInEpoch,1),1);
    % If maximum voltage within the epoch exceeds the maximum
    % threshold, mark as artifact:
    for epoch = 1:size(MaxInEpoch,2)
        if sum(MaxInEpoch(:,epoch)>ExtremeAbsoluteVoltageThreshold(:,1),1)>0
            epochedEEG.RELAXProcessing.Details.AbsoluteAmplitudeExceededThreshold(1,epoch)=1;
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1;
        end
    end
    % If minimum voltage within the epoch exceeds the minimum
    % threshold, mark as artifact:
    for epoch = 1:size(MinInEpoch,2)
        if sum(MinInEpoch(:,epoch)<-ExtremeAbsoluteVoltageThreshold(:,1),1)>0
            epochedEEG.RELAXProcessing.Details.AbsoluteAmplitudeExceededThreshold(1,epoch)=1;
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1;
        end
    end
    
    %% Kurtosis to identify epochs with abnormally peaky or flat distributions of voltage values and mark those epochs in the extreme outlier template:
    epochedEEG = pop_rejkurt(epochedEEG,1,(1:epochedEEG.nbchan),RELAX_cfg.ExtremeSingleChannelKurtosisThreshold,RELAX_cfg.ExtremeAllChannelKurtosisThreshold,0,0);
    % Obtain the highest single electrode kurtosis value within each epoch across all electrodes:
    epochedEEG.RELAXProcessing.Details.MarksForSingleChannelKurt=max(epochedEEG.stats.kurtE);
    % If the highest value in the epoch is below the threshold, mark as 0, if above, mark as 1:
    epochedEEG.RELAXProcessing.Details.MarksForSingleChannelKurt(epochedEEG.RELAXProcessing.Details.MarksForSingleChannelKurt <= RELAX_cfg.ExtremeSingleChannelKurtosisThreshold) = 0;
    epochedEEG.RELAXProcessing.Details.MarksForSingleChannelKurt(epochedEEG.RELAXProcessing.Details.MarksForSingleChannelKurt > RELAX_cfg.ExtremeSingleChannelKurtosisThreshold) = 1;
    % For all channel kurtosis, if the highest value in the epoch is below the threshold, mark as 0, if above, mark as 1:
    epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt=(epochedEEG.stats.kurt)';
    epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt(epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt <= RELAX_cfg.ExtremeAllChannelKurtosisThreshold) = 0;
    epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt(epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt > RELAX_cfg.ExtremeAllChannelKurtosisThreshold) = 1;
    % Add markings from single and all channel kurtosis together in single matrix:
    epochedEEG.RELAXProcessing.Details.MarksForBothSingleAndAllChannelKurt=epochedEEG.RELAXProcessing.Details.MarksForSingleChannelKurt;
    for e=1:size(epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt,2)
        if epochedEEG.RELAXProcessing.Details.MarksForAllChannelKurt(1,e)==1
            epochedEEG.RELAXProcessing.Details.MarksForBothSingleAndAllChannelKurt(1,e)=1;
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1; % Add kurtosis extreme markings to all method marking matrix
        end
    end
    
    %% Marking bad epochs for improbable distributions of voltages within the epoch:
    % Entering into the mask any epochs that show any channel with more than
    % the threshold SD deviation in absolute voltage:
    epochedEEG = pop_jointprob(epochedEEG, 1, [1:epochedEEG.nbchan], RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold, 1000, 0, 0, 0);
    % Obtain the highest single electrode improbable voltage distribution values within each epoch across all electrodes:
    epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold=max(epochedEEG.stats.jpE); 
    % If the highest value in the epoch is below the threshold, mark as 0, if above, mark as 1:
    epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold(epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold <= RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold) = 0;
    epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold(epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold > RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold) = 1;
    for epoch = 1:size(epochedEEG.data,3)
        if epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionExceededThreshold(1,epoch)==1
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1; % Add improbably voltage distribution extreme markings to all method marking matrix
        end
    end
    
    %% Test whether only drift is present for each channel in each epoch based on frequency slope:
    if isfield(epochedEEG.RELAXProcessing.Details, 'DriftRatioEpochsxChannels')==0 % only complete this if it wasn't already computed by RELAX_excluding_channels_and_epoching   
        temp=eeglab2fieldtrip(epochedEEG, 'preprocessing');
        % Set fieldtrip defaults:
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = {'all'};
        cfg.method       = 'mtmfft';
        cfg.taper        = 'hanning';
        cfg.foi          = 1:75;  % compute spectrum from 1-75Hz
        cfg.toi          = 0:0.05:size(temp.time{1,1},2);  % time window "slides" from -0.5 to 1 sec in steps of 0.05 sec (50 ms)
        cfg.pad          = 'nextpow2';
        cfg.keeptrials   = 'yes';

        warning('ignore the following warnings about trial definition, they are not relevant to this function');
        FFTPower = ft_freqanalysis(cfg, temp); % Compute power across the spectrum
        
        % The following removes influence of line noise (or line noise filtering)
        % from the muscle slope computations. It might not be necessary
        % if no online notch filter was applied, and nt_zapline fixes
        % the influence of line noise of the spectrum slope:
        if RELAX_cfg.LineNoiseOrNotchFilterAffectedData==1  
            FreqHz=FFTPower.cfg.foi;
            FreqHz(FreqHz<RELAX_cfg.LineNoiseFrequency-2)=0;
            LowerNotchFreq=min(FreqHz(FreqHz>0));
            LowerNotchFreqLocation=find(FreqHz==LowerNotchFreq);
            FreqHz=FFTPower.cfg.foi;
            FreqHz(FreqHz>RELAX_cfg.LineNoiseFrequency+2)=0;
            UpperNotchFreq=max(FreqHz(FreqHz>0));
            UpperNotchFreqLocation=find(FreqHz==UpperNotchFreq);
            FFTPower.powspctrm(:,:,LowerNotchFreqLocation:UpperNotchFreqLocation)=[];
            FFTPower.cfg.foi(LowerNotchFreqLocation:UpperNotchFreqLocation)=[];
        end

        % Calculate log-frequency log-power slope for 1-75Hz data:
        DriftRatioEpochsxChannels=zeros(size(FFTPower.powspctrm,2),size(FFTPower.powspctrm,1));
        for chan=1:size(FFTPower.powspctrm,2)
            parfor trial=1:size(FFTPower.powspctrm,1)
                powspctrm=squeeze(FFTPower.powspctrm(trial,chan,:))';
                % Fit linear regression to log-log data
                p = polyfit(log(FFTPower.cfg.foi(1,:)),log(powspctrm(1,1:size(FFTPower.cfg.foi,2))),1);
                % Store the slope
                DriftRatioEpochsxChannels(chan,trial) = p(1);         
            end
        end
        epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels=DriftRatioEpochsxChannels;
    end
    % If log-frequency log-power slope exceeds threshold (suggesting
    % only drift), then mark as artifact for each epoch:
    Threshold = RELAX_cfg.ExtremeDriftSlopeThreshold*ones(size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels,1),1);
    for epoch = 1:size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels,2)
        if sum(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels(:,epoch)<Threshold(:,1),1)>0 % if any channel shows drift exceeding the threshold
            epochedEEG.RELAXProcessing.Details.driftslopemaskepochs(1,epoch)=1; % then mark that epoch as an extreme outlier for rejection
            epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,epoch)=1; % Add drift extreme markings to all method marking matrix
        end
    end

    %% Combine the pop masks into the full noise mask:
    % Set templates for marking artifacts, initially just a EEG.data length
    % of NaNs that the MWF template will ignore, which will be added to
    % with 0's (clean periods) and 1's (artifacts) through the RELAX
    % functions:
    epochedEEG.RELAXProcessing.Details.ExtremeOutlierPeriods=epochedEEG.RELAXProcessing.Details.NaNsForNonEvents;
    OneSecondOfNaNs=NaN(1,(round(1000/RELAX_cfg.ms_per_sample)));
    OneSecondOf1s=ones(1,(round(1000/RELAX_cfg.ms_per_sample)));           
    epochedEEG.RELAX.NaNsForExtremeOutlierPeriods=epochedEEG.RELAXProcessing.Details.NaNsForNonEvents;
    % Use the epoch markings from all extreme rejection methods to create
    % template of NaN's to be excluded from MWF cleaning artifact and clean
    % period templates, and to be rejected from data before ICA:
    for e=1:size(epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections,2)
        if epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections(1,e)==1
            epochedEEG.RELAXProcessing.Details.NoiseMaskFullLength(epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOfNaNs;
            epochedEEG.RELAX.NaNsForExtremeOutlierPeriods(epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOfNaNs;
            epochedEEG.RELAXProcessing.Details.ExtremeOutlierPeriods(epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf1s;
        end
    end
    % Take record of the proportion of the data marked as an extreme
    % outlier for rejection, and the epochs that were marked:
    epochedEEG.RELAXProcessingExtremeRejections.ProportionExcludedForExtremeOutlier=mean(epochedEEG.RELAXProcessing.Details.ExtremeOutlierPeriods,'omitnan');
    epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents=epochedEEG.RELAXProcessing.Details.AllMethodsExtremeEpochRejections;
    % Obtain a list of periods of start and end timepoints for the extreme
    % outlier periods. These are kept on record for deletion before the
    % wICA step (or before any analysis)
    [~, ~, bi] = RunLength(epochedEEG.RELAX.NaNsForExtremeOutlierPeriods);
    bi(1,(size(bi,2)+1))=size(epochedEEG.RELAX.NaNsForExtremeOutlierPeriods,2);
    epochedEEG.RELAXProcessingExtremeRejections.ExtremelyBadPeriodsForDeletion = reshape(bi,[2,(size(bi,2))/2])';
    epochedEEG.RELAX.ExtremelyBadPeriodsForDeletion=epochedEEG.RELAXProcessingExtremeRejections.ExtremelyBadPeriodsForDeletion;
    % Transfer details into continuousEEG struct also:
    continuousEEG.RELAXProcessing=epochedEEG.RELAXProcessing;
    continuousEEG.RELAXProcessingExtremeRejections=epochedEEG.RELAXProcessingExtremeRejections;
    if isfield(epochedEEG, 'RELAX')==1
        continuousEEG.RELAX=epochedEEG.RELAX;
    end
end