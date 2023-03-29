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

%% RELAX_excluding_channels_and_epoching:

function [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(continuousEEG, RELAX_cfg)

    %% This function epochs the data for artifact detection, then excludes bad channels that exceed the extreme outlier detection thresholds 

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'ExtremeVoltageShiftThreshold')==0
            RELAX_cfg.ExtremeVoltageShiftThreshold=25; % Threshold MAD from the median all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe voltage shifts within the epoch
        end
        if isfield(RELAX_cfg, 'ExtremeImprobableVoltageDistributionThreshold')==0
            RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; % Threshold SD from the mean of all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe improbable data
        end
        if isfield(RELAX_cfg, 'ms_per_sample')==0
            RELAX_cfg.ms_per_sample=(1000/continuousEEG.srate); 
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
        if isfield(RELAX_cfg, 'ProportionOfExtremeNoiseAboveWhichToRejectChannel')==0
            RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=0.05; % If a channel shows extreme outlying data for more than this proportion of the total data, it is deleted
        end
        if isfield(RELAX_cfg, 'ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel')==0
            RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=0.05; % If the proportion of epochs showing muscle activity from an electrode is higher than this, the electrode is deleted. 
            % Set muscle proportion before deletion to 1 to not delete electrodes based on muscle activity
        end
        if isfield(RELAX_cfg, 'MuscleSlopeThreshold')==0
            RELAX_cfg.MuscleSlopeThreshold=-0.59; % log power log frequency slope threshold above which an epoch is marked as containing muscle
        end
        if isfield(RELAX_cfg, 'MaxProportionOfElectrodesThatCanBeDeleted')==0
            RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted=0.20; % Sets the maximum proportion of electrodes that are allowed to be deleted after PREP's bad electrode deletion step
        end
    elseif exist('RELAX_cfg', 'var')==0     
        RELAX_cfg.ExtremeVoltageShiftThreshold=25; %MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
        RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; %SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
        RELAX_cfg.ms_per_sample=(1000/continuousEEG.srate);
        RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=8; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis 
        RELAX_cfg.ExtremeAllChannelKurtosisThreshold=8; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis  
        RELAX_cfg.ExtremeAbsoluteVoltageThreshold=500; % microvolts max or min above which will be excluded from cleaning and deleted from data
        RELAX_cfg.ExtremeBlinkShiftThreshold=3; % How many MAD from the median of blink affected epochs to exclude as extreme data
        RELAX_cfg.ExtremeDriftSlopeThreshold=-4; % slope of log frequency log power below which to reject as drift without neural activity
        RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=0.05; % Proportion of data where channel shows extreme noise before channel is deleted
        RELAX_cfg.MuscleSlopeThreshold=-0.59; % log power log frequency slope threshold above which an epoch is marked as containing muscle
        RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=0.05; % If the proportion of epochs showing muscle activity from an electrode is higher than this, the electrode is deleted. 
        % Set muscle proportion before deletion to 1 to not delete electrodes based on muscle activity
        RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted=0.20; % Sets the maximum proportion of electrodes that are allowed to be deleted after PREP's bad electrode deletion step
    end
    
    continuousEEG.RELAX.ListOfChannelsAfterRejections={continuousEEG.chanlocs.labels}; % Get list of channels currently present

    % Set templates for marking artifacts, initially just a EEG.data length
    % of NaNs that the MWF template will ignore, which will be added to
    % with 0's (clean periods) and 1's (artifacts) through the RELAX
    % functions:
    continuousEEG.RELAXProcessing.Details.NaNsForNonEventsAllChannels=NaN(size(continuousEEG.data));
    continuousEEG.RELAXProcessing.Details.NaNsForNonEvents=NaN(1,size(continuousEEG.data,2));
    OneSecondOf1s=ones(1,(round(1000/RELAX_cfg.ms_per_sample)));
    OneSecondOf0sAllChannels=zeros((size(continuousEEG.data,1)),round(1000/RELAX_cfg.ms_per_sample)); 
    OneSecondOf0s=zeros(1,(round(1000/RELAX_cfg.ms_per_sample)));

    % Mark triggers into data each 1s with 0.5s overlaps:
    epochedEEG=eeg_regepochs(continuousEEG,'recurrence',0.500,'eventtype','X','extractepochs','off');    
    for e=1:size(epochedEEG.event,2)
        epochedEEG.event(e).originallatency=epochedEEG.event(e).latency;
    end
    epochedEEG = pop_selectevent( epochedEEG, 'type', 'X', 'deleteevents','on');
    % Change X's to Y's if they're the first or last 11 events, to be ignored by RELAX functions, because they're often noisy and will contaminate potential MWF masks:
    for e=1:11
        if (strcmp(epochedEEG.event(e).type, 'X'))
            epochedEEG.event(e).type = 'Y';
        end
    end
    for e=size(epochedEEG.event,2)-11:size(epochedEEG.event,2)
        if (strcmp(epochedEEG.event(e).type, 'X'))
            epochedEEG.event(e).type = 'Y';
        end
    end
    %Epoch data 1s with 0.5s overlaps:
    epochedEEG = pop_epoch( epochedEEG, {'X'}, [0 1.0], 'epochinfo', 'yes');
    epochedEEG = eeg_checkset( epochedEEG );
    epochedEEG = pop_selectevent( epochedEEG, 'type', 'X', 'omitlatency', '1<=1999', 'deleteevents','on'); 
    % Insert 0's into the masking template for all periods except the first
    % and last 5s of the data:
    for e=1:size(epochedEEG.event,2)
        if (strcmp(epochedEEG.event(e).type, 'X'))
            epochedEEG.RELAXProcessing.Details.NaNsForNonEventsAllChannels(:,epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf0sAllChannels;
            epochedEEG.RELAXProcessing.Details.NaNsForNonEvents(:,epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf0s;
        end
    end

    TotalInitialChannels=size(continuousEEG.allchan,2);
    CurrentChannels=size(continuousEEG.chanlocs,2);

    % Calculating how many channels can be deleted in this function based
    % on user setting:
    YouCanRejectThisManyChannelsHere=floor(RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*TotalInitialChannels)-(TotalInitialChannels-CurrentChannels);
    
    epochedEEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject=0;
    epochedEEG.RELAXProcessingExtremeRejections.ExtremeDataBasedChannelToReject={}; % create empty cells for insertion into RELAX outputs table in case the following section is not performed

    %% Detect voltage shift in Epoch to identify outlying channels:
    if YouCanRejectThisManyChannelsHere>0 
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
        % Set blank templates to insert extreme data detections into:
        dimensions=size(epochedEEG.data);
        epochedEEG.RELAXProcessing.Details.ShiftInEpochToRejectChannels=zeros(dimensions(1),dimensions(3)); % extreme voltage shift within epoch artifact rejection type recorded in this template
        epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels=zeros(dimensions(1),dimensions(3)); % all artifact rejection types recorded in this template
        % Mark epochs showing extreme data exceeding thresholds into the
        % template for each electrode:
        for epoch = 1:size(epochedEEG.data,3)
            epochedEEG.RELAXProcessing.Details.ShiftInEpochToRejectChannels(AmplitudeShiftWithinEachEpoch(:,:,epoch)>ShiftUpperBound(:,1),epoch)=1;   
            epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(AmplitudeShiftWithinEachEpoch(:,:,epoch)>ShiftUpperBound(:,1),epoch)=1;   
        end
        % If a channel doesn't show a voltage shift of more than 2 microvolts,
        % mark as extreme outlier (assumed to indicate no neural activity):
        ShiftLowerBound=2*ones(size(epochedEEG.data,1),1);
        for epoch = 1:size(epochedEEG.data,3)
            epochedEEG.RELAXProcessing.Details.FlatRejectChannels(AmplitudeShiftWithinEachEpoch(:,:,epoch)<ShiftLowerBound(:,1),epoch)=1; 
            epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(AmplitudeShiftWithinEachEpoch(:,:,epoch)<ShiftLowerBound(:,1),epoch)=1; % v1.1.3 update - flatline epochs weren't being added to the cumulative bad epoch/electrode detections. They are now added.
        end

        %% Absolute threshold to identify absolute amplitude extreme values:
        MaxInEpoch=squeeze(max(epochedEEG.data,[],2));
        MinInEpoch=squeeze(min(epochedEEG.data,[],2));
        ExtremeAbsoluteVoltageThreshold = RELAX_cfg.ExtremeAbsoluteVoltageThreshold*ones(size(MaxInEpoch,1),1);
        % If maximum voltage within the epoch exceeds the maximum
        % threshold, mark as artifact:
        for epoch = 1:size(MaxInEpoch,2)
            epochedEEG.RELAXProcessing.Details.ExtremeAbsoluteVoltageRejectChannels(MaxInEpoch(:,epoch)>ExtremeAbsoluteVoltageThreshold(:,1),epoch)=1;
            epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(MaxInEpoch(:,epoch)>ExtremeAbsoluteVoltageThreshold(:,1),epoch)=1;
        end
        % If minimum voltage within the epoch exceeds the minimum
        % threshold, mark as artifact:
        for epoch = 1:size(MaxInEpoch,2)
            epochedEEG.RELAXProcessing.Details.ExtremeAbsoluteVoltageRejectChannels(MinInEpoch(:,epoch)<-ExtremeAbsoluteVoltageThreshold(:,1),epoch)=1;
            epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(MinInEpoch(:,epoch)<-ExtremeAbsoluteVoltageThreshold(:,1),epoch)=1;
        end

        %% Kurtosis to identify epochs with abnormally peaky or flat distributions of voltage values and mark those epochs in the extreme outlier template:
        epochedEEG = pop_rejkurt(epochedEEG,1,(1:epochedEEG.nbchan),RELAX_cfg.ExtremeSingleChannelKurtosisThreshold,RELAX_cfg.ExtremeAllChannelKurtosisThreshold,0,0);
        epochedEEG.RELAXProcessing.Details.KurtosisRejectChannels=epochedEEG.reject.rejkurtE; % Extract just single electrode extreme artifact markings
        epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(epochedEEG.RELAXProcessing.Details.KurtosisRejectChannels==1)=1;

        %% Marking bad epochs for improbable distributions of voltages within the epoch:
        % Entering into the mask any epochs that show any channel with more than
        % the threshold SD deviation in absolute voltage:
        epochedEEG = pop_jointprob(epochedEEG, 1, [1:epochedEEG.nbchan], RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold, 1000, 0, 0, 0);
        epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionRejectChannels=epochedEEG.reject.rejjpE; % Extract just single electrode extreme artifact markings
        epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(epochedEEG.RELAXProcessing.Details.ImprobableVoltageDistributionRejectChannels==1)=1;

        %% Test whether only drift is present for each channel in each epoch based on frequency slope:
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
        
        epochedEEG.RELAXProcessing.Details.Muscle_Slopes=FFTPower.powspctrm(:,:,7:75); % Compute power across the spectrum within the frequencies used to examine potential muscle slopes (stored for later)
        foi=FFTPower.cfg.foi(1,7:75); % store the muscle foi for later
        epochedEEG.RELAXProcessing.Details.foi=FFTPower.cfg.foi(1,7:75); % store the muscle foi for later

        % Calculate log-frequency log-power slope for 1-75Hz data:
        DriftRatioEpochsxChannels=zeros(size(FFTPower.powspctrm,2),size(FFTPower.powspctrm,1));
        for chan=1:size(FFTPower.powspctrm,2)
            parfor trial=1:size(FFTPower.powspctrm,1)
                powspctrm=squeeze(FFTPower.powspctrm(trial,chan,:))';
                % Fit linear regression to log-log data
                p = polyfit(log(FFTPower.cfg.foi(1,:)),log(powspctrm(1,1:75)),1);
                % Store the slope
                DriftRatioEpochsxChannels(chan,trial) = p(1);         
            end
        end
        epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels=DriftRatioEpochsxChannels;
        
        % If log-frequency log-power slope exceeds threshold (suggesting
        % only drift), then mark as artifact for each electrode/epoch
        Threshold = RELAX_cfg.ExtremeDriftSlopeThreshold*ones(size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels,1),1);
        for epoch = 1:size(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels,2)
            epochedEEG.RELAXProcessing.Details.driftslopeRejectChannels(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels(:,epoch)<Threshold(:,1),epoch)=1;
            epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels(:,epoch)<Threshold(:,1),epoch)=1;
        end

        % If an epoch shows more extremely bad electrodes than the number of
        % electrodes that the user has set to be deleted, the following
        % lines clear the artifact markings from that epoch, so that the
        % extreme epoch detection function will delete the epoch instead
        % (this ensures epochs where every channel has gone bad do not
        % contribute to the proportion of bad epochs for each single
        % channel, instead opting to just delete that channel later, in
        % order to preserve more electrodes).
        epochedEEG.RELAXProcessing.Details.NumberOfElectrodesBadAtEachTimepoint=sum(epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels,1,'omitnan');
        for e=1:size(epochedEEG.RELAXProcessing.Details.NumberOfElectrodesBadAtEachTimepoint,2)
            if epochedEEG.RELAXProcessing.Details.NumberOfElectrodesBadAtEachTimepoint(1,e)>YouCanRejectThisManyChannelsHere
                epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(1:size(epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels,1),e)=0;
            end
        end

        % Shift extreme artifact detection methods from the epoch space
        % into the full EEG.data length template:
        epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannelsFullLength=zeros(size(continuousEEG.data));
        for c=1:size(epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels,1)
            for e=1:size(epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels,2)
                if epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(c,e)==1
                    epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannelsFullLength(c,epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf1s;
                end
            end
        end
        % Insert NaNs into the full length template
        epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannelsFullLength(isnan(epochedEEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;     
        
        % Calculate the proportion of the full template that has been
        % marked as an extreme outlier for each channel:
        epochedEEG.RELAXProcessing.Details.ProportionExtremeForEachChannel=mean(epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannelsFullLength,2,'omitnan');
        
        ExtremeOutlierChannelRejectionThreshold=RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel;
        epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject={};
        % Obtain a count of the number of channels that show extreme noise above the
        % thresholds in more than the proportion of data you've opted as
        % the threshold before an electrode is deleted. If more electrodes
        % are recommended for deletion than the threshold, then rank the
        % electrodes by the proportion of data showing extreme artifacts,
        % and just reject the worst electrodes up to the user specified
        % threshold:
        epochedEEG.RELAXProcessing.Details.NumberOfExtremeNoiseChannelsRecomendedToDelete=sum(epochedEEG.RELAXProcessing.Details.ProportionExtremeForEachChannel > RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel, 1);
        if epochedEEG.RELAXProcessing.Details.NumberOfExtremeNoiseChannelsRecomendedToDelete > YouCanRejectThisManyChannelsHere
            ExtremeOutlierChannelRejectionThreshold = prctile(epochedEEG.RELAXProcessing.Details.ProportionExtremeForEachChannel,(1-YouCanRejectThisManyChannelsHere/CurrentChannels)*100); % set the proportion of extreme artifact epochs threshold at the % of electrodes available for deletion
        end
        % Record the electrodes to reject:
        if epochedEEG.RELAXProcessing.Details.NumberOfExtremeNoiseChannelsRecomendedToDelete > 0
            for x=1:size(epochedEEG.chanlocs,2)
                if epochedEEG.RELAXProcessing.Details.ProportionExtremeForEachChannel(x,1)>ExtremeOutlierChannelRejectionThreshold
                    epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject{x,1}=epochedEEG.chanlocs(x).labels;
                    % insert NaN where that channel was in the
                    % CumulativeMethodsRejectChannels variable:
                    epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels(x,:)=NaN;
                end
            end
            % Reject that electrode's data in the drift and muscle variables:
            epochedEEG.RELAXProcessing.Details.Muscle_Slopes(:,epochedEEG.RELAXProcessing.Details.ProportionExtremeForEachChannel>ExtremeOutlierChannelRejectionThreshold,:)=[];
            epochedEEG.RELAXProcessing.Details.DriftRatioEpochsxChannels(epochedEEG.RELAXProcessing.Details.ProportionExtremeForEachChannel>ExtremeOutlierChannelRejectionThreshold,:)=[];
        end
        % Delete the channel from the data:
        epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject = epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject(~any(cellfun('isempty', epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject), 2), :);
        epochedEEG=pop_select(epochedEEG,'nochannel',epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject);  
        epochedEEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete=epochedEEG.RELAXProcessing.Details.NumberOfExtremeNoiseChannelsRecomendedToDelete;
        epochedEEG.RELAXProcessingExtremeRejections.ExtremeDataBasedChannelToReject=epochedEEG.RELAXProcessing.Details.ExtremeDataBasedChannelToReject;
        epochedEEG = eeg_checkset( epochedEEG );

        % Record which epochs still show extreme artifacts after the electrode
        % rejections. These epochs will be ignored in the muscle detection
        % step that follows, so potential rejections based on muscle
        % activity do not reflect extreme artifacts:
        epochedEEG.RELAX.ExtremeEpochsToIgnoreInMuscleDetectionStep=sum(epochedEEG.RELAXProcessing.Details.CumulativeMethodsRejectChannels,1,'omitnan');

        %% Muscle slope detection for electrode deletion:
        epochedEEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject={};
        epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject={};
        TotalInitialChannels=size(epochedEEG.allchan,2);
        CurrentChannels=size(epochedEEG.chanlocs,2);

        % Calculating how many electrodes can be deleted based on user
        % settings:
        YouCanRejectThisManyChannelsHere=floor(RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted*TotalInitialChannels)-(TotalInitialChannels-CurrentChannels);
        epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete=0; % v1.1.3 update, so later lines can still use this variable even if the following loop isn't engaged. Thankyou Mana Biabani for the update
        if YouCanRejectThisManyChannelsHere>0
            % Selecting epochs that have slopes that are shallow, suggesting high gamma and muscle artifact:
            MuscleSlopesEpochsxChannels=zeros(size(epochedEEG.RELAXProcessing.Details.Muscle_Slopes,2),size(epochedEEG.RELAXProcessing.Details.Muscle_Slopes,1));
            Muscle_Slopes=epochedEEG.RELAXProcessing.Details.Muscle_Slopes;
            % Compute slope of log-frequency log-power from within the
            % slopes relevant for muscle detection:
            for chan=1:size(Muscle_Slopes,2)
                parfor trial=1:size(Muscle_Slopes,1)
                    powspctrm=squeeze(Muscle_Slopes(trial,chan,:))';
                    % Fit linear regression to log-log data
                    p = polyfit(log(foi(1,:)),log(powspctrm(1,1:69)),1);
                    % Store the slope
                    MuscleSlopesEpochsxChannels(chan,trial) = p(1);         
                end
            end
            epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels=MuscleSlopesEpochsxChannels;

            % Ignore epochs that have already been marked as extreme
            % outliers for deletion:
            for e=1:size(epochedEEG.RELAX.ExtremeEpochsToIgnoreInMuscleDetectionStep,2)
                if epochedEEG.RELAX.ExtremeEpochsToIgnoreInMuscleDetectionStep(1,e)>0
                    epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels(:,e)=NaN;
                end
            end
            
            % Calculate how many epochs show muscle activity from each electrode:
            epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel=sum(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels > RELAX_cfg.MuscleSlopeThreshold, 2,'omitnan');
            % Calculate the proportion of epochs showing muscle activity from each electrode:
            epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel=epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel./size(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels,2);                 

            % Check whether the above steps have recommended more
            % electrodes for rejection than the user's limit, in which case
            % adjust the threshold to allow more muscle activity before
            % electrode rejection:
            ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel;
            epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete=sum(epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel > RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel, 1,'omitnan');
            if epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete > YouCanRejectThisManyChannelsHere
                ProportionOfEpochsShowingMuscleAboveThresholdPerChannelsorted=sort(epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel,1,'descend');
                ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=ProportionOfEpochsShowingMuscleAboveThresholdPerChannelsorted(YouCanRejectThisManyChannelsHere,1); % set the proportion of extreme artifact epochs threshold at the % of electrodes available for deletion
            end
            
            % If there are still electrodes available for deletion, then
            % delete them in the following steps:
            if epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete > 0
                Backup.MuscleSlopesEpochsxChannels=epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels; % Take backup of these variables in case they're needed below
                Backup.Muscle_Slopes=epochedEEG.RELAXProcessing.Details.Muscle_Slopes;
                for x=1:size(epochedEEG.chanlocs,2)
                    if epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel(x,1)>=ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel % If the electrode shows equal to or more muscle contaminated epochs than the threshold...
                        epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject{x,1}=epochedEEG.chanlocs(x).labels; % Mark the electrode for deletion
                    end
                end
                epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels(epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel>=ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel,:)=[]; % If electrode is deleted, remove it's data from this variable
                epochedEEG.RELAXProcessing.Details.Muscle_Slopes(:,epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel>=ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel,:)=[]; % If electrode is deleted, remove it's data from this variable
                epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject = epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject(~any(cellfun('isempty', epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject), 2), :);
                % Using >= the threshold in the above section can sometimes lead to more
                % electrodes being rejected than the threshold. The
                % following lines test if this has happened, and if so,
                % uses > the threshold instead:
                if size(epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject,1)>YouCanRejectThisManyChannelsHere
                    epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject={};
                    epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels=Backup.MuscleSlopesEpochsxChannels;
                    epochedEEG.RELAXProcessing.Details.Muscle_Slopes=Backup.Muscle_Slopes;
                    for x=1:size(epochedEEG.chanlocs,2)
                        if epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel(x,1)>ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel % If the electrode shows more muscle contaminated epochs than the threshold...
                            epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject{x,1}=epochedEEG.chanlocs(x).labels; % Mark the electrode for deletion
                        end
                    end
                    epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels(epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel>ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel,:)=[]; % If electrode is deleted, remove it's data from this variable
                    epochedEEG.RELAXProcessing.Details.Muscle_Slopes(:,epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel>ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel,:)=[]; % If electrode is deleted, remove it's data from this variable
                end
                epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject = epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject(~any(cellfun('isempty', epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject), 2), :);
                epochedEEG=pop_select(epochedEEG,'nochannel',epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject);  % Delete electrodes from epoched data that have been marked as showing more than the threshold number of epochs contaminated by muscle activity
                epochedEEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject=epochedEEG.RELAXProcessing.Details.MuscleBasedElectrodesToReject;
            end
        elseif YouCanRejectThisManyChannelsHere<1
            epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete=0;
            epochedEEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject={};
        end
        epochedEEG.RELAX=rmfield(epochedEEG.RELAX,'ExtremeEpochsToIgnoreInMuscleDetectionStep');
    elseif YouCanRejectThisManyChannelsHere<1
        epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete=0;
        epochedEEG.RELAXProcessingExtremeRejections.NumberOfExtremeNoiseChannelsRecomendedToDelete=0;
        epochedEEG.RELAXProcessingExtremeRejections.MuscleBasedElectrodesToReject={};      
        epochedEEG.RELAXProcessingExtremeRejections.ExtremeDataBasedChannelToReject={};
    end
    epochedEEG.RELAXProcessingExtremeRejections.NumberOfMuscleContaminatedChannelsRecomendedToDelete=epochedEEG.RELAXProcessing.Details.NumberOfMuscleContaminatedChannelsRecomendedToDelete;
    epochedEEG.RELAX.ListOfChannelsAfterRejections={epochedEEG.chanlocs.labels}; % Get list of good channels
    epochedEEG.RELAX.NumberOfElectrodesAfterRejections=epochedEEG.nbchan;
    continuousEEG=pop_select(continuousEEG,'channel',epochedEEG.RELAX.ListOfChannelsAfterRejections); % Delete electrodes from continuous data that have been marked as showing more than the threshold number of epochs contaminated by muscle activity
    continuousEEG.RELAX=epochedEEG.RELAX;
    continuousEEG.RELAXProcessing=epochedEEG.RELAXProcessing;
    continuousEEG.RELAXProcessingExtremeRejections=epochedEEG.RELAXProcessingExtremeRejections;
end    
