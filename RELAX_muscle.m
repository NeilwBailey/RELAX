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

%% RELAX_muscle:

function [continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg)

    %% This function creates a template of clean and muscle artifact affected periods for use in MWF cleaning

    % set some defaults if not specified:
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'MuscleSlopeThreshold')==0
            RELAX_cfg.MuscleSlopeThreshold=-0.59; 
            % log frequency/ log power slope above which to mark as muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or
            % More stringent = -0.72, more negative values removed more muscle, but also some brain activity.
        end
        if isfield(RELAX_cfg, 'MaxProportionOfDataCanBeMarkedAsMuscle')==0
            RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=0.50; % Maximum amount of data periods to be marked as muscle artifact for cleaning by the MWF. You want at least a reasonable amount of both clean and artifact templates for effective cleaning.
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=0.50; % Maximum amount of data periods to be marked as muscle artifact for cleaning by the MWF. You want at least a reasonable amount of both clean and artifact templates for effective cleaning.
        RELAX_cfg.MuscleSlopeThreshold=-0.59;
    end
    
    % Only compute the frequency-power details if they aren't already available from a previous function:
    if isfield(epochedEEG.RELAXProcessing.Details, 'MuscleSlopesEpochsxChannels')==0 || isempty(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels)
        if isfield(epochedEEG.RELAXProcessing.Details, 'Muscle_Slopes')==0 || isempty(epochedEEG.RELAXProcessing.Details.Muscle_Slopes)
            temp=eeglab2fieldtrip(epochedEEG, 'preprocessing');
            %% FFT method:
            % Set fieldtrip cfg:
            cfg              = [];
            cfg.output       = 'pow';
            cfg.channel      = {'all'};
            cfg.method       = 'mtmfft';
            cfg.taper        = 'hanning';
            cfg.foi          = 7:75;  % compute spectrum from 7-75Hz
            cfg.toi          = 0:0.05:size(temp.time{1,1},2);  % time window "slides" from -0.5 to 1 sec in steps of 0.05 sec (50 ms)
            cfg.pad          = 'nextpow2';
            cfg.keeptrials   = 'yes';
            warning('ignore the following warnings about trial definition, they are not relevant to this function');
            FFTPower = ft_freqanalysis(cfg, temp); % Compute power spectrum
            epochedEEG.RELAXProcessing.Details.Muscle_Slopes=FFTPower.powspctrm;
            epochedEEG.RELAXProcessing.Details.foi=FFTPower.cfg.foi;  
        end
        %% Selecting epochs that have slopes that are shallow, suggesting high muscle artifact:     
        MuscleSlopesEpochsxChannels=zeros(size(epochedEEG.RELAXProcessing.Details.Muscle_Slopes,2),size(epochedEEG.RELAXProcessing.Details.Muscle_Slopes,1));
        Muscle_Slopes=epochedEEG.RELAXProcessing.Details.Muscle_Slopes;
        for chan=1:size(Muscle_Slopes,2)
            parfor trial=1:size(Muscle_Slopes,1)
                powspctrm=squeeze(Muscle_Slopes(trial,chan,:))';
                % Fit linear regression to log-log data
                p = polyfit(log(epochedEEG.RELAXProcessing.Details.foi(1,:)),log(powspctrm(1,1:69)),1);
                % Store the slope
                MuscleSlopesEpochsxChannels(chan,trial) = p(1);         
            end
        end
        epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels=MuscleSlopesEpochsxChannels;
    end

    % Mark epoch as contaminated by muscle if muscleSlope exceeds threshold:
    epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel=sum(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels > RELAX_cfg.MuscleSlopeThreshold, 2);
    epochedEEG.RELAXProcessing.Details.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel=epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel./size(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels,2);                 

    %% To create the muscle artifact template:
    % First work out the maximum muscle slope value for each epoch, then
    % replace values in epochs that are below the
    % RELAX_cfg.MuscleSlopeThreshold with 0's, so excluding
    % epochs that aren't showing muscle activity from the muscle mask, then
    % replace values above the RELAX_cfg.MuscleSlopeThreshold with 1's,
    % indicating they should be marked as muscle in the MWF mask
    epochedEEG.RELAXProcessing.Details.HighMuscleMask = max(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels); % Returns max value in each epoch
    epochedEEG.RELAXProcessing.Details.HighMuscleMask(epochedEEG.RELAXProcessing.Details.HighMuscleMask > RELAX_cfg.MuscleSlopeThreshold) = 1; % marks as 1 (artifact) if exceeds the threshold
    epochedEEG.RELAXProcessing.Details.HighMuscleMask(epochedEEG.RELAXProcessing.Details.HighMuscleMask <= RELAX_cfg.MuscleSlopeThreshold) = 0; % marks as 0 (clean) if does not exceed the threshold
    
    % If an epoch has already been marked as "to be excluded" by the
    % extreme artifact detection script, mark as NaN here (to ignore in
    % calculating the proportion of epochs contaminated by muscle):
    if isfield(epochedEEG, 'RELAX')==1
        for e=1:size(epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents,2)
            if epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents(1,e)==1
                epochedEEG.RELAXProcessing.Details.HighMuscleMask(1,e)=NaN;
                epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels(:,e)=NaN;
            end
        end
    end

    % The following replaces all values that aren't above the muscle
    % threshold with NaN, then sums the values that are above the threshold
    % for each epoch to allow identification of the worst epochs:
    epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels;
    epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs < RELAX_cfg.MuscleSlopeThreshold) = NaN;
    % Shift the baseline of the values to the RELAX_cfg.MuscleSlopeThreshold, so that all
    % muscle artifacts can be ranked in severity of EMG starting from 0 (least
    % severe) and moving more positive as more severe:
    epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs-RELAX_cfg.MuscleSlopeThreshold;
    % Sum muscle slopes across all channels that show slopes above the
    % threshold. This gives an indication of how badly each epoch is
    % affected by muscle activity, with more affected electrodes within 
    % the epoch providing higher values: 
    epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=sum(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,1,'omitnan');

    % Take advantage of the 500ms overlap for each epoch to increase the
    % resolution of muscle slopes, by averaging the values from the
    % overlapping periods together:
    epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(1,:)=epochedEEG.RELAX.NaNsForExtremeOutlierPeriods;
    epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(2,:)=epochedEEG.RELAX.NaNsForExtremeOutlierPeriods(1,:);
    epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(2,5500:6000)=NaN;
    epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(1,size(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength,2)-6000:size(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength,2))=NaN;
    epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(2,size(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength,2)-6000:size(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength,2))=NaN;
    % Odd numbered epochs:
    for e=1:2:size(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,2)
        epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(1,epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs(e);
    end
    % Even numbered epochs:
    for e=2:2:size(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,2)
        epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength(2,epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs(e);
    end
    % Average the overlapping periods together:
    epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLengthMean=mean(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLength,'omitnan');    
    MuscleSlopeThresholdAfterAdjustment=0; % Threshold = 0, because all slope values have had the threshold subtracted from them (so the threshold is now 0). 
    epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts=epochedEEG.RELAX.NaNsForExtremeOutlierPeriods; % Insert NaNs from extreme artifact markings to reject (ignored in MWF template)
    epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLengthMean>MuscleSlopeThresholdAfterAdjustment)=1; % Insert markings where muscle slope exceeded threshold
    epochedEEG.RELAXProcessing.ProportionOfDataShowingMuscleActivityTotal = mean(epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts,'omitnan'); % work out proportion of data marked as showing muscle activity
    %% Checks the proportion of muscle artifact periods selected, and reduces the threshold for marking a period if the above section marked more than the RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle threshold
    PercentThresholdAboveWhichToIncludeAsArtifactInMask=100-(RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle*100);
    if epochedEEG.RELAXProcessing.ProportionOfDataShowingMuscleActivityTotal > RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle
        epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts=epochedEEG.RELAX.NaNsForExtremeOutlierPeriods;
        MuscleSlopeThresholdAfterAdjustment=prctile(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLengthMean,PercentThresholdAboveWhichToIncludeAsArtifactInMask);
        epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts(epochedEEG.RELAXProcessing.Details.MuscleStrengthFullLengthMean>MuscleSlopeThresholdAfterAdjustment)=1;
    end
    epochedEEG.RELAXProcessing.Details.NoiseMaskFullLength(epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts==1)=1; % Transfer muscle artifact markings into cleaning template     
    epochedEEG.RELAXProcessing.ProportionMarkedAsMuscle=mean(epochedEEG.RELAXProcessing.Details.TemplateMarkedForMuscleArtifacts,'omitnan'); % calculate the proportion of the data that has been marked as containing muscle artifact to be cleaned by MWF
    continuousEEG.RELAXProcessing=epochedEEG.RELAXProcessing; % copy all details into the continuousEEG struct
    if isfield(epochedEEG, 'RELAX')==1
        continuousEEG.RELAX=epochedEEG.RELAX;
    end  
end