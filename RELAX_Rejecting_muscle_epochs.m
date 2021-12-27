function [EEG] = RELAX_Rejecting_muscle_epochs(EEG, RELAX_epoching_cfg)

    %% This function rejects epochs showing muscle activity (can be implemented after cleaning to ensure analysis is not affected by muscle activity)
    % set some defaults if not specified:
    if exist('RELAX_epoching_cfg', 'var')==1
        if isfield(RELAX_epoching_cfg, 'MuscleSlopeThreshold')==0
            RELAX_epoching_cfg.MuscleSlopeThreshold=-0.31; 
    % log frequency/ log power slope above which to mark as muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or
    % More stringent = -0.72, more negative values removed more muscle, but also some brain activity.
    % -0.31 was the threshold above which a slope was definitely affected by muscle activity, 
    % so this less stringent setting ensures no epochs that are not
    % affected by muscle are rejected
        end
        if isfield(RELAX_epoching_cfg, 'MaxProportionOfDataCanBeMarkedAsMuscle')==0
            RELAX_epoching_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=0.50; % maximum that can be deleted for showing muscle slopes
        end
    elseif exist('RELAX_epoching_cfg', 'var')==0
        RELAX_epoching_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=0.50;  
        RELAX_epoching_cfg.MuscleSlopeThreshold=-0.31;
    end
    
    temp=eeglab2fieldtrip(EEG, 'preprocessing');

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
    FFTPower = ft_freqanalysis(cfg, temp);

    %% Selecting epochs that have slopes that are shallow, suggesting muscle artifact:
    EEG.RELAXProcessing.Details.muscleSlopesEpochsxChannels=zeros(size(FFTPower.powspctrm,2),size(FFTPower.powspctrm,1));
    for chan=1:size(FFTPower.powspctrm,2)
        for trial=1:size(FFTPower.powspctrm,1)
            powspctrm=squeeze(FFTPower.powspctrm(trial,chan,:))';
            % Fit linear regression to log-log data
            p = polyfit(log(FFTPower.cfg.foi(1,:)),log(powspctrm(1,1:69)),1);
            % Store the slope
            EEG.RELAXProcessing.Details.muscleSlopesEpochsxChannels(chan,trial) = p(1);         
        end
    end

    %% For epoch mask (find the max value from all channels for gamma in the epochs:
    % First work out the maximum muscle noise value for each epoch, then 
    % replace values in epochs that are below the RELAX_epoching_cfg.MuscleSlopeThreshold with 0's, so excluding
    % epochs that aren't showing muscle activity from the muscle mask, then
    % replace values above the RELAX_epoching_cfg.MuscleSlopeThreshold with 1's, indicating they
    % should be included in the MWF mask
    EEG.RELAXProcessing.Details.MaxMuscleSlopeInEpochAcrossAllChannels = max(EEG.RELAXProcessing.Details.muscleSlopesEpochsxChannels); % Returns max value in each column
    EEG.RELAXProcessing.Details.HighMuscleMask=EEG.RELAXProcessing.Details.MaxMuscleSlopeInEpochAcrossAllChannels;
    EEG.RELAXProcessing.Details.HighMuscleMask(EEG.RELAXProcessing.Details.HighMuscleMask > RELAX_epoching_cfg.MuscleSlopeThreshold) = 1;
    EEG.RELAXProcessing.Details.HighMuscleMask(EEG.RELAXProcessing.Details.HighMuscleMask <= RELAX_epoching_cfg.MuscleSlopeThreshold) = 0;

    % The following NaNs all values that aren't above the muscle threshold,
    % then sums the values that are above for each epoch to allow
    % identification of the worst epochs
    EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=EEG.RELAXProcessing.Details.muscleSlopesEpochsxChannels;
    % Replace all values below the threshold with NaN's (as they don't reflect
    % muscle activity)
    EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs(EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs < RELAX_epoching_cfg.MuscleSlopeThreshold) = NaN;
    % Shift the baseline of the values to the RELAX_epoching_cfg.MuscleSlopeThreshold, so that all
    % muscle artifacts can be ranked in severity of EMG starting from 0 (least
    % severe) and moving more positive as more severe:
    EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs-RELAX_epoching_cfg.MuscleSlopeThreshold;
    % Sum muscle activity across all channels that show it above the threshold:
    EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=nansum(EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,1);
    % rank epochs from worst muscle contamination at the top to least
    % contamination at the bottom:
    EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsBestToWorst=sort(EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,2,'descend')';
    EEG.RELAXProcessing.ProportionOfEpochsShowingMuscleActivityTotal = mean(EEG.RELAXProcessing.Details.HighMuscleMask,'omitnan');

    %% Checks the proportion of muscle artifact epochs selected, and reduces the threshold for rejection of an epoch if above the threshold were selected for rejection:
    PercentThresholdAboveWhichToIncludeAsArtifactInMask=100-(RELAX_epoching_cfg.MaxProportionOfDataCanBeMarkedAsMuscle*100);

    EEG.RELAXProcessing.Details.MuscleSlopeThresholdBasedOnPercentToReject={};
    if EEG.RELAXProcessing.ProportionOfEpochsShowingMuscleActivityTotal > RELAX_epoching_cfg.MaxProportionOfDataCanBeMarkedAsMuscle
        EEG.RELAXProcessing.Details.MuscleSlopeThresholdBasedOnPercentToReject = prctile(EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsBestToWorst,PercentThresholdAboveWhichToIncludeAsArtifactInMask);
        EEG.RELAXProcessing.Details.HighMuscleMask=EEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs;
        EEG.RELAXProcessing.Details.HighMuscleMask(EEG.RELAXProcessing.Details.HighMuscleMask <= EEG.RELAXProcessing.Details.MuscleSlopeThresholdBasedOnPercentToReject) = 0; 
        EEG.RELAXProcessing.Details.HighMuscleMask(EEG.RELAXProcessing.Details.HighMuscleMask > EEG.RELAXProcessing.Details.MuscleSlopeThresholdBasedOnPercentToReject) = 1;
    end
    
    EEG.RELAXProcessing.NumberOfEpochsToRejectForMuscle=sum(EEG.RELAXProcessing.Details.HighMuscleMask,'omitnan');
    EEG.RELAXProcessing.ProportionRejectedAsMuscle=mean(EEG.RELAXProcessing.Details.HighMuscleMask,'omitnan');
    EEG=pop_rejepoch(EEG,EEG.RELAXProcessing.Details.HighMuscleMask,0);
    
end