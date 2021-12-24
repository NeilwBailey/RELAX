function [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg)

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'MuscleSlopeThreshold')==0
            RELAX_cfg.MuscleSlopeThreshold=-0.59; 
        end
        if isfield(RELAX_cfg, 'LowPassFilter')==0
            RELAX_cfg.LowPassFilter=80;
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.MuscleSlopeThreshold=-0.59;
        RELAX_cfg.LowPassFilter=80;
    end
    
    computefourier=0; computemuscleslope=0;
    
    if continuousEEG.RELAX.Data_has_been_cleaned==1
        epochedEEG=eeg_regepochs(continuousEEG,'recurrence',0.500,'eventtype','X','extractepochs','off');
        %Epoch data
        epochedEEG = pop_epoch( epochedEEG, {'X'}, [0 1.0], 'epochinfo', 'yes');
        epochedEEG = eeg_checkset( epochedEEG );
        epochedEEG = pop_selectevent( epochedEEG, 'type', 'X', 'omitlatency', '1<=1999', 'deleteevents','on');
    end
    
    if RELAX_cfg.LowPassFilter<75
        msg='Data has been low pass filtered below 75Hz so muscle artifact slopes are not computable';
        warning(msg)
    end    
    
    if (continuousEEG.RELAX.Data_has_been_cleaned==1 && RELAX_cfg.LowPassFilter>=75)
        computefourier=1; computemuscleslope=1; % compute fourier power and muscle slope if data has been cleaned (the values may be present for the raw data, but the raw data muscle details will be inappropriate for calculating the clean data metrics)
    end
    
    % If required details aren't already present, set to compute them:
    if isfield(epochedEEG,'RELAXProcessing')~=1 && RELAX_cfg.LowPassFilter>=75
        computefourier=1; computemuscleslope=1;
    end
    if isfield(epochedEEG,'RELAXProcessing') && RELAX_cfg.LowPassFilter>=75
        if isfield(epochedEEG.RELAXProcessing, 'Details')~=1
            computefourier=1; computemuscleslope=1;
        end
    end
    if (RELAX_cfg.LowPassFilter>=75 && continuousEEG.RELAX.Data_has_been_cleaned~=1) 
        if isfield(epochedEEG,'RELAXProcessing')
            if isfield(epochedEEG.RELAXProcessing, 'Details')
                if isfield(epochedEEG.RELAXProcessing.Details, 'muscleRatioEpochsxChannels')==0 || isempty(epochedEEG.RELAXProcessing.Details.muscleRatioEpochsxChannels)
                    computemuscleslope=1;
                    if isfield(epochedEEG.RELAXProcessing.Details, 'musclePower')==0 || isempty(epochedEEG.RELAXProcessing.Details.musclePower)
                        computefourier=1;
                    end
                end
            end
        end
    end
    
    if computefourier==1
        temp=eeglab2fieldtrip(epochedEEG, 'preprocessing');
        %% FFT method:
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = {'all'};
        cfg.method       = 'mtmfft';
        cfg.taper        = 'hanning';
        cfg.foi          = 7:75;  
        cfg.toi          = 0:0.05:1; % time window "slides" from -0.5 to 1 sec in steps of 0.05 sec (50 ms)
        cfg.pad          = 'nextpow2';
        cfg.keeptrials   = 'yes';
        FFTPower = ft_freqanalysis(cfg, temp);
        epochedEEG.RELAXProcessing.Details.musclePower=FFTPower.powspctrm;
    end
    if computemuscleslope==1
        %% Selecting epochs that have slopes that are shallow, suggesting high gamma and muscle artifact:     
        muscleRatioEpochsxChannels=zeros(size(epochedEEG.RELAXProcessing.Details.musclePower,2),size(epochedEEG.RELAXProcessing.Details.musclePower,1));
        musclePower=epochedEEG.RELAXProcessing.Details.musclePower;
        for chan=1:size(musclePower,2)
            parfor trial=1:size(musclePower,1)
                powspctrm=squeeze(musclePower(trial,chan,:))';
                % Fit linear regression to log-log data
                p = polyfit(log(FFTPower.cfg.foi(1,:)),log(powspctrm(1,1:69)),1);
                % Store the slope
                muscleRatioEpochsxChannels(chan,trial) = p(1);         
            end
        end
        epochedEEG.RELAXProcessing.Details.muscleRatioEpochsxChannels=muscleRatioEpochsxChannels;
    end  

    if RELAX_cfg.LowPassFilter>=75
        % The following NaNs all values that aren't above the muscle threshold,
        % then sums the values that are above for each epoch to allow
        % identification of the worst epochs
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=epochedEEG.RELAXProcessing.Details.muscleRatioEpochsxChannels;
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs < RELAX_cfg.MuscleSlopeThreshold) = NaN;
        % Shift the baseline of the values to the RELAX_cfg.MuscleSlopeThreshold, so that all
        % muscle artifacts can be ranked in severity of EMG starting from 0 (least
        % severe) and moving more positive as more severe. This way, 0
        % reflects the muscle slope threshold decided upon, and values
        % closer to 0 mean less muscle activity above the threshold,
        % whereas higher values reflect more muscle activity.
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs-RELAX_cfg.MuscleSlopeThreshold;
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs;
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero(isnan(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero))=0;

        epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel=nansum(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,1);
        epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel(epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel>1)=1;
        epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel=sum(epochedEEG.RELAXProcessing.Details.muscleRatioEpochsxChannels > RELAX_cfg.MuscleSlopeThreshold, 2);

        if continuousEEG.RELAX.Data_has_been_cleaned==0
            epochedEEG.RELAX.Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel=epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel./size(epochedEEG.RELAXProcessing.Details.muscleRatioEpochsxChannels,2); 
            epochedEEG.RELAX.Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel=mean(epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel,'all','omitnan'); 
            epochedEEG.RELAX.Metrics.Raw.MeanMuscleStrengthFromOnlySuperThresholdValues=mean(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,'all','omitnan');
            epochedEEG.RELAX.Metrics.Raw.MeanMuscleStrengthScaledByProportionShowingMuscle=mean(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero,'all','omitnan');
        end
        if continuousEEG.RELAX.Data_has_been_cleaned==1
            epochedEEG.RELAX.Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel=epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel./size(epochedEEG.RELAXProcessing.Details.muscleRatioEpochsxChannels,2); 
            epochedEEG.RELAX.Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel=mean(epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel,'all','omitnan');
            epochedEEG.RELAX.Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues=mean(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,'all','omitnan');
            epochedEEG.RELAX.Metrics.Cleaned.MeanMuscleStrengthScaledByProportionShowingMuscle=mean(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero,'all','omitnan');
        end
        continuousEEG.RELAX.Metrics=epochedEEG.RELAX.Metrics;
    end
end