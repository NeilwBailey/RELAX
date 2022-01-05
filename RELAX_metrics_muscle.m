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

%% RELAX_metrics_muscle:

function [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg)

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'MuscleSlopeThreshold')==0
            RELAX_cfg.MuscleSlopeThreshold=-0.59; 
        end
        if isfield(RELAX_cfg, 'LowPassFilter')==0
            RELAX_cfg.LowPassFilter=80;
        end
        if isfield(RELAX_cfg, 'LineNoiseOrNotchFilterAffectedData')==0
            RELAX_cfg.LineNoiseOrNotchFilterAffectedData=1; % Removes influence of line noise (or line noise filtering) from the muscle slope computations.
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.MuscleSlopeThreshold=-0.59;
        RELAX_cfg.LowPassFilter=80;
        RELAX_cfg.LineNoiseOrNotchFilterAffectedData=1;
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
                if isfield(epochedEEG.RELAXProcessing.Details, 'MuscleSlopesEpochsxChannels')==0 || isempty(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels)
                    computemuscleslope=1;
                    if isfield(epochedEEG.RELAXProcessing.Details, 'Muscle_Slopes')==0 || isempty(epochedEEG.RELAXProcessing.Details.Muscle_Slopes)
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
        warning('ignore the following warnings about trial definition, they are not relevant to this function');
        FFTPower = ft_freqanalysis(cfg, temp);
        epochedEEG.RELAXProcessing.Details.Muscle_Slopes=FFTPower.powspctrm;
        foi=FFTPower.cfg.foi;
    else 
        foi=epochedEEG.RELAXProcessing.Details.foi;
    end
    if computemuscleslope==1
        %% Selecting epochs that have slopes that are shallow, suggesting muscle artifact:     
        MuscleSlopesEpochsxChannels=zeros(size(epochedEEG.RELAXProcessing.Details.Muscle_Slopes,2),size(epochedEEG.RELAXProcessing.Details.Muscle_Slopes,1));
        Muscle_Slopes=epochedEEG.RELAXProcessing.Details.Muscle_Slopes;
        
        % The following removes influence of line noise (or line noise filtering)
        % from the muscle slope computations. It might not be necessary
        % if no online notch filter was applied, and nt_zapline fixes
        % the influence of line noise of the spectrum slope:
        if RELAX_cfg.LineNoiseOrNotchFilterAffectedData==1  
            FreqHz=foi;
            FreqHz(FreqHz<RELAX_cfg.LineNoiseFrequency-2)=0;
            LowerNotchFreq=min(FreqHz(FreqHz>0));
            LowerNotchFreqLocation=find(FreqHz==LowerNotchFreq);
            FreqHz=FFTPower.cfg.foi;
            FreqHz(FreqHz>RELAX_cfg.LineNoiseFrequency+2)=0;
            UpperNotchFreq=max(FreqHz(FreqHz>0));
            UpperNotchFreqLocation=find(FreqHz==UpperNotchFreq);
            Muscle_Slopes(:,:,LowerNotchFreqLocation:UpperNotchFreqLocation)=[];
            foi(LowerNotchFreqLocation:UpperNotchFreqLocation)=[];
        end
        
        for chan=1:size(Muscle_Slopes,2)
            parfor trial=1:size(Muscle_Slopes,1)
                powspctrm=squeeze(Muscle_Slopes(trial,chan,:))';
                % Fit linear regression to log-log data
                p = polyfit(log(foi(1,:)),log(powspctrm(1,1:size(foi,2))),1);
                % Store the slope
                MuscleSlopesEpochsxChannels(chan,trial) = p(1);         
            end
        end
        epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels=MuscleSlopesEpochsxChannels;
    end  

    if RELAX_cfg.LowPassFilter>=75
        % The following NaNs all values that aren't above the muscle threshold,
        % then sums the values that are above for each epoch to allow
        % identification of the worst epochs
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels;
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs < RELAX_cfg.MuscleSlopeThreshold) = NaN;
        epochedEEG.RELAXProcessing.Details.MuscleSlopeSeverity=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs;
        % Shift the baseline of the values to the RELAX_cfg.MuscleSlopeThreshold, so that all
        % muscle artifacts can be ranked in severity of EMG starting from 0 (least
        % severe) and moving more positive as more severe. This way, 0
        % reflects the muscle slope threshold decided upon, and values
        % closer to 0 mean less muscle activity above the threshold,
        % whereas higher values reflect more muscle activity.
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs-RELAX_cfg.MuscleSlopeThreshold;
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero=epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs;
        epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero(isnan(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochsSubthresholdEqualsZero))=0;

        epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel=sum(epochedEEG.RELAXProcessing.Details.SortingOutWorstMuscleEpochs,1,'omitnan');
        epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel(epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel>1)=1;
        epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel=sum(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels > RELAX_cfg.MuscleSlopeThreshold, 2);

        if continuousEEG.RELAX.Data_has_been_cleaned==0
            epochedEEG.RELAX_Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel=epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel./size(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels,2); 
            epochedEEG.RELAX_Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel=mean(epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel,'all','omitnan'); 
            epochedEEG.RELAX_Metrics.Raw.MeanMuscleStrengthFromOnlySuperThresholdValues=mean(epochedEEG.RELAXProcessing.Details.MuscleSlopeSeverity,'all','omitnan');
        end
        if continuousEEG.RELAX.Data_has_been_cleaned==1
            epochedEEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdPerChannel=epochedEEG.RELAXProcessing.Details.NumberOfBadMuscleEpochsFromEachChannel./size(epochedEEG.RELAXProcessing.Details.MuscleSlopesEpochsxChannels,2); 
            epochedEEG.RELAX_Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel=mean(epochedEEG.RELAXProcessing.Details.EpochsContainingMuscleAnyChannel,'all','omitnan');
            epochedEEG.RELAX_Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues=mean(epochedEEG.RELAXProcessing.Details.MuscleSlopeSeverity,'all','omitnan');
        end
        continuousEEG.RELAX_Metrics=epochedEEG.RELAX_Metrics;
    end
end