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

%% RELAX_drift:

function [continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg)

    % This function detects whether a single channel has drifted away from
    % other channels and creates a template for the MWF to clean the
    % artifact.

    % set some defaults for included channels and trials, if not already
    % specified:
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'ProportionWorstEpochsForDrift')==0
            RELAX_cfg.ProportionWorstEpochsForDrift=0.30; % Maximum proportion of data that can be marked as showing drift in the template
        end
        if isfield(RELAX_cfg, 'DriftSeverityThreshold')==0
            RELAX_cfg.DriftSeverityThreshold=10; % MAD from the median of all electrodes. This could be set lower and would catch less severe drift
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.ProportionWorstEpochsForDrift=0.30; % Maximum proportion of data that can be marked as showing drift in the template    
        RELAX_cfg.DriftSeverityThreshold=10; % MAD from the median of all electrodes. This could be set lower and would catch less severe drift 
    end
        
    OneSecondOf1s=ones(1,round(1000/RELAX_cfg.ms_per_sample)); % set a variable filled with 1 second worth of 1's to insert into the artifact template when drift is detected

    % The following seems to detect drift in electrodes within epochs fairly
    % well. It measures whether a single channel shows an amplitude that
    % has drifted away from the other channels, so that it's value is more
    % than the threshold away from the all channel median for >=100 of the
    % timepoints in the epoch.
    
    % It is a modified version of the process used in the FASTER pipeline 
    % (which uses +/- 3SD from the mean, which I find to not be sensitive 
    % as the SD is heavily influenced by outliers, so if a few electrodes 
    % are drifting it won't catch them). So instead, this function uses 
    % deviation from the median by the median absolute deviation (10MAD
    % was a good limit in our testing)

    % this may not be helpful if you're just looking at oscillations above 2Hz,
    % as that data won't be adversely affected by drift below 1Hz

    % The multiplication factor of the median absolute deviation
    % and the number of electrode timepoints above that factor could be
    % experimented with to fine tune.

    epochedEEG.RELAXProcessing.Details.driftmaskepochs=zeros(size(epochedEEG.data,3),1)';
    
    % Filter out the frequencies >5Hz to avoid marking high amplitude
    % alpha as drift (only implemented to detect drift, the output of this
    % script is not filtered differently to the input):
    lowpassfilteredEEG = RELAX_filtbutter( epochedEEG, [], 5, 4, 'lowpass' );
    Message = ['Filtering here only performed to better detect drift, output data will still be bandpass filtered from ', num2str(RELAX_cfg.HighPassFilter), ' to ', num2str(RELAX_cfg.LowPassFilter)];
    disp(Message);
    
    % Robust average re-reference the data first (this creates a more equitable
    % spread of the data from each channel. Prior to this, channels close
    % to the common reference show very little signal, biasing the drift
    % detection threshold towards 0).
    [EEGavgref] = RELAX_average_rereference(lowpassfilteredEEG);
    Message = ['Average re-referencing here only performed to better detect drift, data for the next processing steps will not been average re-referenced yet at this point'];
    disp(Message);

    % Calculate the median across all channels, and the MAD, and the
    % threshold based on these from the number of MAD from the median that
    % you have chosen as the threshold:
    MedianAmplitude=median(EEGavgref.data(:,:,:),1);
    MADAmplitude=mad(EEGavgref.data(:,:,:),1);
    UpperBound=MedianAmplitude+(RELAX_cfg.DriftSeverityThreshold*MADAmplitude);
    LowerBound=MedianAmplitude-(RELAX_cfg.DriftSeverityThreshold*MADAmplitude);

    % If voltage exceeds MAD * your drift severity threshold for a total of
    % 100ms then mark the epoch as showing drift in the MWF mask 
    for epoch = 1:size(epochedEEG.data,3)
        for channel = 1:size(epochedEEG.data,1)
            if sum(EEGavgref.data(channel,:,epoch)>UpperBound(1,:,epoch),2)>round(100/RELAX_cfg.ms_per_sample) % do 100ms worth of data within this epoch exceed the drift threshold for this electrode?
                epochedEEG.RELAXProcessing.Details.driftmaskepochs(1,epoch)=1; % values in first column (electrodes) exceed upper bound, so mark as drift
            end
            if sum(EEGavgref.data(channel,:,epoch)<LowerBound(1,:,epoch),2)>round(100/RELAX_cfg.ms_per_sample) % do 100ms worth of data within this epoch exceed the drift threshold for this electrode?
                epochedEEG.RELAXProcessing.Details.driftmaskepochs(1,epoch)=1; % values in first column (electrodes) exceed lower bound, so mark as drift
            end
        end
    end
    
    % Exclude epochs that were previously marked as to be excluded because
    % they were extreme outliers:
    for e=1:size(epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents,2)
        if epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents(1,e)==1
            epochedEEG.RELAXProcessing.Details.driftmaskepochs(1,e)=NaN;
        end
    end

    epochedEEG.RELAXProcessing.ProportionMarkedAsDrift=mean(epochedEEG.RELAXProcessing.Details.driftmaskepochs,'omitnan');
    
    % If more than the threshold of epochs are marked, restrict number of epochs
    % marked to the worst epochs up to the threshold:
    if epochedEEG.RELAXProcessing.ProportionMarkedAsDrift>RELAX_cfg.ProportionWorstEpochsForDrift
        epochedEEG.RELAXProcessing.Details.driftmaskepochs=zeros(size(epochedEEG.RELAXProcessing.Details.driftmaskepochs));
        % Extract upper and lower bounds for each timepoint across all
        % channels (all channels have the same threshold to be compared
        % against)
        UpperBoundFullMatrix=zeros(size(epochedEEG.data));
        LowerBoundFullMatrix=zeros(size(epochedEEG.data));
        ExceedsUpperBound=nan(size(epochedEEG.data));
        ExceedsLowerBound=nan(size(epochedEEG.data));
        for x=1:size(epochedEEG.data,1)
            UpperBoundFullMatrix(x,:,:)=UpperBound(1,:,:);
            LowerBoundFullMatrix(x,:,:)=LowerBound(1,:,:);
        end
        % If any value at any timepoint for any channel exceeds the adjusted 
        % threshold, keep that value, otherwise the rest are left
        % as NaNs:
        for epoch = 1:size(epochedEEG.data,3)
            for timepoint = 1:size(epochedEEG.data,2)
                for channel = 1:size(epochedEEG.data,1)
                    if EEGavgref.data(channel,timepoint,epoch)>UpperBoundFullMatrix(channel,timepoint,epoch)
                        ExceedsUpperBound(channel,timepoint,epoch)=epochedEEG.data(channel,timepoint,epoch);
                    end
                    if EEGavgref.data(channel,timepoint,epoch)<LowerBoundFullMatrix(channel,timepoint,epoch)
                        ExceedsLowerBound(channel,timepoint,epoch)=epochedEEG.data(channel,timepoint,epoch);
                    end
                end
            end
        end
        
        % Sum all values that exceed the threshold within all channels for
        % each epoch:
        CumulativeExceedsUpperBoundPerEpoch=sum(ExceedsUpperBound,[1 2],'omitnan');
        CumulativeExceedsLowerBoundPerEpoch=sum(ExceedsLowerBound,[1 2],'omitnan');
        for e=1:size(epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents,2)
            if epochedEEG.RELAX.EpochsMarkedAsNaNsForExtremeEvents(1,e)==1
                CumulativeExceedsUpperBoundPerEpoch(1,e)=NaN;
                CumulativeExceedsLowerBoundPerEpoch(1,e)=NaN;
            end
        end
        upper=squeeze(sort(CumulativeExceedsUpperBoundPerEpoch,'descend'));
        lower=squeeze(sort(CumulativeExceedsLowerBoundPerEpoch,'ascend'));
        % Obtain new thresholds based on the percent you've selected to
        % mask as bad:
        DriftThresholdUpper = prctile(upper,100-((RELAX_cfg.ProportionWorstEpochsForDrift/2)*100));
        DriftThresholdLower = prctile(lower,((RELAX_cfg.ProportionWorstEpochsForDrift/2)*100));
        % If it exceeds the threshold, mark it as bad in the drift mask:
        for epoch = 1:size(epochedEEG.data,3)
            if CumulativeExceedsUpperBoundPerEpoch(1,1,epoch)>DriftThresholdUpper
                epochedEEG.RELAXProcessing.Details.driftmaskepochs(1,epoch)=1; 
            end
            if CumulativeExceedsLowerBoundPerEpoch(1,1,epoch)<DriftThresholdLower
                epochedEEG.RELAXProcessing.Details.driftmaskepochs(1,epoch)=1;
            end
        end  
    end
    epochedEEG.RELAXProcessing.Details.justdriftmask=epochedEEG.RELAXProcessing.Details.NaNsForNonEvents;
    % Add drift epoch marks to overall noise mask:
    for e=1:size(epochedEEG.RELAXProcessing.Details.driftmaskepochs,2)
        if epochedEEG.RELAXProcessing.Details.driftmaskepochs(1,e)==1
            epochedEEG.RELAXProcessing.Details.NoiseMaskFullLength(epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf1s;
            epochedEEG.RELAXProcessing.Details.justdriftmask(epochedEEG.event(e).originallatency:epochedEEG.event(e).originallatency+(round(1000/RELAX_cfg.ms_per_sample)-1))=OneSecondOf1s;
        end
    end
    epochedEEG.RELAXProcessing.ProportionMarkedAsDrift=mean(epochedEEG.RELAXProcessing.Details.justdriftmask,'omitnan');
    continuousEEG.RELAXProcessing=epochedEEG.RELAXProcessing;
    if isfield(epochedEEG, 'RELAX')==1
        continuousEEG.RELAX=epochedEEG.RELAX;
    end

end