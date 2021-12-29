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

%% RELAX_EPOCH_CLEAN_DATA_FOR_ANALYSIS:

clear all; close all; clc;
%% DEPENDENCIES (toolboxes you need to install, and cite if you use this script):

% MATLAB signal processing toolbox (from MATLAB website)

% EEGLAB:
% https://sccn.ucsd.edu/eeglab/index.php
% Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of neuroscience methods, 134(1), 9-21.
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_current\eeglab2019_1');
eeglab;

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.
addpath('C:\Program Files\MATLAB\fieldtrip-20180805');

% Specify  RELAX folder location (this toolbox):
addpath('D:\Data_Analysis\RELAX_v0_91\');

% Specify rejection parameters:
RELAX_epoching_cfg.SingleChannelImprobableDataThreshold=5; %MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
RELAX_epoching_cfg.AllChannelImprobableDataThreshold=3; %SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
RELAX_epoching_cfg.SingleChannelKurtosisThreshold=5; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis 
RELAX_epoching_cfg.AllChannelKurtosisThreshold=3; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis 
RELAX_epoching_cfg.reject_amp=60; % Absolute voltage amplitude threshold - if an epoch has voltages that deviate from 0 by more than this value the epoch is rejected

% The following muscle affected epoch rejection may 
% be useful if you want to ensure no epochs included in
% the analysis are still affected by muscle activity:
RELAX_epoching_cfg.MuscleSlopeThreshold=-0.31; % log frequency/ log power slope threshold for muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or more stringent = -0.72, more negative thresholds remove more muscle. 
RELAX_epoching_cfg.MaxProportionOfMuscleEpochsToClean=0.50; % Maximum proportion of muscle epochs to remove (if more epochs than this are affected, only the worst effected up to this proportion will be removed)
RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity='yes'; % 'yes' or 'no'

RELAX_epoching_cfg.DataType='Task'; % 'Task' for cognitive tasks, 'Resting' for data without stimuli presented

% Baseline correct the data?
RELAX_epoching_cfg.Perform_BL_Correction='yes'; % 'yes' or 'no'
RELAX_epoching_cfg.BLperiod=[-200 0]; % Set baseline period for baseline correction of data to reduce potential influence of drift

% Specify the to be processed file locations:
RELAX_epoching_cfg.myPath='D:\DATA_TO_BE_PREPROCESSED\';
RELAX_epoching_cfg.CleanedPath=[RELAX_epoching_cfg.myPath filesep 'RELAXProcessed' filesep 'Cleaned_Data'];

% Load pre-processing statistics file for these participants if it already
% exists:
RELAX_epoching_cfg.OutputPath=[RELAX_epoching_cfg.CleanedPath filesep 'Epoched'];
mkdir(RELAX_epoching_cfg.OutputPath);

%% List all files in directory
cd(RELAX_epoching_cfg.CleanedPath);
RELAX_epoching_cfg.dirList=dir('*.set');
RELAX_epoching_cfg.files={RELAX_epoching_cfg.dirList.name};
if isempty(RELAX_epoching_cfg.files)
    disp('No files found..')
end

%% Loop over each file in your list 
for FileNumber=1:numel(RELAX_epoching_cfg.files)
    RELAX_epoching_cfg.filename=RELAX_epoching_cfg.files{FileNumber};
    clearvars -except 'FilesWithoutConvergence' 'RELAX_epoching_cfg' 'FileNumber' 'FileName' 'Participant_IDs' 'Medianvoltageshiftwithinepoch' 'EpochRejections';
    %% Load data (assuming the data is in EEGLAB .set format):
    EEG = pop_loadset(RELAX_epoching_cfg.filename);
    FileName = extractBefore(RELAX_epoching_cfg.files{FileNumber},".");
    Participant_IDs{1,FileNumber} = extractBefore(RELAX_epoching_cfg.files{FileNumber},".");
    EEG.RELAXProcessing.ChannelsRemaining=EEG.nbchan;

    %% Interpolate channels that were excluded back into the data:
    EEG = pop_interp(EEG, EEG.allchan, 'spherical');

    % Record how many channels had to be interpolated back into the data after cleaning:
    EEG.RELAXProcessing.TotalChannels=EEG.nbchan;
    EEG.RELAXProcessing.ProportionOfChannelsInterpolated=(EEG.RELAXProcessing.TotalChannels-EEG.RELAXProcessing.ChannelsRemaining)/EEG.RELAXProcessing.TotalChannels;

    %% Epoch Data:
    
    % If resting data, you can use this to insert triggers at a specified
    % interval and then epoch data around triggers
    % (this example creates 5 second epochs with 1.5 second overlaps):
 
    if strcmp(RELAX_epoching_cfg.DataType,'Resting')==1
        restingdatatriggerinterval=3.5;
        EEG=eeg_regepochs(EEG,'recurrence',restingdatatriggerinterval,'eventtype','X','extractepochs','off');      
        EEG = pop_epoch( EEG, { 'X'}, [-2.5 2.5], 'epochinfo', 'yes');
    end
     
    % If cognitive data with triggers embedded, this epochs data from -0.5 to 1 seconds and deletes triggers other than the epoched trigger within each epoch. Example below with an emotional Go Nogo task:  
    if strcmp(RELAX_epoching_cfg.DataType,'Task')
        EEG = pop_epoch( EEG, {'HappyGo' 'HappyNogo' 'SadGo' 'SadNogo' }, [-0.5 1.0], 'epochinfo', 'yes');
        % Remove triggers that are in the epoch, but aren't the trigger
        % the data is being epoched around (this can be helpful if epoching
        % data separately by condition following this script):
        EEG = pop_selectevent( EEG, 'omitlatency', '-1001<=-1','type', {'HappyGo' 'HappyNogo' 'SadGo' 'SadNogo' }, 'deleteevents','on');
        EEG = pop_selectevent( EEG, 'omitlatency', '1<=1999','type', {'HappyGo' 'HappyNogo' 'SadGo' 'SadNogo' }, 'deleteevents','on');
    end
    EEG = eeg_checkset( EEG );

    %% Baseline Correct Data:
    
    if strcmp(RELAX_epoching_cfg.Perform_BL_Correction,'yes')

        % Regression based baseline correction method (recommended):
        [EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,'Factor_1_Level_1', {'HappyGo' 'SadGo' }, 'Factor_2_Level_1', {'HappyGo' 'HappyNogo' }); 
        % if a 2 x 2 design, this will code triggers other than 'HappyGo'/'SadGo' as -1, and Go as 1 in the first factor, and triggers other than 'HappyGo'/'HappyNogo' as -1 in the second factor
        
        %[EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,'Factor_1_Level_1',{'Go'}); if a 2 condition design, this will code triggers other than 'Go' as -1, and Go as 1
        %[EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg); % if only 1 stimulus condition present for each participant
        
        % (the benefits of this method are explained in Alday, 2019, and the specific implementation performed in Bailey et al. 2022)
        % Alday, P. M. (2019). How much baseline correction do we need in ERP research? Extended GLM model can replace baseline correction while lifting its limits. Psychophysiology, 56(12), e13451.
        % Bailey et al. (2022). Meditators probably show increased behaviour-monitoring related neural activity. 
        
        % Note that for designs with more than 3 conditions for a single
        % factor, linear mixed modelling to baseline correct the data is more
        % appropriate (not provided here).

        % Also note that the method provided below is effective for 
        % nice clean data with a larger number of epochs per participant,
        % but for data that is not clean or very few epochs per participant (<5 is my guess) 
        % then a method that includes all participants in the single regression 
        % and includes individual as a factor is more appropriate.
        
        % Traditional subtraction based baseline correction (not recommended):
        % EEG = pop_rmbase( EEG, [RELAX_epoching_cfg.BLperiod]);
        
    end
 
    %% THIS SECTION CONTAINS FUNCTIONS WHICH DETECT AND REJECT ARTIFACTS

    % Count initial epochs:
    EEG.RELAXProcessing.InitialEpochs=size(EEG.data,3);

    % Any one of these functions can be commented out to ignore those artifacts
    % when creating the mask:

    % This section uses traditional amplitude, improbable voltage distributions within epochs, and kurtosis to reject epochs:
    ROIidx= 1:EEG.nbchan; 
    EEG = pop_eegthresh(EEG,1,[ROIidx],[-RELAX_epoching_cfg.reject_amp],[RELAX_epoching_cfg.reject_amp],[EEG.xmin],[EEG.xmax],1,0);
    EEG = pop_jointprob(EEG,1,[ROIidx],RELAX_epoching_cfg.SingleChannelImprobableDataThreshold,RELAX_epoching_cfg.AllChannelImprobableDataThreshold,1,0);
    EEG = pop_rejkurt(EEG,1,(1:EEG.nbchan),RELAX_epoching_cfg.SingleChannelKurtosisThreshold,RELAX_epoching_cfg.AllChannelKurtosisThreshold,1,0);
    EEG = eeg_rejsuperpose(EEG, 1, 0, 1, 1, 1, 1, 1, 1);
    EEG = pop_rejepoch(EEG, [EEG.reject.rejglobal] ,0);
        
    %% If you have not filtered out data below 75Hz, you could use an objective muscle slope measure to reject epochs with remaining muscle activity:   
    % Use epoched data and FFT to detect slope of log frequency log
    % power, add periods exceeding muscle threshold to mask. This method is
    % designed for use with data that includes up to 75Hz, so  is not
    % useful if frequencies below 75Hz are filtered out

    if strcmp(RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity,'yes')
        [EEG] = RELAX_Rejecting_muscle_epochs(EEG, RELAX_epoching_cfg);
    end
    
    %% CHECK MEDIAN VOLTAGE SHIFT IN EACH EPOCH FOR EACH PARTICIPANT, AND RECORD ARTIFACT REJECTION DETAILS
    % The following measures the median size of the voltage shift in each epoch
    % for each participant. Used at the end of this script to provide
    % advice on which participants to manually check as potential bad
    % data:        
    voltageshiftwithinepoch=range(EEG.data(:,:,:),2);
    Medianvoltageshiftwithinepoch(:,FileNumber)=median(voltageshiftwithinepoch,3);
    EEG.RELAXProcessing.EpochsRemaining=size(EEG.data,3);
    EEG.RELAXProcessing.ProportionOfEpochsRejected=(EEG.RELAXProcessing.InitialEpochs-EEG.RELAXProcessing.EpochsRemaining)/EEG.RELAXProcessing.InitialEpochs;
    EEG.RELAXProcessing.aFileName=FileName;

     % Order the MWF Processing statistics structure in alphabetical order:
    [~, neworder] = sort(lower(fieldnames(EEG.RELAXProcessing)));
    EEG.EpochRejections = orderfields(EEG.RELAXProcessing, neworder);
    if isfield(EEG.EpochRejections, 'NaNsForExtremeOutlierPeriods')
        EEG.EpochRejections=rmfield(EEG.EpochRejections,'NaNsForExtremeOutlierPeriods');
    end
    EpochRejections(FileNumber,:) = struct2table(EEG.EpochRejections,'AsArray',true);
    
    % If the ICA didn't converge when cleaning the data, a note is made
    % here:
    if isfield(EEG.RELAX, 'NonConvergence')
        FilesWithoutConvergence{FileNumber}=FileName;
    end

    %% Save data:
    SaveSetMWF2 =[RELAX_epoching_cfg.OutputPath filesep FileName '_Epoched.set'];    
    EEG = pop_saveset( EEG, SaveSetMWF2 );  
end

%% The following checks for participants who show outlying values for the median voltage shift within each epoch:
% The following detects outlier files in the median amount of their max-min
% voltage shift within an epoch, after adjusting for the fact that the data
% across all participants is likely to be positively skewed with a log transform.
MedianvoltageshiftwithinepochLogged=log10(Medianvoltageshiftwithinepoch);
InterQuartileRange=iqr(MedianvoltageshiftwithinepochLogged,2);
Upper25 = prctile(MedianvoltageshiftwithinepochLogged,75,2);
Lower25 = prctile(MedianvoltageshiftwithinepochLogged,25,2);

% 75th% and 25th% +/- (1.5 x IQR) is the recommended outlier detection method, 
% so this is used to recommend which participants to manually check
% However, I find this can be a bit too sensitive upon manual inspection,
% and that 1.75, 2, or even 2.5 can be a better threshold.
LowerBound=size(MedianvoltageshiftwithinepochLogged,1);
UpperBound=size(MedianvoltageshiftwithinepochLogged,1);
for x=1:size(MedianvoltageshiftwithinepochLogged,1)
    LowerBound(x,1)=Lower25(x,1)-(2*InterQuartileRange(x,1));
    UpperBound(x,1)=Upper25(x,1)+(2*InterQuartileRange(x,1));
end

VoltageShiftsTooLow=MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooLow=VoltageShiftsTooLow-LowerBound;
VoltageShiftsTooLow(0<VoltageShiftsTooLow)=0;
CumulativeSeverityOfAmplitudesBelowThreshold=sum(VoltageShiftsTooLow,1)';

VoltageShiftsTooHigh=MedianvoltageshiftwithinepochLogged;
VoltageShiftsTooHigh=VoltageShiftsTooHigh-UpperBound;
VoltageShiftsTooHigh(0>VoltageShiftsTooHigh)=0;
CumulativeSeverityOfAmplitudesAboveThreshold=sum(VoltageShiftsTooHigh,1)';

% Plot:
plot(LowerBound); hold on; plot(UpperBound); 
% for c=1:size(EEG.chanlocs,2); electrode{c}=EEG.chanlocs(c).labels;end
hold on; plot(MedianvoltageshiftwithinepochLogged); xticks([1:1:60]);xticklabels({EEG.chanlocs.labels});legend('LowerBound', 'UpperBound');

OutlierParticipantsToManuallyCheck = table(Participant_IDs', CumulativeSeverityOfAmplitudesBelowThreshold,CumulativeSeverityOfAmplitudesAboveThreshold);

savefileone=[RELAX_epoching_cfg.myPath 'RELAXProcessed' filesep 'OutlierParticipantsToManuallyCheck'];
save(savefileone,'OutlierParticipantsToManuallyCheck')
 
savefileone=[RELAX_epoching_cfg.myPath 'RELAXProcessed' filesep 'EpochRejections'];
save(savefileone,'EpochRejections')

LoggedMedianVoltageShiftAcrossEpochs=array2table(MedianvoltageshiftwithinepochLogged);
LoggedMedianVoltageShiftAcrossEpochs.Properties.VariableNames =Participant_IDs';
for c=1:size(EEG.chanlocs,2); chanlist{c}=EEG.chanlocs(c).labels; end
LoggedMedianVoltageShiftAcrossEpochs.Properties.RowNames =chanlist;
savefileone=[RELAX_epoching_cfg.myPath 'RELAXProcessed' filesep 'LoggedMedianVoltageShiftAcrossEpochs'];
save(savefileone,'LoggedMedianVoltageShiftAcrossEpochs')

savefileone=[RELAX_epoching_cfg.myPath 'RELAXProcessed' filesep 'RELAX_epoching_cfg'];
save(savefileone,'RELAX_epoching_cfg')

clearvars -except 'OutlierParticipantsToManuallyCheck' 'RELAX_epoching_cfg' 'EpochRejections' 'EpochRejectionStats' 'FilesWithoutConvergence' 'LoggedMedianVoltageShiftAcrossEpochs';
        