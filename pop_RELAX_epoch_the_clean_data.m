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

%% pop_RELAX:
% Clean data with RELAX via the EEGLAB gui:
function [OutlierParticipantsToManuallyCheck,EpochRejections,RELAX_epoching_cfg] = pop_RELAX_epoch_the_clean_data(RELAX_epoching_cfg)

% Check if input is provided:
if nargin == 0
    RELAX_epoching_cfg = [];
end

%% Location of files for processing

% Specify the to be processed file locations:
if ~isfield(RELAX_epoching_cfg,'CleanedPath')
    RELAX_epoching_cfg.CleanedPath='D:\DATA_TO_BE_PREPROCESSED\RELAXProcessed\Cleaned_Data\';
end

if ~isfield(RELAX_epoching_cfg,'InterpolateRejectedChannels')
    RELAX_epoching_cfg.InterpolateRejectedChannels='yes';
end

% Specify rejection parameters:
if ~isfield(RELAX_epoching_cfg,'SingleChannelImprobableDataThreshold')
    RELAX_epoching_cfg.SingleChannelImprobableDataThreshold=5; %MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
end

if ~isfield(RELAX_epoching_cfg,'AllChannelImprobableDataThreshold')
    RELAX_epoching_cfg.AllChannelImprobableDataThreshold=3; %SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
end

if ~isfield(RELAX_epoching_cfg,'SingleChannelKurtosisThreshold')
    RELAX_epoching_cfg.SingleChannelKurtosisThreshold=5; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis 
end

if ~isfield(RELAX_epoching_cfg,'AllChannelKurtosisThreshold')
    RELAX_epoching_cfg.AllChannelKurtosisThreshold=3; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis 
end

if ~isfield(RELAX_epoching_cfg,'reject_amp')
    RELAX_epoching_cfg.reject_amp=60; % Absolute voltage amplitude threshold - if an epoch has voltages that deviate from 0 by more than this value the epoch is rejected
end

% The following muscle affected epoch rejection may 
% be useful if you want to ensure no epochs included in
% the analysis are still affected by muscle activity:
if ~isfield(RELAX_epoching_cfg,'MuscleSlopeThreshold')
    RELAX_epoching_cfg.MuscleSlopeThreshold=-0.31; % log frequency/ log power slope threshold for muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or more stringent = -0.72, more negative thresholds remove more muscle. 
end

if ~isfield(RELAX_epoching_cfg,'MaxProportionOfMuscleEpochsToClean')
    RELAX_epoching_cfg.MaxProportionOfMuscleEpochsToClean=0.50; % Maximum proportion of muscle epochs to remove (if more epochs than this are affected, only the worst effected up to this proportion will be removed)
end

if ~isfield(RELAX_epoching_cfg,'RemoveEpochsShowingMuscleActivity')
    RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity='no'; % 'yes' or 'no'
end

if ~isfield(RELAX_epoching_cfg,'DataType')
    RELAX_epoching_cfg.DataType='Task'; % 'Task' for cognitive tasks, 'Resting' for data without stimuli presented
end

if ~isfield(RELAX_epoching_cfg,'BLperiod')
    RELAX_epoching_cfg.BLperiod=[-200 0]; % Set baseline period for baseline correction of data to reduce potential influence of drift
end

if ~isfield(RELAX_epoching_cfg,'restingdatatriggerinterval')
    RELAX_epoching_cfg.restingdatatriggerinterval=3.5; % (seconds) sets how often an 'X' trigger is inserted into resting data for epoching. 3.5 with [-2.5 2.5] PeriodToEpoch provides 5s epochs with 1.5s overlaps
end

if ~isfield(RELAX_epoching_cfg,'TriggersToEpoch')
    RELAX_epoching_cfg.TriggersToEpoch={'HappyGo' 'HappyNogo' 'SadGo' 'SadNogo' }; % triggers to use for epoching
end

if ~isfield(RELAX_epoching_cfg,'PeriodToEpoch')
    RELAX_epoching_cfg.PeriodToEpoch=[-0.5 1.0]; % [start end] surrounding trigger
end

if ~isfield(RELAX_epoching_cfg,'RemoveOtherTriggers')
    RELAX_epoching_cfg.RemoveOtherTriggers='yes'; % 'yes'|'no' - removes all other triggers except the trigger used for timelocking the epoch
end

if ~isfield(RELAX_epoching_cfg,'BL_correction_method')
    RELAX_epoching_cfg.BL_correction_method='regression'; % Set type of baseline correction to perform ('regression' = recommended, 'subtraction' = traditional but not recommended)
end

if ~isfield(RELAX_epoching_cfg,'BL_correction_Factor_1_Level_1')
    RELAX_epoching_cfg.BL_correction_Factor_1_Level_1={'HappyGo' 'SadGo' }; % triggers to include in factor 1's level 1 for regression BL correction (all other triggers are included in factor 1's level 2)
end

if ~isfield(RELAX_epoching_cfg,'BL_correction_Factor_2_Level_1')
    RELAX_epoching_cfg.BL_correction_Factor_2_Level_1={'HappyGo' 'HappyNogo' }; % triggers to include in factor 2's level 1 for regression BL correction (all other triggers are included in factor 2's level 2)
end

if ~isfield(RELAX_epoching_cfg,'NumberOfFactors')
    RELAX_epoching_cfg.NumberOfFactors=2;
end

if ~isfield(RELAX_epoching_cfg,'FilesToProcess')
    RELAX_epoching_cfg.FilesToProcess=[1 3];
end

%% GENERATE THE GUI

% Only launch GUI if cfg stucture not provided
if nargin == 0

% Callback command for getting data folder
commandload1 = [ '[filepath] = uigetdir(''*'', ''Select the folder containing cleaned data'');' ...
    'if filepath ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''datadir''), ''string'', filepath );' ...
    'end;' ...
    'clear filepath tagtest;' ];

% ICA options
YesNoOptions = {'yes','no'};
TaskOrRestOptions = {'Task', 'Resting'};
BLCorrectionOptions = {'regression', 'subtraction', 'none'};

% GUI layout
geometry = {[0.3 0.8 0.1 0.5 0.2] ... % setting directory / file to process & interpolate rejected channels
            1 ...
            [0.4 0.4 1 0.2] ... % data type and resting interval for triggers
            1 ...
            [1 1.5 .5 0.3] ... % triggers to epoch around, period to epoch 
            1 ...
            [1 0.1] ... % remove other triggers?
            1 ...
            [1 0.6 1 0.3 1] ... % BL correction method and window
            1 ...
            [1 1 1 2] ... % regression BL correction design
            1 ...
            [1 1 1 1] ... % regression BL correction design
            1 ...
            [1] ...
            1 ...
            [1 0.2 1 0.2 1 0.2]... % epoch rejection thresholds
            1 ...
            [1 0.2 1 0.2]... % kurtosis epoch rejection thresholds
            1 ...
            [1 0.2] ...
            1 ...
            [1 0.2 1 0.2]... % muscle rejection settings 
            1 ...
            [5 1]... % files to process
            };
            
% GUI settings with defaults:        
uilist = {{'style', 'text', 'string', 'Cleaned data folder:','fontweight','bold','fontsize', 9} ...
          { 'style' 'edit'       'string' RELAX_epoching_cfg.CleanedPath 'tag' 'datadir' ,'fontsize', 9} ... 
          { 'style' 'pushbutton' 'string' '...' 'callback' commandload1 ,'fontsize', 9}... 
          {'style', 'text', 'string', 'Interpolate Rejected Channels?','fontsize', 9} ...
          {'style', 'popupmenu', 'string', YesNoOptions, 'tag', 'InterpolateOpts','Value',1,'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Data Type?','fontsize', 9} ...
          {'style', 'popupmenu', 'string', TaskOrRestOptions, 'tag', 'TaskOrRestOpts','Value',1,'fontsize', 9} ...
          {'style', 'text', 'string', 'If resting data, how long between inserted triggers? (s)','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.restingdatatriggerinterval,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'Triggers to epoch around:','fontsize', 9} ...
          {'style', 'edit', 'string', [RELAX_epoching_cfg.TriggersToEpoch{1,1} ' ' RELAX_epoching_cfg.TriggersToEpoch{1,2} ' ' RELAX_epoching_cfg.TriggersToEpoch{1,3}...
          ' ' RELAX_epoching_cfg.TriggersToEpoch{1,4}],'fontsize', 9} ...
          {'style', 'text', 'string', 'Period to epoch (s):','fontsize', 9} ...
          {'style', 'edit', 'string', [num2str(RELAX_epoching_cfg.PeriodToEpoch(1,1)), '  ', num2str(RELAX_epoching_cfg.PeriodToEpoch(1,2))],'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'Remove other triggers? (necessary for regression BL correction method)','fontsize', 9} ...
          {'style', 'popupmenu', 'string', YesNoOptions, 'tag', 'RemoveOtherTriggersOpts','Value',1,'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Baseline Correction Method:','fontsize', 9} ...
          {'style', 'popupmenu', 'string', BLCorrectionOptions, 'tag', 'BLCorrectionOpts','Value',1,'fontsize', 9} ...
          {'style', 'text', 'string', 'Baseline Correction Window (ms):','fontsize', 9} ...
          {'style', 'edit', 'string', [num2str(RELAX_epoching_cfg.BLperiod(1,1)), '  ', num2str(RELAX_epoching_cfg.BLperiod(1,2))],'fontsize', 9}...
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Regression BL correction design:','fontsize', 9} ...
          {'style', 'text', 'string', 'Number of Factors (max = 2):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.NumberOfFactors,'fontsize', 9}...
          {'style', 'text', 'string', '(if triggers are not listed below at Level 1, they are included in Level 2 for that factor)','fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Triggers to include in Factor 1, Level 1:','fontsize', 9} ...
          {'style', 'edit', 'string', [RELAX_epoching_cfg.BL_correction_Factor_1_Level_1{1,1} ' ' RELAX_epoching_cfg.BL_correction_Factor_1_Level_1{1,2}],'fontsize', 9} ...
          {'style', 'text', 'string', 'Triggers to include in Factor 2, Level 1:','fontsize', 9} ...
          {'style', 'edit', 'string', [RELAX_epoching_cfg.BL_correction_Factor_2_Level_1{1,1} ' ' RELAX_epoching_cfg.BL_correction_Factor_2_Level_1{1,2}],'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Epoch Rejection Thresholds:','fontweight','bold','fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Single Channel Improbable Data Threshold (SD):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.SingleChannelImprobableDataThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'All Channel Improbable Data Threshold (SD):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.AllChannelImprobableDataThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'Absolute voltage amplitude threshold (microvolts):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.reject_amp,'fontsize', 9}...
          {}...  
          {'style', 'text', 'string', 'Single Channel Kurtosis Threshold (SD):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.SingleChannelKurtosisThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'All Channel Kurtosis Threshold (SD):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.AllChannelKurtosisThreshold,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'Reject Muscle Contaminated Epochs?','fontsize', 9} ...
          {'style', 'popupmenu', 'string', YesNoOptions, 'tag', 'RejectMuscleOpts','Value',2,'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Log-Frequency Log-Power Slope Threshold for detecting muscle activity:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.MuscleSlopeThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'Maximum Proportion of epochs can be removed because of muscle:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_epoching_cfg.MaxProportionOfMuscleEpochsToClean,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'File numbers to process this session (from start file # to finish file # and files in between):','fontsize', 9} ...
          {'style', 'edit', 'string', ['1', '  ', '3'],'fontsize', 9}...
          };

result = inputgui('geometry', geometry, 'geomvert', [1 .4 1 .4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1],  'uilist', uilist, 'title', 'RELAX Parameter Setting',  'helpcom', 'pophelp(''pop_RELAX'')');

% Replace default values with inputs
RELAX_epoching_cfg.CleanedPath = result{1};
RELAX_epoching_cfg.InterpolateRejectedChannels = YesNoOptions{1,result{2}};
RELAX_epoching_cfg.DataType = TaskOrRestOptions{1,result{3}};
RELAX_epoching_cfg.restingdatatriggerinterval = str2double(result{4});
RELAX_epoching_cfg.TriggersToEpoch = strtrim(strsplit(strrep(result{5},',','')));
PeriodToEpoch=strsplit(strrep(result{6},',',''));
RELAX_epoching_cfg.PeriodToEpoch=[str2double(PeriodToEpoch{1}) str2double(PeriodToEpoch{2})];
RELAX_epoching_cfg.RemoveOtherTriggers= YesNoOptions{1,result{7}};
RELAX_epoching_cfg.BL_correction_method=BLCorrectionOptions{1,result{8}};
BLperiod=strsplit(strrep(result{9},',',''));
RELAX_epoching_cfg.BLperiod=[str2double(BLperiod{1}) str2double(BLperiod{2})];
RELAX_epoching_cfg.NumberOfFactors=str2double(result{10});
RELAX_epoching_cfg.BL_correction_Factor_1_Level_1=strtrim(strsplit(strrep(result{11},',','')));
RELAX_epoching_cfg.BL_correction_Factor_2_Level_1=strtrim(strsplit(strrep(result{12},',','')));
RELAX_epoching_cfg.SingleChannelImprobableDataThreshold= str2double(result{13});
RELAX_epoching_cfg.AllChannelImprobableDataThreshold= str2double(result{14});
RELAX_epoching_cfg.reject_amp= str2double(result{15});
RELAX_epoching_cfg.SingleChannelKurtosisThreshold= str2double(result{16});
RELAX_epoching_cfg.AllChannelKurtosisThreshold= str2double(result{17});
RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity= YesNoOptions{1,result{18}};
RELAX_epoching_cfg.MuscleSlopeThreshold= str2double(result{19});
RELAX_epoching_cfg.MaxProportionOfMuscleEpochsToClean= str2double(result{20});
FilesToProcess=strsplit(strrep(result{21},',',''));
RELAX_epoching_cfg.FilesToProcess=str2double(FilesToProcess(1,1)):str2double(FilesToProcess(1,2));

if ~isfield(RELAX_epoching_cfg,'filename')
    RELAX_epoching_cfg.filename = [];
end

end

%% Check for dependencies:

if (exist('ft_freqanalysis','file')==0) && strcmp(RELAX_epoching_cfg.RemoveEpochsShowingMuscleActivity,'yes')
    warndlg('fieldtrip may not be installed, or the folder path for fieldtrip has not been set in MATLAB. Plugin can be installed via EEGLAB: "File" > "Manage EEGLAB Extensions"','fieldtrip might not be installed');
end

%%
% List all files in directory
cd(RELAX_epoching_cfg.CleanedPath);
RELAX_epoching_cfg.dirList=dir('*.set');
RELAX_epoching_cfg.files={RELAX_epoching_cfg.dirList.name};
if isempty(RELAX_epoching_cfg.files)
    disp('No files found..')
end

if ~isfield(RELAX_epoching_cfg,'FilesToProcess')
    RELAX_epoching_cfg.FilesToProcess=1;
end

[OutlierParticipantsToManuallyCheck,EpochRejections,RELAX_epoching_cfg] = RELAX_epoch_the_clean_data_Wrapper (RELAX_epoching_cfg);

% To enable debugging if necessary:
% 
% OutlierParticipantsToManuallyCheck={}; 
% EpochRejections={};
% RawMetrics={};
% RELAXProcessingRoundOneAllParticipants={};
% RELAXProcessingRoundTwoAllParticipants={};
% RELAXProcessing_wICA_AllParticipants={};
% RELAXProcessingRoundThreeAllParticipants={};
% RELAX_issues_to_check={};
% RELAXProcessingExtremeRejectionsAllParticipants={};