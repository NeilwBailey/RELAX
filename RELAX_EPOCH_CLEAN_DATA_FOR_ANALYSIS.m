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
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_2021_1');
eeglab;

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_2021_1\plugins\Fieldtrip-lite20210601');

% Specify  RELAX folder location (this toolbox):
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_2021_1\plugins\RELAX-1.0.0\');

% Specify rejection parameters:
RELAX_epoching_cfg.InterpolateRejectedChannels='yes';
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
RELAX_epoching_cfg.restingdatatriggerinterval=3.5; % (seconds) sets how often an 'X' trigger is inserted into resting data for epoching. 3.5 with [-2.5 2.5] PeriodToEpoch provides 5s epochs with 1.5s overlaps
RELAX_epoching_cfg.TriggersToEpoch={'HappyGo' 'HappyNogo' 'SadGo' 'SadNogo' }; % triggers to use for epoching
RELAX_epoching_cfg.PeriodToEpoch=[-0.5 1.0]; % [start end] surrounding trigger
RELAX_epoching_cfg.RemoveOtherTriggers='yes'; % 'yes'|'no' - removes all other triggers except the trigger used for timelocking the epoch

% Baseline correct the data?
RELAX_epoching_cfg.BL_correction_method='regression'; % Set type of baseline correction to perform ('regression' = recommended, 'subtraction' = traditional but not recommended)
RELAX_epoching_cfg.BLperiod=[-200 0]; % Set baseline period for baseline correction of data to reduce potential influence of drift
RELAX_epoching_cfg.NumberOfFactors=2;
RELAX_epoching_cfg.BL_correction_Factor_1_Level_1={'HappyGo' 'SadGo' }; % triggers to include in factor 1's level 1 for regression BL correction (all other triggers are included in factor 1's level 2)
RELAX_epoching_cfg.BL_correction_Factor_2_Level_1={'HappyGo' 'HappyNogo' }; % triggers to include in factor 2's level 1 for regression BL correction (all other triggers are included in factor 2's level 2)

% Specify the to be processed file locations:
RELAX_epoching_cfg.CleanedPath=['D:\DATA_TO_BE_PREPROCESSED' filesep 'RELAXProcessed' filesep 'Cleaned_Data'];

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

%% RUN SCRIPT BELOW:
RELAX_epoching_cfg.FilesToProcess=1:numel(RELAX_epoching_cfg.files); % Set which files to process

[OutlierParticipantsToManuallyCheck,EpochRejections,RELAX_epoching_cfg] = RELAX_epoch_clean_data_Wrapper(RELAX_epoching_cfg);
    
        