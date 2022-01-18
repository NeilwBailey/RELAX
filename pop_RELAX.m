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
function [RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingRoundOneAllParticipants, RELAXProcessingRoundTwoAllParticipants, RELAXProcessing_wICA_AllParticipants,...
        RELAXProcessingRoundThreeAllParticipants, RELAX_issues_to_check, RELAXProcessingExtremeRejectionsAllParticipants] = pop_RELAX(varargin)

%% DEPENDENCIES (toolboxes you need to install, and cite if you use this script):
% use fileseparators 'filesep' for increased compatability if necessary (replace the \ with a filesep [outside of quotes]) 

% MATLAB signal processing toolbox (from MATLAB website)
% MATLAB statistics and machine learning toolbox (from MATLAB website)

% PREP pipeline to reject bad electrodes (install plugin to EEGLAB, or via the github into the EEGLAB plugin folder): 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
% Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., & Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for large-scale EEG analysis. Frontiers in neuroinformatics, 9, 16.
% http://vislab.github.io/EEG-Clean-Tools/ or install via EEGLAB extensions

% Specify the MWF path:
% https://github.com/exporl/mwf-artifact-removal
% Somers, B., Francart, T., & Bertrand, A. (2018). A generic EEG artifact removal algorithm based on the multi-channel Wiener filter. Journal of neural engineering, 15(3), 036007.
% addpath(genpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_current\eeglab2019_1\plugins\mwf-artifact-removal-master'));

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.
% addpath('C:\Program Files\MATLAB\fieldtrip-20180805');

% TESA for filtering (install plugin to EEGLAB):
% https://github.com/nigelrogasch/TESA
% Rogasch, N. C., Sullivan, C., Thomson, R. H., Rose, N. S., Bailey, N. W., Fitzgerald, P. B., ... & Hernandez-Pavon, J. C. (2017). Analysing concurrent transcranial magnetic stimulation and electroencephalographic data: A review and introduction to the open-source TESA software. Neuroimage, 147, 934-951.

% fastica:
% http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml 
% addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_current\eeglab2019_1\plugins\FastICA_25\');

% ICLabel in your eeglab folder as a plugin or via the github:
% https://github.com/sccn/ICLabel

% Need to install RunLength: https://au.mathworks.com/matlabcentral/fileexchange/41813-runlength?s_tid=mwa_osa_a
% This can be tricky on a mac: you need to manually compile the 'RunLength' function as it wouldn't auto-install through the Matlab script. 
% Macs have very tight security settings around MEX files. To get around this, install Xcode for Mac and then run the line of code to get it working
% addpath('D:\Data_Analysis\Analysis_Tools_and_Software\RunLength'); 

% Need to install MinGW-w64 if on windows if you haven't already:
% http://mingw-w64.org/doku.php

% Specify  RELAX folder location (this toolbox):
% addpath('D:\Data_Analysis\RELAX_v0_91\');

%% Location of files for processing

% Specify your electrode locations with the correct cap file:
RELAX_cfg.caploc='D:\Data_Analysis\Cap_Location_Files\standard-10-5-cap385.elp'; %path containing electrode positions

% Specify the to be processed file locations:
RELAX_cfg.myPath='D:\DATA_TO_BE_PREPROCESSED\';

% List all files in directory
cd(RELAX_cfg.myPath);
RELAX_cfg.dirList=dir('*.set');
RELAX_cfg.files={RELAX_cfg.dirList.name};
if isempty(RELAX_cfg.files)
    disp('No files found..')
end

%% Parameters that can be specified:

% The following selections determine how many epochs of each type of
% artifact to include in the mask. It's not obvious what the best choices are.
% I have selected as defaults the parameters that work best for our data.
% If only performing one MWF run, perhaps including all artifacts in the mask
% is best (as long as that leaves enough clean data for the clean data
% mask). Somers et al. (2018) suggest that it doesn't hurt to
% include clean data in the artifact mask, and so it's helpful to have wide
% boundaries around artifacts. 

% However, I wonder if some files might include the majority of trials in the artifact
% mask could mean that most of the task related ERPs might be considered
% artifacts, and cleaned out from the data. Including ERP trials in the clean
% mask reduces this issue, as the clean mask defines what the good
% data should look like, and if it has ERPs, then the artifact periods will
% be cleaned of only the difference between the clean mask and the artifact
% mask, leaving the ERP in the data. This relies upon a consistent
% ERP across both the clean and artifact masks however. An alternative fix
% could be to implement NaNs (which the MWF mask ignores) over ERP periods
% instead of marking them as artifacts, so the ERP is not included in the
% artifact mask, but this would mean any artifacts that are present in the
% ERPs but not outside the ERP periods aren't cleaned.

% Lastly, you may want to delete more unused electrodes than specified by
% this script. If so, modify the section titled: "Delete channels that
% are not relevant if present" to include the electrodes you would like to
% delete.

%% DEFAULT SETTINGS:

RELAX_cfg.Do_MWF_Once=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Twice=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Thrice=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no). I think cleaning drift in this is a good idea.
RELAX_cfg.Perform_wICA_on_ICLabel=1; % 1 = Perform wICA on artifact components marked by ICLabel (1 for yes, 0 for no).
RELAX_cfg.ICA_method='fastica_symm';

RELAX_cfg.computerawmetrics=1; % Compute blink and muscle metrics from the raw data?
RELAX_cfg.computecleanedmetrics=1; % Compute SER, ARR, blink and muscle metrics from the cleaned data?

RELAX_cfg.MWFRoundToCleanBlinks=2; % Which round to clean blinks in (1 for the first, 2 for the second...)
RELAX_cfg.ProbabilityDataHasNoBlinks=0; % 0 = data almost certainly has blinks, 1 = data might not have blinks, 2 = data definitely doesn't have blinks.
% 0 = eg. task related data where participants are focused with eyes open, 
% 1 = eg. eyes closed recordings, but with participants who might still open their eyes at times, 
% 2 = eg. eyes closed resting with highly compliant participants and recordings that were strictly made only when participants had their eyes closed.

RELAX_cfg.DriftSeverityThreshold=10; %MAD from the median of all electrodes. This could be set lower and would catch less severe drift 
RELAX_cfg.ProportionWorstEpochsForDrift=0.30; % Maximum proportion of epochs to include in the mask from drift artifact type.

RELAX_cfg.ExtremeVoltageShiftThreshold=20; % Threshold MAD from the median all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe voltage shifts within the epoch
RELAX_cfg.ExtremeAbsoluteVoltageThreshold=500; % microvolts max or min above which will be excluded from cleaning and deleted from data
RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; % Threshold SD from the mean of all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe improbable data
RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=8; % Threshold kurtosis of each electrode against the same electrode in different epochs. This could be set lower and would catch less severe kurtosis 
RELAX_cfg.ExtremeAllChannelKurtosisThreshold=8; % Threshold kurtosis across all electrodes. This could be set lower and would catch less severe kurtosis
RELAX_cfg.ExtremeDriftSlopeThreshold=-4; % slope of log frequency log power below which to reject as drift without neural activity
RELAX_cfg.ExtremeBlinkShiftThreshold=8; % How many MAD from the median across blink affected epochs to exclude as extreme data 
% (applies the higher value out of this value and the
% RELAX_cfg.ExtremeVoltageShiftThreshold above as the
% threshold, which caters for the fact that blinks don't affect
% the median, so without this, if data is clean and blinks are
% large, blinks can get excluded as extreme outliers)
            
% Clean periods that last for a shorter duration than the following value to be marked as artifacts, 
% and pad short artifact periods out into artifact periods of at least the following length when 
% they are shorter than this value to reduce rank deficiency issues in MWF cleaning). 
% Note that it's better to include clean periods in the artifact mask rather than the including artifact in the clean mask.
RELAX_cfg.MinimumArtifactDuration=1200; % in ms. It's better to make this value longer than 1000ms, as doing so will catch diminishing artifacts that aren't detected in a neighbouring 1000ms period, which might still be bad
RELAX_cfg.MinimumBlinkArtifactDuration=800; % blink marking is based on the maximum point of the blink rather than the 1000ms divisions for muscle artifacts, so this can be shorter than the value above (blinks do not typically last >500ms)

RELAX_cfg.BlinkElectrodes={'FP1';'FPZ';'FP2';'AF3';'AF4';'F3';'F1';'FZ';'F2';'F4'}; % sets the electrodes to average for blink detection using the IQR method. These should be frontal electrodes maximally affected by blinks. The order is the order of preference for icablinkmetrics.
% A single HOEG electrode for each side is selected by the script, prioritized in the following order (if the electrode in position 1 isn't present, the script will check for electrode in position 2, and so on...).
RELAX_cfg.HEOGLeftpattern = ["AF7", "F7", "FT7", "F5", "T7", "FC5", "C5", "TP7", "AF3"]; % sets left side electrodes to use for horizontal eye movement detection. These should be lateral electrodes maximally effected by blinks.
RELAX_cfg.HEOGRightpattern = ["AF8", "F8","FT8","F6","T8", "FC6", "C6", "TP8", "AF4"]; % sets right side electrodes to use for horizontal eye movement detection. These should be lateral electrodes maximally effected by blinks.
RELAX_cfg.BlinkMaskFocus=150; % this value decides how much data before and after the right and left base of the eye blink to mark as part of the blink artifact window. 
% I found 100ms on either side of the blink bases works best with a delay of 7 on the MWF. However, it also seemed to create too short artifact masks at times, which may lead to insufficient rank for MWF, so I left the default as 150ms.
RELAX_cfg.HorizontalEyeMovementType=2; % 1 to use the IQR method, 2 to use the MAD method for identifying threshold. IQR method less effective for smaller sample sizes (shorter files).
RELAX_cfg.HorizontalEyeMovementThreshold=2; % MAD deviation from the median that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below).
RELAX_cfg.HorizontalEyeMovementThresholdIQR=1.5; % If IQR method set above, IQR deviation that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below).
RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold=25; % The number of timepoints (ms) that exceed the horizontal eye movement threshold within the test period (set below) before the period is marked as horizontal eye movement.
RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow=(2*RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold)-1; % Window duration to test for horizontal eye movement, set to 2x the value above by default.
RELAX_cfg.HorizontalEyeMovementFocus=200; % Buffer window, masking periods earlier and later than the time where horizontal eye movements exceed the threshold.

RELAX_cfg.HighPassFilter=0.25; % Sets the high pass filter. 1Hz is best for ICA decomposition if you're examining just oscillatory data, 0.25Hz seems to be the highest before ERPs are adversely affected by filtering 
%(lower than 0.2Hz may be better, but I find a minority of my files show drift at 0.3Hz even).
if RELAX_cfg.HighPassFilter>0.25
    Warning='You have high pass filtered above 0.25, which can adversely affect ERP analyses';
end
RELAX_cfg.LowPassFilter=80; % If you filter out data below 75Hz, you can't use the objective muscle detection method
RELAX_cfg.LineNoiseFrequency=50; % Frequencies for bandstop filter in order to address line noise (set to 60 in countries with 60Hz line noise, and 50 in countries with 50Hz line noise).

RELAX_cfg.ElectrodesToDelete={'CB1'; 'CB2'; 'HEOG'; 'IO1'; 'M1'; 'M2'; 'LO1'; 'LO2'; 'E1'; 'E3'; 'ECG'; 'SO1'; 'ECG'; 'SPARE1'; 'SPARE2'; 'SPARE3'; 'BP1'; 'BP2'; 'VEOG'};
% If your EEG recording includes non-scalp electrodes or electrodes that you want to delete before cleaning, you can set them to be deleted here. 
% The RELAX cleaning pipeline does not need eye, heart, or mastoid electrodes for effective cleaning.

RELAX_cfg.KeepAllInfo=0; % setting this value to 1 keeps all the details from the MWF pre-processing and MWF computation. Helpful for debugging if necessary but makes for large file sizes.
RELAX_cfg.saveextremesrejected=0; % setting this value to 1 tells the script to save the data after only filtering, extreme channels have been rejected and extreme periods have been noted
RELAX_cfg.saveround1=0; % setting this value to 1 tells the script to save the first round of MWF pre-processing
RELAX_cfg.saveround2=0; % setting this value to 1 tells the script to save the second round of MWF pre-processing
RELAX_cfg.saveround3=0; % setting this value to 1 tells the script to save the third round of MWF pre-processing

RELAX_cfg.OnlyIncludeTaskRelatedEpochs=0; % If this =1, the MWF clean and artifact templates will only include data within 5 seconds of a task trigger (other periods will be marked as NaN, which the MWF script ignores).

RELAX_cfg.MuscleSlopeThreshold=-0.59; %log-frequency log-power slope threshold for muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or more stringent = -0.72, more negative thresholds remove more muscle.
RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=0.50;  % Maximum amount of data periods to be marked as muscle artifact for cleaning by the MWF. You want at least a reasonable amount of both clean and artifact templates for effective cleaning.
% I set this reasonably high, because otherwise muscle artifacts could considerably influence the clean mask and add noise into the data
RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=0.05; % If the proportion of epochs showing muscle activity from an electrode is higher than this, the electrode is deleted. 
% Set muscle proportion before deletion to 1 to not delete electrodes based on muscle activity
RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=0.05; % If the proportion of all epochs from a single electrode that are marked as containing extreme artifacts is higher than this, the electrode is deleted

RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted=0.20; % Sets the maximum proportion of electrodes that are allowed to be deleted after PREP's bad electrode deletion step

RELAX_cfg.MWFDelayPeriod=8; % The MWF includes both spatial and temporal information when filtering out artifacts. Longer delays apparently improve performance. 

%% Notes on the above parameter settings:
% OnlyIncludeTaskRelatedEpochs: this means the MWF cleaning templates 
% won't be distracted by potential large artifacts outside of task related
% periods. However, it also may mean that event related periods are
% included in the artifact template (to be removed). If artifact masks land
% disproportionately around triggers (or certain types of triggers), my
% sense is that all the EEG activity time locked to that trigger will be
% considered an artifact and filtered out by the MWF (including ERPs),
% potentially reducing your ERP for example. An alternative solution to
% this would be to record a sufficient amount of activity before / after
% the task related activity, which includes all the same artifacts as
% within the task, but also a sufficient amount of clean data. Then use
% this data to create the artifact masks (then perhaps marking the task 
% related data as only either clean, or NaNs if they contain artifacts
% (which will be cleaned but aren't thought of as artifacts for the
% artifact template). This would prevent the ERP activity of interest from
% being included in the artifact tempalte (and potentially cleaned)

% Delay periods >5 can lead to generalised eigenvector rank deficiency
% in some files, and if this occurs cleaning is ineffective. Delay
% period = 5 was used by Somers et al (2018). The rank deficiency is
% likely to be because data filtering creates a temporal
% dependency between consecutive datapoints, reducing their
% independence when including the temporal aspect in the MWF
% computation. To address this, the MWF function attempts MWF cleaning
% at the delay period set above, but if rank deficiency occurs,
% it reduces the delay period by 1 and try again (for 3
% iterations).

% Using robust detrending (which does not create any temporal
% dependence,unlike filtering) may be an alternative which avoids rank
% deficiency (but our initial test suggested this led to worse cleaning
% than filtering)

% MuscleSlopeThreshold: (-0.31 = data from paralysed participants showed no
% independent components with a slope value more positive than this (so
% excluding slopes above this threshold means only excluding data that we
% know must be influenced by EMG). Using -0.31 as the threshold means
% possibly leaving low level EMG data in, and only eliminating the data we
% know is definitely EMG)
% (-0.59 is where the histograms between paralysed ICs and EMG ICs cross,
% so above this value contains a very small amount of the brain data, 
% and over 50% of the EMG data. Above this point, data is more likely to be
% EMG than brain)
% (-0.72 is the maximum of the histogram of the paralysed IC data, so
% excluding more positive values than this will exclude most of the EMG
% data, but also some brain data).

% Fitzgibbon, S. P., DeLosAngeles, D., Lewis, T. W., Powers, D. M. W., Grummett, T. S., Whitham, E. M., ... & Pope, K. J. (2016). Automatic determination of EMG-contaminated components and validation of independent component analysis using EEG during pharmacologic paralysis. Clinical Neurophysiology, 127(3), 1781-1793.

% Find the files to process
RELAX_cfg.FilesToProcess=1:numel(RELAX_cfg.files); % Set which files to process

%% GENERATE THE GUI

% Callback command for getting data folder
commandload1 = [ '[filepath] = uigetdir(''*'', ''Select a text file'');' ...
    'if filepath ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''datadir''), ''string'', filepath );' ...
    'end;' ...
    'clear filepath tagtest;' ];

% ICA options
icaOptions = {'runica','cudaica','fastica_symm','fastica_defl','amica'};
ProbabilityOfBlinksOptions = {'data almost certainly has blinks', 'data might not have blinks', 'data definitely does not have blinks'};

% GUI layout
geometry = {[3.5 5 1]...% title
            1 ...
            [0.58 1.5 0.2 1.78] ... % setting directory
            1 ...
            [0.5 1 2] ... % electrodes to exclude
            1 ...
            [5 .8 .5 2 0.5 1.6] ... % bandpass filter settings
            1 ...
            [1 0.2 0.2 1.2 0.2 0.2 1 0.2] ... % electrode rejection settings
            1 ...
            [1] ... % extreme period rejection setting title
            1 ...
            [0.7 0.2 0.3 1.1 0.2 0.3 0.7 0.2 0.3 0.7 0.2 0.8] ... % extreme period rejection settings
            1 ...
            [0.7 0.2 1 0.2 1.1 0.2 0.6]... % extreme period rejection settings MAD blinks and drift
            1 ...
            [0.2 0.09 0.2 0.2 0.12 0.12 0.18 0.6] ... %Wiener filter settings
            1 ...
            [0.8 0.25 0.3 0.6 0.6 0.9 1 1.6] ... % wICA settings and Blink probability options
            1 ...
            [1] ... % eye movement title
            1 ...
            [0.7 1 0.95 0.9 1 0.9]... % eye movement affected electrodes
            1 ...
            [1.3 0.25 0.3 1.9 0.15]... % muscle marking for MWF
            1 ...
            [1.7 0.2 0.05 2.2 0.2]... % drift marking for MWF
            1 ...
            [2 0.3 4]... % horizontal eye movement marking for MWF
            1 ...
            [0.6 0.3 0.4 1.6]... % computing metrics?
            1 ...
            [1.2 1.7 1 1 1 0.4]... % saving intermediate files?
            1 ...
            [3 0.5 3]... % saving intermediate files?
            };
            
% GUI settings with defaults:        
uilist = {{'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'RELAX Cleaning Pipeline - Parameter Setting','fontweight','bold','fontsize', 12.5} ...
          {'style', 'text', 'string', ' '} ...
          {}...
          {'style', 'text', 'string', 'Folder with data:','fontweight','bold'} ...
          { 'style' 'edit'       'string' RELAX_cfg.myPath 'tag' 'datadir' } ... 
          { 'style' 'pushbutton' 'string' '...' 'callback' commandload1 }... 
          {'style', 'text', 'string', ' '} ...
          {}...
          {'style', 'text', 'string', 'Electrodes to exclude:','fontweight','bold'} ...
          {'style', 'edit', 'string', ''} ...
          {'style', 'text', 'string', ' '} ...
          {}...
          {'style', 'text', 'string', 'Bandpass Filter (Hz) [highpass, lowpass] (0.25Hz highpass is best for ERPs, 1Hz is better for cleaning artifacts):'} ...
          {'style', 'edit', 'string', [num2str(RELAX_cfg.HighPassFilter), '  ', num2str(RELAX_cfg.LowPassFilter)]}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Line Noise Frequency (Hz) [eg. 50 or 60]:'} ...
          {'style', 'edit', 'string', RELAX_cfg.LineNoiseFrequency}...
          {'style', 'text', 'string', ' '} ...
          {}...
          {'style', 'text', 'string', 'Max Proportion of electrodes that can be deleted as bad:'} ...
          {'style', 'edit', 'string', RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Delete electrode if it shows > proportion of extreme noise than:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Delete electrode if it shows > proportion of muscle than:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel}...
          {}...
          {'style', 'text', 'string', 'Extreme outlier detection thresholds, applied to each 1s period:','fontweight','bold'} ...
          {}...
          {'style', 'text', 'string', 'Absolute voltage shift:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeAbsoluteVoltageThreshold}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Improbable voltage value distributions (SD):'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Single channel kurtosis:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeSingleChannelKurtosisThreshold}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'All channel kurtosis:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeAllChannelKurtosisThreshold}...
          {'style', 'text', 'string', ' '} ...
          {}...
          {'style', 'text', 'string', 'MAD from median voltage shift:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeVoltageShiftThreshold}...
          {'style', 'text', 'string', 'MAD from median voltage shift in blink affected epochs:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeBlinkShiftThreshold}...
          {'style', 'text', 'string', 'Log-power log-frequency slope for detecting drift threshold:'} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeDriftSlopeThreshold}...
          {'style', 'text', 'string', ' '} ...
          {}...  
          {'style', 'text', 'string', 'MWF Delay Period:','fontweight','bold'} ...
          {'style', 'edit', 'string', RELAX_cfg.MWFDelayPeriod}...
          {'style', 'text', 'string', ''} ...
          {'style', 'text', 'string', 'Use MWF to clean:','fontweight','bold'} ... 
          {'Style', 'checkbox', 'string' 'Muscle' 'value' RELAX_cfg.Do_MWF_Once 'tag' 'once' } ... % Input 1
          {'Style', 'checkbox', 'string' 'Blinks' 'value' RELAX_cfg.Do_MWF_Twice 'tag' 'twice' } ... % Input 2
          {'Style', 'checkbox', 'string' 'HEOG/drift' 'value' RELAX_cfg.Do_MWF_Thrice 'tag' 'thrice' } ... % Input 3; Line 2
          {'style', 'text', 'string', ' '} ...
          {} ...
          {'style', 'text', 'string', 'Clean with wICA?','fontweight','bold'} ...
          {'Style', 'checkbox', 'string' 'yes' 'value' RELAX_cfg.Perform_wICA_on_ICLabel 'tag' 'wICA' } ... % Input 1
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'ICA method:'} ...
          {'style', 'popupmenu', 'string', icaOptions, 'tag', 'icaOpts','Value',3} ...
          {'style', 'text', 'string', 'Does the data contain blinks?','fontweight','bold'} ...
          {'style', 'popupmenu', 'string', ProbabilityOfBlinksOptions, 'tag', 'blinkOpts','Value',1} ...
          {'style', 'text', 'string', ' '} ...
          {} ...
          {'style', 'text', 'string', 'Electrodes affected by eye movements / blinks (for marking artifacts to be cleaned with MWF, list HEOG electrodes in order of use preference):','fontweight','bold'} ...
          {}...
          {'style', 'text', 'string', 'Blink affected electrodes:'} ...
          {'style', 'edit', 'string', [RELAX_cfg.BlinkElectrodes{1,1} ' ' RELAX_cfg.BlinkElectrodes{2,1} ' ' RELAX_cfg.BlinkElectrodes{3,1}...
          ' ' RELAX_cfg.BlinkElectrodes{4,1} ' ' RELAX_cfg.BlinkElectrodes{5,1} ' ' RELAX_cfg.BlinkElectrodes{6,1} ...
          ' ' RELAX_cfg.BlinkElectrodes{7,1} ' ' RELAX_cfg.BlinkElectrodes{8,1} ' ' RELAX_cfg.BlinkElectrodes{9,1} ' ' RELAX_cfg.BlinkElectrodes{10,1}]} ...
          {'style', 'text', 'string', 'Left sided HEOG affected electrodes:'} ...
          {'style', 'edit', 'string', [RELAX_cfg.HEOGLeftpattern{1,1} ' ' RELAX_cfg.HEOGLeftpattern{1,2} ' ' RELAX_cfg.HEOGLeftpattern{1,3}...
          ' ' RELAX_cfg.HEOGLeftpattern{1,4} ' ' RELAX_cfg.HEOGLeftpattern{1,5} ' ' RELAX_cfg.HEOGLeftpattern{1,6} ...
          ' ' RELAX_cfg.HEOGLeftpattern{1,7} ' ' RELAX_cfg.HEOGLeftpattern{1,8} ' ' RELAX_cfg.HEOGLeftpattern{1,9}]} ...
          {'style', 'text', 'string', 'Right sided HEOG affected electrodes:'} ...
          {'style', 'edit', 'string', [RELAX_cfg.HEOGRightpattern{1,1} ' ' RELAX_cfg.HEOGRightpattern{1,2} ' ' RELAX_cfg.HEOGRightpattern{1,3}...
          ' ' RELAX_cfg.HEOGRightpattern{1,4} ' ' RELAX_cfg.HEOGRightpattern{1,5} ' ' RELAX_cfg.HEOGRightpattern{1,6} ...
          ' ' RELAX_cfg.HEOGRightpattern{1,7} ' ' RELAX_cfg.HEOGRightpattern{1,8} ' ' RELAX_cfg.HEOGRightpattern{1,9}]} ...
          {}...
          {'style', 'text', 'string', 'Log-frequency Log-power slope muscle artifact threshold for MWF cleaning:'} ...
          {'style', 'edit', 'string', RELAX_cfg.MuscleSlopeThreshold}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Max proportion marked as muscle for MWF cleaning (only most severe data up to max is marked if data exceeds threshold):'} ...
          {'style', 'edit', 'string', RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle}...
          {}...
          {'style', 'text', 'string', 'Single electrode voltage MAD from median of all electrodes to mark as drift for MWF cleaning:'} ...
          {'style', 'edit', 'string', RELAX_cfg.DriftSeverityThreshold}...
          {'style', 'text', 'string', ' '} ...
          {'style', 'text', 'string', 'Max proportion marked as drift for MWF cleaning (only most severe data up to max is marked if data exceeds threshold):'} ...
          {'style', 'edit', 'string', RELAX_cfg.ProportionWorstEpochsForDrift}...
          {}...
          {'style', 'text', 'string', 'Horizontal eye movement threshold (MAD from median) for MWF cleaning:'} ...
          {'style', 'edit', 'string', RELAX_cfg.HorizontalEyeMovementThreshold}...
          {'style', 'text', 'string', ''}...
          {}...
          {'style', 'text', 'string', 'Compute Artifact/Cleaning Metrics?','fontweight','bold'} ... 
          {'Style', 'checkbox', 'string' 'For raw data' 'value' RELAX_cfg.computerawmetrics 'tag' 'raw' } ... % Input 1
          {'Style', 'checkbox', 'string' 'For cleaned data' 'value' RELAX_cfg.computecleanedmetrics 'tag' 'clean' } ... % Input 2
          {'style', 'text', 'string', ' '}...
          {}...
          {'style', 'text', 'string', 'Save Intermediate Cleaning Steps?','fontweight','bold'} ... 
          {'Style', 'checkbox', 'string' 'After electrode rejections and extreme period marking only' 'value' RELAX_cfg.saveextremesrejected 'tag' 'saveextremesrejected' } ... % Input 1
          {'Style', 'checkbox', 'string' 'After 1st MWF cleaning step' 'value' RELAX_cfg.saveround1 'tag' 'MWF1' } ... % Input 2
          {'Style', 'checkbox', 'string' 'After 2nd MWF cleaning step' 'value' RELAX_cfg.saveround2 'tag' 'MWF2' } ... % Input 2
          {'Style', 'checkbox', 'string' 'After 3rd MWF cleaning step' 'value' RELAX_cfg.saveround3 'tag' 'MWF3' } ... % Input 2
          {'style', 'text', 'string', ' '}...
          {}...
          {'style', 'text', 'string', 'File numbers to process this session (from [start file number] to [finish file number] and all files in between):'} ...
          {'style', 'edit', 'string', [num2str(RELAX_cfg.FilesToProcess(1,1)), '  ', num2str(numel(RELAX_cfg.files))]}...
          {'style', 'text', 'string', ''}...
          };

result = inputgui('geometry', geometry, 'uilist', uilist, 'title', 'RELAX Parameter Setting',  'helpcom', 'pophelp(''pop_RELAX'')');

% Replace default values with inputs
RELAX_cfg.myPath = result{1};
RELAX_cfg.ElectrodesToDelete = strtrim(strsplit(strrep(result{2},',','')))';
BandPassFilterSettings=strsplit(strrep(result{3},',',''));
RELAX_cfg.HighPassFilter=str2double(BandPassFilterSettings(1,1));
RELAX_cfg.LowPassFilter=str2double(BandPassFilterSettings(1,2));
RELAX_cfg.LineNoiseFrequency=str2double(result{4});
RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted=str2double(result{5});
RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=str2double(result{6});
RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=str2double(result{7});
RELAX_cfg.ExtremeAbsoluteVoltageThreshold=str2double(result{8});
RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=str2double(result{9});
RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=str2double(result{10});
RELAX_cfg.ExtremeAllChannelKurtosisThreshold=str2double(result{11});
RELAX_cfg.ExtremeVoltageShiftThreshold=str2double(result{12});
RELAX_cfg.ExtremeBlinkShiftThreshold=str2double(result{13});
RELAX_cfg.ExtremeDriftSlopeThreshold=str2double(result{14});
RELAX_cfg.MWFDelayPeriod=str2double(result{15});
RELAX_cfg.Do_MWF_Once=result{16}; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Twice=result{17}; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Thrice=result{18};
RELAX_cfg.Perform_wICA_on_ICLabel=result{19};
RELAX_cfg.ICA_method = icaOptions{result{20}};
RELAX_cfg.ProbabilityDataHasNoBlinks=(result{21})-1;
RELAX_cfg.BlinkElectrodes= strtrim(strsplit(strrep(result{22},',','')))';
RELAX_cfg.HEOGLeftpattern = string(strtrim(strsplit(strrep(result{23},',',''))));
RELAX_cfg.HEOGRightpattern= string(strtrim(strsplit(strrep(result{24},',',''))));
RELAX_cfg.MuscleSlopeThreshold=str2double(result{25});
RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=str2double(result{26});
RELAX_cfg.DriftSeverityThreshold=str2double(result{27});
RELAX_cfg.ProportionWorstEpochsForDrift=str2double(result{28});
RELAX_cfg.HorizontalEyeMovementThreshold=str2double(result{29});
RELAX_cfg.computerawmetrics=result{30};
RELAX_cfg.computecleanedmetrics=result{31};
RELAX_cfg.saveextremesrejected=result{32};
RELAX_cfg.saveround1=result{33};
RELAX_cfg.saveround2=result{34};
RELAX_cfg.saveround3=result{35};
FilesToProcess=strsplit(strrep(result{36},',',''));
RELAX_cfg.FilesToProcess=str2double(FilesToProcess(1,1)):str2double(FilesToProcess(1,2));

[RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingRoundOneAllParticipants, RELAXProcessingRoundTwoAllParticipants, RELAXProcessing_wICA_AllParticipants,...
        RELAXProcessingRoundThreeAllParticipants, RELAX_issues_to_check, RELAXProcessingExtremeRejectionsAllParticipants] = RELAX_Wrapper (RELAX_cfg);

