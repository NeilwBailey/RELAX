%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2022) Neil Bailey - visit https://github.com/NeilwBailey/RELAX/wiki for full instructions

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
        RELAXProcessing_ICA_AllParticipants, RELAXProcessingRoundThreeAllParticipants, RELAX_issues_to_check, RELAXProcessingExtremeRejectionsAllParticipants] = pop_RELAX_beta(RELAX_cfg)

%% DEPENDENCIES (toolboxes you need to install, and cite if you use this script):
% use fileseparators 'filesep' for increased compatability if necessary (replace the \ with a filesep [outside of quotes]) 

% MATLAB signal processing toolbox (from MATLAB website)
% MATLAB statistics and machine learning toolbox (from MATLAB website)

% PREP pipeline to reject bad electrodes (install plugin to EEGLAB, or via the github into the EEGLAB plugin folder): 
% Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., & Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for large-scale EEG analysis. Frontiers in neuroinformatics, 9, 16.
% http://vislab.github.io/EEG-Clean-Tools/ or install via EEGLAB extensions

% Specify the MWF path:
% https://github.com/exporl/mwf-artifact-removal
% Somers, B., Francart, T., & Bertrand, A. (2018). A generic EEG artifact removal algorithm based on the multi-channel Wiener filter. Journal of neural engineering, 15(3), 036007.

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.

% fastica:
% http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml 

% ICLabel in your eeglab folder as a plugin or via the github:
% https://github.com/sccn/ICLabel

% Check if input is provided:
if nargin == 0
    RELAX_cfg = [];
end

%% Location of files for processing

% Specify your electrode locations with the correct cap file:
if ~isfield(RELAX_cfg,'caploc')
    RELAX_cfg.caploc='D:\Cap_Location_Files\standard-10-5-cap385.elp'; %path containing electrode positions
end

% Specify the to be processed file locations:
if ~isfield(RELAX_cfg,'myPath')
    RELAX_cfg.myPath='D:\DATA_TO_BE_PREPROCESSED\';
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
if ~isfield(RELAX_cfg,'Do_MWF_Once')
    RELAX_cfg.Do_MWF_Once=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
end
if ~isfield(RELAX_cfg,'Do_MWF_Twice')
    RELAX_cfg.Do_MWF_Twice=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
end
if ~isfield(RELAX_cfg,'Do_MWF_Thrice')
    RELAX_cfg.Do_MWF_Thrice=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no). I think cleaning drift in this is a good idea.
end
if ~isfield(RELAX_cfg,'Perform_wICA_on_ICLabel')
    RELAX_cfg.Perform_wICA_on_ICLabel=1; % 1 = Perform wICA on artifact components marked by ICLabel (1 for yes, 0 for no).
end
if ~isfield(RELAX_cfg,'Perform_ICA_subtract')
    RELAX_cfg.Perform_ICA_subtract=0; % 1 = Perform ICA subtract on artifact components marked by ICLabel (1 for yes, 0 for no) (non-optimal, intended to be optionally used separately to wICA rather than additionally)
end
if ~isfield(RELAX_cfg,'ICA_method')
    RELAX_cfg.ICA_method='fastica_symm';
end
if ~isfield(RELAX_cfg,'Report_all_ICA_info')
    RELAX_cfg.Report_all_ICA_info='no'; % set to yes to provide detailed report of ICLabel artifact information. Runs ~20s slower per file.
end
if ~isfield(RELAX_cfg,'computerawmetrics')
    RELAX_cfg.computerawmetrics=1; % Compute blink and muscle metrics from the raw data?
end
if ~isfield(RELAX_cfg,'computecleanedmetrics')
    RELAX_cfg.computecleanedmetrics=1; % Compute SER, ARR, blink and muscle metrics from the cleaned data?
end
if ~isfield(RELAX_cfg,'MWFRoundToCleanBlinks')
    RELAX_cfg.MWFRoundToCleanBlinks=2; % Which round to clean blinks in (1 for the first, 2 for the second...)
end
if ~isfield(RELAX_cfg,'LowPassFilterAt_6Hz_BeforeDetectingBlinks')
    RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks='no'; % low pass filters the data @ 6Hz prior to blink detection (helps if high power alpha is disrupting blink detection, not necessary in the vast majority of cases, default = 'no')
end
if ~isfield(RELAX_cfg,'ProbabilityDataHasNoBlinks')
    RELAX_cfg.ProbabilityDataHasNoBlinks=0; % 0 = data almost certainly has blinks, 1 = data might not have blinks, 2 = data definitely doesn't have blinks.
end
% 0 = eg. task related data where participants are focused with eyes open,
% 1 = eg. eyes closed recordings, but with participants who might still open their eyes at times, 
% 2 = eg. eyes closed resting with highly compliant participants and recordings that were strictly made only when participants had their eyes closed.
if ~isfield(RELAX_cfg,'DriftSeverityThreshold')
    RELAX_cfg.DriftSeverityThreshold=10; %MAD from the median of all electrodes. This could be set lower and would catch less severe drift 
end
if ~isfield(RELAX_cfg,'ProportionWorstEpochsForDrift')
    RELAX_cfg.ProportionWorstEpochsForDrift=0.30; % Maximum proportion of epochs to include in the mask from drift artifact type.
end
if ~isfield(RELAX_cfg,'ExtremeVoltageShiftThreshold')
    RELAX_cfg.ExtremeVoltageShiftThreshold=20; % Threshold MAD from the median all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe voltage shifts within the epoch
end
if ~isfield(RELAX_cfg,'ExtremeAbsoluteVoltageThreshold')
    RELAX_cfg.ExtremeAbsoluteVoltageThreshold=500; % microvolts max or min above which will be excluded from cleaning and deleted from data
end
if ~isfield(RELAX_cfg,'ExtremeImprobableVoltageDistributionThreshold')
    RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; % Threshold SD from the mean of all epochs for each electrode against the same electrode in different epochs. This could be set lower and would catch less severe improbable data
end
if ~isfield(RELAX_cfg,'ExtremeSingleChannelKurtosisThreshold')
    RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=8; % Threshold kurtosis of each electrode against the same electrode in different epochs. This could be set lower and would catch less severe kurtosis 
end
if ~isfield(RELAX_cfg,'ExtremeAllChannelKurtosisThreshold')
    RELAX_cfg.ExtremeAllChannelKurtosisThreshold=8; % Threshold kurtosis across all electrodes. This could be set lower and would catch less severe kurtosis
end
if ~isfield(RELAX_cfg,'ExtremeDriftSlopeThreshold')
    RELAX_cfg.ExtremeDriftSlopeThreshold=-4; % slope of log frequency log power below which to reject as drift without neural activity
end
if ~isfield(RELAX_cfg,'ExtremeBlinkShiftThreshold')
    RELAX_cfg.ExtremeBlinkShiftThreshold=8; % How many MAD from the median across blink affected epochs to exclude as extreme data 
end
% (applies the higher value out of this value and the
% RELAX_cfg.ExtremeVoltageShiftThreshold above as the
% threshold, which caters for the fact that blinks don't affect
% the median, so without this, if data is clean and blinks are
% large, blinks can get excluded as extreme outliers)
            
% Clean periods that last for a shorter duration than the following value to be marked as artifacts, 
% and pad short artifact periods out into artifact periods of at least the following length when 
% they are shorter than this value to reduce rank deficiency issues in MWF cleaning). 
% Note that it's better to include clean periods in the artifact mask rather than the including artifact in the clean mask.
if ~isfield(RELAX_cfg,'MinimumArtifactDuration')
    RELAX_cfg.MinimumArtifactDuration=1200; % in ms. It's better to make this value longer than 1000ms, as doing so will catch diminishing artifacts that aren't detected in a neighbouring 1000ms period, which might still be bad
end
if ~isfield(RELAX_cfg,'MinimumBlinkArtifactDuration')
    RELAX_cfg.MinimumBlinkArtifactDuration=800; % blink marking is based on the maximum point of the blink rather than the 1000ms divisions for muscle artifacts, so this can be shorter than the value above (blinks do not typically last >500ms)
end
if ~isfield(RELAX_cfg,'BlinkElectrodes')
    RELAX_cfg.BlinkElectrodes={'FP1';'FPZ';'FP2';'AF3';'AF4';'F3';'F1';'FZ';'F2';'F4'}; % sets the electrodes to average for blink detection using the IQR method. These should be frontal electrodes maximally affected by blinks. The order is the order of preference for icablinkmetrics.
end
% A single HOEG electrode for each side is selected by the script, prioritized in the following order (if the electrode in position 1 isn't present, the script will check for electrode in position 2, and so on...).
if ~isfield(RELAX_cfg,'HEOGLeftpattern')
    RELAX_cfg.HEOGLeftpattern = ["AF7", "F7", "FT7", "F5", "T7", "FC5", "C5", "TP7", "AF3"]; % sets left side electrodes to use for horizontal eye movement detection. These should be lateral electrodes maximally effected by blinks.
end
if ~isfield(RELAX_cfg,'HEOGRightpattern')
    RELAX_cfg.HEOGRightpattern = ["AF8", "F8","FT8","F6","T8", "FC6", "C6", "TP8", "AF4"]; % sets right side electrodes to use for horizontal eye movement detection. These should be lateral electrodes maximally effected by blinks.
end
if ~isfield(RELAX_cfg,'BlinkMaskFocus')
    RELAX_cfg.BlinkMaskFocus=150; % this value decides how much data before and after the right and left base of the eye blink to mark as part of the blink artifact window. 
end
% I found 100ms on either side of the blink bases works best with a delay of 7 on the MWF. However, it also seemed to create too short artifact masks at times, which may lead to insufficient rank for MWF, so I left the default as 150ms.
if ~isfield(RELAX_cfg,'HorizontalEyeMovementType')
    RELAX_cfg.HorizontalEyeMovementType=2; % 1 to use the IQR method, 2 to use the MAD method for identifying threshold. IQR method less effective for smaller sample sizes (shorter files).
end
if ~isfield(RELAX_cfg,'HorizontalEyeMovementThreshold')
    RELAX_cfg.HorizontalEyeMovementThreshold=2; % MAD deviation from the median that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below).
end
if ~isfield(RELAX_cfg,'HorizontalEyeMovementThresholdIQR')
    RELAX_cfg.HorizontalEyeMovementThresholdIQR=1.5; % If IQR method set above, IQR deviation that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below).
end
if ~isfield(RELAX_cfg,'HorizontalEyeMovementTimepointsExceedingThreshold')
RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold=25; % The number of timepoints (ms) that exceed the horizontal eye movement threshold within the test period (set below) before the period is marked as horizontal eye movement.
end
if ~isfield(RELAX_cfg,'HorizontalEyeMovementTimepointsTestWindow')
    RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow=(2*RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold)-1; % Window duration to test for horizontal eye movement, set to 2x the value above by default.
end
if ~isfield(RELAX_cfg,'HorizontalEyeMovementFocus')
    RELAX_cfg.HorizontalEyeMovementFocus=200; % Buffer window, masking periods earlier and later than the time where horizontal eye movements exceed the threshold.
end
if ~isfield(RELAX_cfg,'LowPassFilterBeforeMWF')
    RELAX_cfg.LowPassFilterBeforeMWF='no'; % set as no for the updated implementation, avoiding low pass filtering prior to MWF reduces chances of rank deficiencies, increasing potential values for MWF delay period 
end
if ~isfield(RELAX_cfg,'FilterType')
    RELAX_cfg.FilterType='Butterworth'; % set as 'pop_eegfiltnew' to use EEGLAB's filter or 'Butterworth' to use Butterworth filter
end
if ~isfield(RELAX_cfg,'NotchFilterType')
    RELAX_cfg.NotchFilterType='Butterworth'; % set as 'Butterworth' to use Butterworth filter or 'ZaplinePlus' to use ZaplinePlus. ZaplinePlus works best on data sampled at 512Hz or below, consider downsampling if above this.
end
if ~isfield(RELAX_cfg,'HighPassFilter')
    RELAX_cfg.HighPassFilter=0.25; % Sets the high pass filter. 1Hz is best for ICA decomposition if you're examining just oscillatory data, 0.25Hz seems to be the highest before ERPs are adversely affected by filtering 
end
%(lower than 0.2Hz may be better, but I find a minority of my files show drift at 0.3Hz even).
if ~isfield(RELAX_cfg,'LowPassFilter')
    RELAX_cfg.LowPassFilter=80; % If you filter out data below 75Hz, you can't use the objective muscle detection method
end
if ~isfield(RELAX_cfg,'LineNoiseFrequency')
    RELAX_cfg.LineNoiseFrequency=50; % Frequencies for bandstop filter in order to address line noise (set to 60 in countries with 60Hz line noise, and 50 in countries with 50Hz line noise).
end
if ~isfield(RELAX_cfg,'ElectrodesToDelete')
    RELAX_cfg.ElectrodesToDelete={'CB1'; 'CB2'; 'HEOG'; 'IO1'; 'M1'; 'M2'; 'LO1'; 'LO2'; 'E1'; 'E3'; 'ECG'; 'SO1'; 'ECG'; 'SPARE1'; 'SPARE2'; 'SPARE3'; 'BP1'; 'BP2'; 'VEOG'};
end
% If your EEG recording includes non-scalp electrodes or electrodes that you want to delete before cleaning, you can set them to be deleted here. 
% The RELAX cleaning pipeline does not need eye, heart, or mastoid electrodes for effective cleaning.
if ~isfield(RELAX_cfg,'KeepAllInfo')
    RELAX_cfg.KeepAllInfo=0; % setting this value to 1 keeps all the details from the MWF pre-processing and MWF computation. Helpful for debugging if necessary but makes for large file sizes.
end
if ~isfield(RELAX_cfg,'saveextremesrejected')
    RELAX_cfg.saveextremesrejected=0; % setting this value to 1 tells the script to save the data after only filtering, extreme channels have been rejected and extreme periods have been noted
end
if ~isfield(RELAX_cfg,'saveround1')
    RELAX_cfg.saveround1=0; % setting this value to 1 tells the script to save the first round of MWF pre-processing
end
if ~isfield(RELAX_cfg,'saveround2')
    RELAX_cfg.saveround2=0; % setting this value to 1 tells the script to save the second round of MWF pre-processing
end
if ~isfield(RELAX_cfg,'saveround3')
    RELAX_cfg.saveround3=0; % setting this value to 1 tells the script to save the third round of MWF pre-processing
end
if ~isfield(RELAX_cfg,'OnlyIncludeTaskRelatedEpochs')
    RELAX_cfg.OnlyIncludeTaskRelatedEpochs=0; % If this =1, the MWF clean and artifact templates will only include data within 5 seconds of a task trigger (other periods will be marked as NaN, which the MWF script ignores).
end
if ~isfield(RELAX_cfg,'MuscleSlopeThreshold')
    RELAX_cfg.MuscleSlopeThreshold=-0.59; %log-frequency log-power slope threshold for muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or more stringent = -0.72, more negative thresholds remove more muscle.
end
if ~isfield(RELAX_cfg,'MaxProportionOfDataCanBeMarkedAsMuscle')
    RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=0.50;  % Maximum amount of data periods to be marked as muscle artifact for cleaning by the MWF. You want at least a reasonable amount of both clean and artifact templates for effective cleaning.
end
% I set this reasonably high, because otherwise muscle artifacts could considerably influence the clean mask and add noise into the data
if ~isfield(RELAX_cfg,'ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel')
    RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=0.05; % If the proportion of epochs showing muscle activity from an electrode is higher than this, the electrode is deleted. 
end
% Set muscle proportion before deletion to 1 to not delete electrodes based on muscle activity
if ~isfield(RELAX_cfg,'ProportionOfExtremeNoiseAboveWhichToRejectChannel')
    RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=0.05; % If the proportion of all epochs from a single electrode that are marked as containing extreme artifacts is higher than this, the electrode is deleted
end
if ~isfield(RELAX_cfg,'MaxProportionOfElectrodesThatCanBeDeleted')
    RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted=0.20; % Sets the maximum proportion of electrodes that are allowed to be deleted after PREP's bad electrode deletion step
end
if ~isfield(RELAX_cfg,'InterpolateRejectedElectrodesAfterCleaning')
    RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning='no'; % Interpolate rejected electrodes back into the data after each file has been cleaned and before saving the cleaned data?
end
if ~isfield(RELAX_cfg,'MWFDelayPeriod')
    RELAX_cfg.MWFDelayPeriod=16; % The MWF includes both spatial and temporal information when filtering out artifacts. Longer delays apparently improve performance. 
end
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

%% GENERATE THE GUI

% Only launch GUI if cfg stucture not provided
if nargin == 0

% Callback command for getting data folder
commandload1 = [ '[filepath] = uigetdir(''*'', ''Select the folder containing your data'');' ...
    'if filepath ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''datadir''), ''string'', filepath );' ...
    'end;' ...
    'clear filepath tagtest;' ];

% Callback command for getting cap location file:
commandload2 = [ '[capfile cap_path] = uigetfile(''*'', ''Select the cap location file appropriate for your data'');' ...
    'if capfile ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''caploc''), ''string'', strcat(cap_path,capfile) );' ...
    'end;' ...
    'clear capfile tagtest;' ];

% Callback command for getting cap location file:
commandload3 = [ '[filetoprocess file_path] = uigetfile(''*.set'', ''Select the file you want to process (leave blank to process a range)'');' ...
    'if filetoprocess ~=0,' ...
    '   set(findobj(''parent'', gcbf, ''tag'', ''filetoprocess''), ''string'', strcat(file_path,filetoprocess) );' ...
    'end;' ...
    'clear filetoprocess tagtest;' ];

% ICA options
wICA_or_ICA={'Reduce ICA artifacts with wICA','Subtract ICA artifacts','No ICA artifact reduction'};
icaOptions = {'extended_infomax_ICA','cudaica','fastica_symm','fastica_defl','amica','picard'};
ProbabilityOfBlinksOptions = {'data almost certainly has blinks', 'data might not have blinks', 'data definitely does not have blinks'};
YesNoOptions={'no','yes'};
YesNoOptionsLowPassFilterForBlinks={'no - default (works in the vast majority of cases)','yes - apply if RELAX is not detecting blinks because of high power alpha'};
BandPassFilterOptions={'Butterworth','pop_eegfiltnew'};
LineNoiseOptions={'Butterworth','ZaplinePlus'};

% GUI layout
geometry = {[0.5 1.0 0.2 1.0 1.0 0.2] ... % setting directory / file to process
            1 ...
            [0.6 1.5 0.2 0.75 1.2] ... % setting cap location file and electrodes to exclude
            1 ...
            [4 1.1 3 0.7 2.75 .75 2 1.2 2 1.2] ... % bandpass filter settings
            1 ...
            [1.4 0.2 0.2 1.4 0.2 0.2 1.3 0.2] ... % electrode rejection settings
            1 ...
            [1] ...
            1 ...
            [0.7 0.3 1.3 0.3 0.8 0.3 0.7 0.3] ... % extreme period rejection settings
            1 ...
            [0.9 0.2 1.3 0.2 1.3 0.2]... % extreme period rejection settings MAD blinks and drift
            1 ...
            [0.6 0.2 0.6 0.4 0.4 0.6 1]... %Wiener filter settings 
            1 ...
            [0.6 0.6 0.4 0.8]... % wICA settings 
            1 ...
            [1.3 1 0.3 0.9 1.5 1.6 1.5]... % Blink probability and electrode options
            1 ...
            [0.95 0.9 1 0.9]... % eye movement affected electrodes
            1 ...
            [2.5 0.3 0.3 1.9 0.3]... % muscle marking for MWF
            1 ...
            [1.7 0.2 0.05 2.2 0.2]... % drift marking for MWF
            1 ...
            [2.2 0.3 1 2.2 0.5]... % horizontal eye movement marking for MWF
            1 ...
            [2 0.4 0.8 0.7 0.8]... % computing metrics? 
            1 ...
            [1.1 1.9 1.1 1.1 1.1]... % saving intermediate files?
            };
            
% GUI settings with defaults:        
uilist = {{'style', 'text', 'string', 'Raw data folder:','fontweight','bold','fontsize', 9} ...
          { 'style' 'edit'       'string' RELAX_cfg.myPath 'tag' 'datadir' ,'fontsize', 9} ... 
          { 'style' 'pushbutton' 'string' '...' 'callback' commandload1 ,'fontsize', 9}... 
          {'style', 'text', 'string', 'File to clean (leave blank to clean >1):','fontweight','bold','fontsize', 9} ...
          { 'style' 'edit'       'string' '' 'tag' 'filetoprocess' ,'fontsize', 9} ... 
          { 'style' 'pushbutton' 'string' '...' 'callback' commandload3 ,'fontsize', 9}... 
          {}...
          {'style', 'text', 'string', 'Cap Location File:','fontweight','bold','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.caploc 'tag' 'caploc' ,'fontsize', 9} ...
          { 'style' 'pushbutton' 'string' '...' 'callback' commandload2 ,'fontsize', 9}... 
          {'style', 'text', 'string', 'Electrodes to exclude:','fontweight','bold','fontsize', 9} ...
          {'style', 'edit', 'string', ''} ...
          {}...
          {'style', 'text', 'string', 'Bandpass Filter [highpass, lowpass]:','fontsize', 9} ...
          {'style', 'edit', 'string', [num2str(RELAX_cfg.HighPassFilter), '  ', num2str(RELAX_cfg.LowPassFilter)],'fontsize', 9}...
          {'style', 'text', 'string', 'Line Noise Frequency [eg. 50 or 60]:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.LineNoiseFrequency,'fontsize', 9}...
          {'style', 'text', 'string', 'Lowpass filter before MWF?','fontweight','bold','fontsize', 9} ...
          {'style', 'popupmenu', 'string', YesNoOptions, 'tag', 'YesNoOpts','Value',1,'fontsize', 9} ...
          {'style', 'text', 'string', 'Bandpass Filter Type:','fontweight','bold','fontsize', 9} ...
          {'style', 'popupmenu', 'string', BandPassFilterOptions, 'tag', 'BandPassFilterOpts','Value',1,'fontsize', 9} ...
          {'style', 'text', 'string', 'Clean Line Noise With:','fontweight','bold','fontsize', 9} ...
          {'style', 'popupmenu', 'string', LineNoiseOptions, 'tag', 'linenoiseOpts','Value',1,'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Max Proportion of electrodes that can be deleted as bad:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted,'fontsize', 9}...
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {'style', 'text', 'string', 'Extreme noise proportion electrode deletion threshold:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel,'fontsize', 9}...
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {'style', 'text', 'string', 'Muscle noise proportion electrode deletion threshold:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'Extreme outlier detection thresholds, applied to each 1s period:','fontweight','bold','fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Absolute voltage shift:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeAbsoluteVoltageThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'Improbable voltage value distributions (SD):','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'Single channel kurtosis:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeSingleChannelKurtosisThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'All channel kurtosis:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeAllChannelKurtosisThreshold,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'MAD from median voltage shift:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeVoltageShiftThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'MAD from median voltage shift in blink affected epochs:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeBlinkShiftThreshold,'fontsize', 9}...
          {'style', 'text', 'string', 'Log-frequency Log-power slope for detecting drift threshold:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ExtremeDriftSlopeThreshold,'fontsize', 9}...
          {}...  
          {'style', 'text', 'string', 'MWF Delay Period:','fontweight','bold','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.MWFDelayPeriod,'fontsize', 9}...
          {'style', 'text', 'string', 'Use MWF to clean:','fontweight','bold','fontsize', 9} ... 
          {'Style', 'checkbox', 'string' 'Muscle' 'value' RELAX_cfg.Do_MWF_Once 'tag' 'once' ,'fontsize', 9} ... % Input 1
          {'Style', 'checkbox', 'string' 'Blinks' 'value' RELAX_cfg.Do_MWF_Twice 'tag' 'twice' ,'fontsize', 9} ... % Input 2
          {'Style', 'checkbox', 'string' 'HEOG/drift' 'value' RELAX_cfg.Do_MWF_Thrice 'tag' 'thrice' ,'fontsize', 9} ... % Input 3; Line 2
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Clean Artifacts with ICA?','fontsize', 9} ...
          {'style', 'popupmenu', 'string', wICA_or_ICA, 'tag', 'wICA_or_ICAOpts','Value',1,'fontsize', 9} ...
          {'style', 'text', 'string', 'ICA method:','fontsize', 9} ...
          {'style', 'popupmenu', 'string', icaOptions, 'tag', 'icaOpts','Value',3,'fontsize', 9} ...
          {} ...
          {'style', 'text', 'string', 'Does the data contain blinks?','fontweight','bold','fontsize', 9} ...
          {'style', 'popupmenu', 'string', ProbabilityOfBlinksOptions, 'tag', 'blinkOpts','Value',1,'fontsize', 9} ...
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {'style', 'text', 'string', 'Blink affected electrodes:','fontsize', 9} ...
          {'style', 'edit', 'string', [RELAX_cfg.BlinkElectrodes{1,1} ' ' RELAX_cfg.BlinkElectrodes{2,1} ' ' RELAX_cfg.BlinkElectrodes{3,1}...
          ' ' RELAX_cfg.BlinkElectrodes{4,1} ' ' RELAX_cfg.BlinkElectrodes{5,1} ' ' RELAX_cfg.BlinkElectrodes{6,1} ...
          ' ' RELAX_cfg.BlinkElectrodes{7,1} ' ' RELAX_cfg.BlinkElectrodes{8,1} ' ' RELAX_cfg.BlinkElectrodes{9,1} ' ' RELAX_cfg.BlinkElectrodes{10,1}],'fontsize', 9} ...
          {'style', 'text', 'string', '6Hz low pass filter before blink detection?','fontweight','bold','fontsize', 9} ...
          {'style', 'popupmenu', 'string', YesNoOptionsLowPassFilterForBlinks, 'tag', 'YesNoOpts','Value',1,'fontsize', 9} ...
          {} ...
          {'style', 'text', 'string', 'Left sided HEOG affected electrodes:','fontsize', 9} ...
          {'style', 'edit', 'string', [RELAX_cfg.HEOGLeftpattern{1,1} ' ' RELAX_cfg.HEOGLeftpattern{1,2} ' ' RELAX_cfg.HEOGLeftpattern{1,3}...
          ' ' RELAX_cfg.HEOGLeftpattern{1,4} ' ' RELAX_cfg.HEOGLeftpattern{1,5} ' ' RELAX_cfg.HEOGLeftpattern{1,6} ...
          ' ' RELAX_cfg.HEOGLeftpattern{1,7} ' ' RELAX_cfg.HEOGLeftpattern{1,8} ' ' RELAX_cfg.HEOGLeftpattern{1,9}],'fontsize', 9} ...
          {'style', 'text', 'string', 'Right sided HEOG affected electrodes:','fontsize', 9} ...
          {'style', 'edit', 'string', [RELAX_cfg.HEOGRightpattern{1,1} ' ' RELAX_cfg.HEOGRightpattern{1,2} ' ' RELAX_cfg.HEOGRightpattern{1,3}...
          ' ' RELAX_cfg.HEOGRightpattern{1,4} ' ' RELAX_cfg.HEOGRightpattern{1,5} ' ' RELAX_cfg.HEOGRightpattern{1,6} ...
          ' ' RELAX_cfg.HEOGRightpattern{1,7} ' ' RELAX_cfg.HEOGRightpattern{1,8} ' ' RELAX_cfg.HEOGRightpattern{1,9}],'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'Log-frequency Log-power slope muscle artifact threshold for MWF cleaning:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.MuscleSlopeThreshold,'fontsize', 9}...
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {'style', 'text', 'string', 'Max proportion marked as muscle for MWF cleaning:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'Single electrode drift threshold for MWF cleaning:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.DriftSeverityThreshold,'fontsize', 9}...
          {'style', 'text', 'string', ' ','fontsize', 9} ...
          {'style', 'text', 'string', 'Max proportion marked as drift for MWF cleaning:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.ProportionWorstEpochsForDrift,'fontsize', 9}...
          {}...
          {'style', 'text', 'string', 'Horizontal eye movement threshold (MAD from median) for MWF cleaning:','fontsize', 9} ...
          {'style', 'edit', 'string', RELAX_cfg.HorizontalEyeMovementThreshold,'fontsize', 9}...
          {'style', 'text', 'string', '','fontsize', 9}...
          {'style', 'text', 'string', 'Interpolate rejected electrodes back into data after cleaning?','fontsize', 9} ...
          {'style', 'popupmenu', 'string', YesNoOptions, 'tag', 'YesNoOpts','Value',1,'fontsize', 9} ...
          {}...
          {'style', 'text', 'string', 'File numbers to process this session (from start file # to finish file # and files in between):','fontsize', 9} ...
          {'style', 'edit', 'string', ['1', '  ', '2'],'fontsize', 9}...
          {'style', 'text', 'string', 'Compute Metrics?','fontweight','bold','fontsize', 9} ... 
          {'Style', 'checkbox', 'string' 'Raw metrics' 'value' RELAX_cfg.computerawmetrics 'tag' 'raw' ,'fontsize', 9} ... % Input 1
          {'Style', 'checkbox', 'string' 'Cleaned metrics' 'value' RELAX_cfg.computecleanedmetrics 'tag' 'clean' ,'fontsize', 9} ... % Input 2
          {}...
          {'style', 'text', 'string', 'Save Intermediate Steps?','fontweight','bold','fontsize', 9} ... 
          {'Style', 'checkbox', 'string' 'After electrode rejection/extreme period marking' 'value' RELAX_cfg.saveextremesrejected 'tag' 'saveextremesrejected' ,'fontsize', 9} ... % Input 1
          {'Style', 'checkbox', 'string' 'After 1st MWF cleaning' 'value' RELAX_cfg.saveround1 'tag' 'MWF1' ,'fontsize', 9} ... % Input 2
          {'Style', 'checkbox', 'string' 'After 2nd MWF cleaning' 'value' RELAX_cfg.saveround2 'tag' 'MWF2' ,'fontsize', 9} ... % Input 2
          {'Style', 'checkbox', 'string' 'After 3rd MWF cleaning' 'value' RELAX_cfg.saveround3 'tag' 'MWF3' ,'fontsize', 9} ... % Input 2
          };

result = inputgui('geometry', geometry, 'geomvert', [1 .4 1 .4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1 0.4 1],  'uilist', uilist, 'title', 'RELAX Parameter Setting',  'helpcom', 'pophelp(''pop_RELAX_beta'')');

% Replace default values with inputs
RELAX_cfg.myPath = result{1};
RELAX_cfg.filename = result{2};
RELAX_cfg.caploc = result{3};
RELAX_cfg.ElectrodesToDelete = strtrim(strsplit(strrep(result{4},',','')))';
BandPassFilterSettings=strsplit(strrep(result{5},',',''));
RELAX_cfg.HighPassFilter=str2double(BandPassFilterSettings(1,1));
RELAX_cfg.LowPassFilter=str2double(BandPassFilterSettings(1,2));
RELAX_cfg.LineNoiseFrequency=str2double(result{6});
RELAX_cfg.LowPassFilterBeforeMWF=YesNoOptions{result{7}};
RELAX_cfg.FilterType=BandPassFilterOptions{result{8}};
RELAX_cfg.NotchFilterType=LineNoiseOptions{result{9}};
RELAX_cfg.MaxProportionOfElectrodesThatCanBeDeleted=str2double(result{10});
RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=str2double(result{11});
RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=str2double(result{12});
RELAX_cfg.ExtremeAbsoluteVoltageThreshold=str2double(result{13});
RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=str2double(result{14});
RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=str2double(result{15});
RELAX_cfg.ExtremeAllChannelKurtosisThreshold=str2double(result{16});
RELAX_cfg.ExtremeVoltageShiftThreshold=str2double(result{17});
RELAX_cfg.ExtremeBlinkShiftThreshold=str2double(result{18});
RELAX_cfg.ExtremeDriftSlopeThreshold=str2double(result{19});
RELAX_cfg.MWFDelayPeriod=str2double(result{20});
RELAX_cfg.Do_MWF_Once=result{21}; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Twice=result{22}; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Thrice=result{23};
RELAX_cfg.Perform_wICA_on_ICLabel=result{24};
RELAX_cfg.Perform_ICA_subtract=result{24}-1;
RELAX_cfg.ICA_method = icaOptions{result{25}};
RELAX_cfg.ProbabilityDataHasNoBlinks=(result{26})-1;
RELAX_cfg.BlinkElectrodes= strtrim(strsplit(strrep(result{27},',','')))';
RELAX_cfg.LowPassFilterAt_6Hz_BeforeDetectingBlinks=YesNoOptions{result{28}};
RELAX_cfg.HEOGLeftpattern = string(strtrim(strsplit(strrep(result{29},',',''))));
RELAX_cfg.HEOGRightpattern= string(strtrim(strsplit(strrep(result{30},',',''))));
RELAX_cfg.MuscleSlopeThreshold=str2double(result{31});
RELAX_cfg.MaxProportionOfDataCanBeMarkedAsMuscle=str2double(result{32});
RELAX_cfg.DriftSeverityThreshold=str2double(result{33});
RELAX_cfg.ProportionWorstEpochsForDrift=str2double(result{34});
RELAX_cfg.HorizontalEyeMovementThreshold=str2double(result{35});
RELAX_cfg.InterpolateRejectedElectrodesAfterCleaning=YesNoOptions{result{36}};
FilesToProcess=strsplit(strrep(result{37},',',''));
RELAX_cfg.FilesToProcess=str2double(FilesToProcess(1,1)):str2double(FilesToProcess(1,2));
RELAX_cfg.computerawmetrics=result{38};
RELAX_cfg.computecleanedmetrics=result{39};
RELAX_cfg.saveextremesrejected=result{40};
RELAX_cfg.saveround1=result{41};
RELAX_cfg.saveround2=result{42};
RELAX_cfg.saveround3=result{43};

if ~isfield(RELAX_cfg,'filename')
    RELAX_cfg.filename = [];
end

end

%% Check for dependencies:

eeglabPath = fileparts(which('eeglab'));
MWFPluginPath=strcat(eeglabPath,'\plugins\mwf-artifact-removal-master\');
addpath(genpath(MWFPluginPath));
if (exist('mwf_process','file')==0)
    warndlg('MWF toolbox may not be installed in EEGLAB plugins folder. Toolbox can be installed from: "https://github.com/exporl/mwf-artifact-removal"','MWF Cleaning Not Available');
end

toolboxlist=ver;
if isempty(find(strcmp({toolboxlist.Name}, 'Signal Processing Toolbox')==1, 1))
    warndlg('Signal Processing Toolbox may not be installed. Toolbox can be installed through MATLAB "Add-Ons" button','Signal Processing Toolbox not installed');
end
if isempty(find(strcmp({toolboxlist.Name}, 'Wavelet Toolbox')==1, 1))
    warndlg('Wavelet Toolbox may not be installed. Toolbox can be installed through MATLAB "Add-Ons" button','Wavelet Toolbox not installed');
end
if isempty(find(strcmp({toolboxlist.Name}, 'Statistics and Machine Learning Toolbox')==1, 1))
    warndlg('Statistics and Machine Learning Toolbox may not be installed. Toolbox can be installed through MATLAB "Add-Ons" button','Statistics and Machine Learning Toolbox not installed');
end

PluginPath=strcat(eeglabPath,'\plugins\');
PluginList=dir(PluginPath);

if (exist('iclabel','file')==0)
    warndlg('ICLabel may not be installed. Plugin can be installed via EEGLAB: "File" > "Manage EEGLAB Extensions"','ICLabel not installed');
end

if (exist('findNoisyChannels','file')==0)
    warndlg('PrepPipeline may not be installed or the folder path for Prep has not been set in MATLAB. Plugin can be installed via EEGLAB: "File" > "Manage EEGLAB Extensions"','PrepPipeline might not be installed');
end

if (exist('ft_freqanalysis','file')==0)
    warndlg('fieldtrip may not be installed, or the folder path for fieldtrip has not been set in MATLAB. Plugin can be installed via EEGLAB: "File" > "Manage EEGLAB Extensions"','fieldtrip might not be installed');
end

if (strcmp(RELAX_cfg.ICA_method,'fastica_symm')||strcmp(RELAX_cfg.ICA_method,'fastica_defl')) && (exist('fastica','file')==0) 
    warndlg('FastICA may not be installed, or the folder path for FastICA has not been set in MATLAB. Plugin can be installed from: "http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml"','FastICA might not be installed');
end

if (strcmp(RELAX_cfg.ICA_method,'cudaica')) && (exist('cudaica','file')==0)
    warndlg('CudaICA may not be installed. Instructions for installation can be found at: "https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code". Warning - installation can be difficult','CudaICA not installed');
end

if (strcmp(RELAX_cfg.ICA_method,'amica')) && (exist('runamica15','file')==0)
    warndlg('AMICA may not be installed. Plugin can be installed via EEGLAB: "File" > "Manage EEGLAB Extensions"','AMICA not installed');
end

if (strcmp(RELAX_cfg.NotchFilterType,'ZaplinePlus')) && (exist('clean_data_with_zapline_plus_eeglab_wrapper','file')==0)
    warndlg('ZaplinePlus may not be installed. Plugin can be installed via EEGLAB: "File" > "Manage EEGLAB Extensions"','ZaplinePlus not installed');
end

%%

if RELAX_cfg.HighPassFilter>0.25
    Warning='You have high pass filtered above 0.25, which can adversely affect ERP analyses';
end

% List all files in directory
cd(RELAX_cfg.myPath);
RELAX_cfg.dirList=dir('*.set');
RELAX_cfg.files={RELAX_cfg.dirList.name};
if isempty(RELAX_cfg.files)
    disp('No files found..')
end

if ~isfield(RELAX_cfg,'FilesToProcess')
    RELAX_cfg.FilesToProcess=1;
end

[RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingRoundOneAllParticipants, RELAXProcessingRoundTwoAllParticipants, RELAXProcessing_wICA_AllParticipants,...
        RELAXProcessing_ICA_AllParticipants, RELAXProcessingRoundThreeAllParticipants, RELAX_issues_to_check, RELAXProcessingExtremeRejectionsAllParticipants] = RELAX_Wrapper_beta (RELAX_cfg);
   
if RELAX_cfg.HighPassFilter>0.25
    Warning='You have high pass filtered above 0.25, which can adversely affect ERP analyses';
end

% To enable debugging if necessary:

% FileNumber={}; 
% CleanedMetrics={};
% RawMetrics={};
% RELAXProcessingRoundOneAllParticipants={};
% RELAXProcessingRoundTwoAllParticipants={};
% RELAXProcessing_wICA_AllParticipants={};
% RELAXProcessing_ICA_AllParticipants={};
% RELAXProcessingRoundThreeAllParticipants={};
% RELAX_issues_to_check={};
% RELAXProcessingExtremeRejectionsAllParticipants={};
