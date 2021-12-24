clear all; close all; clc;

%% RELAX EEG CLEANING PIPELINE, COPYRIGHT NEIL BAILEY (2022)

%% This script firstly filters your data, then deletes bad electrodes and marks extremely bad periods for exclusion from further processing.
%% It then creates a mask of clean and artifact data specific to the length of your data file, which is submitted to a multiple weiner filter (MWF) in order to clean your data.
%% The MWF cleaned data is then submitted to an independent component analysis, and artifactual components detected by ICLabel are cleaned by wavelet enhanced ICA.
%% The resulting saved files are continuous, very well cleaned of artifacts, while still preserving the neural signal.

%% NOTE THAT THE PIPELINE WILL ONLY WORK ON CONTINUOUS DATA

% RELAX is intended to be fully automated, using the most empirical approaches
% we could find in order to identify periods of EEG data containing
% artifacts, as well as removing channels that are either bad, or contain
% above a specified threshold for the proportion of epochs containing
% muscle artifacts.

% In order to use this script, do the following:

% 1) Install the dependencies listed below, and define the folder where
% they are installed (where required)

% 2) Load the "to be processed" files into a single folder, EEGLAB
% formatted (.set format), and specify this folder location in RELAX_cfg.myPath
% (these files should have triggers recoded to reflect participant accuracy
% if desirable, and not have been re-referenced or have any other
% pre-processing completed yet).

% 3) Specify the parameters in RELAX_cfg. These parameters let you determine
% a range of factors that influence the mask creation, eg. "how many 
% epochs can show muscle activity in a channel before that channel  is
% deleted" (note that a maximum of 20% of channel deletions are allowed 
% after the PREP electrode deletion, "what is the maximum proportion of 
% epochs that can be masked as an artifact for each type of artifact", 
% "should I run the MWF twice or just once", and "what is the threshold 
% of artifact severity,  above which sections are marked as an artifact, 
% below which they are not?". 
% Additionally, particular artifact detection and masking functions
% could be excluded if you're not concerned about those artifacts by
% commenting that function out (adding a % symbol at the start of the
% line). For example, you might only want to include muscle, and
% drift artifacts in the MWF cleaning, and not blink and horizontal 
% eye movements if you are analysing eyes closed resting data.

% 4)Click "Run" in the menu up the top of this script
% The script will then process all participants, taking ~10 min per participant. 
% You can run just a subset of the participants first by altering the
% following line: "for Subjects=1:numel(RELAX_cfg.files)"

% 5) The script will produce files that contain cleaning quality values for
% relevant processes in this script. Means, SD's (or medians and MADs) 
% and ranges across all participants can be examined to determine 
% the success of the script for the study as a whole, and to identify
% potential issues with specific files (eg. the artifact to residue ratio
% is ~1, suggesting artifacts weren't removed successfully, or >75% of the
% epochs contained muscle activity, suggesting bad data throughout). These
% statistics should be reported in publications.

% Note that it's particularly important to check for files that show rank
% deficiency in the MWF processing, as this could prevent the data from
% being cleaned. If this happens, I've found that adjusting the RELAX_cfg
% parameters can fix the problem. Particularly the amount of data that can
% be included in the mask (the mask needs to contain both enough artifact
% and enough clean data to create artifact and clean templates), and the
% minimum length of a marked patch of clean or artifact data (if
% patches are too short, the MWF can't create an adequate mask. The amount
% of data marked as artifact around eye blinks can particularly affect
% this).

% 6) The script will produce output EEG files that have been cleaned into a
% folder labelled 1xMWF for the first run, 2xMWF, and 3xMWF (if
% you select to save these). Cleaned data will be saved into the Cleaned_data 
% folder. The cleaned files will have electrodes removed, and are referenced to the 
% average reference. From this point, I would interpolate missing electrodes, 
% re-code the task related triggers as to whether the participant
% was correct or not, epoch the data, reject epochs that still show artifacts 
% (the epoch rejection wrapper could be adapted to achieve this)
% and separate different conditions into different files in preparation 
% for analysis.

% Note that if you have very large files (or your computer doesn't have much RAM),
% the MWF cleaning may run out of RAM (and suggest a created array is too large).
% This could be addressed by downsampling the data, or reducing 
% RELAX_cfg.MWFDelayPeriod=8; to 4-5. This will reduce the number of samples 
% that the MWF cleaning takes into account when characterising artifacts, 
% but still cleans very effectively.

%% DEPENDENCIES (toolboxes you need to install, and cite if you use this script):
% (use fileseparators 'filesep' for increased compatability if necessary) 

% MATLAB signal processing toolbox (from MATLAB website)
% MATLAB statistics and machine learning toolbox (from MATLAB website)

% EEGLAB:
% https://sccn.ucsd.edu/eeglab/index.php
% Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of neuroscience methods, 134(1), 9-21.
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_current\eeglab2019_1');
eeglab;

% PREP pipeline to reject bad electrodes (install plugin to EEGLAB, or via the github into the EEGLAB plugin folder): 
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
% Bigdely-Shamlo, N., Mullen, T., Kothe, C., Su, K. M., & Robbins, K. A. (2015). The PREP pipeline: standardized preprocessing for large-scale EEG analysis. Frontiers in neuroinformatics, 9, 16.
% http://vislab.github.io/EEG-Clean-Tools/ or install via EEGLAB extensions

% Specify the MWF path:
% https://github.com/exporl/mwf-artifact-removal
% Somers, B., Francart, T., & Bertrand, A. (2018). A generic EEG artifact removal algorithm based on the multi-channel Wiener filter. Journal of neural engineering, 15(3), 036007.
addpath(genpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_current\eeglab2019_1\plugins\mwf-artifact-removal-master'));

% Fieldtrip:
% http://www.fieldtriptoolbox.org/
% Robert Oostenveld, Pascal Fries, Eric Maris, and Jan-Mathijs Schoffelen. FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data. Computational Intelligence and Neuroscience, vol. 2011, Article ID 156869, 9 pages, 2011. doi:10.1155/2011/156869.
addpath('C:\Program Files\MATLAB\fieldtrip-20180805');

% BLINKER to detect eyeblinks (install plugin to EEGLAB):
% http://vislab.github.io/EEG-Blinks/
% Kleifges, K., Bigdely-Shamlo, N., Kerick, S. E., & Robbins, K. A. (2017). BLINKER: Automated extraction of ocular indices from EEG enabling large-scale analysis. Frontiers in neuroscience, 11, 12.
% Sometimes Blinker needs to be installed and unpacked from github rather than via EEGLAB if it doesn't work.

% TESA for filtering (install plugin to EEGLAB):
% https://github.com/nigelrogasch/TESA
% Rogasch, N. C., Sullivan, C., Thomson, R. H., Rose, N. S., Bailey, N. W., Fitzgerald, P. B., ... & Hernandez-Pavon, J. C. (2017). Analysing concurrent transcranial magnetic stimulation and electroencephalographic data: A review and introduction to the open-source TESA software. Neuroimage, 147, 934-951.

% fastica:
% http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml 
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\eeglab_current\eeglab2019_1\plugins\FastICA_25\');

% ICLabel in your eeglab folder as a plugin or via the github:
% https://github.com/sccn/ICLabel

% Need to install RunLength: https://au.mathworks.com/matlabcentral/fileexchange/41813-runlength?s_tid=mwa_osa_a
% This can be tricky on a mac: you need to manually compile the 'RunLength' function as it wouldn't auto-install through the Matlab script. 
% Macs have very tight security settings around MEX files. To get around this, install Xcode for Mac and then run the line of code to get it working
addpath('D:\Data_Analysis\Analysis_Tools_and_Software\RunLength'); 

% Need to install MinGW-w64 if on windows if you haven't already:
% http://mingw-w64.org/doku.php

% Specify  RELAX folder location (this toolbox):
addpath('D:\Data_Analysis\RELAX_v0_9\');

% Cite muscle activity detection methods (the manuscript applied to ICA,
% whereas this method applies them to single electrodes rather than ICA):
% Fitzgibbon, S. P., DeLosAngeles, D., Lewis, T. W., Powers, D. M. W., Grummett, T. S., Whitham, E. M., ... & Pope, K. J. (2016). Automatic determination of EMG-contaminated components and validation of independent component analysis using EEG during pharmacologic paralysis. Clinical Neurophysiology, 127(3), 1781-1793.

% Specify your electrode locations with the correct cap file:
RELAX_cfg.caploc='D:\Data_Analysis\Cap_Location_Files\standard-10-5-cap385.elp'; %path containing electrode positions

% Specify the to be processed file locations:
RELAX_cfg.myPath='D:\DATA_TO_BE_PREPROCESSED\';

% Lastly, you may want to delete more unused electrodes than specified by
% this script. If so, search for and modify the section titled: "Delete channels that
% are not relevant if present" to include the electrodes you would like to
% delete.

%% Parameters that can be specified:

% The following selections determine how many epochs of each type of
% artifact to include in the mask. It's not clear what the best choice is.
% I have selected as defaults the parameters that seem to work for my data.
% If only performing one run, perhaps including all artifacts in the mask
% is best (as long as that leaves enough clean data for the clean data
% mask). The authors of the MWF manuscript suggest that it doesn't hurt to
% include clean data in the artifact mask, and so are happy to have wide
% boundaries around artifacts. 

% However, I wonder if including the majority of trials in the artifact
% mask could mean that most of the task related ERPs might be considered
% artifacts, and cleaned out from the data (including them in the clean
% mask would reduce this issue, as the clean mask defines what the good
% data should look like, and if it has ERPs, then the artifact periods will
% be cleaned of only the difference between the clean mask and the artifact
% mask, leaving the ERP in the data. I guess this relies upon a consistent
% ERP across both the clean and artifact masks however.

RELAX_cfg.Do_MWF_Once=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Twice=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no).
RELAX_cfg.Do_MWF_Thrice=1; % 1 = Perform the MWF cleaning a second time (1 for yes, 0 for no). I think cleaning drift in this is a good idea.
RELAX_cfg.Perform_wICA_on_ICLabel=1; % 1 = Perform wICA on artifact components marked by ICLabel (1 for yes, 0 for no).

RELAX_cfg.MWFRoundToCleanBlinks=2; % Which round to clean blinks in (1 for the first, 2 for the second...)
RELAX_cfg.ProbabilityDataHasNoBlinks=0; % 0 = data almost certainly has blinks, 1 = data might not have blinks, 2 = data definitely doesn't have blinks.
% 0 = eg. task related data where participants are focused with eyes open, 
% 1 = eg. eyes closed recordings, but with participants who might still open their eyes at times, 
% 2 = eg. eyes closed resting with highly compliant participants and recordings that were strictly made only when participants had their eyes closed.

RELAX_cfg.DriftSeverityThreshold=10; %MAD from the median of all electrodes. This could be set lower and would catch less severe drift 
RELAX_cfg.PercentWorstEpochsForDrift=0.30; % Maximum percent of epochs to include in the mask from drift artifact type.

RELAX_cfg.MWFDelayPeriod=8; % The MWF includes both spatial and temporal information when filtering out artifacts. Longer delays apparently improve performance. 
% Delay periods >5 can lead to generalised eigenvector rank deficiency in some files, 
% and if this occurs cleaning is ineffective. 5 was used by Somers et al (2018). 
% This is likely to be because data filtering creates a temporal
% dependency between consecutive datapoints, reducing their independence
% when including the temporal aspect in the MWF computation.
% To address this, I have set the MWF function to attempt MWF cleaning
% at the delay period set above, but if rank deficiency occurs, to reduce
% the delay period by 1 and try again (for 3 iterations).
% Using robust detrending (which does not create any temporal dependence,
% unlike filtering) may be an alternative which avoids rank deficiency 
% (but our initial test suggested this led to worse cleaning than filtering).
  
RELAX_cfg.ExtremeVoltageShiftThreshold=20; % How many MAD from the median of all epochs for each electrode against itself. This could be set lower and would catch less severe pops
RELAX_cfg.ExtremeAbsoluteVoltageThreshold=500; % microvolts max or min above which will be excluded from cleaning and deleted from data
RELAX_cfg.ExtremeImprobableVoltageDistributionThreshold=8; % SD from the mean of all epochs for each electrode against itself. This could be set lower and would catch less severe improbable data
RELAX_cfg.ExtremeSingleChannelKurtosisThreshold=8; % SD from the mean of the single electrodes. This could be set lower and would catch less severe kurtosis 
RELAX_cfg.ExtremeAllChannelKurtosisThreshold=8; % SD from the mean of all electrodes. This could be set lower and would catch less severe kurtosis 
RELAX_cfg.ExtremeBlinkShiftThreshold=8; % How many MAD from the median of blink affected epochs to exclude as extreme data (uses the higher value out of this value and RELAX_cfg.PopSeverityShiftThreshold above, which caters for the fact that blinks don't affect the median, so without this, if data is clean and blinks are large, blinks can get excluded as extreme
RELAX_cfg.ExtremeDriftSlopeThreshold=-4; % slope of log frequency log power below which to reject as drift without neural activity

% Clean periods that last for a shorter duration than the following value to be marked as artifacts, 
% and pad short artifact periods out into artifact periods of at least the following length when 
% they are shorter than this value to reduce rank deficiency issues in MWF cleaning). 
% Note that it's better to include clean periods in the artifact mask rather than the including artifact in the clean mask.
RELAX_cfg.MinimumArtifactDuration=1200; % in ms. It's better to make this value longer than 1000ms, as doing so will catch diminishing artifacts that aren't detected in a neighbouring 1000ms period, which might still be bad
RELAX_cfg.MinimumBlinkArtifactDuration=800; % blink marking is based on the maximum point of the blink rather than the 1000ms divisions for muscle artifacts, so this can be shorter than the value above (blinks do not typically last >500ms)

RELAX_cfg.BlinkElectrodes={'FP1';'FPZ';'FP2';'AF3';'AF4';'F3';'F1';'FZ';'F2';'F4'}; % sets the electrodes to average for blink detection using the IQR method. These should be frontal electrodes maximally affected by blinks
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
RELAX_cfg.PercentWorstEpochsForMuscle=0.50;  % Maximum amount of data periods to be marked as muscle artifact for cleaning by the MWF. You want at least a reasonable amount of both clean and artifact templates for effective cleaning.
% I set this reasonably high, because otherwise muscle artifacts could considerably influence the clean mask and add noise into the data
RELAX_cfg.ProportionOfMuscleContaminatedEpochsAboveWhichToRejectChannel=0.05; % If the proportion of all epochs from a single electrode that are marked as containing muscle activity is higher than this, the electrode is deleted
RELAX_cfg.ProportionOfExtremeNoiseAboveWhichToRejectChannel=0.05; % If the proportion of all epochs from a single electrode that are marked as containing extreme artifacts is higher than this, the electrode is deleted

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

%% RUN SCRIPT BELOW:

% Load pre-processing statistics file for these participants if it already
% exists (note that this can cause errors if the number of variables
% inserted into the output table differs between participants, which can be
% caused by using different parameters in the preceding section):

tic;

RELAX_cfg.OutputPath=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep];   % use fileseparators for increased compatability 
if ~exist([RELAX_cfg.OutputPath, 'dir']); mkdir(RELAX_cfg.OutputPath); end % make dir if non-existant 

cd(RELAX_cfg.OutputPath);
dirList=dir('*.mat');
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatisticsRoundOne.mat')==1
        load('ProcessingStatisticsRoundOne.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatisticsRoundTwo.mat')==1
        load('ProcessingStatisticsRoundTwo.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'ProcessingStatisticsRoundThree.mat')==1
        load('ProcessingStatisticsRoundThree.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'RawMetrics.mat')==1
        load('RawMetrics.mat');
    end
end
for x=1:numel(dirList)
    if  strcmp(dirList(x).name,'CleanedMetrics.mat')==1
        load('CleanedMetrics.mat');
    end
end

%% List all files in directory
cd(RELAX_cfg.myPath);
RELAX_cfg.dirList=dir('*.set');
RELAX_cfg.files={RELAX_cfg.dirList.name};
if isempty(RELAX_cfg.files)
    disp('No files found..')
end
%% Loop over each file in the directory list 
for Subjects=1:numel(RELAX_cfg.files)
    RELAX_cfg.filename=RELAX_cfg.files{Subjects};
    clearvars -except 'RELAX_cfg' 'Subjects' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingRoundOneAllParticipants' 'RELAXProcessingRoundTwoAllParticipants'...
        'RELAXProcessingRoundThreeAllParticipants' 'FilesWithRankDeficiencyRoundOne' 'FilesWithRankDeficiencyRoundTwo' 'FilesWithRankDeficiencyRoundThree' 'NoBlinksDetected' 'Warning';
    %% Load data (assuming the data is in EEGLAB .set format):
    
    cd(RELAX_cfg.myPath);
    EEG = pop_loadset(RELAX_cfg.filename);

    ParticipantID = extractBefore(RELAX_cfg.filename,".");
    EEG.RELAXProcessing.aParticipantID=cellstr(ParticipantID);
    
    EEG.RELAX.Data_has_been_averagerereferenced=0;
    EEG.RELAX.Data_has_been_cleaned=0;
    RELAX_cfg.ms_per_sample=(1000/EEG.srate);

    %% Select channels     
    EEG=pop_chanedit(EEG,  'lookup', RELAX_cfg.caploc);
    
    %% Delete channels that are not relevant if present       
    EEG=pop_select(EEG,'nochannel',RELAX_cfg.ElectrodesToDelete);
    EEG = eeg_checkset( EEG );
    EEG.allchan=EEG.chanlocs; % take list of all included channels before any rejections

    %% Band Pass filter data from 0.25 to 80Hz

    % Fourth order butterworth filtering at 0.3Hz and above distorts ERPs, whereas 0.2Hz
    % has a negligible effect on ERPs. Some of our data showed drift at
    % 0.3Hz, so we've gotten as close as possible to filtering that out
    % without falling in to the range of having a negative effect on the
    % ERPs. Cleaner data might be better to be filtered at 0.1Hz.
        
    % https://www.frontiersin.org/articles/10.3389/fpsyg.2012.00131/full
    
    % We also experimented with EEGLAB's cleanline to reconstruct the notch
    % filtered data (which reconstructs the 50Hz data rather than removing
    % it. However, when this was used, it introduced strange
    % high frequency artifacts at the MWF cleaning steps, which seemed to
    % depend on the duration setting of the MWF filter (artifacts were
    % longer when the filter was set for long). Similarly, using a higher 
    % order than 2 for the notch led to common rank insufficiency in the MWF step, 
    % as did any use of EEGLAB's pop_eegfiltnew. I think this is due to increased
    % temporal dependencies created with these filters, which adversely
    % affect the independence of timepoints required for the temporal
    % filter cleaning by the MWF. Some good alternatives might be robust
    % detrending and zapline reported by De Cheveigne (2018) and available
    % in the noisetools package (nt_detrend and nt_zapline). However,
    % automatic parameter setting is not available for this yet, and our
    % initial testing of nt_detrend resulted in worse performance (despite
    % enabling MWF cleaning delay periods up to 30 with no rank
    % deficiencies, demonstrating the benefits of using these methods to
    % avoid adding temporal dependencies to the data, unlike filtering
    % approaches).
    
    % de Cheveigné, A., & Arzounian, D. (2018). Robust detrending, rereferencing, outlier detection, and inpainting for multichannel data. NeuroImage, 172, 903-912.
    
    % Use TESA to apply butterworth filter: 
    EEG = pop_tesa_filtbutter( EEG, RELAX_cfg.LineNoiseFrequency-3, RELAX_cfg.LineNoiseFrequency+3, 2, 'bandstop' );
    EEG = pop_tesa_filtbutter( EEG, RELAX_cfg.HighPassFilter, RELAX_cfg.LowPassFilter, 4, 'bandpass' );

    %% Clean flat channels and bad channels showing improbable data:
    % PREP pipeline: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4471356/
    noisyOut = findNoisyChannels(EEG);  
    EEG.RELAXProcessing.PREPBasedChannelToReject={};
    for x=1:size(noisyOut.noisyChannels.all,2) % loop through output of PREP's findNoisyChannels and take a record of noisy electrodes for deletion:
        PREPBasedChannelToReject{x}=EEG.chanlocs(noisyOut.noisyChannels.all(x)).labels;
        EEG.RELAXProcessing.PREPBasedChannelToReject = PREPBasedChannelToReject';
    end
    EEG=pop_select(EEG,'nochannel',noisyOut.noisyChannels.all); % delete noisy electrodes detected by PREP
    EEG.RELAX.ListOfChannelsAfterRejections={EEG.chanlocs.labels}; % Get list of good channels
    continuousEEG=EEG;

    [continuousEEG, epochedEEG] = RELAX_excluding_channels_and_epoching(continuousEEG, RELAX_cfg); % Epoch data, detect extremely bad data, delete channels if over the set threshold for proportion of data affected by extreme outlier for each electrode
    [continuousEEG, epochedEEG] = RELAX_excluding_extreme_values(continuousEEG, epochedEEG, RELAX_cfg); % Mark extreme periods for exclusion from MWF cleaning, and deletion before wICA cleaning

    % Use the continuous data to detect eye blinks and mark
    % these in the EEG.event as well as in the mask. The output is
    % continuous data but includes all the previous extreme period 
    % markings from the epoched data.
    if RELAX_cfg.ProbabilityDataHasNoBlinks<2
        [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg); % use an IQR threshold method to detect and mark blinks
        if continuousEEG.RELAX.IQRmethodDetectedBlinks(1,1)==0 % If a participants doesn't show any blinks, make a note
            NoBlinksDetected{Subjects,1}=ParticipantID; 
        end
        [continuousEEG, epochedEEG] = RELAX_metrics_blinks(continuousEEG, epochedEEG); % record blink amplitude ratio from raw data for comparison.
    end
    
    rawEEG=continuousEEG; % Take a copy of the not yet cleaned data for calculation of all cleaning SER and ARR at the end
    
    if RELAX_cfg.saveextremesrejected==1
        mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Extremes_Rejected'])
        SaveSetExtremes_Rejected =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep 'Extremes_Rejected', filesep ParticipantID '_Extremes_Rejected.set'];    
        EEG = pop_saveset( rawEEG, SaveSetExtremes_Rejected ); % If desired, save data here with bad channels deleted, filtering applied, extreme outlying data periods marked
    end

    %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK AND CLEAN MUSCLE ARTIFACTS
    % Any one of these functions can be commented out to ignore those artifacts
    % when creating the mask    
    if RELAX_cfg.Do_MWF_Once==1

        % Use epoched data and FFT to detect slope of log frequency log
        % power, add periods exceeding muscle threshold to mask:
        [continuousEEG, epochedEEG] = RELAX_muscle(continuousEEG, epochedEEG, RELAX_cfg);        
        [continuousEEG, epochedEEG] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg); % record muscle contamination metrics from raw data for comparison.

        EEG=continuousEEG; % Return continuousEEG to the "EEG" variable for MWF processing

        % If including eye blink cleaning in first round MWF, then insert
        % eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==1
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=nanmean(EEG.RELAX.eyeblinkmask);
        end

        % Combine the extreme period mask NaNs into the full noise mask, so
        % these periods are ignored by the MWF cleaning template when
        % constructing an artifact and clean period mask.
        for e=1:size(EEG.RELAX.NaNsForExtremeOutlierPeriods,2)
            if isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods(1,e))
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(1,e)=NaN;
            end
        end
        
        % The following eliminates very brief lengths of marked periods 
        % during the mask for MWF cleaning (without doing this, very short periods can
        % lead to rank deficiency)    
        % If period has been marked as shorter than
        % RELAX_cfg.MinimumArtifactDuration, then pad it out:
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks');
        
        EEG.RELAX.NoiseMaskFullLengthR1=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedAllArtifacts=nanmean(EEG.RELAXProcessing.Details.NoiseMaskFullLength);
              
        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);          
          
        EEG.RELAXProcessingRoundOne=EEG.RELAXProcessing;          
        RELAXProcessingRoundOne=EEG.RELAXProcessingRoundOne;
        
        if isfield(RELAXProcessingRoundOne,'Details')
            RELAXProcessingRoundOne=rmfield(RELAXProcessingRoundOne,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundOne,'Details')
                EEG.RELAXProcessingRoundOne=rmfield(EEG.RELAXProcessingRoundOne,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingRoundOneAllParticipants(Subjects,:) = struct2table(RELAXProcessingRoundOne,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        % Save round 1 MWF pre-processing:
        if RELAX_cfg.saveround1==1
            mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF'])
            SaveSetMWF1 =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '1xMWF', filesep ParticipantID '_MWF1.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF1 ); 
        end
    end

    %% PERFORM A SECOND ROUND OF MWF. THIS IS HELPFUL IF THE FIRST ROUND DOESN'T SUFFICIENTLY CLEAN ARTIFACTS. 

    % It has been suggested to be useful by Somers et al (2018)
    % (particularly when used in a cascading fashion). 

    % However, I can see risks. If artifact masks fall on task relevant
    % activity in both rounds of the MWF, it may be that the task relevant data
    % is just cleaned right out of the signal.
    
    if RELAX_cfg.Do_MWF_Twice==1

        ParticipantID = extractBefore(RELAX_cfg.filename,".");
        EEG.RELAXProcessing.aParticipantID=cellstr(ParticipantID);
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        
        % If blinks weren't initially detected, detect them here
        % (this happens in <1/200 cases, but is a good back up.
        if RELAX_cfg.ProbabilityDataHasNoBlinks==0
            if EEG.RELAX.IQRmethodDetectedBlinks(1,1)==0
                continuousEEG=EEG;
                [continuousEEG, epochedEEG] = RELAX_blinks_IQR_method(continuousEEG, epochedEEG, RELAX_cfg);
                EEG=continuousEEG;
            end
        end
        
        % If including eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask:
        if isfield(EEG.RELAX, 'eyeblinkmask')
            if RELAX_cfg.MWFRoundToCleanBlinks==2
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                EEG.RELAXProcessing.ProportionMarkedBlinks=nanmean(EEG.RELAX.eyeblinkmask);
            end
        end
  
        % Combine the extreme masks NaNs into the full noise mask:
        for e=1:size(EEG.RELAX.NaNsForExtremeOutlierPeriods,2)
            if isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods(1,e))
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(1,e)=NaN;
            end
        end
        % The following eliminates very brief lengths of currently clean marked
        % periods during the mask (without doing this, very short periods can
        % lead to rank deficiency)  
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'blinks');
        
        EEG.RELAX.NoiseMaskFullLengthR2=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedAllArtifacts=nanmean(EEG.RELAXProcessing.Details.NoiseMaskFullLength);

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);           
        
        EEG.RELAXProcessingRoundTwo=EEG.RELAXProcessing;
        EEG.RELAX.ProportionMarkedAllArtifactsR2=EEG.RELAXProcessing.ProportionMarkedAllArtifacts;        
        RELAXProcessingRoundTwo=EEG.RELAXProcessingRoundTwo;
        if isfield(RELAXProcessingRoundTwo,'Details')
            RELAXProcessingRoundTwo=rmfield(RELAXProcessingRoundTwo,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundTwo,'Details')
                EEG.RELAXProcessingRoundTwo=rmfield(EEG.RELAXProcessingRoundTwo,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingRoundTwoAllParticipants(Subjects,:) = struct2table(RELAXProcessingRoundTwo,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');
        % Save round 1 MWF pre-processing:
        if RELAX_cfg.saveround2==1
            mkdir([RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF'])
            SaveSetMWF2 =[RELAX_cfg.myPath, filesep 'RELAXProcessed' filesep '2xMWF', filesep ParticipantID '_MWF2.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF2 ); 
        end     
    end
    
    %% PERFORM A THIRD ROUND OF MWF. 

    % It has been suggested to be useful by Somers et al (2018)
    % (particularly when used in a cascading fashion). 

    % However, I can see risks. If artifact masks fall on task relevant
    % activity in both rounds of the MWF, it may be that the task relevant data
    % is just cleaned right out of the signal.
    
    if RELAX_cfg.Do_MWF_Thrice==1
        
        ParticipantID = extractBefore(RELAX_cfg.filename,".");
        EEG.RELAXProcessing.aParticipantID=cellstr(ParticipantID);
        
        EEG.RELAXProcessing.ProportionMarkedBlinks=0;
        % If less than 5% of data was masked as eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask in round 3:
        if EEG.RELAX.ProportionMarkedAllArtifactsR2<0.05
            if isfield(EEG.RELAX, 'eyeblinkmask')
                EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
                EEG.RELAX.eyeblinkmask(isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods))=NaN;
                EEG.RELAXProcessing.ProportionMarkedBlinks=nanmean(EEG.RELAX.eyeblinkmask);
            end
        end

        % Epoch the data into 1 second epochs with a 500ms overlap. Outputs
        % both the ContinuousEEG (which has been filtered above by this
        % point) and the epoched data as EEG.
        [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
        
        %% THIS SECTION CONTAINS FUNCTIONS WHICH MARK ARTIFACTS

        % Any one of these functions can be commented out to ignore those artifacts
        % when creating the mask
        
        % Use epoched data to add periods showing excessive amplitudes at
        % any time, excessive amplitude shifts across the epoch, and epochs
        % with excessive deviation from the all channel median to the mask
        % (based on median absolute deviation):
        
        [continuousEEG, epochedEEG] = RELAX_drift(continuousEEG, epochedEEG, RELAX_cfg);
        
        % Use the filtered continuous data to detect horizontal eye
        % movements and mark these in the EEG.event as well as in the mask.
        % You may want to simply reject horizontal eye movements at a later
        % stage if your task requires participants to look straight ahead
        % for the entire task. Alternatively, if your task requires
        % participants to complete horizontal eye movements time locked to
        % a stimuli, this section will mark every event with these
        % horizontal eye movements as an artifact, and should not be
        % implemented.
        
        % The output is continuous data:
        [continuousEEG] = RELAX_horizontaleye(continuousEEG, RELAX_cfg);
    
        %% Return to the "EEG" variable for MWF processing:
        EEG=continuousEEG;
        
        % If including eye blink cleaning in second round MWF, then insert
        % eye blink mask into noise mask:
        if RELAX_cfg.MWFRoundToCleanBlinks==3
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(EEG.RELAX.eyeblinkmask==1)=1;
            EEG.RELAX.eyeblinkmask(isnan(EEG.RELAXProcessing.Details.NaNsForNonEvents))=NaN;
            EEG.RELAXProcessing.ProportionMarkedBlinks=nanmean(EEG.RELAX.eyeblinkmask);
        end
        
        [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, 'notblinks');
        
        EEG.RELAX.NoiseMaskFullLengthR3=EEG.RELAXProcessing.Details.NoiseMaskFullLength;
        EEG.RELAXProcessing.ProportionMarkedAllArtifacts=nanmean(EEG.RELAXProcessing.Details.NoiseMaskFullLength);

        %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED ABOVE:
        [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg);               
        
        EEG.RELAXProcessingRoundThree=EEG.RELAXProcessing; 
        RELAXProcessingRoundThree=EEG.RELAXProcessing;
        
        if isfield(RELAXProcessingRoundThree,'Details')
            RELAXProcessingRoundThree=rmfield(RELAXProcessingRoundThree,'Details');
        end
        if RELAX_cfg.KeepAllInfo==0
            if isfield(EEG.RELAXProcessingRoundThree,'Details')
                EEG.RELAXProcessingRoundThree=rmfield(EEG.RELAXProcessingRoundThree,'Details');
            end
        end
        % Record processing statistics for all participants in single table:
        RELAXProcessingRoundThreeAllParticipants(Subjects,:) = struct2table(RELAXProcessingRoundThree,'AsArray',true);
        EEG = rmfield(EEG,'RELAXProcessing');

        if RELAX_cfg.saveround3==1
            mkdir([RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep '3xMWF'])
            SaveSetMWF3 =[RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep '3xMWF', filesep ParticipantID '_MWF3.set'];    
            EEG = pop_saveset( EEG, SaveSetMWF3 ); 
        end         
    end
    
    %% Perform robust average re-referencing of the data, reject periods marked as extreme outliers
    
    if RELAX_cfg.Do_MWF_Once==0
        EEG=continuousEEG;
    end

    [EEG] = RELAX_average_rereference(EEG);

    EEG = eeg_checkset( EEG );        
    EEG.NumberOfChannelsAfterRejections=EEG.nbchan;

    % Reject periods that were marked as NaNs in the MWF masks because they 
    % showed extreme shift within the epoch or extremely improbable data:
    EEG = eeg_eegrej( EEG, EEG.RELAX.ExtremelyBadPeriodsForDeletion);

    %% Perform wICA on ICLabel identified artifacts that remain:
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        % The following performs wICA, implemented on only the components
        % marked as artifact by ICLabel.
        
        % fastica_symm is repeated up to 3 times in the case of non-convergence
        % to ensure non-convergence doesn't impair cleaning. fastica
        % performs quickly on continuous data and doesn't seem to be
        % inferior at cleaning compared to AMICA (which is much slower). It
        % also seems to be comparable (or only slightly worse) than
        % extended infomax (run via cudaICA for speed).
        EEG.RELAXProcessing_wICA.aParticipantID=cellstr(ParticipantID);
        [EEG,~, ~, ~, ~] = RELAX_wICA_on_ICLabel_artifacts(EEG,'fastica_symm', 1, 0, EEG.srate, 5,'coif5'); % add: 'Report_all_wICA_info' to the end = optionally report proportion of ICs categorized as each category, and variance explained by ICs from each category (~20s slower if on)
        EEG = eeg_checkset( EEG );

        RELAXProcessing_wICA=EEG.RELAXProcessing_wICA;
        % Record processing statistics for all participants in single table:
        RELAXProcessing_wICA_AllParticipants(Subjects,:) = struct2table(RELAXProcessing_wICA,'AsArray',true);
    end
    
    EEG.RELAX.Data_has_been_cleaned=1;
    
    %% COMPUTE CLEANED METRICS:
        
    [continuousEEG, epochedEEG] = RELAX_epoching(EEG, RELAX_cfg);
    [continuousEEG, ~] = RELAX_metrics_blinks(continuousEEG, epochedEEG);
    [continuousEEG, ~] = RELAX_metrics_muscle(continuousEEG, epochedEEG, RELAX_cfg);
    [continuousEEG] = RELAX_metrics_final_SER_and_ARR(rawEEG, continuousEEG); % this is only a good metric for testing only the cleaning of artifacts marked for cleaning by MWF, see notes in function.
    EEG=continuousEEG;
    EEG = rmfield(EEG,'RELAXProcessing');

    if isfield(EEG.RELAX, 'Metrics')
        if isfield(EEG.RELAX.Metrics, 'Cleaned')
            if isfield(EEG.RELAX.Metrics.Cleaned,'BlinkAmplitudeRatio')
                CleanedMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX.Metrics.Cleaned.BlinkAmplitudeRatio,1),Subjects)=EEG.RELAX.Metrics.Cleaned.BlinkAmplitudeRatio;
                CleanedMetrics.BlinkAmplitudeRatio(CleanedMetrics.BlinkAmplitudeRatio==0)=NaN;
            end
            if isfield(EEG.RELAX.Metrics.Cleaned,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(Subjects)=EEG.RELAX.Metrics.Cleaned.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                CleanedMetrics.MeanMuscleStrengthScaledByProportionShowingMuscle(Subjects)=EEG.RELAX.Metrics.Cleaned.MeanMuscleStrengthScaledByProportionShowingMuscle;
                CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(Subjects)=EEG.RELAX.Metrics.Cleaned.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
            end
            if isfield(EEG.RELAX.Metrics.Cleaned,'All_SER')
                CleanedMetrics.All_SER.(pipeline{pipe})(Subjects,1)=EEG.RELAX.Metrics.Cleaned.All_SER;
                CleanedMetrics.All_ARR.(pipeline{pipe})(Subjects,1)=EEG.RELAX.Metrics.Cleaned.All_ARR;
            end
        end
        if isfield(EEG.RELAX.Metrics, 'Raw')
            if isfield(EEG.RELAX.Metrics.Raw,'BlinkAmplitudeRatio')
                RawMetrics.BlinkAmplitudeRatio(1:size(EEG.RELAX.Metrics.Raw.BlinkAmplitudeRatio,1),Subjects)=EEG.RELAX.Metrics.Raw.BlinkAmplitudeRatio;
                RawMetrics.BlinkAmplitudeRatio(RawMetrics.BlinkAmplitudeRatio==0)=NaN;
            end
            if isfield(EEG.RELAX.Metrics.Raw,'MeanMuscleStrengthFromOnlySuperThresholdValues')
                RawMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues(Subjects)=EEG.RELAX.Metrics.Raw.MeanMuscleStrengthFromOnlySuperThresholdValues; 
                RawMetrics.MeanMuscleStrengthScaledByProportionShowingMuscle(Subjects)=EEG.RELAX.Metrics.Raw.MeanMuscleStrengthScaledByProportionShowingMuscle;
                RawMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel(Subjects)=EEG.RELAX.Metrics.Raw.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel;
            end
        end   
    end
    
    %% SAVE FILE:
    
    mkdir([RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep 'Cleaned_Data'])
    SaveSetMWF2 =[RELAX_cfg.myPath,filesep 'RELAXProcessed' filesep 'Cleaned_Data', filesep ParticipantID '_RELAX.set'];    
    EEG = pop_saveset( EEG, SaveSetMWF2 );  
    
    %% Save statistics for each participant and across participants, graph cleaning metrics:

    if RELAX_cfg.Do_MWF_Once==1
        savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundOne'];
        save(savefileone,'RELAXProcessingRoundOneAllParticipants')
    end
    if RELAX_cfg.Do_MWF_Twice==1
        savefiletwo=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundTwo'];
        save(savefiletwo,'RELAXProcessingRoundTwoAllParticipants')
    end
    if RELAX_cfg.Do_MWF_Thrice==1
        savefilethree=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatisticsRoundThree'];
        save(savefilethree,'RELAXProcessingRoundThreeAllParticipants')
    end
    if RELAX_cfg.Perform_wICA_on_ICLabel==1
        savefilefour=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'ProcessingStatistics_wICA'];
        save(savefilefour,'RELAXProcessing_wICA_AllParticipants')
    end
    if exist('CleanedMetrics','var')
        savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'CleanedMetrics'];
        save(savemetrics,'CleanedMetrics')
    end
    if exist('RawMetrics','var')
        savemetrics=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RawMetrics'];
        save(savemetrics,'RawMetrics')
    end
    savefileone=[RELAX_cfg.myPath filesep 'RELAXProcessed' filesep 'RELAX_cfg'];
    save(savefileone,'RELAX_cfg')
end

clearvars -except 'RELAX_cfg' 'Subjects' 'CleanedMetrics' 'RawMetrics' 'RELAXProcessingRoundOneAllParticipants' 'RELAXProcessingRoundTwoAllParticipants'...
        'RELAXProcessingRoundThreeAllParticipants' 'FilesWithRankDeficiencyRoundOne' 'FilesWithRankDeficiencyRoundTwo' 'FilesWithRankDeficiencyRoundThree' 'NoBlinksDetected' 'Warning';
    
if exist('CleanedMetrics','var')
    try
        figure('Name','BlinkAmplitudeRatio');
        boxplot(CleanedMetrics.BlinkAmplitudeRatio);
    catch
    end
    try
        figure('Name','MeanMuscleStrengthFromOnlySuperThresholdValues');
        plot(CleanedMetrics.MeanMuscleStrengthFromOnlySuperThresholdValues);
    catch
    end
    try
        figure('Name','MeanMuscleStrengthScaledByProportionShowingMuscle');
        plot(CleanedMetrics.MeanMuscleStrengthScaledByProportionShowingMuscle);
    catch
    end
    try
        figure('Name','ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel');
        plot(CleanedMetrics.ProportionOfEpochsShowingMuscleAboveThresholdAnyChannel);
    catch
    end
end

toc
%% POTENTIAL IMPROVEMENTS THAT COULD BE MADE:
% 1) work out a way to threshold horizontal eye movements so the script
% catches the onset and offset, rather than just the absolute +/-2MAD on
% opposite sides of the head?
% 2) an adaptive threshold for the wICA cleaning
% 3) adding a requirement that the IQR blink detection method detects that
% positive amplitude shifts are biased towards frontal electrodes?
% 4) instead of using typical ICLabel, adding objective muscle slope and eye
% movement metrics to detect these artifactual components
% 5) Independent Vector Analysis might also improve component separation,
% although I couldn't get this to work with ICLabel.
