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

%% RELAX_perform_MWF_cleaning:

function [EEG] = RELAX_perform_MWF_cleaning (EEG, RELAX_cfg)

    %% This function uses the MWF cleaning toolbox (Somers et al. 2018) and the artifact and clean data period template created by the RELAX functions to clean EEG data
    % Somers, B., Francart, T. and Bertrand, A. (2018). A generic EEG artifact removal algorithm based on the multi-channel Wiener filter. Journal of Neural Engineering, 15(3), 036007. DOI: 10.1088/1741-2552/aaac92

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'MWFDelayPeriod')==0
            RELAX_cfg.MWFDelayPeriod=8; % The MWF includes both spatial and temporal information when filtering out artifacts. Longer delays apparently improve performance. 
        end
        if isfield(RELAX_cfg, 'KeepAllInfo')==0
            RELAX_cfg.KeepAllInfo=0; 
        end
    elseif exist('RELAX_cfg', 'var')==0  
        RELAX_cfg.MWFDelayPeriod=8; 
        RELAX_cfg.KeepAllInfo=0; 
    end
    
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
    
    % First, order the MWF Processing statistics structure in alphabetical order:
    [~, neworder] = sort(lower(fieldnames(EEG.RELAXProcessing)));
    EEG.RELAXProcessing = orderfields(EEG.RELAXProcessing, neworder);
    EEG.RELAXProcessing.SignalToErrorRatio=NaN;
    EEG.RELAXProcessing.ArtifactToResidueRatio=NaN;
    EEG.RELAXProcessing.RankDeficiency={};
    EEG.RELAXProcessing.DelayPeriod=NaN;
    
    %% RUN MWF TO CLEAN DATA BASED ON MASKS CREATED IN RELAX FUNCTIONS:
    if EEG.RELAXProcessing.ProportionMarkedInMWFArtifactMaskTotal>0.05    
        [cleanEEG, d, W, SER, ARR] = mwf_process (EEG.data, EEG.RELAXProcessing.Details.NoiseMaskFullLength, RELAX_cfg.MWFDelayPeriod);
        EEG.RELAXProcessing.DelayPeriod=RELAX_cfg.MWFDelayPeriod;
        % If Generalized eigenvectors are not scaled as assumed, try again
        % with a shorter delay period, up to 3 times:
        [warnmsg] = lastwarn;
        pattern = "eigenvectors";
        if contains(warnmsg,pattern)
            warning(['Trying again with a shorter delay period.',[]]);
            clear warnmsg
            lastwarn('')
            [cleanEEG, d, W, SER, ARR] = mwf_process (EEG.data, EEG.RELAXProcessing.Details.NoiseMaskFullLength, (RELAX_cfg.MWFDelayPeriod-1));
            EEG.RELAXProcessing.DelayPeriod=RELAX_cfg.MWFDelayPeriod-1;
        end
        [warnmsg] = lastwarn;
        if contains(warnmsg,pattern)
            warning(['Trying again with a shorter delay period.',[]]);
            clear warnmsg
            lastwarn('')
            [cleanEEG, d, W, SER, ARR] = mwf_process (EEG.data, EEG.RELAXProcessing.Details.NoiseMaskFullLength, (RELAX_cfg.MWFDelayPeriod-2));
            EEG.RELAXProcessing.DelayPeriod=RELAX_cfg.MWFDelayPeriod-2;
        end
        [warnmsg] = lastwarn;
        if contains(warnmsg,pattern)
            warning(['Trying again with a shorter delay period.',[]]);
            clear warnmsg
            lastwarn('')
            [cleanEEG, d, W, SER, ARR] = mwf_process (EEG.data, EEG.RELAXProcessing.Details.NoiseMaskFullLength, (RELAX_cfg.MWFDelayPeriod-3));
            EEG.RELAXProcessing.DelayPeriod=RELAX_cfg.MWFDelayPeriod-3;
        end

        EEG.data=cleanEEG; % Copy cleaned EEG data to the EEG struct

        EEG.RELAXProcessing.Details.EstimatedArtifactInEachChannel=d;
        EEG.RELAXProcessing.Details.MatrixUsedToEstimateArtifacts=W;
        EEG.RELAXProcessing.SignalToErrorRatio=SER;
        EEG.RELAXProcessing.ArtifactToResidueRatio=ARR;

        % Check for rank deficiency and note the issue if there is an issue:
        [warnmsg] = lastwarn;
        pattern = "eigenvectors";
        if contains(warnmsg,pattern)
            EEG.RELAXProcessing.RankDeficiency=warnmsg;
        end
        EEG = eeg_checkset( EEG ); 
        clear warnmsg
        lastwarn('')
        if RELAX_cfg.KeepAllInfo==0
            EEG.RELAXProcessing=rmfield(EEG.RELAXProcessing,'Details');
        end
    end
end