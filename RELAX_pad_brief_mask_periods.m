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

%% RELAX_pad_brief_mask_periods:

function [EEG] = RELAX_pad_brief_mask_periods (EEG, RELAX_cfg, type)

    %Check inputs
    if ~(strcmp(type,'blinks') || strcmp(type,'notblinks'))
        error('Input ''type'' needs to be either ''blinks'' or ''notblinks''.')
    end

    % set some defaults for included channels and trials, if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'ms_per_sample')==0
            RELAX_cfg.ms_per_sample=(1000/EEG.srate); 
        end
        if isfield(RELAX_cfg, 'MinimumArtifactDuration')==0
            RELAX_cfg.MinimumArtifactDuration=1200; 
        end
        if isfield(RELAX_cfg, 'MinimumBlinkArtifactDuration')==0
            RELAX_cfg.MinimumBlinkArtifactDuration=800; 
        end
    elseif exist('RELAX_cfg', 'var')==0     
        RELAX_cfg.MinimumArtifactDuration=1200;
        RELAX_cfg.MinimumBlinkArtifactDuration=800;
        RELAX_cfg.ms_per_sample=(1000/EEG.srate);
    end
    
    if type(strcmp(type,'blinks'))
        MinimumArtifactDuration=RELAX_cfg.MinimumBlinkArtifactDuration;
    elseif type(strcmp(type,'notblinks'))
        MinimumArtifactDuration=RELAX_cfg.MinimumArtifactDuration;
    end
    
    %% The following eliminates very brief lengths of currently clean marked periods during the mask (without doing this, very short periods can lead to rank deficiency).    
    % Combine the extreme period mask NaNs into the full noise mask, so
    % these periods are ignored by the MWF cleaning template when
    % constructing an artifact and clean period mask. This is best done at
    % this stage, just in case the padding steps above encroach into
    % periods marked as NaN for extreme outlier artifacts, combine the pop
    % masks NaNs into the full noise mask:
    for e=1:size(EEG.RELAX.NaNsForExtremeOutlierPeriods,2)
        if isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods(1,e))
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(1,e)=NaN;
        end
    end

    % If an artifact has been marked as shorter than
    % MinimumArtifactDuration, then pad it out:
    clear shortOneRunIndex; clear b; clear n; clear bi;
    [b, n, bi] = RunLength(EEG.RELAXProcessing.Details.NoiseMaskFullLength);
    shortOneRunIndex = find(b==1 & n<round(MinimumArtifactDuration/RELAX_cfg.ms_per_sample));
    for ns = 1:size(shortOneRunIndex,2)
        Start=(round((bi(shortOneRunIndex(ns))+bi(shortOneRunIndex(ns)+1))/2)-(MinimumArtifactDuration/2));
        End=(round((bi(shortOneRunIndex(ns))+bi(shortOneRunIndex(ns)+1))/2)+(MinimumArtifactDuration/2));
        EEG.RELAXProcessing.Details.NoiseMaskFullLength(Start:End)=1;
    end
    
    % If a clean patch lasts less the MinimumArtifactDuration and either side is marked as an artifact, 
    % close the two artifact markings on either side together:
    clear shortZeroRunIndex; clear b; clear n; clear bi; clear Start; clear End;
    [b, n, bi] = RunLength(EEG.RELAXProcessing.Details.NoiseMaskFullLength);
    shortZeroRunIndex = find(b==0 & n<round(MinimumArtifactDuration/RELAX_cfg.ms_per_sample));
    for ns = 1:size(shortZeroRunIndex,2)
        if b(shortZeroRunIndex(ns)-1)==1 && b(shortZeroRunIndex(ns)+1)==1
            Start=(round((bi(shortZeroRunIndex(ns))+bi(shortZeroRunIndex(ns)+1))/2)-(MinimumArtifactDuration/2));
            End=(round((bi(shortZeroRunIndex(ns))+bi(shortZeroRunIndex(ns)+1))/2)+(MinimumArtifactDuration/2));
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(Start:End)=1;
        end
    end

    % Combine the extreme period mask NaNs into the full noise mask again
    % in case the preceding steps encroached into the periods already
    % marked as NaN:
    for e=1:size(EEG.RELAX.NaNsForExtremeOutlierPeriods,2)
        if isnan(EEG.RELAX.NaNsForExtremeOutlierPeriods(1,e))
            EEG.RELAXProcessing.Details.NoiseMaskFullLength(1,e)=NaN;
        end
    end

    % If a patch still lasts less than minimum artifact duration, mark as NaNs:
    clear shortRunIndex; clear b; clear n; clear bi;
    [b, n, bi] = RunLength(EEG.RELAXProcessing.Details.NoiseMaskFullLength);
    shortRunIndex = find(b<2 & n<round(MinimumArtifactDuration/RELAX_cfg.ms_per_sample));
    for ns = 1:size(shortRunIndex,2)
        EEG.RELAXProcessing.Details.NoiseMaskFullLength(bi(shortRunIndex(ns)):bi(shortRunIndex(ns)+1)-1) = NaN(1,bi(shortRunIndex(ns)+1)-(bi(shortRunIndex(ns))));
    end
end