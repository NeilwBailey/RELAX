%% RELAX EEG CLEANING PIPELINE, Copyright (C) (2022) Neil Bailey

% Adapted using code from: tesa_compselect.m:

% Copyright (C) 2016  Nigel Rogasch, Monash University,
% nigel.rogasch@monash.edu

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

%% RELAX_muscle:

function [ICsMostLikelyMuscle] = RELAX_muscle_IC_detection(EEG,RELAX_cfg)

    % set some defaults if not specified:
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'MuscleSlopeThresholdFor_wICA_cleaning')==0
            RELAX_cfg.MuscleSlopeThresholdFor_wICA_cleaning=-0.59;
            % log frequency/ log power slope above which to mark ICs as muscle artifact. Less stringent = -0.31, Middle Stringency = -0.59 or
            % More stringent = -0.72, more negative values clean more muscle, but also perhaps some brain activity.
        end
        if isfield(RELAX_cfg, 'LineNoiseOrNotchFilterAffectedData')==0
            RELAX_cfg.LineNoiseOrNotchFilterAffectedData=1; % Removes influence of line noise (or line noise filtering) from the muscle slope computations.
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.MuscleSlopeThresholdFor_wICA_cleaning=-0.59;
        RELAX_cfg.LineNoiseOrNotchFilterAffectedData=1;
    end

    % Calculate pwelch
    [pxx,fp] = pwelch(EEG.icaact',size(EEG.icaact,2),[],size(EEG.icaact,2),EEG.srate);
    FFTout = pxx';
    fp = fp';
    
    % Calculate FFT bins
    freq=7:0.5:75;
    fftBins = zeros(size(FFTout,1),size(freq,2)); %preallocate
    for a=1:size(freq,2)
        [~, index1]=min(abs(fp-((freq(1,a)-0.25))));
        [~, index2]=min(abs(fp-((freq(1,a)+0.25))));
        fftBins(:,a)=mean(FFTout(:,index1:index2),2); %creates bins for 0.5 Hz in width centred around whole frequencies (i.e. 0.5, 1, 1.5 Hz etc)
    end
    
    % Define frequencies to include in the analysis
    muscleSlopesFromICs=zeros(1,size(EEG.icaact,1));
    for compNum =1:size(EEG.icaweights,1)
        [~,fin1] = min(abs(7 - freq));
        [~,fin2] = min(abs(75 - freq));

        freqHz = freq(1,fin1:fin2);
        freqPow = fftBins(compNum,fin1:fin2);

        if RELAX_cfg.LineNoiseOrNotchFilterAffectedData==1 
            [~,fex1] = min(abs(RELAX_cfg.LineNoiseFrequency-2 - freqHz));
            [~,fex2] = min(abs(RELAX_cfg.LineNoiseFrequency+2 - freqHz));
            freqHz(fex1:fex2) = [];
            freqPow(fex1:fex2) = [];
        end

        % Fit linear regression to log-log data
        p = polyfit(log(freqHz),log(freqPow),1);

        % Store the slope
        muscleSlopesFromICs(compNum) = p(1);
    end
    ICsMostLikelyMuscle=(muscleSlopesFromICs>-0.59);
end
