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

%% RELAX_RegressionBL_Correction:

function [Regression_BL_Corrected_EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,varargin)

    % This script performs a regression based baseline correction method
    % (the benefits of this over subtraction BL correction are explained in Alday, 2019, 
    % and the specific implementation performed in Bailey et al. 2022)
    
    % Alday, P. M. (2019). How much baseline correction do we need in ERP research? Extended GLM model can replace baseline correction while lifting its limits. Psychophysiology, 56(12), e13451.
    % Bailey et al. (2022). Meditators probably show increased behaviour-monitoring related neural activity. 
    
    %% Use as: 
    
    % [EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,'Factor_1_Level_1', {'HappyGo' 'SadGo' }, 'Factor_2_Level_1', {'HappyGo' 'HappyNogo' }); 
    % if a 2 x 2 design, this will code triggers other than 'HappyGo'/'SadGo' as -1, and Go as 1 in the first factor, and triggers other than 'HappyGo'/'HappyNogo' as -1 in the second factor

    %[EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg,'Factor_1_Level_1',{'Go'}); 
    % if a 2 condition design, this will code triggers other than 'Go' as -1, and Go as 1

    %[EEG]=RELAX_RegressionBL_Correction(EEG,RELAX_epoching_cfg); 
    % if only 1 stimulus condition present for each participant
    
    %% check inputs
    Factor1_present=0;
    for k = 1:length(varargin)
        if strcmpi(varargin{k},'Factor_1_Level_1')
            Factor1_present=1;
            Factor_1_Level_1_Codes = varargin{k+1}; 
            varargin{k+1}=[]; 
            varargin{k}=[];
        end
    end
    Factor2_present=0;
    for k = 1:length(varargin)
        if strcmpi(varargin{k},'Factor_2_Level_1')
            Factor2_present=1;
            Factor_2_Level_1_Codes = varargin{k+1}; 
            varargin{k+1}=[]; 
            varargin{k}=[];
        end
    end
    
    %% Set levels for each factor present:

    % Obtain a list of conditions, using sum coding to preserve the main
    % effect between conditions, which is included in the "confounds"
    % variable but not removed by the regression. Note that this design is
    % for a 2 x 2 interaction. If there is only 1 factor with 2 conditions
    % (eg. Go vs Nogo trials) then simply use the same sum coding (with
    % condition 1 coded as 1, condition 2 coded as -1), but reduce the 
    % number of columns in the "confounds" variable by 1.

    % 2nd column codes for emotion, 3rd for Go/Nogo, 1st column will include mean baseline voltage for regressing this out of the data
    for e=1:size(EEG.data,3) % for each event:
        confounds(e,1)=0;
        if Factor1_present==1
            if contains(EEG.event(e).type,Factor_1_Level_1_Codes)==1 
                confounds(e,2)=1; % code HappyGo as [1 1] 
            elseif contains(EEG.event(e).type,Factor_1_Level_1_Codes)==0 
                confounds(e,2)=-1; % code HappyNogo as [1 -1]
            end
        end
        if Factor2_present==1
            if contains(EEG.event(e).type,Factor_2_Level_1_Codes)==1 
                confounds(e,3)=1; % code SadGo as [-1 1]
            elseif contains(EEG.event(e).type,Factor_2_Level_1_Codes)==0 
                confounds(e,3)=-1; % code SadNogo as [-1 -1]
            end
        end
    end
    confounds(:,size(confounds,2)+1)=1; % add a constant

    % Note that for designs with more than 3 conditions for a single
    % factor, linear mixed modelling to baseline correct the data is more
    % appropriate (not provided here).

    % Also note that the method provided below is effective for 
    % nice clean data with a larger number of epochs per participant,
    % but for data that is not clean or few epochs per participant, 
    % a method that includes all participants in the single regression 
    % and includes individual as a factor is more appropriate.
    
    %% Determine mean amplitude in baseline period and remove with regression:

    BLperiodstart=find(EEG.times==RELAX_epoching_cfg.BLperiod(1,1)); % set time window for baseline correction by reference to timepoints in the epoch (in seconds)
    BLperiodend=find(EEG.times==RELAX_epoching_cfg.BLperiod(1,2)); % this window will be established in the following line based on sample points in the epoch (rather than timepoints)

    for c=1:size(EEG.data,1) % c = channels, for each channel
        confounds(:,1) = squeeze(mean(EEG.data(c,BLperiodstart:BLperiodend,:),2)); % obtain a list of the mean amplitude in the BL period from each trial to regress this out from each timepoint for that channel
        for s = 1:size(EEG.data, 2) % s = samples, for each sample
            Beta(:,s) =  confounds \  squeeze(EEG.data(c, s, :)); % B = X\Y = confounds \ data
        end
        for s=1:size(Beta,2)
            model(1,s,:) = confounds(:,1) .* squeeze(Beta(1,s)); % model = confounds * weights = X * X\Y 
            % (because only column 1 is included in the model, only the influence of the mean BL period is included in the model)
        end
        regression_BL_corrected_data(c, :, :) = EEG.data(c, :, :) - model(:,:,:); % Yclean = data - model = Y - X * X\Y  (the model of the confounding mean BL data is subtracted from the data)
    end
    EEG.data=regression_BL_corrected_data;
    Regression_BL_Corrected_EEG = eeg_checkset( EEG );
    
end