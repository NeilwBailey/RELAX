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

%% RELAX_horizontaleye:

function [continuousEEG] = RELAX_horizontaleye(continuousEEG, RELAX_cfg)

    % This function marks periods that show lateral electrode amplitude shifts
    % of more than the MAD threshold in the opposite direction on the opposite side of the
    % head as a horizontal eye movement (adapted from TESA [Rogasch et al.
    % 2018], which uses SD instead of MAD). You may not want to include this in
    % your mask if your task requires participants to focus straight ahead (and
    % may instead like to reject epochs showing horizontal eye movements altogether)

    % set some defaults if not specified
    if exist('RELAX_cfg', 'var')==1
        if isfield(RELAX_cfg, 'OnlyIncludeTaskRelatedEpochs')==0
            RELAX_cfg.OnlyIncludeTaskRelatedEpochs=0; % If this =1, the clean and artifact templates will only include data within 5 seconds of a task trigger, periods outside this will be marked as NaN.
        end
        if isfield(RELAX_cfg, 'HorizontalEyeMovementType')==0
            RELAX_cfg.HorizontalEyeMovementType=2; % 1 to use the IQR method, 2 to use the MAD method for identifying threshold. MAD method works better.
        end
        if isfield(RELAX_cfg, 'HorizontalEyeMovementThreshold')==0
            RELAX_cfg.HorizontalEyeMovementThreshold=2; % MAD deviation from the median that will be marked as horizontal eye movement if both lateral electrodes show activity above this threshold and in the opposite direction for a certain duration (duration set below).
        end
        if isfield(RELAX_cfg, 'HorizontalEyeMovementThresholdIQR')==0
            RELAX_cfg.HorizontalEyeMovementThresholdIQR=1.5; % IQR deviation that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration (duration set below). MAD method works better.
        end
        if isfield(RELAX_cfg, 'HorizontalEyeMovementTimepointsExceedingThreshold')==0
            RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold=25; % Timepoints exceeding horizontal threshold within the test period set below before marked as horizontal eye movement.
        end
        if isfield(RELAX_cfg, 'HorizontalEyeMovementTimepointsExceedingThreshold')==0
            RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow=(2*RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold)-1; % Window duration to test for horizontal eye movement, set to 2x the value above by default.
        end
        if isfield(RELAX_cfg, 'HorizontalEyeMovementFocus')==0
            RELAX_cfg.HorizontalEyeMovementFocus=400; % Buffer window, masking periods earlier and later than the time where horizontal eye movements start exceeding the threshold.
        end       
        if isfield(RELAX_cfg, 'HEOGLeftpattern')==0
            RELAX_cfg.HEOGLeftpattern = ["AF7", "F7", "FT7", "F5", "T7", "FC5", "C5", "TP7", "AF3"]; % sets electrodes to use for horizontal eye movement detection in priority order (if the electrode in position 1 isn't present, the script will check for electrode in position 2, and so on...).
        end
        if isfield(RELAX_cfg, 'HEOGRightpattern')==0
            RELAX_cfg.HEOGRightpattern = ["AF8", "F8","FT8","F6","T8", "FC6", "C6", "TP8", "AF4"]; % sets electrodes to use for horizontal eye movement detection in priority order (if the electrode in position 1 isn't present, the script will check for electrode in position 2, and so on...).
        end
    elseif exist('RELAX_cfg', 'var')==0
        RELAX_cfg.OnlyIncludeTaskRelatedEpochs=0; % If this =1, the clean and artifact templates will only include data within 5 seconds of a task trigger, periods outside this will be marked as NaN.
        RELAX_cfg.HorizontalEyeMovementType=2; % 1 to use the IQR method, 2 to use the MAD method for identifying threshold
        RELAX_cfg.HorizontalEyeMovementThreshold=2; % MAD deviation from the median that will be marked as horizontal eye movement if both lateral electrodes show activity above this threshold and in the opposite direction for a certain duration (duration set below).
        RELAX_cfg.HorizontalEyeMovementThresholdIQR=1.5; % IQR deviation that will be marked as horizontal eye movement if both lateral electrodes show activity above this for a certain duration set below (duration set below). MAD method works better.
        RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold=25; % Datapoints exceeding horizontal threshold within the test period set below before marked as horizontal eye movement.
        RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow=(2*RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold)-1; % Window duration to test for horizontal eye movement, set to 2x the value above by default.
        RELAX_cfg.HorizontalEyeMovementFocus=400; % Buffer window, masking periods earlier and later than the time where horizontal eye movements exceed the threshold.
        RELAX_cfg.HEOGLeftpattern = ["AF7", "F7", "FT7", "F5", "T7", "FC5", "C5", "TP7", "AF3"];
        RELAX_cfg.HEOGRightpattern = ["AF8", "F8","FT8","F6","T8", "FC6", "C6", "TP8", "AF4"]; 
    end
    
%% HORIZONTAL EYE MOVEMENT MARKING BY LATERAL ELECTRODES ON OPPOSITE SIDES OF THE HEAD SHOWING VALUES OF THE OPPOSITE POLARITY, BOTH LARGER THAN MEDIAN +/- X MAD:

        continuousEEG.RELAXProcessing.Details.horizontaleyemask=continuousEEG.RELAXProcessing.Details.NaNsForNonEvents;
        
        % Selecting electrodes to detect horizontal eye movement. These
        % should be lateral electrodes, and preferentially towards the
        % front if available:
        
        % Left side:
        Leftpattern = RELAX_cfg.HEOGLeftpattern;
        for x = 1:size(continuousEEG.chanlocs,2)
            if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,1))==1
                LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,1));  
            end
        end
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,2))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,2));
                end
            end
        end           
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,3))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,3));
                end
            end
        end      
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,4))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,4));
                end
            end
        end       
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,5))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,5));
                end
            end
        end  
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,6))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,6));
                end
            end
        end
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,7))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,7));
                end
            end
        end
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,8))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,8));
                end
            end
        end
        if exist('LeftHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,9))==1
                    LeftHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Leftpattern(1,9));
                end
            end
        end
        
        % Right side:
        Rightpattern = RELAX_cfg.HEOGRightpattern;
        for x = 1:size(continuousEEG.chanlocs,2)
            if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,1))==1
                RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,1));  
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,2))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,2));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,3))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,3));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,4))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,4));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,5))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,5));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,6))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,6));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,7))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,7));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,8))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,8));
                end
            end
        end
        if exist('RightHEyeChannelLocationList','var') ~= 1
            for x = 1:size(continuousEEG.chanlocs,2)
                if strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,9))==1
                    RightHEyeChannelLocationList(1,x) = strcmpi(continuousEEG.chanlocs(x).labels, Rightpattern(1,9));
                end
            end
        end

        % If HEOG electrodes were present, then detect HEOG movements:
        if exist('RightHEyeChannelLocationList','var')==1 && exist('LeftHEyeChannelLocationList','var')==1
            LeftHEyeLocation=find(LeftHEyeChannelLocationList==1);
            RightHEyeLocation=find(RightHEyeChannelLocationList==1);
            if RELAX_cfg.HorizontalEyeMovementType==1% 1 to use the IQR method, 2 to use the MAD method for identifying threshold
                % IQR method:        
                continuousEEG.RELAXProcessing.Details.LeftHEOGIQR = iqr(continuousEEG.data(LeftHEyeLocation,:),2);
                continuousEEG.RELAXProcessing.Details.LeftUpper25 = prctile(continuousEEG.data(LeftHEyeLocation,:),75,2);  
                continuousEEG.RELAXProcessing.Details.LeftLower25 = prctile(continuousEEG.data(LeftHEyeLocation,:),25,2);  
                continuousEEG.RELAXProcessing.Details.HEOGLeftUpperBoundIQR=continuousEEG.RELAXProcessing.Details.LeftUpper25+(continuousEEG.RELAXProcessing.Details.LeftHEOGIQR*RELAX_cfg.HorizontalEyeMovementThresholdIQR);
                continuousEEG.RELAXProcessing.Details.HEOGLeftLowerBoundIQR=continuousEEG.RELAXProcessing.Details.LeftLower25-(continuousEEG.RELAXProcessing.Details.LeftHEOGIQR*RELAX_cfg.HorizontalEyeMovementThresholdIQR);
                continuousEEG.RELAXProcessing.Details.RightHEOGIQR = iqr(continuousEEG.data(RightHEyeLocation,:),2);
                continuousEEG.RELAXProcessing.Details.RightUpper25 = prctile(continuousEEG.data(RightHEyeLocation,:),75,2);  
                continuousEEG.RELAXProcessing.Details.RightLower25 = prctile(continuousEEG.data(RightHEyeLocation,:),25,2);  
                continuousEEG.RELAXProcessing.Details.HEOGRightUpperBoundIQR=continuousEEG.RELAXProcessing.Details.RightUpper25+(continuousEEG.RELAXProcessing.Details.RightHEOGIQR*RELAX_cfg.HorizontalEyeMovementThresholdIQR);
                continuousEEG.RELAXProcessing.Details.HEOGRightLowerBoundIQR=continuousEEG.RELAXProcessing.Details.RightLower25-(continuousEEG.RELAXProcessing.Details.RightHEOGIQR*RELAX_cfg.HorizontalEyeMovementThresholdIQR);
                % If more than the specified number of consecutive datapoints show HEOG affected
                % electrodes on the opposite sides of the head with voltage values
                % above the IQR threshold in opposite directions in those
                % opposite electrodes, mark as HEOG artifact for buffer window period on either
                % side of the datapoint:
                for x=round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample)+1:size(continuousEEG.data,2)-round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample)
                    if ((sum(continuousEEG.data(RightHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            >continuousEEG.RELAXProcessing.Details.HEOGRightUpperBoundIQR)>round(RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold/RELAX_cfg.ms_per_sample))...
                            &&(sum(continuousEEG.data(LeftHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            <continuousEEG.RELAXProcessing.Details.HEOGLeftLowerBoundIQR)>round(RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold/RELAX_cfg.ms_per_sample)))==1
                        continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x-round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample):x+round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample))=1;
                    elseif ((sum(continuousEEG.data(RightHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            <continuousEEG.RELAXProcessing.Details.HEOGRightLowerBoundIQR)>round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            &&(sum(continuousEEG.data(LeftHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            >continuousEEG.RELAXProcessing.Details.HEOGRightUpperBoundIQR)>round(RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold/RELAX_cfg.ms_per_sample)))==1
                        continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x-round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample):x+round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample))=1;
                    end
                end
            end

            % MAD method:
            if RELAX_cfg.HorizontalEyeMovementType==2% 1 to use the IQR method, 2 to use the MAD method for identifying threshold
                % Work out the median voltage for the HEOG affected electrodes, and
                % calculate +/- X MAD from the median
                continuousEEG.RELAXProcessing.Details.HEOGLeftMedian = median(continuousEEG.data(LeftHEyeLocation,:));
                continuousEEG.RELAXProcessing.Details.LeftMAD = mad(continuousEEG.data(LeftHEyeLocation,:));   
                continuousEEG.RELAXProcessing.Details.HEOGLeftMedianPlus2MAD=continuousEEG.RELAXProcessing.Details.HEOGLeftMedian+(continuousEEG.RELAXProcessing.Details.LeftMAD*RELAX_cfg.HorizontalEyeMovementThreshold);
                continuousEEG.RELAXProcessing.Details.HEOGLeftMedianMinus2MAD=continuousEEG.RELAXProcessing.Details.HEOGLeftMedian-(continuousEEG.RELAXProcessing.Details.LeftMAD*RELAX_cfg.HorizontalEyeMovementThreshold);
                continuousEEG.RELAXProcessing.Details.HEOGRightMedian = median(continuousEEG.data(RightHEyeLocation,:));
                continuousEEG.RELAXProcessing.Details.HEOGRightMAD = mad(continuousEEG.data(RightHEyeLocation,:));   
                continuousEEG.RELAXProcessing.Details.HEOGRightMedianPlus2MAD=continuousEEG.RELAXProcessing.Details.HEOGRightMedian+(continuousEEG.RELAXProcessing.Details.HEOGRightMAD*RELAX_cfg.HorizontalEyeMovementThreshold);
                continuousEEG.RELAXProcessing.Details.HEOGRightMedianMinus2MAD=continuousEEG.RELAXProcessing.Details.HEOGRightMedian-(continuousEEG.RELAXProcessing.Details.HEOGRightMAD*RELAX_cfg.HorizontalEyeMovementThreshold);
                % If more than the specified number of consecutive datapoints show HEOG affected
                % electrodes on the opposite sides of the head with voltage values
                % above the MAD threshold from the median in opposite directions in those
                % opposite electrodes, mark as HEOG artifact for buffer window period on either
                % side of the datapoint:
                for x=round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample)+1:size(continuousEEG.data,2)-round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample)
                    if ((sum(continuousEEG.data(RightHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            >continuousEEG.RELAXProcessing.Details.HEOGRightMedianPlus2MAD)>round(RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold/RELAX_cfg.ms_per_sample))...
                            &&(sum(continuousEEG.data(LeftHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            <continuousEEG.RELAXProcessing.Details.HEOGLeftMedianMinus2MAD)>round(RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold/RELAX_cfg.ms_per_sample)))==1
                        continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x-round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample):x+round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample))=1;
                    elseif ((sum(continuousEEG.data(RightHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            <continuousEEG.RELAXProcessing.Details.HEOGRightMedianMinus2MAD)>round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            &&(sum(continuousEEG.data(LeftHEyeLocation,x:x+round(RELAX_cfg.HorizontalEyeMovementTimepointsTestWindow/RELAX_cfg.ms_per_sample))...
                            >continuousEEG.RELAXProcessing.Details.HEOGRightMedianPlus2MAD)>round(RELAX_cfg.HorizontalEyeMovementTimepointsExceedingThreshold/RELAX_cfg.ms_per_sample)))==1
                        continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x-round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample):x+round(RELAX_cfg.HorizontalEyeMovementFocus/RELAX_cfg.ms_per_sample))=1;
                    end
                end
            end

            % Insert NaNs into masks for non-task periods if that option is
            % selected
            for x=1:size(continuousEEG.RELAXProcessing.Details.NaNsForNonEvents,2)
                if  isnan(continuousEEG.RELAXProcessing.Details.NaNsForNonEvents(1,x))==1
                    continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x)=NaN;
                end
            end 
            % Insert hotizontal eye movement events into EEG trace:
            EEGeventlength=length(continuousEEG.event);
            CumulatingNumberOfHorizontalEyeMovements=0;
            for x=1:size(continuousEEG.RELAXProcessing.Details.horizontaleyemask,2)-1
                if continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x)-continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x+1)==-1
                    CumulatingNumberOfHorizontalEyeMovements=CumulatingNumberOfHorizontalEyeMovements+1;                
                    continuousEEG.event(CumulatingNumberOfHorizontalEyeMovements+EEGeventlength).type='StartHorizontalEyeMovement';
                    continuousEEG.event(CumulatingNumberOfHorizontalEyeMovements+EEGeventlength).latency=x;
                elseif continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x)-continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,x+1)==1
                    CumulatingNumberOfHorizontalEyeMovements=CumulatingNumberOfHorizontalEyeMovements+1;
                    continuousEEG.event(CumulatingNumberOfHorizontalEyeMovements+EEGeventlength).type='EndHorizontalEyeMovement';
                    continuousEEG.event(CumulatingNumberOfHorizontalEyeMovements+EEGeventlength).latency=x;
                end
            end
            % Add HEOG arifacts to overall noise mask
            for e=1:size(continuousEEG.RELAXProcessing.Details.NoiseMaskFullLength,2)
                if (continuousEEG.RELAXProcessing.Details.horizontaleyemask(1,e)==1)
                    continuousEEG.RELAXProcessing.Details.NoiseMaskFullLength(1,e)=1;
                end
            end
        end
        % Calculate proportion of HEOG affected data:
        continuousEEG.RELAXProcessing.ProportionMarkedFromHEOG=mean(continuousEEG.RELAXProcessing.Details.horizontaleyemask,'omitnan');
end