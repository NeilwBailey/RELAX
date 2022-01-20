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

%% eegplugin_RELAX:
% Add RELAX to the EEGLAB gui:
function vers = eegplugin_RELAX(fig, try_strings, catch_strings)

    vers = 'RELAX 0.91';
    if nargin < 2
        error('eegplugin_RELAX requires 2 arguments');
    end

    if ispc      % windows
            wfactor1 = 1.20;
            wfactor2 = 1.21;
    elseif ismac % Mac OSX
            wfactor1 = 1.45;
            wfactor2 = 1.46;
    else
            wfactor1 = 1.30;
            wfactor2 = 1.31;
    end
    posmainfig = get(gcf,'Position');
    hframe     = findobj('parent', gcf,'tag','Frame1');
    posframe   = get(hframe,'position');
    set(gcf,'position', [posmainfig(1:2) posmainfig(3)*wfactor1 posmainfig(4)]);
    set(hframe,'position', [posframe(1:2) posframe(3)*wfactor2 posframe(4)]);

    menuRELAX = findobj(fig,'tag','EEGLAB');   % At EEGLAB Main Menu

    submenu = uimenu( menuRELAX,'Label','RELAX','separator','on','tag','RELAX','userdata','startup:on;continuous:on;epoch:on;study:on;erpset:on');
    
    % menu callbacks
        % --------------
    comProcessData = [try_strings.no_check...
        '[RELAX_cfg, FileNumber, CleanedMetrics, RawMetrics, RELAXProcessingRoundOneAllParticipants, RELAXProcessingRoundTwoAllParticipants, RELAXProcessing_wICA_AllParticipants, RELAXProcessingRoundThreeAllParticipants, RELAX_issues_to_check, RELAXProcessingExtremeRejectionsAllParticipants] = pop_RELAX();'...
        catch_strings.add_to_hist];
     % create menus
        % -------------------------   
    uimenu( submenu, 'Label', 'Preprocess EEG Data'  , 'CallBack', comProcessData);
    
end
