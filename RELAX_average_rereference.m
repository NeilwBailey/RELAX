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

%% RELAX_average_rereference:

function [EEG] = RELAX_average_rereference(EEG)
    EEGaveragerereference=0;
    if exist('EEG','var')
        if isfield(EEG, 'RELAX')==0
            EEGaveragerereference=1;
        end
        if isfield(EEG, 'RELAX')==1
            if isfield(EEG.RELAX, 'Data_has_been_averagerereferenced')==0 || EEG.RELAX.Data_has_been_averagerereferenced~=1 
                EEGaveragerereference=1;
            end
        end
        if EEGaveragerereference==1
            %% Re-reference EEG data to the common average, then reject bad channels again as interpolated channels lead to rank deficiency:
            EEG = pop_interp(EEG, EEG.allchan, 'spherical');
            % Apply average reference after adding initial reference electrode
            % back into the data
            EEG.nbchan = EEG.nbchan+1;
            EEG.data(end+1,:,:) = zeros(1, EEG.pnts, EEG.trials);
            EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
            EEG = pop_reref(EEG, []);
            EEG = pop_select( EEG,'nochannel',{'initialReference'});
            EEG.RELAX.Data_has_been_averagerereferenced=1;
            EEG=pop_select(EEG,'channel',EEG.RELAX.ListOfChannelsAfterRejections);
        end
    end
end
