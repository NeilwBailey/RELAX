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
