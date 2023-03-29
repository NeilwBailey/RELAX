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

%% RELAX_ICA_subtract:

function [EEG] = RELAX_ICA_subtract(EEG,RELAX_cfg)
    fastica_symm_Didnt_Converge=[0 0 0];
    % run ICA:
    if strcmp(RELAX_cfg.ICA_method,'runica')
        [OUTEEG, ~] = pop_runica_nwb(EEG, 'extended',1,'interupt','on'); %runica for parametric, default extended for finding subgaussian distributions
        W = OUTEEG.icaweights*OUTEEG.icasphere;
        A = inv(W);
        OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
        IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'cudaica')
        [OUTEEG, ~] = pop_runica_nwb(EEG, 'cudaica', 'extended',1); %runica for parametric, default extended for finding subgaussian distributions
        W = OUTEEG.icaweights*OUTEEG.icasphere;
        A = inv(W);
        OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
        IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'fastica_symm')
        % NWB - the following lines repeat fastica_symm up to 3 times in
        % the case of non-convergence, then switches to fastica_defl to
        % ensure ICA convergence (as cleaning as adversely affected by
        % non-convergence issues).
         [OUTEEG, ~, NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
         fastica_symm_Didnt_Converge(1,1)=NonConvergence;
         if NonConvergence==1
             [OUTEEG, ~, NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
             fastica_symm_Didnt_Converge(1,2)=NonConvergence;
         end
         if NonConvergence==1
             [OUTEEG, ~, NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
             fastica_symm_Didnt_Converge(1,3)=NonConvergence;
         end
         if NonConvergence==1
             OUTEEG = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'defl', 'g', 'tanh', 'stabilization', 'on');
         end
         W = OUTEEG.icaweights*OUTEEG.icasphere;
         A = inv(W);
         OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
         IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'fastica_defl')
         OUTEEG = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'defl', 'g', 'tanh', 'stabilization', 'on');
         W = OUTEEG.icaweights*OUTEEG.icasphere;
         A = inv(W);
         OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
         IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'amica')
        OUTEEG=EEG;
        % You'll need to install amica15 first, and in the folder that you
        % specify in the line below (with no spaces in any part of the folder or subfolders):
        % You can download amica15 via EEGLAB    
        cd('D:\Data\AMICA1.5.1');    
        % define parameters
        numprocs = 1;       % # of nodes (default = 1)
        max_threads = 4;    % # of threads per node
        num_models = 1;     % # of models of mixture ICA 
        max_iter = 2000;    % max number of learning steps 
        mkdir('D:\Data\AMICAtmp');
        outdir = 'D:\Data\AMICAtmp\';
        % Run AMICA:    
        [OUTEEG.icaweights, OUTEEG.icasphere, ~] = runamica15(OUTEEG.data, 'num_chans', EEG.nbchan, 'num_models',num_models,'outdir',outdir,'numprocs', numprocs, 'max_threads', max_threads, 'max_iter',max_iter,'pcakeep', EEG.nbchan, 'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1);
        W = OUTEEG.icaweights*OUTEEG.icasphere;
        A = inv(W);
        OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
        IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    elseif strcmp(RELAX_cfg.ICA_method,'picard') % RELAX 1.1.3 - NWB added to allow PICARD-O to be run using default settings
        [OUTEEG, ~] = pop_runica_nwb(EEG, 'picard', 'mode','ortho','tol',1e-6,'maxiter',500); % run picard
        W = OUTEEG.icaweights*OUTEEG.icasphere;
        A = inv(W);
        OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
        IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
    end
    % Use ICLabel to identify artifactual components, so that wICA can be
    % performed on them only:
    % (https://github.com/sccn/ICLabel)
    OUTEEG = iclabel(OUTEEG);
    [~, I]=max(OUTEEG.etc.ic_classification.ICLabel.classifications, [], 2);
    ICsMostLikelyNotBrain=(I>1)';
    
    ArtifactICList=find(ICsMostLikelyNotBrain==1);
    EEG = pop_subcomp( OUTEEG,ArtifactICList, 0); % removes artifactual components by subtraction
    EEG = eeg_checkset( EEG );
    
    EEG.RELAXProcessing_ICA.fastica_symm_Didnt_Converge=fastica_symm_Didnt_Converge; % Tracks whether fastica_symm showed convergence issues (1) or not (0), and how many non-convergences. If 3 non-convergences, then fastica_defl was implemented.
    
    % Check if data might have been too short for effective ICA, using Makoto's rule
    % of thumb that ICA requires data length of ((number of channels)^2)*30
    % if data were sampled at 250 Hz (assuming that higher sampling
    % rates require the same time duration of data as low sampling rates,
    % so 1000Hz sampling rates require ((number of channels)^2)*120)
    % (https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code)
    ms_per_sample=(1000/EEG.srate);
    if ((EEG.nbchan^2)*(120/ms_per_sample))>EEG.pnts
        EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA='yes';
    else
        EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA='no';
    end
    
    if strcmp (EEG.RELAXProcessing_ICA.DataMaybeTooShortForValidICA,'yes')
        warning('Data may have been shorter than recommended for effective ICA decomposition')
    end
    
    EEG.RELAXProcessing_ICA.Proportion_artifactICs_reduced_by_ICA=mean(ICsMostLikelyNotBrain);
    
    if strcmp(RELAX_cfg.Report_all_ICA_info,'yes')
    
        EEG.RELAXProcessing_ICA.ProportionICs_was_Brain=sum(I==1)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_ICA.ProportionICs_was_Muscle=sum(I==2)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_ICA.ProportionICs_was_Eye=sum(I==3)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_ICA.ProportionICs_was_Heart=sum(I==4)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_ICA.ProportionICs_was_LineNoise=sum(I==5)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_ICA.ProportionICs_was_ChannelNoise=sum(I==6)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);
        EEG.RELAXProcessing_ICA.ProportionICs_was_Other=sum(I==7)/size(OUTEEG.etc.ic_classification.ICLabel.classifications,1);    

        ICsMostLikelyBrain=(I==1)';
        ICsMostLikelyMuscle=(I==2)';
        ICsMostLikelyEye=(I==3)';
        ICsMostLikelyHeart=(I==4)';
        ICsMostLikelyLineNoise=(I==5)';
        ICsMostLikelyChannelNoise=(I==6)';
        ICsMostLikelyOther=(I==7)';

        for x=1:size(OUTEEG.etc.ic_classification.ICLabel.classifications,1)
            [~, varianceWav(x)] =compvar(OUTEEG.data, OUTEEG.icaact, OUTEEG.icawinv, x);
        end

        BrainVariance=sum(abs(varianceWav(ICsMostLikelyBrain)));
        ArtifactVariance=sum(abs(varianceWav(~ICsMostLikelyBrain)));
        EEG.RELAXProcessing_wICA.ProportionVariance_was_BrainICs=(BrainVariance/(BrainVariance+ArtifactVariance));

        MuscleVariance=sum(abs(varianceWav(ICsMostLikelyMuscle)));
        EyeVariance=sum(abs(varianceWav(ICsMostLikelyEye)));
        HeartVariance=sum(abs(varianceWav(ICsMostLikelyHeart)));
        LineNoiseVariance=sum(abs(varianceWav(ICsMostLikelyLineNoise)));
        ChannelNoiseVariance=sum(abs(varianceWav(ICsMostLikelyChannelNoise)));
        OtherVariance=sum(abs(varianceWav(ICsMostLikelyOther)));

        EEG.RELAXProcessing_ICA.ProportionVariance_was_MuscleICs=(MuscleVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_ICA.ProportionVariance_was_EyeICs=(EyeVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_ICA.ProportionVariance_was_HeartICs=(HeartVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_ICA.ProportionVariance_was_LineNoiseICs=(LineNoiseVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_ICA.ProportionVariance_was_ChannelNoiseICs=(ChannelNoiseVariance/(BrainVariance+ArtifactVariance));
        EEG.RELAXProcessing_ICA.ProportionVariance_was_OtherICs=(OtherVariance/(BrainVariance+ArtifactVariance));
    
    end
end