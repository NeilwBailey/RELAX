function [EEG,wIC,A,W,IC] = RELAX_wICA_on_ICLabel_artifacts(EEG,varargin) % NWB altered to output EEGLAB struct after wICA removal of artifacts identified by ICLabel
    %--------------- function [wIC,A,W] = wICA(data,varargin) -----------------
    %
    % Performs ICA on data matrix (row vector) and subsequent wavelet
    % thresholding to remove low-amplitude activity from the computed ICs.
    % This is useful for extracting artifact-only ICs in EEG (for example), and
    % then subtracting the artifact-reconstruction from the original data. 
    %
    %               >>> INPUTS >>>
    % Required: 
    %   data = data matrix in row format
    % Optional:
    %   type = "fastica" or "radical"...two different ICA algorithms based on
    %       entropy. "fastica" (default) is parametric, "radical" is nonparametric.
    %   mult = threshold multiplier...multiplies the computed threshold from
    %       "ddencmp" by this number. Higher thresh multipliers = less
    %       "background" (or low amp. signal) is kept in the wICs.
    %   plotting = 1 or 0. If 1, plots wIC vs. non-wavelet thresholded ICs
    %   Fs = sampling rate, (for plotting...default = 1);
    %   L = level set for stationary wavelet transform. Higher levels give
    %       better frequency resolution, but less temporal resolution. 
    %       Default = 5
    %   wavename = wavelet family to use. type "wavenames" to see a list of
    %       possible wavelets. (default = "coif5");
    %
    %               <<< OUTPUTS <<<
    %   wIC = wavelet-thresholded ICs
    %   A = mixing matrix (inv(W)) (optional)
    %   W = demixing matrix (inv(A)) (optional)
    %   IC = non-wavelet ICs (optional)
    %   
    %       * you can reconstruct the artifact-only signals as:
    %               artifacts = A*wIC;
    %       - upon reconstruction, you can then subtract the artifacts from your
    %       original data set to remove artifacts, for instance.
    %
    % Example:
    %  n = rand(10,1000);
    %  a = [zeros(1,400),[.5,.8,1,2,2.4,2.5,3.5,5,6.3,6,4,3.2,3,1.7,1,-.6,-2.2,-4,-3.6,-3,-1,0],zeros(1,578)];
    %  data = n + linspace(0,2,10)'*a;
    %  [wIC,A] = wICA(data,[],5,1);
    %  ahat = A*wIC;
    %  nhat = data-ahat;
    %  err = sum(sqrt((nhat-n).^2));

    % By JMS, 11/10/2015
    %---------------------------------------------------------------------------------------

    % check inputs
    if nargin>1 && ~isempty(varargin{1})
    type=varargin{1}; else type='runica';end
    if nargin>2 && ~isempty(varargin{2})
    mult=varargin{2};else mult=1;end
    if nargin>3 && ~isempty(varargin{3})
    plotting=varargin{3}; else plotting=0;end
    if nargin>4 && ~isempty(varargin{4})
    Fs=varargin{4};else Fs=1;end
    if nargin>5 && ~isempty(varargin{5})
    L=varargin{5}; else L=5;end
    if nargin>6 && ~isempty(varargin{6})
    wavename=varargin{6}; else wavename='coif5';end

    fastica_symm_Didnt_Converge=[0 0 0]; % NWB addition to track whether fastica_symm doesn't converge

    % run ICA using "runica" or "radical"
    if strcmp(type,'runica')
        [OUTEEG, com] = pop_runica(EEG, 'extended',1,'interupt','on'); %runica for parametric, default extended for finding subgaussian distributions
        W = OUTEEG.icaweights*OUTEEG.icasphere;
        A = inv(W);
        %% NWB added section to ensure ICA details are updated in EEGLAB struct
        OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
        %% NWB addition ended
        IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
        %com = pop_export(OUTEEG,'ICactivationmatrix','ica','on','elec','off','time','off','precision',4);
        %IC = ICactivationmatrix;
    elseif strcmp(type,'radical')
        [IC,W] = radical(data); % radical ICA for non-parametric
        A = inv(W);
        %% NWB added section to enable cudaICA / fastICA to be used as an ICA option to enable quicker ICA computation, or AMICA for theoretically optimal performance
    elseif strcmp(type,'cudaica')
        [OUTEEG, com] = pop_runica(EEG, 'cudaica', 'extended',1); %runica for parametric, default extended for finding subgaussian distributions
        W = OUTEEG.icaweights*OUTEEG.icasphere;
        A = inv(W);
        OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
        IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
        %com = pop_export(OUTEEG,'ICactivationmatrix','ica','on','elec','off','time','off','precision',4);
        %IC = ICactivationmatrix;
    elseif strcmp(type,'fastica_symm')
        % NWB - the following lines repeat fastica_symm up to 3 times in
        % the case of non-convergence, then switches to fastica_defl to
        % ensure ICA convergence (as cleaning as adversely affected by
        % non-convergence issues).
         NonConvergence=0;
         [OUTEEG com NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
         fastica_symm_Didnt_Converge(1,1)=NonConvergence;
         if NonConvergence==1
             [OUTEEG com NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
             fastica_symm_Didnt_Converge(1,2)=NonConvergence;
         end
         if NonConvergence==1
             [OUTEEG com NonConvergence] = pop_runica_nwb( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'symm', 'g', 'tanh', 'stabilization', 'on');
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
        %com = pop_export(OUTEEG,'ICactivationmatrix','ica','on','elec','off','time','off','precision',4);
        %IC = ICactivationmatrix;
    elseif strcmp(type,'fastica_defl')
         OUTEEG = pop_runica( EEG, 'icatype', 'fastica','numOfIC', EEG.nbchan, 'approach', 'defl', 'g', 'tanh', 'stabilization', 'on');
         W = OUTEEG.icaweights*OUTEEG.icasphere;
         A = inv(W);
         OUTEEG = eeg_checkset(OUTEEG, 'ica'); 
        if isempty(OUTEEG.icaact)==1
            OUTEEG.icaact = (OUTEEG.icaweights*OUTEEG.icasphere)*OUTEEG.data(OUTEEG.icachansind,:);      
            OUTEEG.icaact = reshape( OUTEEG.icaact, size(OUTEEG.icaact,1), OUTEEG.pnts, OUTEEG.trials);
        end
         IC=reshape(OUTEEG.icaact, size(OUTEEG.icaact,1), []);
        %com = pop_export(OUTEEG,'ICactivationmatrix','ica','on','elec','off','time','off','precision',4);
        %IC = ICactivationmatrix;
    elseif strcmp(type,'amica')
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
        outdir = [ 'D:\Data\AMICAtmp\' ];
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
    end
    %% NWB added section to identify artifactual ICs with ICLabel:
    % Use ICLabel to identify artifactual components, so that wICA can be
    % performed on them only:
    % (https://github.com/sccn/ICLabel)
    OUTEEG = iclabel(OUTEEG);
    [~, I]=max(OUTEEG.etc.ic_classification.ICLabel.classifications, [], 2);
    ICsMostLikelyNotBrain=(I>1)';
    %% NWB addition ended

    % padding data for proper wavelet transform...data must be divisible by
    % 2^L, where L = level set for the stationary wavelet transform
    modulus = mod(size(IC,2),2^L); %2^level (level for wavelet)
    if modulus ~=0
        extra = zeros(1,(2^L)-modulus);
    else
        extra = [];
    end

    % loop through ICs and perform wavelet thresholding
    %% Neil Bailey (2020) edited - only do so on artifact components identified by ICLabel
    disp('Performing wavelet thresholding');
    for s = 1:size(IC,1)
        if ICsMostLikelyNotBrain(s)==1 % NWB added to perform this only on artifacts identified by ICLabel
            if ~isempty(extra)
                sig = [IC(s,:),extra]; % pad with zeros
            else
                sig = IC(s,:);
            end
            [thresh,sorh,~] = ddencmp('den','wv',sig); % get automatic threshold value
            thresh = thresh*mult; % multiply threshold by scalar
            swc = swt(sig,L,wavename); % use stationary wavelet transform (SWT) to wavelet transform the ICs
            Y = wthresh(swc,sorh,thresh); % threshold the wavelet to remove small values
            wIC(s,:) = iswt(Y,wavename); % perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
            clear y sig thresh sorh swc 
        end
    end
    %% NWB added section to pad non-artifact components with 0s in the same way that the artifact components were padded:
    if sum(ICsMostLikelyNotBrain)==0
        wIC(1,:)=zeros(1,size(EEG.data,2));
        wIC = [wIC(:,:),extra]; % pad with zeros
    end
    for s = 1:size(IC,1)
        if ICsMostLikelyNotBrain(s)==0
            wIC(s,:)=zeros(1,size(wIC,2));
        end
    end
    %% NWB addition ended

    % remove extra padding
    if ~isempty(extra)
        wIC = wIC(:,1:end-numel(extra));
    end

    % plot the ICs vs. wICs
    if plotting>0
        disp('Plotting');
        subplot(3,1,1);
            multisignalplot(IC,Fs,'r');
            title('ICs');
        subplot(3,1,2);
            multisignalplot(wIC,Fs,'r');
            title('wICs')
        subplot(3,1,3);
            multisignalplot(IC-wIC,Fs,'r');
            title('Difference (IC - wIC)');
    end
    %% NWB added section to remove wICA artifact and reconstruct data within this function (rather than in the main script):
    artifacts = A*wIC;
    %reshape EEG signal from EEGlab format to channelsxsamples format
    EEG2D=reshape(EEG.data, size(EEG.data,1), []);
    %subtract out wavelet artifact signal from EEG signal
    wavcleanEEG=EEG2D-artifacts;
    EEG.data = wavcleanEEG;
    EEG.RELAXProcessing_wICA.fastica_symm_Didnt_Converge=fastica_symm_Didnt_Converge; % Tracks whether fastica_symm showed convergence issues (1) or not (0), and how many non-convergences. If 3 non-convergences, then fastica_defl was implemented.
    %% NWB addition ended

end