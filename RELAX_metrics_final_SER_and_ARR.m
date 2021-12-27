function [continuousEEG, d] = RELAX_metrics_final_SER_and_ARR(rawEEG, continuousEEG)

    %% This function uses all artifact templates from the MWF cleaning steps to compute the SER and ARR cleaning efficacy metrics

    % Remove the same periods in the raw EEG as have been removed in the
    % cleaned EEG:
    rawEEG = eeg_eegrej( rawEEG, continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion);
    % And also from the masks, firstly by setting up a full EEG.data length
    % template with 1's to denote periods that were excluded as extreme
    % outliers:
    reject = zeros(1,size(rawEEG.data,2));
    for i=1:size(continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion,1)
        reject(continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion(i,1):continuousEEG.RELAX.ExtremelyBadPeriodsForDeletion(i,2)) = 1;
    end
    % Then by using that template to remove periods from the masks, and
    % combining the masks together so all artifacts from all masks are
    % present in a single template:
    if isfield(continuousEEG.RELAX,'NoiseMaskFullLengthR1')
        NoiseMaskFullLengthR1 = continuousEEG.RELAX.NoiseMaskFullLengthR1;
        NoiseMaskFullLengthR1(:,reject == 1) = [];
        NoiseMaskFullLengthAll=NoiseMaskFullLengthR1;
    end
    if isfield(continuousEEG.RELAX,'NoiseMaskFullLengthR2')
        NoiseMaskFullLengthR2 = continuousEEG.RELAX.NoiseMaskFullLengthR2;
        NoiseMaskFullLengthR2(:,reject == 1) = [];
        if isfield(continuousEEG.RELAX,'NoiseMaskFullLengthR1')==0
            NoiseMaskFullLengthAll=NoiseMaskFullLengthR2;
        end
        NoiseMaskFullLengthAll(NoiseMaskFullLengthR2==1)=1;
    end
    if isfield(continuousEEG.RELAX,'NoiseMaskFullLengthR3')
        NoiseMaskFullLengthR3 = continuousEEG.RELAX.NoiseMaskFullLengthR3;
        NoiseMaskFullLengthR3(:,reject == 1) = [];
        if exist('NoiseMaskFullLengthAll','var')==0
            NoiseMaskFullLengthAll=NoiseMaskFullLengthR3;
        end
        NoiseMaskFullLengthAll(NoiseMaskFullLengthR3==1)=1;
    end

    % Average re-reference the raw data first so the SER and ARR are
    % computed from average re-referenced raw and cleaned data:
    [averageref_rawEEG] = RELAX_average_rereference(rawEEG);

    % Calculate the artifact that has been removed:
    d = averageref_rawEEG.data-continuousEEG.data;

    % Calculate SER and ARR for each type of artifact in the MWF masks:
    [continuousEEG.RELAX_Metrics.Cleaned.All_SER, continuousEEG.RELAX_Metrics.Cleaned.All_ARR] = mwf_performance(averageref_rawEEG.data, d, NoiseMaskFullLengthAll);

end