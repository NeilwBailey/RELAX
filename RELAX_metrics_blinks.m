function [continuousEEG, epochedEEG, BlinkMetricProblem] = RELAX_metrics_blinks(continuousEEG, epochedEEG)
    % Metric suggested by: Robbins, K. A., Touryan, J., Mullen, T., Kothe, C., & Bigdely-Shamlo, N. (2020). How Sensitive are EEG Results to Preprocessing Methods: A Benchmarking Study. IEEE Transactions on Neural Systems and Rehabilitation Engineering, 28(5), 1081-1090.
    BlinkMetricProblem=[];
    try
        [EEG] = RELAX_average_rereference(continuousEEG); % This is only implemented if the data hasn't already been average re-referenced by RELAX
        EEG = pop_interp(EEG, EEG.allchan, 'spherical');
        % Epoch around the blink max:
        blinkEEG = pop_epoch( EEG, { 'EyeBlinkMax'}, [-2 2.001], 'epochinfo', 'yes');
        % Exclude epochs that contain more than one blink (as this can reduce
        % the blink amplitude ratio since the baseline period could contain
        % increased blink related amplitudes:
        for e=1:size(blinkEEG.epoch,2)
            EpochsWithMultipleBlinks(e)=sum(contains(string(blinkEEG.epoch(e).eventtype),"EyeBlinkMax"));     
        end
        blinkEEG = pop_select(blinkEEG, 'notrial', find(EpochsWithMultipleBlinks>1));
        % baseline correct data:
        blinkEEG.data=blinkEEG.data-mean(blinkEEG.data(:,[1:501,3501:4001],:),2);
        % convert to absolute values:
        absolutevaluesblink=abs(blinkEEG.data);
        % calculate 
        for e=1:size(blinkEEG.epoch,2)
            BlinkAmplitudeRatioAllEpochs(1:size(absolutevaluesblink,1),1,e)=mean(absolutevaluesblink(:,1500:2500,e),2)./mean(absolutevaluesblink(:,[1:501,3501:4001],e),2);
        end
        BlinkAmplitudeRatio(1:size(BlinkAmplitudeRatioAllEpochs,1),1)=mean(BlinkAmplitudeRatioAllEpochs,3);
        if continuousEEG.RELAX.Data_has_been_cleaned==0
            epochedEEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio=BlinkAmplitudeRatio;
            continuousEEG.RELAX_Metrics.Raw.BlinkAmplitudeRatio=BlinkAmplitudeRatio;
        end
        if continuousEEG.RELAX.Data_has_been_cleaned==1
            epochedEEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio=BlinkAmplitudeRatio;
            continuousEEG.RELAX_Metrics.Cleaned.BlinkAmplitudeRatio=BlinkAmplitudeRatio;
        end
    catch
        BlinkMetricProblem='no blink free baseline periods available';
    end
end