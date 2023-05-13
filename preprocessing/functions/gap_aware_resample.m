% gap_aware_resample resample a nonuniform signal while filling gaps with autoregression 
%   
%   fs: sampling rate to resample the nonuniform signal at
%   
%   threshold: the number of missed samplesat which to switch to
%   autoregression. If threshold is 2 and fs*2 seconds pass in the
%   nonuniform signal without another sample, autorgression is used to fill
%   the gaps before switching back to resamling when the next sample in the
%   nonuniform signal arrives.
function X = gap_aware_resample(X, t, fs, threshold)
    [segment_values, segment_times] = segment_signal(X, t, fs, threshold);
    X = segmented_resampling(segment_values, segment_times, fs);
end

function [segment_values, segment_times] = segment_signal(X, t, fs, threshold)
    segment_ends = [find(diff(t) > threshold/fs); 0];
    segment_ends(end) = length(t);

    i_start = 1; % first segment starts with the first index
    segment_times = cell(length(segment_ends),1);
    segment_values = cell(length(segment_ends),1);
    for index = 1:numel(segment_ends)
        i_end = segment_ends(index);
        % TODO: correct treatment of t_0
        times = t(i_start:i_end);
        times(1) = floor(fs*times(1))/fs;
        [Y, t_Y] = resample(X(i_start:i_end, :), times, fs, 'pchip');
        
        segment_times{index} = t_Y;
        segment_values{index} = Y;
    
        i_start = i_end + 1;
    end
end
function X = segmented_resampling(segment_values, segment_times, fs)
    % fill gaps in segments with NaN
    X = nan(round(fs*(segment_times{end}(end) - segment_times{1}(1))) + 1, 3);
    
    for index = 1:numel(segment_times)
        i_start = round(segment_times{index}(1)*fs) + 1;
        i_end = round(segment_times{index}(end)*fs) + 1;
        X(i_start:i_end, :) =   segment_values{index};
    end
    
    % use autoreg to replace NaN values
    autoreg_len = 150;
    autoreg_order = 150;
    X = fillgaps(X, autoreg_len, autoreg_order);
end