% TODO: current, only handles overlap times of multiples of 1 second
function [s2, t] = nustft(x, t, fs, window_time, window_time_shift, max_window_error)
    nufft_length = 100;
    [starts, ends, targets] = partition_time_sequence(t, window_time, window_time * max_window_error, nufft_length);
    
    [~, x_dim] = size(x);
    window_indices = 1:window_time_shift:numel(starts);
    s2 = zeros(length(targets), nufft_length/2, x_dim);
    parfor i = window_indices
        n = ends(i) - starts(i) + 1;
        f = double(0:(n/2-1))/double(n)*fs;
        Y = nufft(x(starts(i):ends(i), :), t(starts(i):ends(i)), f, 1);
        bin_size = idivide(ends(i)-starts(i)+1, nufft_length);
        if(bin_size == 1)
            s2(i,:,:) = abs(Y .* conj(Y));
        else
            s2(i,:,:) = squeeze(sum(reshape(abs(Y.*conj(Y)),bin_size, [], 3), 1))./double(bin_size);
        end
    end
    t = targets;
end

function [starts, ends, targets] = partition_time_sequence(t, t_win, t_err, nufft_length)
    t_0 = floor(t(1));
    t_f = floor(t(end));
    start_of_sec = zeros(uint32(t_f - t_0 + 1), 1, 'uint32');
    starts = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    ends = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    targets = (t_0+5):t_f;
    % get initial value for start times
    for i = 1:numel(start_of_sec)
        start_of_sec(i) = find(t>=t_0+i-1, 1, 'first');
    end
    % fix sequences to be a multiple of 100 in length
    parfor i = (t_win+1):(numel(start_of_sec))
        i_start = start_of_sec(i-t_win);
        i_end = start_of_sec(i) - 1;
        seq_len = i_end - i_start + 1;
        mod_nufft_len = mod(seq_len, nufft_length);
        error = nufft_length-mod_nufft_len;
        if(mod_nufft_len == 0)
            % the window lined up perfectly
        elseif(seq_len < nufft_length)
            % the window was too short, extend it
            if(i_start - error <= 0)
                % cannot extend backwards, extend forwards. 
                % 
                % This complicates analysis, but this corner case can happen 
                % only once per activity so it should hopefully not be major.
                % TODO: if the data towards the end gets discard, might as well
                % shift everything forward such that there is enough data
                % before the first window that backwards extension works
                i_end = i_end + error;
            elseif(i_start > i_end)
                % This case occurs when an entire time window passes
                % without a new measurement being acquired. This is handled
                % by just leaving the output steady until new input is
                % received.
                % TODO: is there a more clever way to handle this?
                i_end = i_start;
                i_start = i_end - 100 + 1;
            else
                % can extend backwards
                i_start = i_start - error;
            end
        elseif(t(i_end) - t(i_start+mod_nufft_len) > t_win - t_err) % if this guard evaluated, must be too long
            % Window can be shortened without making the time covered by the
            % window too short.
            i_start = i_start + mod_nufft_len;
        else
            % the window is too long, but must be further lengthened to cover
            % the minimum desired timespan while having the samples be a factor
            % of 100
            if(i_start - error <= 0)
                % cannot extend backwards, shrink forwards. Similar drawbacks as
                % before hold.
                i_end = i_end - mod_nufft_len;
            else
                % can extend backwards
                i_start = i_start - error;
            end
        end
        starts(i-t_win) = i_start;
        ends(i-t_win) = i_end;
    end
end