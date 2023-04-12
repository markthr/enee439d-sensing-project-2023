%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";

window_time = 5; % how often to process a new window in seconds
max_window_error = 0.1; % the maximum percent a window can be too short by 
window_time_shift = 1;
nufft_length = 100;
% Setting window_memory = 0 indicates that only the current window should be included for processing.

% get letter for the 18 activities, A-S but without N
activities = char([1:13, 15:19] + 'A' - 1);

time_scale = int64(1E9); % seconds to nanos
fs = 20;
%% load directories
watch_accel = dir(path + "/mat/watch/accel/*.mat");
watch_gyro = dir(path + "/mat/watch/gyro/*.mat");
phone_accel = dir(path + "/mat/phone/accel/*.mat");
phone_gyro = dir(path + "/mat/phone/gyro/*.mat");

%% check data to ensure subjects line up
for i = 1:numel(watch_accel)
    get_subj = "^data_(\d{4})_[a-z]+_[a-z]+\.mat$";
    subject = regexp(watch_accel(i).name, get_subj, 'tokens');
    match_subj = "^data_" + subject{1}{1} + "_[a-z]+_[a-z]+\.mat$";
    assert(regexp(watch_gyro(i).name, match_subj) == 1, "Subject IDs do not match");
    assert(regexp(phone_accel(i).name, match_subj) == 1, "Subject IDs do not match");
    assert(regexp(phone_gyro(i).name, match_subj) == 1, "Subject IDs do not match");
end

%%
fsamp = 20;
fcuts = [9 9.5];
mags = [1 0];
devs = [0.1 0.05];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
%%
% iterate over subjects
for file_index = 1:numel(watch_accel)
    watch_accel_data = load_data_struct(watch_accel(file_index));
    watch_gyro_data = load_data_struct(watch_gyro(file_index));
    phone_accel_data = load_data_struct(phone_accel(file_index));
    phone_gyro_data = load_data_struct(phone_gyro(file_index));
    
    % iterate over activities
    for activity_index = 1:numel(activities)
        % load activity data for each sensor
        activity = activities(activity_index);
        ds = struct; % "data struct" for holding each sensors data per activity 
        ds.w_acc = watch_accel_data.activity_data.(activity);
        ds.w_gyr = watch_gyro_data.activity_data.(activity);
        ds.p_acc = phone_accel_data.activity_data.(activity);
        ds.p_gyr = phone_gyro_data.activity_data.(activity);
        fn = fieldnames(ds);
        
        % collate data from sensors

        % ASSUMPTION: any 1+ second differences between initial timestamps
        % across sensors is a result of the clocks not lining up and not a
        % result of 1+ second lag between sensor startup time. A lag on the
        % order of miliseconds is considered reasonable and not removed.
        for i = 1:numel(fn)
            % remove the part of offset that is greater than 1 second
            ds.(fn{i}).TimeStampNanos = remove_offset_from_nanos(ds.(fn{i}).TimeStampNanos, time_scale);
        end
        beginnings = [ds.w_acc.TimeStampNanos(1), ds.w_gyr.TimeStampNanos(1), ds.p_acc.TimeStampNanos(1), ds.p_gyr.TimeStampNanos(1)];
        nano_alignments = get_alignment(beginnings, time_scale);
        for i = 1:numel(fn)
            % remove the part of offset that is greater than 1 second
            ds.(fn{i}).TimeStampNanos = ds.(fn{i}).TimeStampNanos + nano_alignments(i);
        end
        
        for i = 1:numel(fn)
            X = xyz_to_mat(ds.(fn{i}));
            [s2, t] = nustft(X, double(ds.(fn{i}).TimeStampNanos)*1E-9, fs, window_time, window_time_shift, max_window_error);
            f = (0:49)/100*20;
            subplot(2,2,i)
            surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
            axis xy; axis tight; view(0,90); c = colorbar;
            c.Label.String = 'Energy (dB)';
            xlabel('Time (s)')
            ylabel('Frequency (Hz)')
            title([fn{i} '_X'], 'Interpreter', 'none')

        end
        sgtitle(['Subject: ', int2str(ds.(fn{1}).SubjectID(1)), ', Activity: ', activities(activity_index)])
        
        % get windowed data for each sensor
        % for i = 1:numel(fn)
        %     % get first 51 XYZ values as a matrix, flip direction, then
        %     % apply FIR filter to get initial values for delays/taps in
        %     % order to avoid border effects
        %     [~, taps] = filter(hh, 1, flip(xyz_to_mat(ds.(fn{i}), 1:n)), [], 2);
        %     ds.(fn{i}).hh_vec = filter(hh, 1, xyz_to_mat(ds.(fn{i})), taps, 2);
        % end
        
        % apply nufft to each set of windowed data
        for i = 1:numel(fn)
            % TODO: this is incredibly inefficient. This should be updated
            % to only evaluate non uniform STFTs where needed
        end
        break
        
    end
    break
end

%% function declarations
function data_struct = load_data_struct(file_struct)
    data_struct = load([file_struct.folder '\' file_struct.name]);
end

function new_nanos = remove_offset_from_nanos(nanos, scale)
    new_nanos = nanos - idivide(nanos(1), scale)*scale;
end


function alignments = get_alignment(beginnings, time_scale)
    % heuristic: choose whatever ordering minimizes the average delay 
    % between beginning samples
    n = length(beginnings);
    % calculate resulting delays for all choices of t_0
    delays = zeros(n,n, 'int64');
    for i = 1:n
        delays(i,:) = beginnings - beginnings(i);
    end
    % add 1 second to negative values to get delay value
    mask = delays<0;
    % select minimum average delay
    delays(mask) = time_scale - delays(mask);
    [~, index] = min(mean(delays,2));
    alignments = zeros(n,n, 'int64') - beginnings(index);
    alignments(mask(index, :)) = alignments(mask(index, :)) + time_scale;
    
end

function matrix = xyz_to_mat(struct, sequence)
    if(~exist('sequence','var'))
        matrix = [struct.X; struct.Y; struct.Z];
    else
        matrix = [struct.X(sequence); struct.Y(sequence); struct.Z(sequence)];
    end
    
end

function [starts, ends] = partition_time_sequence(t, t_win, t_err, nufft_length)
    start_of_sec = zeros(uint32(floor(t(end))), 1, 'uint32');
    starts = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    ends = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    
    % get initial value for start times
    for i = 1:numel(start_of_sec)
        start_of_sec(i) = find(t>=(i-1), 1, 'first');
    end
    % fix sequences to be a multiple of 100 in length
    for i = (t_win+1):(numel(start_of_sec))
        i_start = start_of_sec(i-t_win);
        i_end = start_of_sec(i) - 1;
        seq_len = i_end - i_start + 1;
        mod_nufft_len = mod(seq_len, nufft_length);
        if(mod_nufft_len == 0)
            % the window lined up perfectly
        elseif(seq_len < nufft_length)
            % the window was too short, extend it
            if(i == t_win+1)
                % cannot extend backwards, extend forwards. 
                % 
                % This complicates analysis, but this corner case can happen 
                % only once per activity so it should hopefully not be major.
                % TODO: if the data towards the end gets discard, might as well
                % shift everything forward such that there is enough data
                % before the first window that backwards extension works
                i_end = i_end + (nufft_length-mod_nufft_len);
            else
                % can extend backwards
                i_start = i_start - (nufft_length-mod_nufft_len);
            end
        % if execution reaches here, sequence must be too long
        elseif(t(i_end) - t(i_start+mod_nufft_len) > t_win - t_err)
            % Window can be shortened without making the time covered by the
            % window too short.
            i_start = i_start + mod_nufft_len;
        else
            % the window is too long, but must be further lengthened to cover
            % the minimum desired timespan while having the samples be a factor
            % of 100
            if(i == t_win+1)
                % cannot extend backwards, extend forwards. Same drawbacks as
                % above hold.
                i_end = i_end + (nufft_length-mod_nufft_len);
            else
                % can extend backwards
                i_start = i_start - (nufft_length-mod_nufft_len);
            end
        end
        starts(i-t_win) = i_start;
        ends(i-t_win) = i_end;
        if(i_end - i_start + 1 ~= nufft_length)
            disp('[Start, End]')
            disp([i_start, i_end])
            disp(t(i_end) - t(i_start))
        end
    end
end

% TODO: current, only handles overlap times of multiples of 1 second
function [s2, t] = nustft(x, t, fs, window_time, window_time_shift, max_window_error)
    nufft_length = 50;
    [starts, ends] = partition_time_sequence(t, window_time, window_time * max_window_error, nufft_length*2);
    
    [x_dim, ~] = size(x);
    window_indices = 1:window_time_shift:numel(starts);
    s2 = zeros(nufft_length, x_dim, nufft_length);
    for i = window_indices
        n = ends(i) - starts(i) + 1;
        f = double(0:(n/2-1))/double(n)*fs;
        Y = nufft(x(:,starts(i):ends(i)), t(starts(i):ends(i)), f, 2);
        bin_size = idivide(ends(i)-starts(i)+1, nufft_length*2);
        if(bin_size == 1)
            s2(i,:,:) = abs(Y .* conj(Y));
        else
            disp('Binning!')
            disp(bin_size)
            s2(i,:,:) = sum(reshape(abs(Y.*conj(Y)),bin_size,[], x_dim))./double(bin_size);
        end
    end
    t = t(ends);
end