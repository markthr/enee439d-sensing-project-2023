%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";
output_path = path + "/preprocessed/wisdm_";

window_time = 5; % how often to process a new window in seconds
max_window_error = 0.1; % the maximum percent a window can be too short by 
window_time_shift = 1;
nufft_length = 100;
% Setting window_memory = 0 indicates that only the current window should be included for processing.

% get letter for the 18 activities, A-S but without N
activities = char([1:13, 15:19] + 'A' - 1);
n_act= length(activities);

time_scale = int64(1E9); % seconds to nanos
fs = 20;
%% load directories
sensor_paths = struct;
sensor_paths.w_acc = dir(path + "/mat/watch/accel/*.mat");
sensor_paths.w_gyr  = dir(path + "/mat/watch/gyro/*.mat");
sensor_paths.p_acc = dir(path + "/mat/phone/accel/*.mat");
sensor_paths.p_gyr = dir(path + "/mat/phone/gyro/*.mat");

sensor_names = ["Watch Acceleration", "Watch Gyro", "Phone Acceleration", "Phone Gyro"];
fn = fieldnames(sensor_paths);
n_subj = length(sensor_paths.(fn{1}));
% check data to ensure subjects line up
get_subj = "^data_(\d{4})_[a-z]+_[a-z]+\.mat$";
for i = 1:numel(sensor_paths.(fn{1}))
    subject = regexp(sensor_paths.(fn{1})(i).name, get_subj, 'tokens');

    for j = 2:numel(fn)
        subject_o = regexp(sensor_paths.(fn{1})(i).name, get_subj, 'tokens');
        assert(strcmp(cell2mat(subject{1}), cell2mat(subject_o{1})), "Subject IDs do not match");
    end
end

%%
fsamp = 20;
fcuts = [9 9.5];
mags = [1 0];
devs = [0.1 0.05];

[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp);
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
%%
frequencies = (0:49)/100*20;
% table of 175 rows, 4 sensors which each have 3 readings
% sensor readings are then passed into a 50 point NUFFT
% add extra columns for timestamp, subject, activity
freq_names = 0:(50-1);
sens_names = {'w_acc', 'w_gyr', 'p_acc', 'p_gyr'};
components = {'x', 'y', 'z'};
n_com = length(components);
n_sens = length(sens_names);
n_freq = length(freq_names);
n_cols = 3 + n_com*n_sens*n_freq;

col_types = cell(1, n_cols);
col_names = cell(1, n_cols);

col_names{1} = 'SubjectID';
col_types{1} = 'int32';
col_names{2} = 'Activity';
col_types{2} = 'string';
col_names{3} = 'Time';
col_types{3} = 'double';

% column order: iterate over sensor, then field, and then frequency
for i_sens = 0:(n_sens-1)
    for i_com = 0:(n_com-1)
        for i_f = 0:(n_freq-1)
            index = 3 + 1 + i_sens*n_com*n_freq + i_com*n_freq + i_f;
            col_names{index} = [sens_names{i_sens+1}, '_', components{i_com+1}, '_', num2str(freq_names(i_f+1))];
            col_types{index} = 'double';
        end
    end
end
%%
% iterate over subjects

f = waitbar(0);
for file_index = 1:n_subj
    % load data for subject
    subject_data = load_subject(sensor_paths, file_index);
    fn = subject_data.sensors;

    waitbar(0, f, "Processing data for subject: " + subject_data.SubjectID);
    % table for collecting data after preprocessing
    output = table('Size',[180 * length(activities), n_cols],'VariableTypes', col_types, 'VariableNames', col_names);

    % iterate over activities
    row_offset = 0;
    for activity_index = 1:n_act
        waitbar(activity_index/n_act, f);
        
        % check that activity data was recorded for all sensors
        missing = false;
        for i = 1:numel(fn)
            missing = missing || ~isfield(subject_data.(fn{i}).activity_data, activities(activity_index));
        end
        if(missing)
            disp("Activity: " + activities(activity_index) + " was recorded on all sensors for subject: " ... 
                + subject_data.SubjectID + ", skipping")
            continue
        end

        % load activity data for subject
        ds = load_activity(subject_data, activities(activity_index));
        % collate data from sensors
        % ASSUMPTION: any 1+ second differences between initial timestamps
        % across sensors is a result of the clocks not lining up and not a
        % result of 1+ second lag between sensor startup time. A lag on the
        % order of miliseconds is considered reasonable and not removed.
        ds = align_sensor_times(ds, time_scale);
        
        % apply nonuniform STFT
        time_lengths = zeros(1, 4);
        for i_sens = 1:numel(fn)
            [s2, t] = nustft(xyz_to_mat(ds.(fn{i_sens})), single(ds.(fn{i_sens}).TimeStampNanos)*1E-9, fs, window_time, window_time_shift, max_window_error);
            % save results to table
            t_len = numel(t);
            output.Time(row_offset + (1:t_len)) = t;
            output.SubjectID(row_offset + (1:t_len)) = subject_data.SubjectID;
            output.Activity(row_offset + (1:t_len)) = ds.Activity;
            first_col = 3 + (i_sens - 1)*length(freq_names)*length(components) + 1;
            last_col = 3 + i_sens*length(freq_names)*length(components);
            
            % flatten frequency matrix at every time and calculate db
            % permute is used to get the desired interleaving of the
            % frequency matrix. Without permute, matlab flattens the matrix
            % in column major order which would order x_f0, y_f0, z_f0,
            % x_f1... and so on but it is desired to group dimensions
            % together so permuting the dimensions of the 3d array gives
            % row major ordering which keeps components together.
            s2_db = log10(reshape(permute(s2, [1 3 2]), [], 150));
            % save result matrix into output, no factor of 20 to improve
            % training
            output{row_offset + (1:size(s2_db, 1)), first_col:last_col} = s2_db;
            % save time length for determining next row offset
            time_lengths(i_sens) = t_len;
        end
        row_offset = row_offset + min(time_lengths);
    end
    % remove empty rows
    writetable(output(1:row_offset,:), output_path ...
        + ds.SubjectID + ".csv");
end
close(f);

%% function declarations
function subject_data_struct = load_subject(sensor_paths_struct, subject_id)
    fn = fieldnames(sensor_paths_struct);
    subject_data_struct = struct;
    for i = 1:numel(fn)
        sensor_path = sensor_paths_struct.(fn{i});
        file_struct = sensor_path(subject_id);
        subject_data_struct.(fn{i}) = load([file_struct.folder '\' file_struct.name]);
    end

    subject_data_struct.SubjectID = subject_data_struct.(fn{1}).subject;
    subject_data_struct.sensors = fn;
end

function activity_data_struct = load_activity(subject_data_struct, activity)
    fn = subject_data_struct.sensors;
    activity_data_struct = struct;
    for i = 1:numel(fn)
        sensor = subject_data_struct.(fn{i});
        activity_data_struct.(fn{i}) = sensor.activity_data.(activity);
    end
    activity_data_struct.SubjectID = subject_data_struct.(fn{1}).subject;
    activity_data_struct.Activity = activity;
    activity_data_struct.sensors = fn;
end

% simple helper function that removes any component of t_0 greater than 1
% second from a series of times in nanoseconds
function new_nanos = remove_offset_from_nanos(nanos, scale)
    new_nanos = nanos - idivide(nanos(1), scale)*scale;
end

function activity_data_struct = align_sensor_times(activity_data_struct, time_scale)
    % heuristic: assume remove any offset greater than 1 second in inital
    % times
    % heuristic: choose whatever ordering minimizes the average delay 
    % between beginning samples
    fn = activity_data_struct.sensors;
    n = numel(fn);

    % save t_0 of each sensor for determing which to place first
    beginnings = zeros(1,4, 'int64');
    % remove the part of offset that is greater than 1 second
    for i = 1:n
        activity_data_struct.(fn{i}).TimeStampNanos = ...
            remove_offset_from_nanos(activity_data_struct.(fn{i}).TimeStampNanos, time_scale);
        beginnings(i) = activity_data_struct.(fn{i}).TimeStampNanos(1);
    end
    
    % calculate resulting delays for all choices of t_0
    delays = zeros(n,n, 'int64');
    for i = 1:n
        delays(i,:) = beginnings - beginnings(i);
    end
    % add 1 second to negative values to get delay value
    mask = delays<0;
    delays(mask) = time_scale - delays(mask);
    % select minimum average delay
    [~, index] = min(mean(delays,2));
    % generate offset to align the data
    alignments = zeros(1, n, 'int64') - beginnings(index);
    alignments(mask(index, :)) = alignments(mask(index, :)) + time_scale;

    % update time for each sensor
    for i = 1:n
        activity_data_struct.(fn{i}).TimeStampNanos = ...
            activity_data_struct.(fn{i}).TimeStampNanos - alignments(i);
    end
end

function matrix = xyz_to_mat(struct, sequence)
    if(~exist('sequence','var'))
        matrix = [struct.X, struct.Y, struct.Z];
    else
        matrix = [struct.X(sequence), struct.Y(sequence), struct.Z(sequence)];
    end
    
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