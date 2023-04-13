%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";
output_path = path + "/preprocessed/wisdm.csv";

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
for file_index = 1:numel(watch_accel)
    disp("Loading")
    watch_accel_data = load_data_struct(watch_accel(file_index));
    watch_gyro_data = load_data_struct(watch_gyro(file_index));
    phone_accel_data = load_data_struct(phone_accel(file_index));
    phone_gyro_data = load_data_struct(phone_gyro(file_index));

    output = table('Size',[175 * length(activities), n_cols],'VariableTypes', col_types, 'VariableNames', col_names);
    % iterate over activities
    disp("Processing")
    for activity_index = 1:numel(activities)
        % load activity data for each sensor
        activity = activities(activity_index);
        
        ds = struct; % "data struct" for holding each sensors data per activity 

        % check that activity data was recorded
        if(~isfield(watch_accel_data.activity_data, activity) ...
                || ~isfield(watch_gyro_data.activity_data, activity) ...
                || ~isfield(phone_accel_data.activity_data, activity) ...
                || ~isfield(phone_gyro_data.activity_data, activity))

            disp("Activity: " + activity + " was not fully record for subject: " + watch_accel_data.subject ...
                + ", skipping")
            break
        end
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
        
        % apply nonuniform STFT
        for i_sens = 1:numel(fn)
            [s2, t] = nustft(xyz_to_mat(ds.(fn{i_sens})), double(ds.(fn{i_sens}).TimeStampNanos)*1E-9, fs, window_time, window_time_shift, max_window_error);
            % save results to table
            
            
            for i_win = 1:numel(t)
                row_index = (activity_index-1)*175 + i_win;
                if(i_sens == 1)
                    % don't bother overwriting if already set
                    output.Time(row_index) = t(i_win);
                    output.SubjectID(row_index) = watch_accel_data.subject;
                    output.Activity(row_index) = activity;
                end
                first_col = 3 + (i_sens - 1)*length(freq_names)*length(components);
                last_col = 3 + i_sens*length(freq_names)*length(components) - 1;
                output{row_index, first_col:last_col} = 20*log10(reshape(squeeze(s2(i_win,:,:)).', 1, n_freq*n_com));
            end
        end      
    end
    % remove empty rows
    disp("Writing")
    if(file_index == 1)
        writetable(output(~ismissing(output.Activity),:), output_path);
    else
        writetable(output(~ismissing(output.Activity),:), output_path, 'WriteMode','Append');
    end
    disp("Processed Subject: " + phone_gyro_data.activity_data.A.SubjectID(1))
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

function [starts, ends, targets] = partition_time_sequence(t, t_win, t_err, nufft_length)
    start_of_sec = zeros(uint32(floor(t(end))), 1, 'uint32');
    starts = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    ends = zeros(length(start_of_sec)- t_win, 1, 'uint32');
    targets = 5:(numel(start_of_sec)-1);
    % get initial value for start times
    for i = 1:numel(start_of_sec)
        start_of_sec(i) = find(t>=(i-1), 1, 'first');
    end
    % fix sequences to be a multiple of 100 in length
    parfor i = (t_win+1):(numel(start_of_sec))
        i_start = start_of_sec(i-t_win);
        i_end = start_of_sec(i) - 1;
        seq_len = i_end - i_start + 1;
        mod_nufft_len = mod(seq_len, nufft_length);
        if(mod_nufft_len == 0)
            % the window lined up perfectly
        elseif(seq_len < nufft_length)
            % the window was too short, extend it
            if(i_start - (nufft_length-mod_nufft_len) <= 0)
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
            if(i_start - (nufft_length-mod_nufft_len) <= 0)
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
    end
end

% TODO: current, only handles overlap times of multiples of 1 second
function [s2, t] = nustft(x, t, fs, window_time, window_time_shift, max_window_error)
    nufft_length = 100;
    [starts, ends, targets] = partition_time_sequence(t, window_time, window_time * max_window_error, nufft_length);
    
    [x_dim, ~] = size(x);
    window_indices = 1:window_time_shift:numel(starts);
    s2 = zeros(nufft_length/2, x_dim, nufft_length/2);
    parfor i = window_indices
        n = ends(i) - starts(i) + 1;
        f = double(0:(n/2-1))/double(n)*fs;
        Y = nufft(x(:,starts(i):ends(i)), t(starts(i):ends(i)), f, 2);
        bin_size = idivide(ends(i)-starts(i)+1, nufft_length);
        if(bin_size == 1)
            s2(i,:,:) = abs(Y .* conj(Y));
        else
            s2(i,:,:) = squeeze(sum(reshape(abs(Y.*conj(Y)),bin_size, x_dim, [])))./double(bin_size);
        end
    end
    t = targets;
end