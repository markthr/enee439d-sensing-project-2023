%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";
output_path = path + "/preprocessed/resamp_with_time/wisdm_";

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
freq_names = 0:(40-1);
t_names = 0:19;
sens_names = {'w_acc', 'w_gyr', 'p_acc', 'p_gyr'};
components = {'x', 'y', 'z'};
n_com = length(components);
n_sens = length(sens_names);
n_freq = length(freq_names);
n_time = length(t_names);
n_cols = 3 + n_com*n_sens*(n_freq + n_time);

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
        offset = 3 + 1 + i_sens*n_com*(n_freq+n_time) + i_com*(n_freq);
        % frequency components
        for i_f = 0:(n_freq-1)
            index = offset + i_f;
            col_names{index} = ['f_', sens_names{i_sens+1}, '_', components{i_com+1}, '_', num2str(freq_names(i_f+1))];
            col_types{index} = 'double';
        end
       
    end
    for i_com = 0:(n_com-1)
        offset = 3 + 1 + i_sens*n_com*(n_freq+n_time) + n_com * n_freq +  i_com*(n_time);
        % time components
        for i_t = 0:(n_time-1)
            index = offset + i_t;
            col_names{index} = ['t_', sens_names{i_sens+1}, '_', components{i_com+1}, '_', num2str(t_names(i_t+1))];
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
            disp("Activity: " + activities(activity_index) + " was not recorded on all sensors for subject: " ... 
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
            col_offset = 3 + (i_sens - 1) * (n_freq+n_time) * n_com;
            % indicates start of frequency components
            first_col = col_offset + 1;
            % indicates end of frequency components, next is start of time components
            mid_col = col_offset + n_freq * n_com;
            % indicates end of time components
            last_col = col_offset + (n_freq+n_time) * n_com;
            
            X = xyz_to_mat(ds.(fn{i_sens}));
            X = gap_aware_resample(X, double(ds.(fn{i_sens}).TimeStampNanos)*1E-9, fs, 2);
            [s, ~,  t] = stft(X, fs, Window=chebwin(100, 35),  FFTLength=100, ...
                OverlapLength=80, FrequencyRange='onesided', OutputTimeDimension = "downrows");
            
            % trim off high frequencies
            s  = s(:, 1:n_freq, :);
            s2 = abs(s).^2; % take magnitude
            % have t point to end of segment instead of middle
            t = t+2.5;


            % save results to table
            t_len = numel(t);
            output.Time(row_offset + (1:t_len)) = t;
            output.SubjectID(row_offset + (1:t_len)) = subject_data.SubjectID;
            output.Activity(row_offset + (1:t_len)) = ds.Activity;
            
            
            % flatten frequency matrix at every time and calculate db
            % permute is used to get the desired interleaving of the
            % frequency matrix. Without permute, matlab flattens the matrix
            % in column major order which would order x_f0, y_f0, z_f0,
            % x_f1... and so on but it is desired to group dimensions
            % together so permuting the dimensions of the 3d array gives
            % row major ordering which keeps components together.
            s2_db = log10(reshape(permute(s2, [1, 3, 2]), [], n_freq*n_com));
            % save result matrix into output, no factor of 20 to improve
            % training

            % saving frequencies
            output{row_offset + (1:size(s2_db, 1)), first_col:mid_col} = s2_db;
            % save last second of samples
            last_seconds = reshape(X(((t(1)-1)*fs+1):(fs*t(end)),:), n_time, [], 3);
            last_seconds = reshape(permute(last_seconds, [2 1 3]), [], n_time*n_com);
            output{row_offset + (1:size(s2_db, 1)), (mid_col+1):last_col} = last_seconds;
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