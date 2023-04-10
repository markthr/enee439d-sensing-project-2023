%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";

window_period = 5; % how often to process a new window in seconds
window_memory = 1; % how many previous windows to include, can be a decimal.
% Setting window_memory = 0 indicates that only the current window should be included for processing.

% get letter for the 18 activities, A-S but without N
activities = char([1:13, 15:19] + 'A' - 1);

time_scale = int64(1E9); % seconds to nanos
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
ds = load_data_struct(watch_accel(1));
%% compare FSST with and without HPF
data_a = ds.activity_data.A;
t = double(data_a.TimeStampNanos) * 1E-9;
ts = mean(diff(t));
fs = 1/ts;
X_hp = highpass(data_a.X, 2E-4);
%%
[s,f,t] = fsst(data_a.X,fs);
ax = subplot(2, 1, 1);
waterfall(ax, f,t,abs(s)'.^2)
% apply LPF to remove gravity

[s,f,t] = fsst(X_hp,fs);
ax = subplot(2, 1, 2);
waterfall(ax, f,t,abs(s)'.^2)

%% fsst vs stft
subplot(2, 1, 1);
stft(X_hp,fs, 'FrequencyRange', "onesided")
% apply LPF to remove gravity
subplot(2, 1, 2);
fsst(X_hp,fs, 'yaxis');

%% iterate over subjects
transform_fields = ['X', 'Y', 'Z'];
disp(ds.subject)
for i = 1:numel(activities)
    activity = activities{i};
    disp(activity)
    activity_data = ds.activity_data.(activity);
    
    % calculate average sampling time
    ts_avg = mean(diff(activity_data.TimeStampNanos));
    print(ts)
    w = kaiser(200, 10);
    for j = 1:numel(transform_fields)
        ts_avg = mean(diff(activity_data.TimeStampNanos));
        transform_field = transform_fields(j);
        disp(transform_field)
        
        ts_avg = mean(diff(activity_data.TimeStampNanos));
        w = kaiser(ceil(10E9/ts_avg),10);
        raw = activity_data.(transform_field);
        [s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
    end
    
    
end
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
        
        disp('Begin')
        beginnings = [ds.w_acc.TimeStampNanos(1), ds.w_gyr.TimeStampNanos(1), ds.p_acc.TimeStampNanos(1), ds.p_gyr.TimeStampNanos(1)];
        disp(beginnings)
        endings = [ds.w_acc.TimeStampNanos(end), ds.w_gyr.TimeStampNanos(end), ds.p_acc.TimeStampNanos(end), ds.p_gyr.TimeStampNanos(end)];
        disp('Endings')
        disp(endings)
        disp('Duration')
        disp(double((endings - beginnings))/1000000000)
    end
    break
end
%%
activity_data = watch_accel_data.activity_data.A;
ts_avg = mean(diff(activity_data.TimeStampNanos));
w = kaiser(ceil(10E9/ts_avg),10);
raw = activity_data.Y;
[s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
subplot(2,1,1)
surf(t, f, 20*log10(abs(s).^2), 'EdgeColor', 'none');
activity_data = watch_accel_data.activity_data.B;
ts_avg = mean(diff(activity_data.TimeStampNanos));
raw = activity_data.Y;
[s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
subplot(2,1,2)
surf(t, f, 20*log10(abs(s).^2), 'EdgeColor', 'none');
%%
fs = 1000;
t = 0:1/fs:2;
ts = duration(0,0,1/fs);

x = chirp(t,100,1,200,'quadratic');

fsst(x,ts,'yaxis')

title('Quadratic Chirp')

%%
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