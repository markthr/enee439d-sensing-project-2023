%% Setup
path = "J:/enee439d/datasets/wisdm-dataset";

window_period = 5; % how often to process a new window in seconds
window_memory = 1; % how many previous windows to include, can be a decimal.
% Setting window_memory = 0 indicates that only the current window should be
% included for processing.

% get letter for the 18 activities, A-S but without N
activities = char([1:13, 15:19] + 'A' - 1);

time_scale = 1E9; % seconds to nanos
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
data_struct = load_data_struct(watch_accel(1));
%% compare FSST with and without HPF
data_a = data_struct.activity_data.A;
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
disp(data_struct.subject)
for i = 1:numel(activities)
    activity = activities{i};
    disp(activity)
    activity_data = data_struct.activity_data.(activity);
    
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
        activity = activities(activity_index);
        w_acc = watch_accel_data.activity_data.(activity);
        w_gyr = watch_gyro_data.activity_data.(activity);
        p_acc = phone_accel_data.activity_data.(activity);
        p_gyr = phone_gyro_data.activity_data.(activity);
        sensors = [w_acc, w_gyr, p_acc, p_gyr];
        disp('Begin')
        beginnings = [w_acc.TimeStampNanos(1), w_gyr.TimeStampNanos(1), p_acc.TimeStampNanos(1), p_gyr.TimeStampNanos(1)];
        t_0 = min(beginnings);
        for i = 1:numel(sensors)
            sensors(i).TimeStampNanos = sensors(i).TimeStampNanos - t_0;
            disp(sensors(i).TimeStampNanos(1))
        end
        disp(beginnings)
        disp('End')
        beginnings = [w_acc.TimeStampNanos(1), w_gyr.TimeStampNanos(1), p_acc.TimeStampNanos(1), p_gyr.TimeStampNanos(1)];
        endings = [w_acc.TimeStampNanos(end), w_gyr.TimeStampNanos(end), p_acc.TimeStampNanos(end), p_gyr.TimeStampNanos(end)];
        disp(endings)
        disp('Duration')
        disp(double((endings - beginnings))/1000000000)
    end
    break
end
%%
transform_field = 'Y';
ts_avg = mean(diff(activity_data.TimeStampNanos));
w = kaiser(ceil(10E9/ts_avg),10);
raw = activity_data.(transform_field);
[s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
surf(t, f, 20*log10(abs(s).^2), 'EdgeColor', 'none');
%%
fs = 1000;
t = 0:1/fs:2;
ts = duration(0,0,1/fs);

x = chirp(t,100,1,200,'quadratic');

fsst(x,ts,'yaxis')

title('Quadratic Chirp')

function data_struct = load_data_struct(file_struct)
    data_struct = load([file_struct.folder '\' file_struct.name]);
end