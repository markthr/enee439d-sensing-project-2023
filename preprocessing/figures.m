%% Vars
path = "J:/enee439d/datasets/wisdm-dataset";

window_time = 5; % how often to process a new window in seconds
max_window_error = 0.1; % the maximum percent a window can be too short by 
window_time_shift = 1;
nufft_length = 100;
fig_pos = [400 200 800 600];
% Setting window_memory = 0 indicates that only the current window should be included for processing.

% get letter for the 18 activities, A-S but without N
activities = char([1:13, 15:19] + 'A' - 1);

time_scale = int64(1E9); % seconds to nanos
fs = 20;
%% Load directories and check sorting
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

%% Plot Time Data for 29.A
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'A');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;
sensor = ds.w_acc;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(1:1000, i))
    xlim([0, 50])
    ylabel('m/s^2')
    subtitle(components(i) + " Component of Watch Acceleration")
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')

%% Plot Time Data for 29.S
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;
sensor = ds.w_acc;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(1:1000, i))
    xlim([0, 50])
    ylabel('m/s^2')
    subtitle(components(i) + " Component of Watch Acceleration")
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')
%% Plot Time Data for 49.M
subject_data = load_subject(sensor_paths, 50);
ds = load_activity(subject_data, 'M');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;
sensor = ds.w_acc;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t, X(:, i))
    ylabel('m/s^2')
    subtitle(components(i) + " Component of Watch Acceleration")
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')
%% Plot Time Data for 29.A
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'A');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;
sensor = ds.p_gyr;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(1:1000, i))
    xlim([0, 50])
    ylabel('rad/s')
    subtitle(components(i) + " Component of " + sensor_names(4))
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')
%% Plot Time Data for 29.S
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;
sensor = ds.p_gyr;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(1:1000, i))
    xlim([0, 50])
    ylabel('rad/s^2')
    subtitle(components(i) + " Component of " + sensor_names(4))
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')

%% Plot Frequency Time Data for 29.A
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'A');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%% Plot Frequency Time Data for 29.S
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);

f = figure;
f.Position = fig_pos;

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%% Plot Frequency Time Data for 30.A
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'A');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%% Plot Frequency Time Data for 30.S
subject_data = load_subject(sensor_paths, 30);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%% Plot Frequency Time Data for 49.M
subject_data = load_subject(sensor_paths, 50);
ds = load_activity(subject_data, 'M');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = fig_pos;
f = (0:49)/100*20;
for i = 1:numel(fn)
    sensor = ds.(fn{i});
    X = xyz_to_mat(sensor);
    t = double(sensor.TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);

    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')
%% plot time delta
subject_data = load_subject(sensor_paths, 1);
ds = load_activity(subject_data, 'S');
seconds_w_acc = double(ds.w_acc.TimeStampNanos(400:900))*1E-9;
plot(401:900, diff(seconds_w_acc))
ylim([0, 0.4])
xlabel('Sample')
ylabel('Time Difference')
%% NUSTFT for nonuniform data
subject_data = load_subject(sensor_paths, 1);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);
f = figure;
f.Position = fig_pos;

freqs = (0:49)/100*20;
for i = 1:2
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos(400:900))*1E-9;
    [s2, t] = nustft(X(400:900, :), t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,1,i)
    surf(t,freqs,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("NUSTFT: Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%% Plot STFT for nonuniform data
subject_data = load_subject(sensor_paths, 1);
ds = load_activity(subject_data, 'S');
ds = align_sensor_times(ds, time_scale);
f = figure;
f.Position = fig_pos;

freqs = (0:49)/100*20;
for i = 1:2
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos(400:900))*1E-9;
    
    subplot(2,1,i)
    fs= 1/mean(diff(t));
    [s, f, tau] = stft(X(400:900,1), fs, Window=rectwin(100), FrequencyRange="onesided", overlap=80, FFTLength = 100);
    s2 = abs(s).^2;
    surf(t(int32(tau*fs)),f,20*log10(s2),'EdgeColor','none');
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("STFT: Subject: " + ds.SubjectID +", Activity: " + ds.Activity)

%%
subject_data = load_subject(sensor_paths, 15);
ds = load_activity(subject_data, 'R');
ds = align_sensor_times(ds, time_scale);
fig = figure;
fig.Position = [400 200 1200 600];

X = xyz_to_mat(ds.w_acc);
seconds_w_acc = double(ds.w_acc.TimeStampNanos(1:end))*1E-9;
subplot(2,2,1)
plot(2:length(seconds_w_acc), diff(seconds_w_acc) * 1.0E3 )
ylim([0, 0.4])
xlabel('Sample')
ylabel('Time Difference (ms)')
xlim([0 length(seconds_w_acc)])
ylim([49 50])
title("Sampling Time")


subplot(2,2,3)
Y = fft(X(1:100, 3));
L = length(X(1:100, 3));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
fs = 1/(mean(diff(ds.w_acc.TimeStampNanos))*1E-9);
f =  fs*(0:(L/2))/L;
plot(f,20*log10(P1.^2)) 
title("Single-Sided FFT of Z(t) for First 5 Seconds")
xlabel("frequency (Hz)")
ylabel("Energy (dB)")
axis tight

subplot(2,2,[2 4])
stft(X(:,3), fs, Window=rectwin(100), OverlapLength=20,FFTLength=100, FrequencyRange="onesided")
title("Single-Sided Amplitude Spectrum of Z(t)")
axis tight


sgtitle("Z Component of Watch Acceleration for Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%%
for i = 1:numel(fn)
    X = xyz_to_mat(ds.(fn{i}));
    [s2, t] = nustft(X, double(ds.(fn{i}).TimeStampNanos)*1E-9, fs, window_time, window_time_shift, max_window_error);
    f = (0:49)/100*20;
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,:,1)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')

end
sgtitle(['Subject: ', int2str(ds.(fn{1}).SubjectID(1)), ', Activity: ', ds.Activity])

%% function declarations
function subject_data_struct = load_subject(sensor_paths_struct, subject_id)
    fn = fieldnames(sensor_paths_struct);
    subject_data_struct = struct;
    for i = 1:numel(fn)
        sensor_path = sensor_paths_struct.(fn{i});
        file_struct = sensor_path(subject_id);
        subject_data_struct.(fn{i}) = load([file_struct.folder '\' file_struct.name]);
    end
end

function activity_data_struct = load_activity(subject_data_struct, activity)
    fn = fieldnames(subject_data_struct);
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



