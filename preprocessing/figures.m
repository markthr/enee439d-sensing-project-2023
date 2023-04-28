%% Vars
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
ds_a = load_activity(subject_data, 'A');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];
sensor = ds_a.w_acc;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(i, 1:1000))
    xlim([0, 50])
    ylabel('m/s^2')
    subtitle(components(i) + " Component of Watch Acceleration")
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
xlabel('Time (s)')

%% Plot Time Data for 29.S
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'S');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];
sensor = ds_a.w_acc;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(i, 1:1000))
    xlim([0, 50])
    ylabel('m/s^2')
    subtitle(components(i) + " Component of Watch Acceleration")
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
xlabel('Time (s)')
%% Plot Time Data for 49.M
subject_data = load_subject(sensor_paths, 50);
ds = load_activity(subject_data, 'M');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = [800 400 800 600];
sensor = ds.w_acc;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t, X(i, :))
    ylabel('m/s^2')
    subtitle(components(i) + " Component of Watch Acceleration")
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')
%% Plot Time Data for 29.A
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'A');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];
sensor = ds_a.p_gyr;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(i, 1:1000))
    xlim([0, 50])
    ylabel('rad/s')
    subtitle(components(i) + " Component of " + sensor_names(4))
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
xlabel('Time (s)')
%% Plot Time Data for 29.S
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'S');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];
sensor = ds_a.p_gyr;
X = xyz_to_mat(sensor);
t = double(sensor.TimeStampNanos)*1E-9;
components = ["X", "Y", "Z"];
for i = 1:3
    subplot(3,1,i)
    plot(t(1:1000), X(i, 1:1000))
    xlim([0, 50])
    ylabel('rad/s^2')
    subtitle(components(i) + " Component of " + sensor_names(4))
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
xlabel('Time (s)')
%% fsst vs stft
subplot(2, 1, 1);
stft(X_hp,fs, 'FrequencyRange', "onesided")
% apply LPF to remove gravity
subplot(2, 1, 2);
fsst(X_hp,fs, 'yaxis');

%% Plot Frequency Time Data for 29.A
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'A');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds_a.(fn{i}));
    t = double(ds_a.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
%% Plot Frequency Time Data for 29.S
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'S');
ds_a = align_sensor_times(ds_a, time_scale);

f = figure;
f.Position = [800 400 800 600];

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds_a.(fn{i}));
    t = double(ds_a.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
%% Plot Frequency Time Data for 30.A
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'A');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds_a.(fn{i}));
    t = double(ds_a.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
%% Plot Frequency Time Data for 30.S
subject_data = load_subject(sensor_paths, 30);
ds_a = load_activity(subject_data, 'S');
ds_a = align_sensor_times(ds_a, time_scale);

fig = figure;
fig.Position = [800 400 800 600];

f = (0:49)/100*20;
for i = 1:numel(fn)
    X = xyz_to_mat(ds_a.(fn{i}));
    t = double(ds_a.(fn{i}).TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds_a.SubjectID +", Activity: " + ds_a.Activity)
%% Plot Frequency Time Data for 49.M
subject_data = load_subject(sensor_paths, 50);
ds = load_activity(subject_data, 'M');
ds = align_sensor_times(ds, time_scale);

fig = figure;
fig.Position = [800 400 800 600];
f = (0:49)/100*20;
for i = 1:numel(fn)
    sensor = ds.(fn{i});
    X = xyz_to_mat(sensor);
    t = double(sensor.TimeStampNanos)*1E-9;
    [s2, t] = nustft(X, t, fs, window_time, window_time_shift, max_window_error);

    subplot(2,2,i)
    surf(t,f,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
    axis xy; axis tight; view(0,90); c = colorbar;
    c.Label.String = 'Energy (dB)';
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    title([fn{i} '_X'], 'Interpreter', 'none')
end
sgtitle("Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
xlabel('Time (s)')
%% make table
tabulate(ds.p_acc.Group)
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
f.Position = [800 400 800 600];

freqs = (0:49)/100*20;
for i = 1:2
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos(400:900))*1E-9;
    [s2, t] = nustft(X(:, 400:900), t, fs, window_time, window_time_shift, max_window_error);
    
    subplot(2,1,i)
    surf(t,freqs,20*log10(squeeze(s2(:,1,:)).'),'EdgeColor','none');   
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
f.Position = [800 400 800 600];

freqs = (0:49)/100*20;
for i = 1:2
    X = xyz_to_mat(ds.(fn{i}));
    t = double(ds.(fn{i}).TimeStampNanos(400:900))*1E-9;
    
    subplot(2,1,i)
    fs= 1/mean(diff(t));
    [s, f, tau] = stft(X(1,400:900), fs, Window=rectwin(100), FrequencyRange="onesided", overlap=80, FFTLength = 100);
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
Y = fft(X(3, 1:100));
L = length(X(3, 1:100));
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
stft(X(3,:), fs, Window=rectwin(100), OverlapLength=20,FFTLength=100, FrequencyRange="onesided")
title("Single-Sided Amplitude Spectrum of Z(t)")
axis tight


sgtitle("Z Component of Watch Acceleration for Subject: " + ds.SubjectID +", Activity: " + ds.Activity)
%%
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
        matrix = [struct.X; struct.Y; struct.Z];
    else
        matrix = [struct.X(sequence); struct.Y(sequence); struct.Z(sequence)];
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