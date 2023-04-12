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

%% make plots using nustft
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
    title
end
%% make table
tabulate(ds.p_acc.Group)
%% plot time delta
seconds_w_acc = double(ds.w_acc.TimeStampNanos(400:900))*1E-9;
plot(401:900, diff(seconds_w_acc))
ylim([0, 0.4])
xlabel('Sample')
ylabel('Time Difference')
%% a
Y = nufft(double(ds.p_acc.hh_vec), double(ds.p_acc.TimeStampNanos), [], 2);
disp(size(Y))
for group_index = 1:numel(ds.p_acc.Group(end))
    unfft(ds.p_gyr.hh_vec(i,:))
    
end

%%

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
