path = "J:/enee439d/datasets/wisdm-dataset";

watch_accel = dir(path + "/mat/watch/accel/*.mat");

%%
data_struct = load([watch_accel(1).folder '\' watch_accel(1).name]);
%% compare FSST with and without HPF
data_a = data_struct.activity_data.A;
t = double(data_a.TimeStampNanos) * 1E-9;
ts = mean(diff(t));
fs = 1/ts;
[s,f,t] = fsst(data_a.X,fs);
ax = subplot(2, 1, 1);
waterfall(ax, f,t,abs(s)'.^2)
% apply LPF to remove gravity
X_hp = highpass(data_a.X, 2E-4);
[s,f,t] = fsst(X_hp,fs);
ax = subplot(2, 1, 2);
waterfall(ax, f,t,abs(s)'.^2)

%% fsst vs stft
subplot(2, 1, 1);
stft(X_hp,fs, 'FrequencyRange', "onesided")
% apply LPF to remove gravity
subplot(2, 1, 2);
[s, f, t] = fsst(X_hp,fs, 'yaxis');

%% iterate over fields and apply FFT
transform_fields = ['X', 'Y', 'Z'];
activities = fieldnames(data_struct.activity_data);

disp(data_struct.subject)
for i = 1:numel(activities)
    activity = activities{i};
    disp(activity)
    activity_data = data_struct.activity_data.(activity);
    % calculate average sampling time
    ts_avg = mean(diff(activity_data.TimeStampNanos));
    w = kaiaser(200, 10);
    for j = 1:numel(transform_fields)
        ts_avg = mean(diff(activity_data.TimeStampNanos));
        transform_field = transform_fields(j);
        disp(transform_field)
        
        ts_avg = mean(diff(activity_data.TimeStampNanos));
        w = kaiser(ceil(10E9/ts_avg),10);
        raw = activity_data.(transform_field);
        [s, f, t] = stft(raw, 1E9/ts_avg, Window=w, FrequencyRange="onesided", overlap=199);
        data
    end
    
    
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